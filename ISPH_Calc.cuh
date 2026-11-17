//-------------------------------------------------------------------------------------------------
// Implicit-PPE file-scope state (only used when Parameters.cuh IMPLICIT_PPE == 1).
// Lazily allocated inside SOPHIA_single_ISPH on first call (count==0).
// Lives across all SOPHIA_single_ISPH calls; freed implicitly at process exit.
//-------------------------------------------------------------------------------------------------
#if IMPLICIT_PPE
static Real*  g_dev_pres_old      = NULL;   // snapshot of P1[i].pres before each Jacobi sweep
static Real*  g_dev_pres_delta    = NULL;   // per-particle |p_new - p_old|, input to Max reduce
static Real*  g_dev_bi            = NULL;   // per-particle |bi| = |rho*div(u*)/dt|, written by PPE_prep
static Real*  g_dev_ppe_resid     = NULL;   // per-particle |A*p - b| = |PPE2*pres - PPE1|, true residual
static Real*  g_d_pres_delta_max  = NULL;   // device scalar: max Jacobi iter delta
static Real*  g_d_bi_max          = NULL;   // device scalar: max |bi|
static Real*  g_d_resid_max       = NULL;   // device scalar: max |A*p - b|
static void*  g_dev_resid_storage = NULL;   // CUB Max reduce temp storage (shared by all 3 reduces)
static size_t g_resid_storage_bytes = 0;
static FILE*  g_fp_ppe_log        = NULL;   // ./<CASE_NAME>/PPE_iter_log.txt
#endif

// Mean-field one-shot dump guard (Run A). Set true after the snapshot is
// written at time >= AVG_END_TIME so later timesteps don't overwrite it.
static bool g_mean_field_dumped = false;

// Host shadows of the device __device__ symbols avg_start_dyn / avg_end_dyn.
// tid==0 owns writes; updates are mirrored to device via cudaMemcpyToSymbol
// and read by other tids after pthread_barrier_wait.
static volatile Real h_avg_start = (Real)1e30;
static volatile Real h_avg_end   = (Real)1e30;

// Quasi-steady detection ring buffer (Σ½|u|² history, owned by tid==0).
static Real g_qs_uu_buf[QS_DETECT_N];
static int  g_qs_uu_fill = 0;
static int  g_qs_uu_head = 0;

// count value at the moment h_avg_start was latched. -1 = not yet.
// Used to derive count_eff = count - count_avg so that output cadence
// (freq_output, LDM 100-step) can be reset relative to the averaging start.
static int_t count_avg = -1;

#if PPE_DIAG
// PPE residual diagnostic state (Parameters.cuh PPE_DIAG). Lazily allocated by
// tid==0 on the first logged step; device-0 arrays + CSV owned by tid==0.
static Real*  g_diag_bi            = NULL;   // per-particle |rho*div(u*)/dt|
static Real*  g_diag_resid         = NULL;   // per-particle |A*p - b|
static Real*  g_d_diag_bi_max      = NULL;   // device scalar: max|bi|
static Real*  g_d_diag_resid_max   = NULL;   // device scalar: max|A*p - b|
static Real*  g_d_diag_sum         = NULL;   // device scalar: Sum (for div(u^{n+1}) average)
static void*  g_diag_storage       = NULL;   // CUB temp storage (sized for Max & Sum)
static size_t g_diag_storage_bytes = 0;
static int_t  g_n_fluid            = -1;     // cached fluid-particle count (for averages)
static FILE*  g_fp_diag            = NULL;   // ./<OUTPUT_CASE_NAME>/PPE_diag_{explicit,implicit}.csv
#endif

#if CONSISTENT_PPE
// IISPH solver state (Parameters.cuh CONSISTENT_PPE). Single-GPU (file-scope static).
static Real*  g_iisph_resid        = NULL;   // per-particle |s_i - Ap_i| (density residual)
static Real*  g_iisph_d_sum        = NULL;   // device scalar: sum of residual
static void*  g_iisph_storage      = NULL;   // CUB Sum temp storage
static size_t g_iisph_storage_bytes= 0;
static FILE*  g_fp_iisph           = NULL;   // ./<OUTPUT_CASE_NAME>/IISPH_log.txt
static FILE*  g_fp_conv            = NULL;   // ./<OUTPUT_CASE_NAME>/IISPH_conv.txt (per-step iter curve)
#endif

// true when this process resumed from a checkpoint: per-step logs are opened
// in append mode so the pre-crash rows are kept (header only if the file is new).
static bool g_restarted = false;
static bool ckpt_log_is_new(FILE*fp){ fseek(fp, 0, SEEK_END); return ftell(fp)==0; }

//-------------------------------------------------------------------------------------------------
// Checkpoint / restart (Parameters.cuh CKPT_EVERY, CKPT_PATH; resume with env SOPHIA_RESTART=<path>)
//
// Everything carried across steps (single GPU):
//   dev_P1, dev_P2               particle state (P3/SP1/SP2 are rebuilt at every step start)
//   g_str/g_end                  fluid cell lists — built ONCE at count==0 and never again
//   dev_LP1, dev_LP2             LDM particles (incl. per-particle RNG seed)
//   g_str_LDM/g_end_LDM          LDM cell lists (rebuilt every freq_LDM steps; saved anyway)
//   device globals               us_update, OB_length_global, avg_count, avg_start_dyn, avg_end_dyn
//   host statics                 time, count, h_avg_start/end, count_avg, g_mean_field_dumped, QS buffer
// File: [ckpt_header][P1][P2][g_str][g_end][LP1][LP2][g_str_LDM][g_end_LDM][qs_buf][end magic]
// Written to <path>.tmp, fsync'd, previous file kept as <path>.prev, then renamed into place.
//-------------------------------------------------------------------------------------------------
#include <unistd.h>

typedef struct ckpt_header {
	char         magic[16];            // "SOPHIA_CKPT_v1"
	int          sz_part1, sz_part2, sz_lpart1, sz_lpart2;
	int          n_part2, n_part_LDM, n_cells, qs_n;
	int          step;                 // count after the saved step (next step to run)
	Real         sim_time, step_dt;
	float        us_update, OB_length_global;
	unsigned int avg_count;
	Real         avg_start_dyn, avg_end_dyn;
	Real         h_avg_start, h_avg_end;
	int          count_avg;
	int          mean_field_dumped;
	int          qs_fill, qs_head;
	int          pad[8];
} ckpt_header;

static const char CKPT_MAGIC[16]     = "SOPHIA_CKPT_v1";
static const char CKPT_END_MAGIC[16] = "SOPHIA_CKPT_END";

static bool ckpt_write(FILE*fp, const void*buf, size_t n){ return fwrite(buf,1,n,fp)==n; }
static bool ckpt_read (FILE*fp,       void*buf, size_t n){ return fread (buf,1,n,fp)==n; }

// Host buffers file_P1/file_P2/file_LP1/file_LP2 are reused as staging (outputs recopy before use).
static bool save_checkpoint(const char*path,
	part1*dev_P1, part2*dev_P2, int_t*g_str, int_t*g_end,
	L_part1*dev_LP1, L_part2*dev_LP2, int_t*g_str_LDM, int_t*g_end_LDM,
	part1*file_P1, part2*file_P2, L_part1*file_LP1, L_part2*file_LP2)
{
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);

	size_t n_cell_bytes = sizeof(int_t)*num_cells;
	int_t*h_gs  =(int_t*)malloc(n_cell_bytes), *h_ge  =(int_t*)malloc(n_cell_bytes);
	int_t*h_gsL =(int_t*)malloc(n_cell_bytes), *h_geL =(int_t*)malloc(n_cell_bytes);

	cudaMemcpy(file_P1,  dev_P1,  sizeof(part1)*num_part2,      cudaMemcpyDeviceToHost);
	cudaMemcpy(file_P2,  dev_P2,  sizeof(part2)*num_part2,      cudaMemcpyDeviceToHost);
	cudaMemcpy(file_LP1, dev_LP1, sizeof(L_part1)*num_part_LDM, cudaMemcpyDeviceToHost);
	cudaMemcpy(file_LP2, dev_LP2, sizeof(L_part2)*num_part_LDM, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_gs,  g_str,     n_cell_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_ge,  g_end,     n_cell_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_gsL, g_str_LDM, n_cell_bytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_geL, g_end_LDM, n_cell_bytes, cudaMemcpyDeviceToHost);

	ckpt_header hdr;
	memset(&hdr, 0, sizeof(hdr));
	memcpy(hdr.magic, CKPT_MAGIC, sizeof(CKPT_MAGIC));
	hdr.sz_part1 = sizeof(part1);   hdr.sz_part2  = sizeof(part2);
	hdr.sz_lpart1= sizeof(L_part1); hdr.sz_lpart2 = sizeof(L_part2);
	hdr.n_part2  = num_part2; hdr.n_part_LDM = num_part_LDM; hdr.n_cells = num_cells; hdr.qs_n = QS_DETECT_N;
	hdr.step     = count; hdr.sim_time = time; hdr.step_dt = dt;
	cudaMemcpyFromSymbol(&hdr.us_update,        us_update,        sizeof(float));
	cudaMemcpyFromSymbol(&hdr.OB_length_global, OB_length_global, sizeof(float));
	cudaMemcpyFromSymbol(&hdr.avg_count,        avg_count,        sizeof(unsigned int));
	cudaMemcpyFromSymbol(&hdr.avg_start_dyn,    avg_start_dyn,    sizeof(Real));
	cudaMemcpyFromSymbol(&hdr.avg_end_dyn,      avg_end_dyn,      sizeof(Real));
	hdr.h_avg_start = h_avg_start; hdr.h_avg_end = h_avg_end;
	hdr.count_avg = count_avg; hdr.mean_field_dumped = g_mean_field_dumped ? 1 : 0;
	hdr.qs_fill = g_qs_uu_fill; hdr.qs_head = g_qs_uu_head;

	char tmp[512], prev[512];
	snprintf(tmp,  sizeof(tmp),  "%s.tmp",  path);
	snprintf(prev, sizeof(prev), "%s.prev", path);

	bool ok = false;
	FILE*fp = fopen(tmp, "wb");
	if(fp){
		ok = ckpt_write(fp, &hdr, sizeof(hdr))
		  && ckpt_write(fp, file_P1,  sizeof(part1)*num_part2)
		  && ckpt_write(fp, file_P2,  sizeof(part2)*num_part2)
		  && ckpt_write(fp, h_gs,  n_cell_bytes) && ckpt_write(fp, h_ge,  n_cell_bytes)
		  && ckpt_write(fp, file_LP1, sizeof(L_part1)*num_part_LDM)
		  && ckpt_write(fp, file_LP2, sizeof(L_part2)*num_part_LDM)
		  && ckpt_write(fp, h_gsL, n_cell_bytes) && ckpt_write(fp, h_geL, n_cell_bytes)
		  && ckpt_write(fp, g_qs_uu_buf, sizeof(Real)*QS_DETECT_N)
		  && ckpt_write(fp, CKPT_END_MAGIC, sizeof(CKPT_END_MAGIC));
		ok = (fflush(fp)==0) && ok;
		ok = (fsync(fileno(fp))==0) && ok;
		ok = (fclose(fp)==0) && ok;
	}
	free(h_gs); free(h_ge); free(h_gsL); free(h_geL);

	if(ok){
		rename(path, prev);              // keep one previous generation (fails harmlessly on the first save)
		ok = (rename(tmp, path)==0);
	}
	clock_gettime(CLOCK_MONOTONIC, &t1);
	double wall = (t1.tv_sec-t0.tv_sec) + 1e-9*(t1.tv_nsec-t0.tv_nsec);
	if(ok) printf("[ckpt] saved count=%d t=%.4f avg_count=%u -> %s (%.1f s)\n",
	              (int)count, (double)time, hdr.avg_count, path, wall);
	else   printf("[ckpt] ERROR: checkpoint write failed (%s); previous checkpoint left untouched\n", tmp);
	fflush(stdout);
	return ok;
}

// Called after device arrays are allocated/initialised and the averaging window is set,
// right before the main loop. Overwrites the initial state; returns false on any mismatch.
static bool load_checkpoint(const char*path,
	part1*dev_P1, part2*dev_P2, int_t*g_str, int_t*g_end,
	L_part1*dev_LP1, L_part2*dev_LP2, int_t*g_str_LDM, int_t*g_end_LDM,
	part1*file_P1, part2*file_P2, L_part1*file_LP1, L_part2*file_LP2)
{
	FILE*fp = fopen(path, "rb");
	if(!fp){ printf("[ckpt] FATAL: cannot open %s\n", path); return false; }

	ckpt_header hdr;
	if(!ckpt_read(fp, &hdr, sizeof(hdr)) || memcmp(hdr.magic, CKPT_MAGIC, sizeof(CKPT_MAGIC))!=0){
		printf("[ckpt] FATAL: %s is not a SOPHIA checkpoint\n", path); fclose(fp); return false;
	}
	if(hdr.sz_part1!=(int)sizeof(part1) || hdr.sz_part2!=(int)sizeof(part2) ||
	   hdr.sz_lpart1!=(int)sizeof(L_part1) || hdr.sz_lpart2!=(int)sizeof(L_part2) ||
	   hdr.n_part2!=num_part2 || hdr.n_part_LDM!=num_part_LDM || hdr.n_cells!=num_cells ||
	   hdr.qs_n!=QS_DETECT_N || hdr.step_dt!=dt){
		printf("[ckpt] FATAL: checkpoint does not match this build/input "
		       "(part1 %d/%d part2 %d/%d L_part1 %d/%d L_part2 %d/%d n_part2 %d/%d n_LDM %d/%d cells %d/%d qs %d/%d dt %g/%g)\n",
		       hdr.sz_part1,(int)sizeof(part1), hdr.sz_part2,(int)sizeof(part2),
		       hdr.sz_lpart1,(int)sizeof(L_part1), hdr.sz_lpart2,(int)sizeof(L_part2),
		       hdr.n_part2,(int)num_part2, hdr.n_part_LDM,(int)num_part_LDM, hdr.n_cells,(int)num_cells,
		       hdr.qs_n,(int)QS_DETECT_N, (double)hdr.step_dt,(double)dt);
		fclose(fp); return false;
	}

	size_t n_cell_bytes = sizeof(int_t)*num_cells;
	int_t*h_gs  =(int_t*)malloc(n_cell_bytes), *h_ge  =(int_t*)malloc(n_cell_bytes);
	int_t*h_gsL =(int_t*)malloc(n_cell_bytes), *h_geL =(int_t*)malloc(n_cell_bytes);
	char end_magic[16];
	bool ok = ckpt_read(fp, file_P1,  sizeof(part1)*num_part2)
	       && ckpt_read(fp, file_P2,  sizeof(part2)*num_part2)
	       && ckpt_read(fp, h_gs,  n_cell_bytes) && ckpt_read(fp, h_ge,  n_cell_bytes)
	       && ckpt_read(fp, file_LP1, sizeof(L_part1)*num_part_LDM)
	       && ckpt_read(fp, file_LP2, sizeof(L_part2)*num_part_LDM)
	       && ckpt_read(fp, h_gsL, n_cell_bytes) && ckpt_read(fp, h_geL, n_cell_bytes)
	       && ckpt_read(fp, g_qs_uu_buf, sizeof(Real)*QS_DETECT_N)
	       && ckpt_read(fp, end_magic, sizeof(end_magic))
	       && memcmp(end_magic, CKPT_END_MAGIC, sizeof(CKPT_END_MAGIC))==0;
	fclose(fp);
	if(!ok){
		printf("[ckpt] FATAL: %s is truncated or corrupt\n", path);
		free(h_gs); free(h_ge); free(h_gsL); free(h_geL);
		return false;
	}

	cudaMemcpy(dev_P1,  file_P1,  sizeof(part1)*num_part2,      cudaMemcpyHostToDevice);
	cudaMemcpy(dev_P2,  file_P2,  sizeof(part2)*num_part2,      cudaMemcpyHostToDevice);
	cudaMemcpy(dev_LP1, file_LP1, sizeof(L_part1)*num_part_LDM, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_LP2, file_LP2, sizeof(L_part2)*num_part_LDM, cudaMemcpyHostToDevice);
	cudaMemcpy(g_str,     h_gs,  n_cell_bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(g_end,     h_ge,  n_cell_bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(g_str_LDM, h_gsL, n_cell_bytes, cudaMemcpyHostToDevice);
	cudaMemcpy(g_end_LDM, h_geL, n_cell_bytes, cudaMemcpyHostToDevice);
	free(h_gs); free(h_ge); free(h_gsL); free(h_geL);

	cudaMemcpyToSymbol(us_update,        &hdr.us_update,        sizeof(float));
	cudaMemcpyToSymbol(OB_length_global, &hdr.OB_length_global, sizeof(float));
	cudaMemcpyToSymbol(avg_count,        &hdr.avg_count,        sizeof(unsigned int));
	cudaMemcpyToSymbol(avg_start_dyn,    &hdr.avg_start_dyn,    sizeof(Real));
	cudaMemcpyToSymbol(avg_end_dyn,      &hdr.avg_end_dyn,      sizeof(Real));
	time = hdr.sim_time; count = hdr.step;
	h_avg_start = hdr.h_avg_start; h_avg_end = hdr.h_avg_end;
	count_avg = hdr.count_avg; g_mean_field_dumped = (hdr.mean_field_dumped!=0);
	g_qs_uu_fill = hdr.qs_fill; g_qs_uu_head = hdr.qs_head;
	g_restarted = true;

	printf("[ckpt] RESTART from %s: count=%d t=%.4f avg_count=%u h_avg=[%.3f, %.3f) count_avg=%d mean_dumped=%d\n",
	       path, (int)count, (double)time, hdr.avg_count, (double)h_avg_start, (double)h_avg_end,
	       (int)count_avg, hdr.mean_field_dumped);
	fflush(stdout);
	return true;
}

//-------------------------------------------------------------------------------------------------
// 모든 입자 계산을 수행하는 주함수 + LDM
//-------------------------------------------------------------------------------------------------
void SOPHIA_single_ISPH(int_t*g_idx,int_t*p_idx,int_t*g_idx_in,int_t*p_idx_in,int_t*g_str,int_t*g_end,
	int_t*g_idx_LDM,int_t*p_idx_LDM,int_t*g_idx_in_LDM,int_t*p_idx_in_LDM,int_t*g_str_LDM,int_t*g_end_LDM,
	int_t*b_idx,int_t*b_idx_in,int_t*b_str,int_t*b_end,
	part1*dev_P1,part1*dev_SP1,part2*dev_P2,part2*dev_SP2,part3*dev_P3,
	L_part1*dev_LP1,L_part1*dev_SLP1,L_part2*dev_LP2,L_part2*dev_SLP2,
	int_t*p2p_af_in,int_t*p2p_idx_in,int_t*p2p_af,int_t*p2p_idx,
	void*dev_sort_storage,size_t*sort_storage_bytes,
	void*dev_sort_storage_LDM,size_t*sort_storage_bytes_LDM,
	part1*file_P1,part2*file_P2,part3*file_P3,
	L_part1*file_LP1,
	int tid)
	//int*apr_num_part,int_t*count_buffer,int_t*dev_count_buffer,int_t*num_buffer_temp,int_t*APR_cell,int_t*dev_APR_cell,Real*computing_time)
{
	dim3 b,t;
	t.x=128;
	b.x=(num_part2-1)/t.x+1;
	int s=sizeof(int)*(t.x+1);

//-------------------------------------------------------------------------------------------------
// 주변입자 검색
//-------------------------------------------------------------------------------------------------

	KERNEL_set_ncell<<<b,t>>>(dev_P1);
	cudaDeviceSynchronize();

	if(count==0){
		// g_str을 리셋
		cudaMemset(g_str,cu_memset,sizeof(int_t)*num_cells);

		// 입자의 셀번호 계산
		b.x=(num_part2-1)/t.x+1;
		KERNEL_index_particle_to_cell<<<b,t>>>(g_idx_in,p_idx_in,dev_P1);
		cudaDeviceSynchronize();

		// 셀번호를 바탕으로 정렬
		cub::DeviceRadixSort::SortPairs(dev_sort_storage,*sort_storage_bytes,g_idx_in,g_idx,p_idx_in,p_idx,num_part2);
		cudaDeviceSynchronize();

		// 정렬한 입자를 재배치
		b.x=(num_part2-1)/t.x+1;
		KERNEL_reorder<<<b,t,s>>>(g_idx,p_idx,g_str,g_end,dev_P1,dev_P2,dev_SP1,dev_SP2);
		cudaDeviceSynchronize();

		// 일부 입자정보 리셋
		cudaMemset(dev_P3,0.0,sizeof(part3)*num_part2);

		// 입자정보를 P1 에 복사
		cudaMemcpy(dev_P1,dev_SP1,sizeof(part1)*num_part2,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);

		// 완화거리 가변 모델 사용
		// if(dim==3) KERNEL_variable_smoothing_length3D<<<b,t>>>(g_str,g_end,dev_SP1);
		// cudaDeviceSynchronize();

#if LDM_MEAN_MODE
		// Run B: load mean field AFTER the first NNPS sort so the record
		// index in mean_field.bin (saved in Run A's sorted order) matches
		// the current sorted dev_P1 indexing. Also keep dev_SP1 in sync.
		{
			int ok = load_mean_field(dev_P1, num_part2, MEAN_FIELD_PATH);
			if(!ok){
				printf("[mean] FATAL: load_mean_field failed (tid=%d). Aborting.\n", tid);
			}
			cudaMemcpy(dev_SP1, dev_P1, sizeof(part1)*num_part2, cudaMemcpyDeviceToDevice);
		}
		pthread_barrier_wait(&barrier);
#endif

	}
	else{
#if !LDM_MEAN_MODE
		// 0스텝 이후 단순 복사 >> 입자정보를 SP1, SP2 에 복사
		cudaMemcpy(dev_SP1,dev_P1,sizeof(part1)*num_part2,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
		cudaMemcpy(dev_SP2,dev_P2,sizeof(part2)*num_part2,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);

		// 일부 입자정보 리셋
		cudaMemset(dev_P3,0.0,sizeof(part3)*num_part2);
		pthread_barrier_wait(&barrier);
#endif
		// 평균장 모드(LDM_MEAN_MODE)에서는 SPH 블록을 전부 건너뛰므로 P1 이 변하지 않는다.
		// 위 복사·초기화는 스텝당 ~20 GB 의 순수 낭비라 생략한다 (SP1/SP2 는 count==0 값 유지).
	}

	// LDM 주변입자 검색 (농도 계산)
	// Gate: LDM dormant until quasi-steady reached (h_avg_start latched).
	// release / advection / NNPS all skipped while h_avg_start == ∞.
	if(time >= h_avg_start){
	if(count%freq_LDM==0){

		b.x=(num_part_LDM-1)/t.x+1;
		particle_releasing_Dispersion_Kernel<<<b,t>>>(dt,freq_LDM,time,time_end,dev_LP1);

		// g_str을 리셋
		cudaMemset(g_str_LDM,cu_memset,sizeof(int_t)*num_cells);

		// 입자의 셀번호 계산
		b.x=(num_part_LDM-1)/t.x+1;
		KERNEL_index_particle_to_cell_LDM<<<b,t>>>(g_idx_in_LDM,p_idx_in_LDM,dev_LP1);
		cudaDeviceSynchronize();

		// 셀번호를 바탕으로 정렬
		cub::DeviceRadixSort::SortPairs(dev_sort_storage_LDM,*sort_storage_bytes_LDM,g_idx_in_LDM,g_idx_LDM,p_idx_in_LDM,p_idx_LDM,num_part_LDM);
		cudaDeviceSynchronize();

		// 정렬한 입자를 재배치
		b.x=(num_part_LDM-1)/t.x+1;
		KERNEL_reorder_LDM<<<b,t,s>>>(g_idx_LDM,p_idx_LDM,g_str_LDM,g_end_LDM,dev_LP1,dev_LP2,dev_SLP1,dev_SLP2);
		cudaDeviceSynchronize();

		// 입자정보를 P1 에 복사
		cudaMemcpy(dev_LP1,dev_SLP1,sizeof(L_part1)*num_part_LDM,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
	}
	else{
		// LDM 입자 배열 업데이트
		cudaMemcpy(dev_SLP1,dev_LP1,sizeof(L_part1)*num_part_LDM,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
		cudaMemcpy(dev_SLP2,dev_LP2,sizeof(L_part2)*num_part_LDM,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
	}
	}  // end of h_avg_start gate (LDM NNPS block)

//-------------------------------------------------------------------------------------------------
// 계산 준비: gradient correction, filter, reference density, p_type switch, penetration, density gradient
//-------------------------------------------------------------------------------------------------
#if !LDM_MEAN_MODE
	b.x=(num_part2-1)/t.x+1;
#if INFLOW_TURB
	KERNEL_inflow_turb3D<<<b,t>>>(dev_SP1,time);     // 유입 버퍼 = U(z) + 합성 변동
	cudaDeviceSynchronize();
#endif
	// 미분보정필터 계산 >> 진우 코드랑 비교 필요
	gradient_correction(g_str,g_end,dev_SP1,dev_P3);
	cudaDeviceSynchronize();
	KERNEL_clc_gradient_correction_3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3);
	cudaDeviceSynchronize();

	// filter, reference density, p_type switch, penetration, normal gradient etc
	KERNEL_clc_prep3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3,count,dt);
	cudaDeviceSynchronize();

	// standard k-epsilon (Launder-Spalding) — used for validation against DDF
	if(turbulence_model==K_E){
		KERNEL_clc_ke_turb3D<<<b,t>>>(1,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
		KERNEL_update_turbulence<<<b,t>>>(dt,dev_SP1,dev_P3);
		cudaDeviceSynchronize();
		KERNEL_turb_viscosity<<<b,t>>>(dev_SP1,dev_P3);
		cudaDeviceSynchronize();
	}

	// Velocity condition for wall
	b.x=(num_part2-1)/t.x+1;
	if(noslip_bc==1){
		KERNEL_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1);
		cudaDeviceSynchronize();
	}

	if(MOST_boundary){
		b.x=(num_part2-1)/t.x+1;
		KERNEL_MOST_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3,count,dt,1);
		cudaDeviceSynchronize();
		KERNEL_MOST_band_reduce<<<1,1>>>();
		cudaDeviceSynchronize();
	}

	if(count==0){
		KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}

//-------------------------------------------------------------------------------------------------
// 압력힘을 제외한 힘 계산
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	KERNEL_advection_force3D<<<b,t>>>(1,g_str,g_end,dev_SP1,dev_SP2,dev_P3,count);
	cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// PREDICTOR (Optional)
//-------------------------------------------------------------------------------------------------
	// Predictor(for projection method)
	if(time_type==Pre_Cor){
		b.x=(num_part2-1)/t.x+1;
		KERNEL_clc_projection<<<b,t>>>(count,dt,time,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}

//-------------------------------------------------------------------------------------------------
// Immersed Boundary Method (IBM)
//-------------------------------------------------------------------------------------------------
	//IBM force interpolation (fb for IB)
	b.x=(num_part2-1)/t.x+1;
	IBM_force_interpolation3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// IBM spreading interpolation (fb for fluid)
	IBM_spreading_interpolation3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	IBM_corrector<<<b,t>>>(dt_structure,time,dev_SP1,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// Direct forcing: u** = u* + fb [m/s], PPE는 u**의 div로 풀이.
	IBM_corrector_fluid<<<b,t>>>(dt,time,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// Zero velocity/temp for interior fluid particles (inside buildings) before PPE
	// KERNEL_zero_interior_velocity<<<b,t>>>(dev_SP1);
	// cudaDeviceSynchronize();


//-------------------------------------------------------------------------------------------------
// 압력 계산
//-------------------------------------------------------------------------------------------------
	//For initial condition(Neumann boundary)

#if CONSISTENT_PPE
	//--- IISPH consistent-operator PPE: prep once, relaxed-Jacobi sweeps with early-exit ---
	// A == div(applied symmetric pressure force) so converging this is a true projection.
	// Stop when avg density error (mean|s_i-Ap_i|/rho0) < IISPH_ETA, capped at maxIteration.
	b.x=(num_part2-1)/t.x+1;

	// Lazy init: residual buffer + Sum-reduce storage + log (once, on first call)
	if(g_iisph_resid==NULL){
		cudaMalloc((void**)&g_iisph_resid, sizeof(Real)*num_part2);
		cudaMalloc((void**)&g_iisph_d_sum, sizeof(Real));
		cudaMemset(g_iisph_resid, 0, sizeof(Real)*num_part2);
		cub::DeviceReduce::Sum(g_iisph_storage, g_iisph_storage_bytes, g_iisph_resid, g_iisph_d_sum, num_part2);
		cudaMalloc(&g_iisph_storage, g_iisph_storage_bytes);
		if(tid==0){
			char fn[512]; sprintf(fn, "./%s/IISPH_log.txt", OUTPUT_CASE_NAME);
			g_fp_iisph=fopen(fn,g_restarted?"a":"w");
			if(g_fp_iisph && ckpt_log_is_new(g_fp_iisph)){ fprintf(g_fp_iisph,"# step\ttime\titer_used\tavg_dens_err\tconverged\n"); fflush(g_fp_iisph); }
			char fc[512]; sprintf(fc, "./%s/IISPH_conv.txt", OUTPUT_CASE_NAME);
			g_fp_conv=fopen(fc,g_restarted?"a":"w");
			if(g_fp_conv && ckpt_log_is_new(g_fp_conv)){ fprintf(g_fp_conv,"# step\titer\tavg_div\n"); fflush(g_fp_conv); }
		}
	}

	KERNEL_iisph_prep<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	{
		int_t nsweep=(maxIteration>0)?maxIteration:1;
		int_t minit =(minIteration>0)?minIteration:1;
		// On conv-log steps run the FULL sweep count (no early-exit) and record the per-iteration
		// div curve, so we can see div(u^{n+1}) decreasing within a single step.
		bool conv_log_step = (tid==0) && (IISPH_CONV_FREQ>0) && (count % (IISPH_CONV_FREQ>0?IISPH_CONV_FREQ:1) == 0) && (g_fp_conv!=NULL);
		int_t used=0; Real avg_err=0.0f; Real prev_err=0.0f; int conv=0;
		for(int_t it=0; it<nsweep; it++){
			cudaMemsetAsync(g_iisph_resid, 0, sizeof(Real)*num_part2);
			KERNEL_iisph_pass1<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
			cudaDeviceSynchronize();
			KERNEL_iisph_pass2<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_iisph_resid);
			cudaDeviceSynchronize();
			KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
			cudaDeviceSynchronize();

			// Early-exit gate in DIVERGENCE units (not density): the IISPH residual is
			// r_i = s_i-Ap_i = -dt*rho*div(u^{n+1}), so |rho*div(u^{n+1})/dt| = |r_i|/dt^2.
			// This all-particle mean (non-fluid=0) is the SAME quantity as div_un1_avg, so the
			// gate now measures the actual incompressibility, not the misleading density error.
			Real rsum=0.0f;
			cub::DeviceReduce::Sum(g_iisph_storage, g_iisph_storage_bytes, g_iisph_resid, g_iisph_d_sum, num_part2);
			cudaDeviceSynchronize();
			cudaMemcpy(&rsum, g_iisph_d_sum, sizeof(Real), cudaMemcpyDeviceToHost);
			avg_err = rsum/((Real)num_part2*dt*dt);
			used=it+1;
			if(conv_log_step) fprintf(g_fp_conv,"%d\t%d\t%.6e\n",(int)count,(int)used,(double)avg_err);
			// Convergence: absolute floor (avg_err<ETA) OR stagnation (this step's floor reached,
			// per-sweep relative drop < STAG_TOL). Stagnation prevents the 200-cap spin when the
			// developed-flow floor rises above IISPH_ETA.
			Real rel_drop = (it>0)? (prev_err-avg_err)/(prev_err+1e-20f) : 1.0f;
			if(!conv_log_step && used>=minit && (avg_err<IISPH_ETA || rel_drop<IISPH_STAG_TOL)){ conv=1; break; }
			prev_err=avg_err;
		}
		if(conv_log_step){ conv=(avg_err<IISPH_ETA); fflush(g_fp_conv); }
		if(tid==0){
			if(g_fp_iisph){ fprintf(g_fp_iisph,"%d\t%.6g\t%d\t%.6e\t%d\n",(int)count,(double)time,(int)used,(double)avg_err,conv); fflush(g_fp_iisph); }
			if((count%50)==0) printf("IISPH: iter=%d/%d  avg_dens_err=%.3e  conv=%d\n",(int)used,(int)nsweep,(double)avg_err,conv);
		}
	}
#else
#if IMPLICIT_PPE
	//--- Implicit Jacobi PPE: iterate up to maxIteration, early-exit on max|delta p| < tol ---
	// Tracks 3 metrics per step:
	//   - Jacobi delta:    max|p^{k+1} - p^k|  (iteration convergence indicator, early-exit gate)
	//   - max|bi|:         max|rho*div(u*)/dt| (input incompressibility error from predictor)
	//   - max|A*p - b|:    true PPE residual   (proportional to |div(u^{n+1})|, the real metric)
	// NOTE: single-GPU only. The g_dev_* buffers are file-scope static (shared across pthreads),
	// so multi-GPU would race/cross-device. Run the implicit comparison with 1 GPU.
	b.x=(num_part2-1)/t.x+1;

	// Lazy init (once, on first SOPHIA_single_ISPH call)
	if(g_dev_pres_old==NULL){
		cudaMalloc((void**)&g_dev_pres_old,    sizeof(Real)*num_part2);
		cudaMalloc((void**)&g_dev_pres_delta,  sizeof(Real)*num_part2);
		cudaMalloc((void**)&g_dev_bi,          sizeof(Real)*num_part2);
		cudaMalloc((void**)&g_dev_ppe_resid,   sizeof(Real)*num_part2);
		cudaMalloc((void**)&g_d_pres_delta_max,sizeof(Real));
		cudaMalloc((void**)&g_d_bi_max,        sizeof(Real));
		cudaMalloc((void**)&g_d_resid_max,     sizeof(Real));
		cudaMemset(g_dev_pres_old,   0, sizeof(Real)*num_part2);
		cudaMemset(g_dev_pres_delta, 0, sizeof(Real)*num_part2);
		cudaMemset(g_dev_bi,         0, sizeof(Real)*num_part2);
		cudaMemset(g_dev_ppe_resid,  0, sizeof(Real)*num_part2);
		// CUB Max-reduce temp storage (sized via dummy call — same array length for all 3 reduces)
		cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
		                       g_dev_pres_delta, g_d_pres_delta_max, num_part2);
		cudaMalloc(&g_dev_resid_storage, g_resid_storage_bytes);
		if(tid==0){
			char fn[512];
			sprintf(fn, "./%s/PPE_iter_log.txt", OUTPUT_CASE_NAME);
			g_fp_ppe_log = fopen(fn, g_restarted?"a":"w");
			if(g_fp_ppe_log && ckpt_log_is_new(g_fp_ppe_log)){
				fprintf(g_fp_ppe_log, "# IMPLICIT_PPE Jacobi solver log\n");
				fprintf(g_fp_ppe_log, "# tolerance=%g  minIter=%d  maxIter=%d\n",
				        (double)PPE_TOLERANCE, (int)minIteration, (int)maxIteration);
				fprintf(g_fp_ppe_log, "# bi_max     = |rho*div(u*)/dt|_inf   (input divergence; same every iter)\n");
				fprintf(g_fp_ppe_log, "# resid_max  = |A*p_final - b|_inf    (true PPE residual; ~|div(u^{n+1})|*rho/dt)\n");
				fprintf(g_fp_ppe_log, "# delta_p    = max|p^{k+1} - p^k|     (Jacobi iteration delta; early-exit gate)\n");
				fprintf(g_fp_ppe_log, "# step\ttime\titer_used\tdelta_p\tbi_max\tresid_max\tconverged\n");
				fflush(g_fp_ppe_log);
			}
		}
	}

	// Snapshot pres -> pres_old so the first delta is well-defined
	KERNEL_snapshot_pres<<<b,t>>>(dev_SP1, g_dev_pres_old);
	cudaDeviceSynchronize();

	int_t ppe_iter_used = 0;
	Real  ppe_delta_p   = 0.0f;
	Real  ppe_bi_max    = 0.0f;
	Real  ppe_resid_max = 0.0f;
	int   ppe_converged = 0;
	int_t maxit_safe = (maxIteration>0) ? maxIteration : 1;
	int_t minit_safe = (minIteration>0) ? minIteration : 1;

	for(int_t it=0; it<maxit_safe; it++){
		// Zero bi array so non-fluid slots don't pollute Max reduce (PPE_prep skips them via early return)
		cudaMemsetAsync(g_dev_bi, 0, sizeof(Real)*num_part2);
		KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_dev_bi);
		cudaDeviceSynchronize();
		KERNEL_PPE_calc3D<<<b,t>>>(dev_SP1);
		cudaDeviceSynchronize();
		KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();

		// Compute per-particle |p_new - p_old|, snapshot p -> pres_old
		KERNEL_compute_pres_delta<<<b,t>>>(dev_SP1, g_dev_pres_old, g_dev_pres_delta);
		cudaDeviceSynchronize();

		// Max reduce -> device scalar, copy to host
		cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
		                       g_dev_pres_delta, g_d_pres_delta_max, num_part2);
		cudaDeviceSynchronize();
		cudaMemcpy(&ppe_delta_p, g_d_pres_delta_max, sizeof(Real), cudaMemcpyDeviceToHost);

		ppe_iter_used = it+1;
		if(ppe_iter_used >= minit_safe && ppe_delta_p < PPE_TOLERANCE){
			ppe_converged = 1;
			break;
		}
	}

	// After loop: one extra PPE_prep with converged pres to refresh PPE1/PPE2,
	// then compute true residual |A*p - b|. bi is also refreshed (same value as last iter — u* unchanged).
	cudaMemsetAsync(g_dev_bi, 0, sizeof(Real)*num_part2);
	KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_dev_bi);
	cudaDeviceSynchronize();
	KERNEL_compute_ppe_residual<<<b,t>>>(dev_SP1, g_dev_ppe_resid);
	cudaDeviceSynchronize();

	cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
	                       g_dev_bi, g_d_bi_max, num_part2);
	cudaDeviceSynchronize();
	cudaMemcpy(&ppe_bi_max, g_d_bi_max, sizeof(Real), cudaMemcpyDeviceToHost);

	cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
	                       g_dev_ppe_resid, g_d_resid_max, num_part2);
	cudaDeviceSynchronize();
	cudaMemcpy(&ppe_resid_max, g_d_resid_max, sizeof(Real), cudaMemcpyDeviceToHost);

	// Per-step logging (only tid==0)
	if(tid==0){
		if(g_fp_ppe_log){
			fprintf(g_fp_ppe_log, "%d\t%.6g\t%d\t%.6e\t%.6e\t%.6e\t%d\n",
			        (int)count, (double)time, (int)ppe_iter_used,
			        (double)ppe_delta_p, (double)ppe_bi_max,
			        (double)ppe_resid_max, ppe_converged);
			fflush(g_fp_ppe_log);
		}
		if((count%50)==0){
			printf("PPE(implicit): iter=%d / %d, dp=%.3e, |bi|=%.3e, |A*p-b|=%.3e, conv=%d\n",
			       (int)ppe_iter_used, (int)maxit_safe,
			       (double)ppe_delta_p, (double)ppe_bi_max,
			       (double)ppe_resid_max, ppe_converged);
		}
	}
#else
	//--- EISPH (Explicit ISPH, 1-sweep Jacobi: current production) ---
	for(int_t it=0;it<=0;it++){   // For iterative(implicit) pressure solver
		// PPE.cuh
		b.x=(num_part2-1)/t.x+1;
		KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, NULL);
		cudaDeviceSynchronize();
		KERNEL_PPE_calc3D<<<b,t>>>(dev_SP1);
		cudaDeviceSynchronize();
		// Pressure condition for wall(Neumann boundary)
		KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}
#endif
#endif  // CONSISTENT_PPE

	// (PPE_DIAG moved AFTER the corrector so div(u^{n+1}) can be measured — see below.)

//-------------------------------------------------------------------------------------------------
// 압력힘 계산
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	KERNEL_pressureforce3D<<<b,t>>>(1,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// 시간 적분 (Time Integration)
//-------------------------------------------------------------------------------------------------
	b.x=(num_part2-1)/t.x+1;
	KERNEL_time_update_projection<<<b,t>>>(dt,dev_SP1,dev_P1,dev_SP2,dev_P2,dev_P3,count);
	cudaDeviceSynchronize();

	// Interior particles skip time update (buffer_type==99 early return in KERNEL_time_update_projection)
	// → dev_P1 retains zeroed state; no separate zeroing needed here

#if PPE_DIAG
	//--- Incompressibility diagnostic (post-corrector, so u^{n+1} exists) every PPE_DIAG_FREQ steps ---
	// bi_max      = max|rho*div(u*)/dt|       (INPUT divergence, before pressure; from dev_SP1)
	// resid_max   = max|A_brookshaw*p - b|     (Brookshaw residual; meaningful only for Brookshaw runs)
	// div_un1_max = max|rho*div(u^{n+1})/dt|   (OUTPUT divergence; solver-agnostic real metric)
	// div_un1_avg = fluid-mean of the same     (geometry dev_SP1, corrected velocity dev_P1)
	if(tid==0 && (count % PPE_DIAG_FREQ == 0)){
		b.x=(num_part2-1)/t.x+1;

		if(g_diag_bi==NULL){
			cudaMalloc((void**)&g_diag_bi,          sizeof(Real)*num_part2);
			cudaMalloc((void**)&g_diag_resid,       sizeof(Real)*num_part2);
			cudaMalloc((void**)&g_d_diag_bi_max,    sizeof(Real));
			cudaMalloc((void**)&g_d_diag_resid_max, sizeof(Real));
			cudaMalloc((void**)&g_d_diag_sum,       sizeof(Real));
			size_t bmax=0,bsum=0;
			cub::DeviceReduce::Max(NULL, bmax, g_diag_bi, g_d_diag_bi_max, num_part2);
			cub::DeviceReduce::Sum(NULL, bsum, g_diag_bi, g_d_diag_sum,    num_part2);
			g_diag_storage_bytes=(bmax>bsum)?bmax:bsum;
			cudaMalloc(&g_diag_storage, g_diag_storage_bytes);

			char fn[512];
#if IMPLICIT_PPE
			sprintf(fn, "./%s/PPE_diag_implicit.csv", OUTPUT_CASE_NAME);
#else
			sprintf(fn, "./%s/PPE_diag_explicit.csv", OUTPUT_CASE_NAME);
#endif
			g_fp_diag=fopen(fn,g_restarted?"a":"w");
			if(g_fp_diag && ckpt_log_is_new(g_fp_diag)){ fprintf(g_fp_diag,"step,time,bi_max,resid_max,div_un1_max,div_un1_avg,div_un1_plain_max,div_un1_plain_avg\n"); fflush(g_fp_diag); }

			// Cache fluid count for the average
			cudaMemsetAsync(g_diag_resid, 0, sizeof(Real)*num_part2);
			KERNEL_fluid_mask<<<b,t>>>(dev_SP1, g_diag_resid);
			cudaDeviceSynchronize();
			Real fsum=0.0f;
			cub::DeviceReduce::Sum(g_diag_storage, g_diag_storage_bytes, g_diag_resid, g_d_diag_sum, num_part2);
			cudaDeviceSynchronize();
			cudaMemcpy(&fsum, g_d_diag_sum, sizeof(Real), cudaMemcpyDeviceToHost);
			g_n_fluid=(int_t)(fsum+0.5f); if(g_n_fluid<1) g_n_fluid=1;
		}

		// bi_max (input divergence, u* in dev_SP1). Also refreshes dev_SP1 PPE1/PPE2 for resid.
		cudaMemsetAsync(g_diag_bi, 0, sizeof(Real)*num_part2);
		KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_diag_bi);
		cudaDeviceSynchronize();
		Real diag_bi_max=0.0f;
		cub::DeviceReduce::Max(g_diag_storage, g_diag_storage_bytes, g_diag_bi, g_d_diag_bi_max, num_part2);
		cudaDeviceSynchronize();
		cudaMemcpy(&diag_bi_max, g_d_diag_bi_max, sizeof(Real), cudaMemcpyDeviceToHost);

		// resid_max (Brookshaw residual)
		cudaMemsetAsync(g_diag_resid, 0, sizeof(Real)*num_part2);
		KERNEL_compute_ppe_residual<<<b,t>>>(dev_SP1, g_diag_resid);
		cudaDeviceSynchronize();
		Real diag_resid_max=0.0f;
		cub::DeviceReduce::Max(g_diag_storage, g_diag_storage_bytes, g_diag_resid, g_d_diag_resid_max, num_part2);
		cudaDeviceSynchronize();
		cudaMemcpy(&diag_resid_max, g_d_diag_resid_max, sizeof(Real), cudaMemcpyDeviceToHost);

		// div(u^{n+1}) measured TWO ways for a fair comparison: KGC (explicit operator) and
		// plain gradW (IISPH operator). Geometry dev_SP1, corrected velocity dev_P1.
		Real dmax=0.0f,dsum=0.0f, dmax_p=0.0f,dsum_p=0.0f;
		// (a) KGC
		cudaMemsetAsync(g_diag_bi, 0, sizeof(Real)*num_part2);
		KERNEL_compute_div_un1<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_P1,dev_SP2,dev_P3, g_diag_bi, 1);
		cudaDeviceSynchronize();
		cub::DeviceReduce::Max(g_diag_storage, g_diag_storage_bytes, g_diag_bi, g_d_diag_bi_max, num_part2);
		cudaDeviceSynchronize(); cudaMemcpy(&dmax, g_d_diag_bi_max, sizeof(Real), cudaMemcpyDeviceToHost);
		cub::DeviceReduce::Sum(g_diag_storage, g_diag_storage_bytes, g_diag_bi, g_d_diag_sum, num_part2);
		cudaDeviceSynchronize(); cudaMemcpy(&dsum, g_d_diag_sum, sizeof(Real), cudaMemcpyDeviceToHost);
		// (b) plain gradW
		cudaMemsetAsync(g_diag_bi, 0, sizeof(Real)*num_part2);
		KERNEL_compute_div_un1<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_P1,dev_SP2,dev_P3, g_diag_bi, 0);
		cudaDeviceSynchronize();
		cub::DeviceReduce::Max(g_diag_storage, g_diag_storage_bytes, g_diag_bi, g_d_diag_bi_max, num_part2);
		cudaDeviceSynchronize(); cudaMemcpy(&dmax_p, g_d_diag_bi_max, sizeof(Real), cudaMemcpyDeviceToHost);
		cub::DeviceReduce::Sum(g_diag_storage, g_diag_storage_bytes, g_diag_bi, g_d_diag_sum, num_part2);
		cudaDeviceSynchronize(); cudaMemcpy(&dsum_p, g_d_diag_sum, sizeof(Real), cudaMemcpyDeviceToHost);
		Real div_un1_avg=dsum/(Real)g_n_fluid, div_un1_plain_avg=dsum_p/(Real)g_n_fluid;

		if(g_fp_diag){
			fprintf(g_fp_diag,"%d,%.6g,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
			        (int)count,(double)time,(double)diag_bi_max,(double)diag_resid_max,
			        (double)dmax,(double)div_un1_avg,(double)dmax_p,(double)div_un1_plain_avg);
			fflush(g_fp_diag);
		}
	}
#endif

	// // Monitor target particle (once per step, after all forces applied)
	// KERNEL_monitor_particle<<<b,t>>>(dev_P1, dev_P3, count, 4);
	// cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// Open Boudnary && APR
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
 	// KERNEL_Neumann_boundary_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	//KERNEL_donothing_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	// KERNEL_open_boundary_extrapolation_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	// Outlet (buffer_type==2): Orlanski convective outflow + pres=0 anchor.
	KERNEL_convective_outflow3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3,dt);
	cudaDeviceSynchronize();
	KERNEL_periodic_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	KERNEL_freeslip_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// Rayleigh damping sponge near the ceiling (uz only). Absorbs grid-scale
	// vertical oscillations trapped under reflective freeslip ceiling.
	if((Real)SPONGE_GAMMA_MAX > (Real)0.0){
		KERNEL_sponge_ceiling<<<b,t>>>(dev_P1, dt);
		cudaDeviceSynchronize();
	}

// #if CASE_ID >= 500 && CASE_ID < 600
// 	// Case CS: sustain inlet buffer at log-law profile (uses us_update; with MOST off it stays at us_ini)
// 	if(dim==3) KERNEL_Recycling3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
// 	cudaDeviceSynchronize();
// #endif

	KERNEL_compute_vorticity3D<<<b,t>>>(g_str,g_end,dev_P1,dev_P3);
	cudaDeviceSynchronize();
#endif  // !LDM_MEAN_MODE  (end of SPH solver block)

#if !LDM_MEAN_MODE
	// Time-averaging accumulation, only within [h_avg_start, h_avg_end).
	// Both bounds default to ∞ until QS detection latches them.
	// Time-averaging accumulation: [h_avg_start, h_avg_end) 구간에서만 누적.
	// AVG_START_TIME<=0 이면 h_avg_start=0, h_avg_end=inf 라 매 스텝 누적되지만
	// h_avg_end 가 무한이라 mean field dump 는 일어나지 않는다(=기존 transient 동작).
	if(time >= h_avg_start && time < h_avg_end){
		b.x=(num_part2-1)/t.x+1;
		KERNEL_accumulate_avg<<<b,t>>>(dev_P1);
		cudaDeviceSynchronize();
	}

	// One-shot mean field dump when the latched window closes (Run A →
	// produces MEAN_FIELD_PATH consumed by Run B with LDM_MEAN_MODE=1).
	if(h_avg_end < (Real)1e29 && time >= h_avg_end && !g_mean_field_dumped){
		if(tid==0){
			cudaMemcpy(file_P1, dev_P1, num_part2*sizeof(part1), cudaMemcpyDeviceToHost);
			save_mean_field(file_P1, num_part2, MEAN_FIELD_PATH);
		}
		pthread_barrier_wait(&barrier);
		g_mean_field_dumped = true;
	}
#endif  // !LDM_MEAN_MODE

//-------------------------------------------------------------------------------------------------
// LDM coupling — gated by h_avg_start (skip predictor/advection/deposition/conc until quasi-steady)
//-------------------------------------------------------------------------------------------------

	if(time >= h_avg_start && (count%freq_LDM==0)){
		// b.x=(num_part_LDM-1)/t.x+1;
		// particle_releasing_Dispersion<<<b,t>>>(dt,freq_LDM,time,time_end,dev_LP1);
		if(time_type==Pre_Cor){
			b.x=(num_part_LDM-1)/t.x+1;
			KERNEL_clc_predictor_Dispersion<<<b,t>>>(freq_LDM*dt,time,dev_LP1,dev_LP2);
			cudaDeviceSynchronize();
		}
		// Main kernel now includes all enhanced features when USE_LDM_ENHANCED_TURBULENCE is defined
		b.x=(num_part_LDM-1)/t.x+1;
		if(CASE_ID<2000)KERNEL_time_update_LDM<<<b,t>>>(freq_LDM*dt,time,time_end,g_str,g_end,dev_LP1,dev_LP2,dev_P1);
		else if(CASE_ID>2000)KERNEL_time_update_LDM_direct<<<b,t>>>(freq_LDM*dt,time,g_str,g_end,dev_LP1,dev_LP2,dev_P1);
		cudaDeviceSynchronize();

		// Dry deposition and resuspension (functions_DD.cuh)
#if !DD_DISABLE
		KERNEL_dry_deposition<<<b,t>>>(freq_LDM*dt, time, dev_LP1, (Real)MOST_z0,
		                               g_str, g_end, dev_P1);
		cudaDeviceSynchronize();
#endif

		// LDM verification: U_mean = 0.3 m/s, D = 0.01 m^2/s
		// Case 1 (pure diffusion):     U_mean=0.0, D=0.01
		// Case 2 (pure advection):     U_mean=0.3, D=0.0
		// Case 3 (advection+diffusion): U_mean=0.3, D=0.01
		if(CASE_ID>2000){
			Real LDM_U_mean = 0.3;   // <-- 변경: 0.0 for Case 1
			Real LDM_D      = 0.01;  // <-- 변경: 0.0 for Case 2
			KERNEL_time_update_LDM_test<<<b,t>>>(freq_LDM*dt,time,dev_LP1,dev_LP2,LDM_U_mean,LDM_D);
		}
		KERNEL_clc_concentration<<<b,t>>>(g_str_LDM,g_end_LDM,dev_LP1);
		cudaDeviceSynchronize();
	}

	if (CASE_ID>2000){

        b.x=(num_part2-1)/t.x+1;
        Real src_x=5.0, src_y=0.0, src_z=3.0, src_r=0.183, src_val=1e5;  // src_r = σ, matched to LDM σ₀

        if(count==0)KERNEL_apply_source_concn_Gaussian<<<b,t>>>(dev_P1,dt,src_x,src_y,src_z,src_r,src_val);
        cudaDeviceSynchronize();

        Real diffusion_coeff=Scalar_Diffusion_coeff; // 확산 계수

        // Scalar transport sub-cycling: freq_scalar 스텝에 1번, dt를 그만큼 키워서 호출
        int freq_scalar = 10;
        if(count % freq_scalar == 0){
            Real sdt = freq_scalar*dt;
            b.x=(num_part2-1)/t.x+1;

            // BiMacCormack 2-pass: forward SL → correct with backward error
            KERNEL_BiMC_forward<<<b,t>>>(g_str, g_end, dev_P1, dev_P3, sdt);
            cudaDeviceSynchronize();
            KERNEL_BiMC_correct<<<b,t>>>(g_str, g_end, dev_P1, dev_P3, diffusion_coeff, sdt);
            cudaDeviceSynchronize();

            b.x=(num_part2-1)/t.x+1;
            KERNEL_copy_concn_adv<<<b,t>>>(dev_P1,dev_P3);
            cudaDeviceSynchronize();
        }

		b.x=(num_part2-1)/t.x+1;
        KERNEL_clc_concn_LDM<<<b,t>>>(g_str_LDM,g_end_LDM,dev_P1,dev_LP1);
        cudaDeviceSynchronize();
    }

//-------------------------------------------------------------------------------------------------
// 출력
//-------------------------------------------------------------------------------------------------

	// Output cadence: relative to the averaging-start latch if it has fired,
	// otherwise relative to the simulation start. Geometry (boundary,
	// buildings) is still one-shot at count==0.
	int_t count_eff = count;

#if LDM_MEAN_MODE
	// 평균장 모드: 유동장이 얼어 있으므로 t=0 에 한 번만 내리고, 이후에는 LDM 만 복사한다.
	// (종전에는 25 step 마다 P1/P2/P3 19 GB 를 D2H 로 옮기고 같은 유동장을 601 번 다시 썼다.)
	if(count == 0){
		cudaMemcpy(file_P1, dev_P1, num_part2*sizeof(part1), cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P2, dev_SP2, num_part2*sizeof(part2), cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P3, dev_P3, num_part2*sizeof(part3), cudaMemcpyDeviceToHost);
		save_vtk_bin_single_boundary(file_P1, file_P2, file_P3);
		save_vtk_bin_single_buildings(file_P1, file_P2, file_P3);
		save_vtk_bin_single_fluid_condition(file_P1, file_P2, file_P3);
		save_vtk_bin_single_fluid(file_P1, file_P2, file_P3);
	}
	if(count_eff % LDM_DUMP_EVERY == 0){
		cudaMemcpy(file_LP1, dev_LP1, num_part_LDM*sizeof(L_part1), cudaMemcpyDeviceToHost);
		save_vtk_bin_single_LDM(file_LP1);
	}
#else
	if(count_eff % 25 == 0){
    	cudaMemcpy(file_P1, dev_P1, num_part2*sizeof(part1), cudaMemcpyDeviceToHost);
    	cudaMemcpy(file_P2, dev_SP2, num_part2*sizeof(part2), cudaMemcpyDeviceToHost);
    	cudaMemcpy(file_P3, dev_P3, num_part2*sizeof(part3), cudaMemcpyDeviceToHost);

    	if(count == 0){
        	save_vtk_bin_single_boundary(file_P1, file_P2, file_P3);
        	save_vtk_bin_single_buildings(file_P1, file_P2, file_P3);
    	}
    	if(TRANSIENT_DUMP_EVERY>0 && count_eff%TRANSIENT_DUMP_EVERY==0 && time>=(Real)TRANSIENT_DUMP_START)
    		save_transient_velocity(file_P1, time);

    	if(count_avg < 0){
        	if(count_eff%SLICE_DUMP_EVERY==0) save_vtk_bin_single_fluid_condition(file_P1, file_P2, file_P3);

    	} 
		else {
        	if(count_eff%SLICE_DUMP_EVERY==0) save_vtk_bin_single_fluid_condition(file_P1, file_P2, file_P3);
        	if(count_eff%FULL_DUMP_EVERY==0) save_vtk_bin_single_fluid(file_P1, file_P2, file_P3);
        	if(count_eff%VORT_DUMP_EVERY==0) save_vtk_bin_single_vorticity(file_P1, file_P2, file_P3);
        	if(count_eff%LDM_DUMP_EVERY==0){
        		cudaMemcpy(file_LP1, dev_LP1, num_part_LDM*sizeof(L_part1), cudaMemcpyDeviceToHost);
        		save_vtk_bin_single_LDM(file_LP1);
        	}
    	}
	}
#endif  // LDM_MEAN_MODE
}

//-------------------------------------------------------------------------------------------------
// SOPHIA 메인 코드
//-------------------------------------------------------------------------------------------------
void*ISPH_Calc(void*arg){

	// 함수의 인자를 받아서 tid 에 저장 (tid = gpu 번호)
	int*idPtr,tid;
	idPtr=(int*)arg;
	tid=*idPtr;

	// 계산 시간 측정
	clock_t start, end;
  	double cpu_time_used;

	start=clock();

	// timestep control 을 위한 변수 설정
	// Real dt_CFL,V_MAX,K_stiff,eta;
	// dt_CFL=V_MAX=K_stiff=eta=0.0;

	num_cells=clc_num_cells();
	count=floor(time/dt+0.5);

	//-------------------------------------------------------------------------------------------------
	// Device 입자 생성
	//-------------------------------------------------------------------------------------------------

	// 계산할 GPU 정의: tid=GPU number
	cudaSetDevice(tid);

	// GPU 내에 분기 생성 : stream is for run the kernel and memcpy peer to peer at the same time.
	cudaStreamCreate(&str1[tid]);
	cudaStreamCreate(&str2[tid]);

	// 출력할 변수 선언 및 메모리 할당
	part1*file_P1;
	file_P1=(part1*)malloc(sizeof(part1)*num_part2);
	memset(file_P1,0,sizeof(part1)*num_part2);

	part2*file_P2;
	part3*file_P3;

	file_P2=(part2*)malloc(sizeof(part2)*num_part2);
	memset(file_P2,0,sizeof(part2)*num_part2);

	file_P3=(part3*)malloc(sizeof(part3)*num_part2);
	memset(file_P3,0,sizeof(part3)*num_part2);

	// LDM 입자 배열
	L_part1*file_LP1;
	file_LP1=(L_part1*)malloc(sizeof(L_part1)*num_part_LDM);
	memset(file_LP1,0,sizeof(L_part1)*num_part_LDM);

	L_part2*file_LP2;

	file_LP2=(L_part2*)malloc(sizeof(L_part2)*num_part_LDM);
	memset(file_LP2,0,sizeof(L_part2)*num_part_LDM);

	//-------------------------------------------------------------------------------------------------
	// Device/GPU 변수 선언 및 메모리 할당
	//-------------------------------------------------------------------------------------------------

	// NNPS 관련 변수
	int_t*g_idx,*p_idx,*g_idx_in,*p_idx_in,*g_str,*g_end;
	int_t*g_idx_LDM,*p_idx_LDM,*g_idx_in_LDM,*p_idx_in_LDM,*g_str_LDM,*g_end_LDM;
	int_t*b_idx,*b_idx_in,*b_str,*b_end;

	num_blocks=NBI*NBJ*NBK;

	// 주요 입자 변수
	part1*dev_P1,*dev_SP1;
	part2*dev_P2,*dev_SP2;
	part3*dev_SP3;

	// 주요 LDM 변수
	L_part1*dev_LP1,*dev_SLP1;
	L_part2*dev_LP2,*dev_SLP2;

	// P2P 데이터 변수 선언
	int*p2p_af_in,*p2p_idx_in,*p2p_af,*p2p_idx;

	// for Open boundary & Adaptive Particle Definement
	int*apr_num_part;
	cudaMalloc((void**)&apr_num_part,sizeof(int));
	cudaMemset(apr_num_part,k_num_part,sizeof(int));

	Real*computing_time;
	int_t*count_buffer,*dev_count_buffer,*num_buffer_temp;
	int_t*APR_cell,*dev_APR_cell;

	count_buffer=(int_t*)malloc(4*sizeof(int_t));
	APR_cell=(int_t*)malloc(num_cells*sizeof(int_t));
	computing_time=(Real*)malloc(20*sizeof(Real));

	cudaMalloc((void**)&dev_count_buffer,4*sizeof(int_t));
	cudaMalloc((void**)&num_buffer_temp,sizeof(int_t));
	cudaMalloc((void**)&dev_APR_cell,num_cells*sizeof(int_t));

	memset(count_buffer,0,4*sizeof(int_t));
	memset(APR_cell,0,num_cells*sizeof(int_t));
	memset(computing_time,0,20*sizeof(Real));

	cudaMemset(dev_count_buffer,0,4*sizeof(int_t));
	cudaMemset(num_buffer_temp,0,sizeof(int_t));
	cudaMemset(dev_APR_cell,0,num_cells*sizeof(int_t));

	// NNPS 입자 메모리 할당
	cudaMalloc((void**)&g_idx,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&p_idx,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&g_idx_in,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&p_idx_in,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&g_str,sizeof(int_t)*num_cells);
	cudaMalloc((void**)&g_end,sizeof(int_t)*num_cells);

	cudaMalloc((void**)&g_idx_LDM,sizeof(int_t)*num_part_LDM);
	cudaMalloc((void**)&p_idx_LDM,sizeof(int_t)*num_part_LDM);
	cudaMalloc((void**)&g_idx_in_LDM,sizeof(int_t)*num_part_LDM);
	cudaMalloc((void**)&p_idx_in_LDM,sizeof(int_t)*num_part_LDM);
	cudaMalloc((void**)&g_str_LDM,sizeof(int_t)*num_cells);
	cudaMalloc((void**)&g_end_LDM,sizeof(int_t)*num_cells);

	cudaMalloc((void**)&b_idx,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&b_idx_in,sizeof(int_t)*num_part2);
	cudaMalloc((void**)&b_str,sizeof(int_t)*num_blocks);
	cudaMalloc((void**)&b_end,sizeof(int_t)*num_blocks);

	// Device 입자 데이터 메모리 할당
	cudaMalloc((void**)&dev_P1,sizeof(part1)*num_part2);
	cudaMalloc((void**)&dev_SP1,sizeof(part1)*num_part2);
	cudaMalloc((void**)&dev_P2,sizeof(part2)*num_part2);
	cudaMalloc((void**)&dev_SP2,sizeof(part2)*num_part2);
	cudaMalloc((void**)&dev_SP3,sizeof(part3)*num_part2);

	// Device LDM 데이터 메모리 할당
	cudaMalloc((void**)&dev_LP1,sizeof(L_part1)*num_part_LDM);
	cudaMalloc((void**)&dev_SLP1,sizeof(L_part1)*num_part_LDM);
	cudaMalloc((void**)&dev_LP2,sizeof(L_part2)*num_part_LDM);
	cudaMalloc((void**)&dev_SLP2,sizeof(L_part2)*num_part_LDM);

	// NNPS 메모리 초기화
	cudaMemset(g_idx_in,0,sizeof(int_t)*num_part2);
	cudaMemset(p_idx_in,0,sizeof(int_t)*num_part2);
	cudaMemset(g_idx,0,sizeof(int_t)*num_part2);
	cudaMemset(p_idx,0,sizeof(int_t)*num_part2);
	cudaMemset(g_str,cu_memset,sizeof(int_t)*num_cells);
	cudaMemset(g_end,0,sizeof(int_t)*num_cells);

	cudaMemset(g_idx_in_LDM,0,sizeof(int_t)*num_part_LDM);
	cudaMemset(p_idx_in_LDM,0,sizeof(int_t)*num_part_LDM);
	cudaMemset(g_idx_LDM,0,sizeof(int_t)*num_part_LDM);
	cudaMemset(p_idx_LDM,0,sizeof(int_t)*num_part_LDM);
	cudaMemset(g_str_LDM,cu_memset,sizeof(int_t)*num_cells);
	cudaMemset(g_end_LDM,0,sizeof(int_t)*num_cells);

	// Device 입자 메모리 초기화
	cudaMemset(dev_P1,0,sizeof(part1)*num_part2);
	cudaMemset(dev_SP1,0,sizeof(part1)*num_part2);
	cudaMemset(dev_P2,0,sizeof(part2)*num_part2);
	cudaMemset(dev_SP2,0,sizeof(part2)*num_part2);
	cudaMemset(dev_SP3,0,sizeof(part3)*num_part2);

	// Device LDM 입자 메모리 초기화
	cudaMemset(dev_LP1,0,sizeof(L_part1)*num_part_LDM);
	cudaMemset(dev_SLP1,0,sizeof(L_part1)*num_part_LDM);
	cudaMemset(dev_LP2,0,sizeof(L_part2)*num_part_LDM);
	cudaMemset(dev_SLP2,0,sizeof(L_part2)*num_part_LDM);

	//-------------------------------------------------------------------------------------------------
	// Device/GPU로 데이터 복사
	//-------------------------------------------------------------------------------------------------

	// Sovler 전역변수 Device로 복사
	cudaMemcpyToSymbol(k_vii,vii,sizeof(int_t)*vii_size);
	cudaMemcpyToSymbol(k_vif,vif,sizeof(Real)*vif_size);

	// 물성 Table 데이터 Device로 복사
	cudaMemcpyToSymbol(k_Tab_T,host_Tab_T,sizeof(Real)*table_size);
	cudaMemcpyToSymbol(k_Tab_h,host_Tab_h,sizeof(Real)*table_size);
	cudaMemcpyToSymbol(k_Tab_k,host_Tab_k,sizeof(Real)*table_size);
	cudaMemcpyToSymbol(k_Tab_cp,host_Tab_cp,sizeof(Real)*table_size);
	cudaMemcpyToSymbol(k_Tab_vis,host_Tab_vis,sizeof(Real)*table_size);

	cudaMemcpyToSymbol(k_table_index,host_table_index,sizeof(int)*10);
	cudaMemcpyToSymbol(k_table_size,host_table_size,sizeof(int)*10);
#if INFLOW_TURB
	init_inflow_turb();     // 유입 합성 난류 모드 (고정 시드, 체크포인트 불필요)
#endif

	// Host 입자정보(HP1)를 분할하여(DHP1) Device로 복사(dev_P1)
	DHP1[tid]=(part1*)malloc(num_part2*sizeof(part1));
	memset(DHP1[tid],0,sizeof(part1)*num_part2);

	// 모든 입자를 더미로 초기화
	for(int i=0;i<num_part2;i++) DHP1[tid][i].i_type=3;

	c_initial_inner_outer_particle_single(HP1,DHP1[tid],tid);											// (CAUTION)
	cudaMemcpy(dev_P1,DHP1[tid],num_part2*sizeof(part1),cudaMemcpyHostToDevice);	// single gpu 이면 그냥 HP1을 device에 복사
	pthread_barrier_wait(&barrier);

#if DEACTIVATE_IBM
	// Building-free run: deactivate IBM (p_type==1000) particles (i_type=3).
	{
		dim3 t_ibm(128), b_ibm((num_part2-1)/128+1);
		KERNEL_deactivate_IBM<<<b_ibm,t_ibm>>>(dev_P1);
		cudaDeviceSynchronize();
		if(tid==0) printf("[DEACTIVATE_IBM] IBM (p_type=1000) particles deactivated -> i_type=3\n");
	}
#endif

	// Host 입자정보(HP1)를 분할하여(DHP1) Device로 복사(dev_P1) (LDM)
	DHLP1[tid]=(L_part1*)malloc(num_part_LDM*sizeof(L_part1));
	memset(DHLP1[tid],0,sizeof(L_part1)*num_part_LDM);

	// 모든 입자를 더미로 초기화
	for(int i=0;i<num_part_LDM;i++) DHLP1[tid][i].i_type=3;

	c_initial_inner_outer_particle_single_LDM(HLP1,DHLP1[tid],tid);											// (CAUTION)
	cudaMemcpy(dev_LP1,DHLP1[tid],num_part_LDM*sizeof(L_part1),cudaMemcpyHostToDevice);	// single gpu 이면 그냥 HP1을 device에 복사
	pthread_barrier_wait(&barrier);

	if(tid==0){
		printf("\n-----------------------------------------------------------\n");
		printf("GPU Domain Division Success\n");
		printf("-----------------------------------------------------------\n\n");
	}

	//-------------------------------------------------------------------------------------------------
	// 화면 출력용 기타 변수들 정의 및 메모리 할당 (최대속도, 최대힘 등)
	//-------------------------------------------------------------------------------------------------

	// host
	Real *max_umag0,*max_rho0,*max_ftotal0,*max_temp0,*min_temp0,*max_fbz0,*max_uz0;
	max_umag0=(Real*)malloc(sizeof(Real));
	max_rho0=(Real*)malloc(sizeof(Real));
	max_ftotal0=(Real*)malloc(sizeof(Real));
	max_temp0=(Real*)malloc(sizeof(Real));
	min_temp0=(Real*)malloc(sizeof(Real));
	max_fbz0=(Real*)malloc(sizeof(Real));
	max_uz0=(Real*)malloc(sizeof(Real));
	max_umag0[0]=max_ftotal0[0]=max_rho0[0]=max_temp0[0]=0.0;
	min_temp0[0]=0.0; max_fbz0[0]=0.0; max_uz0[0]=0.0;

	// Cold-spot diagnostic: host/device buffers for ArgMin (storage allocated below near other CUB temp storage)
	cub::KeyValuePair<int,Real> *h_min_kv=(cub::KeyValuePair<int,Real>*)malloc(sizeof(cub::KeyValuePair<int,Real>));
	cub::KeyValuePair<int,Real> *d_min_kv;
	cudaMalloc((void**)&d_min_kv,sizeof(cub::KeyValuePair<int,Real>));
	void* dev_argmin_storage=NULL;
	size_t argmin_storage_bytes=0;

	// Hot-u diagnostic: find particle with max |u| magnitude (= velocity-limit clipping candidate)
	cub::KeyValuePair<int,Real> *h_max_u_kv=(cub::KeyValuePair<int,Real>*)malloc(sizeof(cub::KeyValuePair<int,Real>));
	cub::KeyValuePair<int,Real> *d_max_u_kv;
	cudaMalloc((void**)&d_max_u_kv,sizeof(cub::KeyValuePair<int,Real>));
	void* dev_argmax_storage=NULL;
	size_t argmax_storage_bytes=0;
	part1 cold_p1;

	// device
	Real*max_rho,*max_umag,*max_ft,*max_temp,*min_temp,*max_fbz,*max_uz;
	Real*d_max_umag0,*d_max_rho0,*d_max_ftotal0,*d_max_temp0,*d_min_temp0,*d_max_fbz0,*d_max_uz0;
	cudaMalloc((void**)&max_rho,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_umag,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_ft,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_temp,sizeof(Real)*num_part2);
	cudaMalloc((void**)&min_temp,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_fbz,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_uz,sizeof(Real)*num_part2);
	cudaMalloc((void**)&d_max_umag0,sizeof(Real));
	cudaMalloc((void**)&d_max_rho0,sizeof(Real));
	cudaMalloc((void**)&d_max_ftotal0,sizeof(Real));
	cudaMalloc((void**)&d_max_temp0,sizeof(Real));
	cudaMalloc((void**)&d_min_temp0,sizeof(Real));
	cudaMalloc((void**)&d_max_fbz0,sizeof(Real));
	cudaMalloc((void**)&d_max_uz0,sizeof(Real));
	cudaMemset(max_umag,0,sizeof(Real)*num_part2);
	cudaMemset(max_rho,0,sizeof(Real)*num_part2);
	cudaMemset(max_ft,0,sizeof(Real)*num_part2);
	cudaMemset(max_temp,0,sizeof(Real)*num_part2);
	cudaMemset(min_temp,0,sizeof(Real)*num_part2);
	cudaMemset(max_fbz,0,sizeof(Real)*num_part2);
	cudaMemset(max_uz,0,sizeof(Real)*num_part2);
	cudaMemset(d_max_umag0,0,sizeof(Real));
	cudaMemset(d_max_rho0,0,sizeof(Real));
	cudaMemset(d_max_ftotal0,0,sizeof(Real));
	cudaMemset(d_max_temp0,0,sizeof(Real));
	cudaMemset(d_min_temp0,0,sizeof(Real));
	cudaMemset(d_max_fbz0,0,sizeof(Real));
	cudaMemset(d_max_uz0,0,sizeof(Real));

	//-------------------------------------------------------------------------------------------------
	// 정렬(Sorting)을 위한 CUB 라이브러리 변수 준비
	//-------------------------------------------------------------------------------------------------

	// Sorting & Max variable to use CUB Library
	void*dev_sort_storage=NULL;
	void*dev_sort_storage_LDM=NULL;
	void*dev_max_storage=NULL;
	size_t sort_storage_bytes=0;
	size_t sort_storage_bytes_LDM=0;
	size_t max_storage_bytes=0;

	// Determine Sorting & Maximum Value Setting for Total Particle Data
	cub::DeviceRadixSort::SortPairs(dev_sort_storage,sort_storage_bytes,g_idx_in,g_idx,p_idx_in,p_idx,num_part2);
	cub::DeviceRadixSort::SortPairs(dev_sort_storage_LDM,sort_storage_bytes_LDM,g_idx_in_LDM,g_idx_LDM,p_idx_in_LDM,p_idx_LDM,num_part_LDM);
	cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_umag,d_max_umag0,num_part2);
	cub::DeviceReduce::ArgMin(dev_argmin_storage,argmin_storage_bytes,min_temp,d_min_kv,num_part2);
	cub::DeviceReduce::ArgMax(dev_argmax_storage,argmax_storage_bytes,max_umag,d_max_u_kv,num_part2);
	cudaDeviceSynchronize();
	cudaMalloc((void**)&dev_sort_storage,sort_storage_bytes);
	cudaMalloc((void**)&dev_sort_storage_LDM,sort_storage_bytes_LDM);
	cudaMalloc((void**)&dev_max_storage,max_storage_bytes);
	cudaMalloc((void**)&dev_argmin_storage,argmin_storage_bytes);
	cudaMalloc((void**)&dev_argmax_storage,argmax_storage_bytes);

	// Domain-integrated diagnostics buffers (Σk_turb, Σ½|u|² over fluid)
	Real *diag_kt, *diag_uu, *d_sum_kt, *d_sum_uu;
	Real sum_kt0=0.0, sum_uu0=0.0;
	cudaMalloc((void**)&diag_kt, sizeof(Real)*num_part2);
	cudaMalloc((void**)&diag_uu, sizeof(Real)*num_part2);
	cudaMalloc((void**)&d_sum_kt, sizeof(Real));
	cudaMalloc((void**)&d_sum_uu, sizeof(Real));
	cudaMemset(diag_kt,0,sizeof(Real)*num_part2);
	cudaMemset(diag_uu,0,sizeof(Real)*num_part2);
	void *dev_sum_storage=NULL;
	size_t sum_storage_bytes=0;
	cub::DeviceReduce::Sum(dev_sum_storage, sum_storage_bytes, diag_kt, d_sum_kt, num_part2);
	cudaMalloc((void**)&dev_sum_storage, sum_storage_bytes);

	pthread_barrier_wait(&barrier);

	//-------------------------------------------------------------------------------------------------
	// 코드 메인
	//-------------------------------------------------------------------------------------------------

	// 초기상태 및 설정 출력
	if(tid==0){
		printf("-----------------------------------------------------------\n");
		printf("Input Summary: \n");
		printf("-----------------------------------------------------------\n");
		printf("	Total number of particles=%d\n",num_part);
		printf("	Device number of particles=%d\n",num_part2);
		printf("	P2P number of particles=%d\n",num_p2p);
		printf("	NI=%d,	NJ=%d,	NK=%d\n",NI,NJ,NK);
		printf("-----------------------------------------------------------\n\n");
		// Input Check
		printf("-----------------------------------------------------------\n");
		printf("Input Check: \n");
		printf("-----------------------------------------------------------\n");
		// check Domain Status
		printf("x min, max : %f %f\n",x_min,x_max);
		printf("y min, max : %f %f\n",y_min,y_max);
		printf("z min, max : %f %f\n",z_min,z_max);
		printf("Cell Size(dcell) %f\n",dcell);
		printf("Number of Cells Per a GPU in x-direction(calc_area) %d\n",calc_area);
		printf("-----------------------------------------------------------\n\n");
		// print out status
		printf("\n");
		printf("-----------------------------\n");
		printf("Start Simultion!!\n");
		printf("-----------------------------\n");
		printf("\n");
	}
	pthread_barrier_wait(&barrier);

	//-------------------------------------------------------------------------------------------------
	// 코드 메인
	//-------------------------------------------------------------------------------------------------

	// t.x=128;
	// b.x=(num_part2-1)/t.x+1;
	// KERNEL_clc_TemptoEnthalpy<<<b,t>>>(dev_P1,dev_P2);
	// cudaDeviceSynchronize();

	// Initialise averaging-window bounds. AVG_START_TIME == 0 → auto-detect
	// via QS criterion in the 50-step diagnostic block below. >0 → fixed.
	// Time-averaging DISABLED (plain transient run): open the LDM gate from
	// t=0 so each LDM particle releases at its own absolute t_release. No QS
	// auto-detect, no accumulation, no mean-field freeze, no early stop.
	if(tid==0){
		if((Real)AVG_START_TIME > (Real)0.0){
			// 고정 평균 구간. AVG_END_TIME>0 이면 그 값, 아니면 start+AVG_WINDOW.
			h_avg_start = (Real)AVG_START_TIME;
			if((Real)AVG_END_TIME > (Real)0.0)     h_avg_end = (Real)AVG_END_TIME;
			else if((Real)AVG_WINDOW > (Real)0.0)  h_avg_end = h_avg_start + (Real)AVG_WINDOW;
			else                                   h_avg_end = (Real)1e30;
			printf("[avg] ENABLED (fixed): window [%.3f, %.3f), LDM gate opens at %.3f\n",
			       (double)h_avg_start, (double)h_avg_end, (double)h_avg_start);
		} else {
			h_avg_start = (Real)0.0;
			h_avg_end   = (Real)1e30;
			printf("[avg] DISABLED: LDM gate open from t=0; releases use absolute t_release.\n");
		}
		Real hs = h_avg_start, he = h_avg_end;
		cudaMemcpyToSymbol(avg_start_dyn, &hs, sizeof(Real));
		cudaMemcpyToSymbol(avg_end_dyn,   &he, sizeof(Real));
	}
	pthread_barrier_wait(&barrier);

#if LDM_MEAN_MODE
	// Mean-flow LDM mode: load_mean_field is now called inside SOPHIA_single_ISPH
	// (count==0 path, AFTER the first NNPS sort) so the mean_field.bin record
	// indices (saved in Run A's sorted order) align with the current sorting.
	// Only the LDM gate is forced open from t=0 here.
	if(tid==0){
		h_avg_start = (Real)0.0;
		h_avg_end   = (Real)1e30;
		Real hs = h_avg_start, he = h_avg_end;
		cudaMemcpyToSymbol(avg_start_dyn, &hs, sizeof(Real));
		cudaMemcpyToSymbol(avg_end_dyn,   &he, sizeof(Real));
		count_avg = 0;
		printf("[mean] LDM gate forced open: h_avg_start=0, h_avg_end=inf, count_avg=0\n");
	}
	pthread_barrier_wait(&barrier);
#endif

	// Resume from checkpoint (env SOPHIA_RESTART=<path>). Any failure aborts instead of
	// silently starting from t=0 and overwriting the run's outputs.
	{
		const char*restart_path = getenv("SOPHIA_RESTART");
		if(restart_path && restart_path[0]){
			bool ok = (ngpu==1);
			if(!ok) printf("[ckpt] FATAL: restart supports ngpu==1 only (ngpu=%d)\n", (int)ngpu);
#if LDM_MEAN_MODE
			ok = false;
			printf("[ckpt] FATAL: restart is not supported with LDM_MEAN_MODE (mean field loads at count==0)\n");
#endif
			if(ok) ok = load_checkpoint(restart_path, dev_P1, dev_P2, g_str, g_end, dev_LP1, dev_LP2,
			                            g_str_LDM, g_end_LDM, file_P1, file_P2, file_LP1, file_LP2);
			if(!ok || count<=0){
				printf("[ckpt] aborting (restart requested but not possible)\n");
				fflush(stdout);
				exit(1);
			}
		}
	}

	while(1){

		if(ngpu==1){
		SOPHIA_single_ISPH(g_idx,p_idx,g_idx_in,p_idx_in,g_str,g_end,g_idx_LDM,p_idx_LDM,g_idx_in_LDM,p_idx_in_LDM,g_str_LDM,g_end_LDM,b_idx,b_idx_in,b_str,b_end,dev_P1,dev_SP1,dev_P2,dev_SP2,dev_SP3,
					dev_LP1,dev_SLP1,dev_LP2,dev_SLP2,p2p_af_in,p2p_idx_in,p2p_af,p2p_idx,dev_sort_storage,&sort_storage_bytes,dev_sort_storage_LDM,&sort_storage_bytes_LDM,file_P1,file_P2,file_P3,file_LP1,tid);
		//SOPHIA_single_ISPH(g_idx,p_idx,g_idx_in,p_idx_in,g_str,g_end,b_idx,b_idx_in,b_str,b_end,dev_P1,dev_SP1,dev_P2,dev_SP2,dev_SP3,
		//			p2p_af_in,p2p_idx_in,p2p_af,p2p_idx,dev_sort_storage,&sort_storage_bytes,file_P1,file_P2,file_P3,tid);
		}

		//-------------------------------------------------------------------------------------------------
		// Time-step Control
		//-------------------------------------------------------------------------------------------------
		if(tid==0){
			time+=dt;
			count++;

			// Unified count_avg latch: covers the manual AVG_START_TIME case
			// where h_avg_start was set at startup (no QS event fires).
			// QS/force-latch paths set count_avg directly inside the 50-step
			// block; the `count_avg<0` guard prevents double-latching.
			if(count_avg < 0 && h_avg_start < (Real)1e29 && time >= h_avg_start){
				count_avg = count;
			}

			//timestep is updated every 10 steps ------------ estimate new timestep (Goswami & Pajarola(2011))
			if((count%50)==0){
				dim3 t,b;
				t.x=128;
				b.x=(num_part2-1)/t.x+1;
				kernel_copy_max<<<b,t>>>(dev_P1,dev_SP2,dev_SP3,max_rho,max_ft,max_umag,max_temp,min_temp,max_fbz,max_uz);
				cudaDeviceSynchronize();

				// Find Max Velocity & Force using CUB - TID=0
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_umag,d_max_umag0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_rho,d_max_rho0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_ft,d_max_ftotal0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_temp,d_max_temp0,num_part2);
				cub::DeviceReduce::Min(dev_max_storage,max_storage_bytes,min_temp,d_min_temp0,num_part2);
				cub::DeviceReduce::ArgMin(dev_argmin_storage,argmin_storage_bytes,min_temp,d_min_kv,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_fbz,d_max_fbz0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_uz,d_max_uz0,num_part2);
				cudaDeviceSynchronize();
				cudaMemcpy(max_umag0,d_max_umag0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_rho0,d_max_rho0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_ftotal0,d_max_ftotal0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_temp0,d_max_temp0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(min_temp0,d_min_temp0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_fbz0,d_max_fbz0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_uz0,d_max_uz0,sizeof(Real),cudaMemcpyDeviceToHost);

				end=clock();
				cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
				printf("Simulation time: %5.2f sec\n", time);
				printf("Computing time: %5.2f sec\n\n", cpu_time_used);
				printf("%d\t tu_max=%5.2f\ttt_max=%5.2f\ttt_min=%5.2f\tfbz_max=%.3e\tuz_max=%5.2f\tftotal_max=%.3e\n\n",
					count,max_umag0[0],max_temp0[0],min_temp0[0],max_fbz0[0],max_uz0[0],max_ftotal0[0]);

				// Cold-spot diagnostic
				cudaMemcpy(h_min_kv,d_min_kv,sizeof(cub::KeyValuePair<int,Real>),cudaMemcpyDeviceToHost);
				int cold_idx=h_min_kv->key;
				if(cold_idx>=0 && cold_idx<num_part2){
					cudaMemcpy(&cold_p1,dev_P1+cold_idx,sizeof(part1),cudaMemcpyDeviceToHost);
					printf("    [cold-spot] idx=%d  T=%6.3f  pos=(%7.2f,%7.2f,%7.2f)  p_type=%d  buffer_type=%d\n\n",
						cold_idx,cold_p1.temp,cold_p1.x,cold_p1.y,cold_p1.z,(int)cold_p1.p_type,(int)cold_p1.buffer_type);
				}

				// Hot-u diagnostic: locate particle with max |u| magnitude (velocity-limit clipping)
				cub::DeviceReduce::ArgMax(dev_argmax_storage,argmax_storage_bytes,max_umag,d_max_u_kv,num_part2);
				cudaDeviceSynchronize();
				cudaMemcpy(h_max_u_kv,d_max_u_kv,sizeof(cub::KeyValuePair<int,Real>),cudaMemcpyDeviceToHost);
				int hot_idx=h_max_u_kv->key;
				if(hot_idx>=0 && hot_idx<num_part2){
					part1 hot_p1;
					cudaMemcpy(&hot_p1,dev_P1+hot_idx,sizeof(part1),cudaMemcpyDeviceToHost);
					Real umag_hot=sqrt(hot_p1.ux*hot_p1.ux+hot_p1.uy*hot_p1.uy+hot_p1.uz*hot_p1.uz);
					printf("    [hot-u]     idx=%d  |u|=%6.3f  u=(%7.3f,%7.3f,%7.3f)  pos=(%7.2f,%7.2f,%7.2f)  p_type=%d  i_type=%d  buffer_type=%d\n",
						hot_idx,umag_hot,hot_p1.ux,hot_p1.uy,hot_p1.uz,
						hot_p1.x,hot_p1.y,hot_p1.z,(int)hot_p1.p_type,(int)hot_p1.i_type,(int)hot_p1.buffer_type);
				}

				// Domain-integrated diagnostics (fluid only): plateau → quasi-steady state
				{
					dim3 td,bd; td.x=128; bd.x=(num_part2-1)/td.x+1;
					kernel_copy_diag_sum<<<bd,td>>>(dev_P1, diag_kt, diag_uu);
					cudaDeviceSynchronize();
					cub::DeviceReduce::Sum(dev_sum_storage, sum_storage_bytes, diag_kt, d_sum_kt, num_part2);
					cub::DeviceReduce::Sum(dev_sum_storage, sum_storage_bytes, diag_uu, d_sum_uu, num_part2);
					cudaMemcpy(&sum_kt0, d_sum_kt, sizeof(Real), cudaMemcpyDeviceToHost);
					cudaMemcpy(&sum_uu0, d_sum_uu, sizeof(Real), cudaMemcpyDeviceToHost);
					unsigned int h_ac=0;
					cudaMemcpyFromSymbol(&h_ac, avg_count, sizeof(unsigned int));
					printf("    [domain] Σk_turb=%.4e   Σ½|u|²=%.4e   avg_count=%u (h_avg_start=%.3f)\n\n",
						sum_kt0, sum_uu0, h_ac, (double)h_avg_start);

					// Quasi-steady auto-detect: relative std of Σ½|u|² over the
					// last QS_DETECT_N samples. Latch h_avg_start once.
					// Sample only every QS_SAMPLE_STRIDE'th 50-step block so the
					// effective window = N × STRIDE × 50 simulation steps.
					if(h_avg_start > (Real)1e29 &&
					   ((count / 50) % QS_SAMPLE_STRIDE == 0)){
						g_qs_uu_buf[g_qs_uu_head] = sum_uu0;
						g_qs_uu_head = (g_qs_uu_head + 1) % QS_DETECT_N;
						if(g_qs_uu_fill < QS_DETECT_N) g_qs_uu_fill++;

						if(g_qs_uu_fill == QS_DETECT_N){
							Real mean = 0.0, var = 0.0;
							for(int q=0; q<QS_DETECT_N; q++) mean += g_qs_uu_buf[q];
							mean /= (Real)QS_DETECT_N;
							for(int q=0; q<QS_DETECT_N; q++){
								Real d = g_qs_uu_buf[q] - mean;
								var += d*d;
							}
							var /= (Real)QS_DETECT_N;
							Real rel_std = (mean > (Real)1e-12) ? sqrt(var) / mean : (Real)1e30;
							printf("    [QS] uu sliding rel_std=%.4e (eps=%.4e, N=%d, mean=%.4e)\n",
								(double)rel_std, (double)QS_DETECT_EPS, (int)QS_DETECT_N, (double)mean);

							if(rel_std < (Real)QS_DETECT_EPS && time >= (Real)QS_MIN_TIME){
								h_avg_start = time;
								if((Real)AVG_END_TIME > (Real)0.0)        h_avg_end = (Real)AVG_END_TIME;
								else if((Real)AVG_WINDOW > (Real)0.0)     h_avg_end = h_avg_start + (Real)AVG_WINDOW;
								else                                      h_avg_end = (Real)1e30;
								Real hs = h_avg_start, he = h_avg_end;
								cudaMemcpyToSymbol(avg_start_dyn, &hs, sizeof(Real));
								cudaMemcpyToSymbol(avg_end_dyn,   &he, sizeof(Real));
								count_avg = count;
								printf("    [QS] >>> quasi-steady DETECTED at t=%.3f (count=%d). avg window=[%.3f, %.3f]\n",
									(double)time, (int)count, (double)h_avg_start, (double)h_avg_end);
							}
						}

						// Force-latch fallback: if QS still not detected by
						// AVG_FORCE_START_TIME, start averaging anyway.
						if(h_avg_start > (Real)1e29 &&
						   (Real)AVG_FORCE_START_TIME > (Real)0.0 &&
						   time >= (Real)AVG_FORCE_START_TIME){
							h_avg_start = time;
							if((Real)AVG_END_TIME > (Real)0.0)        h_avg_end = (Real)AVG_END_TIME;
							else if((Real)AVG_WINDOW > (Real)0.0)     h_avg_end = h_avg_start + (Real)AVG_WINDOW;
							else                                      h_avg_end = (Real)1e30;
							Real hs = h_avg_start, he = h_avg_end;
							cudaMemcpyToSymbol(avg_start_dyn, &hs, sizeof(Real));
							cudaMemcpyToSymbol(avg_end_dyn,   &he, sizeof(Real));
							count_avg = count;
							printf("    [QS] >>> force-latched at t=%.3f (count=%d, no QS by AVG_FORCE_START_TIME=%.3f). avg window=[%.3f, %.3f]\n",
								(double)time, (int)count, (double)AVG_FORCE_START_TIME,
								(double)h_avg_start, (double)h_avg_end);
						}
					}
				}
			}
		}

		// Periodic checkpoint. count has already been advanced, so the file resumes at the next step.
		if(tid==0 && ngpu==1 && (int_t)CKPT_EVERY > 0 && (count % (int_t)CKPT_EVERY)==0 && time < time_end){
			save_checkpoint(CKPT_PATH, dev_P1, dev_P2, g_str, g_end, dev_LP1, dev_LP2,
			                g_str_LDM, g_end_LDM, file_P1, file_P2, file_LP1, file_LP2);
		}

		pthread_barrier_wait(&barrier);
		if(time>=time_end) break;
		// Early stop: LDM_DURATION_AFTER_AVG seconds past h_avg_start latch.
		if((Real)LDM_DURATION_AFTER_AVG > (Real)0.0 &&
		   h_avg_start < (Real)1e29 &&
		   time >= h_avg_start + (Real)LDM_DURATION_AFTER_AVG){
			if(tid==0) printf("[avg] LDM_DURATION_AFTER_AVG reached (h_avg_start=%.3f + %.2f); stopping.\n",
			                 (double)h_avg_start, (double)LDM_DURATION_AFTER_AVG);
			break;
		}

	}


	//-------------------------------------------------------------------------------------------------
	// ##. Save Restart File — periodic checkpoints are written inside the loop (CKPT_EVERY, save_checkpoint)
	//-------------------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------------------
	// ##. Memory Free
	//-------------------------------------------------------------------------------------------------
	free(file_P1);
	free(max_umag0);
	free(max_rho0);
	free(max_ftotal0);
	free(max_temp0);
	free(min_temp0);
	free(max_fbz0);
	free(max_uz0);
	cudaFree(g_idx);
	cudaFree(p_idx);
	cudaFree(g_idx_in);
	cudaFree(p_idx_in);
	cudaFree(g_str);
	cudaFree(g_end);
	cudaFree(g_idx_LDM);
	cudaFree(p_idx_LDM);
	cudaFree(g_idx_in_LDM);
	cudaFree(p_idx_in_LDM);
	cudaFree(g_str_LDM);
	cudaFree(g_end_LDM);
	cudaFree(dev_P1);
	cudaFree(dev_SP1);
	cudaFree(dev_P2);
	cudaFree(dev_SP2);
	cudaFree(dev_SP3);
	cudaFree(dev_LP1);
	cudaFree(dev_SLP1);
	cudaFree(dev_LP2);
	cudaFree(dev_SLP2);
	cudaFree(max_umag);
	cudaFree(max_rho);
	cudaFree(max_ft);
	cudaFree(max_temp);
	cudaFree(min_temp);
	cudaFree(max_fbz);
	cudaFree(max_uz);
	cudaFree(d_max_umag0);
	cudaFree(d_max_rho0);
	cudaFree(d_max_ftotal0);
	cudaFree(d_max_temp0);
	cudaFree(d_min_temp0);
	cudaFree(d_max_fbz0);
	cudaFree(d_max_uz0);
	cudaFree(dev_sort_storage);
	cudaFree(dev_max_storage);
	cudaFree(d_min_kv);
	cudaFree(dev_argmin_storage);
	free(h_min_kv);
	cudaFree(d_max_u_kv);
	cudaFree(dev_argmax_storage);
	free(h_max_u_kv);
	cudaStreamDestroy(str1[tid]);
	cudaStreamDestroy(str2[tid]);
	pthread_barrier_wait(&barrier);

	return 0;
}
