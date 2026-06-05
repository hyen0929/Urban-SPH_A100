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
	
	}
	else{
		// 0스텝 이후 단순 복사 >> 입자정보를 SP1, SP2 에 복사
		cudaMemcpy(dev_SP1,dev_P1,sizeof(part1)*num_part2,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
		cudaMemcpy(dev_SP2,dev_P2,sizeof(part2)*num_part2,cudaMemcpyDeviceToDevice);
		pthread_barrier_wait(&barrier);
		
		// 일부 입자정보 리셋
		cudaMemset(dev_P3,0.0,sizeof(part3)*num_part2);
		pthread_barrier_wait(&barrier);
	}

	// LDM 주변입자 검색 (농도 계산)
	// Gate: LDM dormant until time >= AVG_START_TIME (release / advection / NNPS all skipped)
	if(time >= (Real)AVG_START_TIME){
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
	}  // end of AVG_START_TIME gate (LDM NNPS block)

//-------------------------------------------------------------------------------------------------
// 계산 준비: gradient correction, filter, reference density, p_type switch, penetration, density gradient
//-------------------------------------------------------------------------------------------------
	
	b.x=(num_part2-1)/t.x+1;
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

// #if IMPLICIT_PPE
// 	//--- Implicit Jacobi PPE: iterate up to maxIteration, early-exit on max|delta p| < tol ---
// 	// Tracks 3 metrics per step:
// 	//   - Jacobi delta:    max|p^{k+1} - p^k|  (iteration convergence indicator, early-exit gate)
// 	//   - max|bi|:         max|rho*div(u*)/dt| (input incompressibility error from predictor)
// 	//   - max|A*p - b|:    true PPE residual   (proportional to |div(u^{n+1})|, the real metric)
// 	b.x=(num_part2-1)/t.x+1;

// 	// Lazy init (once, on first SOPHIA_single_ISPH call)
// 	if(g_dev_pres_old==NULL){
// 		cudaMalloc((void**)&g_dev_pres_old,    sizeof(Real)*num_part2);
// 		cudaMalloc((void**)&g_dev_pres_delta,  sizeof(Real)*num_part2);
// 		cudaMalloc((void**)&g_dev_bi,          sizeof(Real)*num_part2);
// 		cudaMalloc((void**)&g_dev_ppe_resid,   sizeof(Real)*num_part2);
// 		cudaMalloc((void**)&g_d_pres_delta_max,sizeof(Real));
// 		cudaMalloc((void**)&g_d_bi_max,        sizeof(Real));
// 		cudaMalloc((void**)&g_d_resid_max,     sizeof(Real));
// 		cudaMemset(g_dev_pres_old,   0, sizeof(Real)*num_part2);
// 		cudaMemset(g_dev_pres_delta, 0, sizeof(Real)*num_part2);
// 		cudaMemset(g_dev_bi,         0, sizeof(Real)*num_part2);
// 		cudaMemset(g_dev_ppe_resid,  0, sizeof(Real)*num_part2);
// 		// CUB Max-reduce temp storage (sized via dummy call — same array length for all 3 reduces)
// 		cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
// 		                       g_dev_pres_delta, g_d_pres_delta_max, num_part2);
// 		cudaMalloc(&g_dev_resid_storage, g_resid_storage_bytes);
// 		if(tid==0){
// 			char fn[512];
// 			sprintf(fn, "./%s/PPE_iter_log.txt", OUTPUT_CASE_NAME);
// 			g_fp_ppe_log = fopen(fn, "w");
// 			if(g_fp_ppe_log){
// 				fprintf(g_fp_ppe_log, "# IMPLICIT_PPE Jacobi solver log\n");
// 				fprintf(g_fp_ppe_log, "# tolerance=%g  minIter=%d  maxIter=%d\n",
// 				        (double)PPE_TOLERANCE, (int)minIteration, (int)maxIteration);
// 				fprintf(g_fp_ppe_log, "# bi_max     = |rho*div(u*)/dt|_inf   (input divergence; same every iter)\n");
// 				fprintf(g_fp_ppe_log, "# resid_max  = |A*p_final - b|_inf    (true PPE residual; ~|div(u^{n+1})|*rho/dt)\n");
// 				fprintf(g_fp_ppe_log, "# delta_p    = max|p^{k+1} - p^k|     (Jacobi iteration delta; early-exit gate)\n");
// 				fprintf(g_fp_ppe_log, "# step\ttime\titer_used\tdelta_p\tbi_max\tresid_max\tconverged\n");
// 				fflush(g_fp_ppe_log);
// 			}
// 		}
// 	}

// 	// Snapshot pres -> pres_old so the first delta is well-defined
// 	KERNEL_snapshot_pres<<<b,t>>>(dev_SP1, g_dev_pres_old);
// 	cudaDeviceSynchronize();

// 	int_t ppe_iter_used = 0;
// 	Real  ppe_delta_p   = 0.0f;
// 	Real  ppe_bi_max    = 0.0f;
// 	Real  ppe_resid_max = 0.0f;
// 	int   ppe_converged = 0;
// 	int_t maxit_safe = (maxIteration>0) ? maxIteration : 1;
// 	int_t minit_safe = (minIteration>0) ? minIteration : 1;

// 	for(int_t it=0; it<maxit_safe; it++){
// 		// Zero bi array so non-fluid slots don't pollute Max reduce (PPE_prep skips them via early return)
// 		cudaMemsetAsync(g_dev_bi, 0, sizeof(Real)*num_part2);
// 		KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_dev_bi);
// 		cudaDeviceSynchronize();
// 		KERNEL_PPE_calc3D<<<b,t>>>(dev_SP1);
// 		cudaDeviceSynchronize();
// 		KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
// 		cudaDeviceSynchronize();

// 		// Compute per-particle |p_new - p_old|, snapshot p -> pres_old
// 		KERNEL_compute_pres_delta<<<b,t>>>(dev_SP1, g_dev_pres_old, g_dev_pres_delta);
// 		cudaDeviceSynchronize();

// 		// Max reduce -> device scalar, copy to host
// 		cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
// 		                       g_dev_pres_delta, g_d_pres_delta_max, num_part2);
// 		cudaDeviceSynchronize();
// 		cudaMemcpy(&ppe_delta_p, g_d_pres_delta_max, sizeof(Real), cudaMemcpyDeviceToHost);

// 		ppe_iter_used = it+1;
// 		if(ppe_iter_used >= minit_safe && ppe_delta_p < PPE_TOLERANCE){
// 			ppe_converged = 1;
// 			break;
// 		}
// 	}

// 	// After loop: one extra PPE_prep with converged pres to refresh PPE1/PPE2,
// 	// then compute true residual |A*p - b|. bi is also refreshed (same value as last iter — u* unchanged).
// 	cudaMemsetAsync(g_dev_bi, 0, sizeof(Real)*num_part2);
// 	KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3, g_dev_bi);
// 	cudaDeviceSynchronize();
// 	KERNEL_compute_ppe_residual<<<b,t>>>(dev_SP1, g_dev_ppe_resid);
// 	cudaDeviceSynchronize();

// 	cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
// 	                       g_dev_bi, g_d_bi_max, num_part2);
// 	cudaDeviceSynchronize();
// 	cudaMemcpy(&ppe_bi_max, g_d_bi_max, sizeof(Real), cudaMemcpyDeviceToHost);

// 	cub::DeviceReduce::Max(g_dev_resid_storage, g_resid_storage_bytes,
// 	                       g_dev_ppe_resid, g_d_resid_max, num_part2);
// 	cudaDeviceSynchronize();
// 	cudaMemcpy(&ppe_resid_max, g_d_resid_max, sizeof(Real), cudaMemcpyDeviceToHost);

// 	// Per-step logging (only tid==0)
// 	if(tid==0){
// 		if(g_fp_ppe_log){
// 			fprintf(g_fp_ppe_log, "%d\t%.6g\t%d\t%.6e\t%.6e\t%.6e\t%d\n",
// 			        (int)count, (double)time, (int)ppe_iter_used,
// 			        (double)ppe_delta_p, (double)ppe_bi_max,
// 			        (double)ppe_resid_max, ppe_converged);
// 			fflush(g_fp_ppe_log);
// 		}
// 		if((count%50)==0){
// 			printf("PPE(implicit): iter=%d / %d, dp=%.3e, |bi|=%.3e, |A*p-b|=%.3e, conv=%d\n",
// 			       (int)ppe_iter_used, (int)maxit_safe,
// 			       (double)ppe_delta_p, (double)ppe_bi_max,
// 			       (double)ppe_resid_max, ppe_converged);
// 		}
// 	}
// #else
	//--- EISPH (Explicit ISPH, 1-sweep Jacobi: current production) ---
	for(int_t it=0;it<=0;it++){   // For iterative(implicit) pressure solver
		// PPE.cuh
		b.x=(num_part2-1)/t.x+1;
		KERNEL_PPE_prep3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
		KERNEL_PPE_calc3D<<<b,t>>>(dev_SP1);
		cudaDeviceSynchronize();
		// Pressure condition for wall(Neumann boundary)
		KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}
//#endif

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
	//cudaDeviceSynchronize();
	KERNEL_periodic_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	KERNEL_freeslip_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

// #if CASE_ID >= 500 && CASE_ID < 600
// 	// Case CS: sustain inlet buffer at log-law profile (uses us_update; with MOST off it stays at us_ini)
// 	if(dim==3) KERNEL_Recycling3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
// 	cudaDeviceSynchronize();
// #endif

	KERNEL_compute_vorticity3D<<<b,t>>>(g_str,g_end,dev_P1,dev_P3);
	cudaDeviceSynchronize();

	// Time-averaging accumulation (gated by AVG_START_TIME defined in Parameters.cuh / case_config.cuh)
	if(time >= (Real)AVG_START_TIME){
		KERNEL_accumulate_avg<<<b,t>>>(dev_P1);
		cudaDeviceSynchronize();
	}

//-------------------------------------------------------------------------------------------------
// LDM coupling — gated by AVG_START_TIME (skip predictor/advection/deposition/conc until steady state)
//-------------------------------------------------------------------------------------------------

	if(time >= (Real)AVG_START_TIME && (count%freq_LDM==0)){
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
		KERNEL_dry_deposition<<<b,t>>>(freq_LDM*dt, time, dev_LP1, (Real)MOST_z0);
		cudaDeviceSynchronize();

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

	if((count%freq_output)==0){
		// printf("save plot...........................\n");
		cudaMemcpy(file_P1,dev_P1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P2,dev_SP2,num_part2*sizeof(part2),cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P3,dev_P3,num_part2*sizeof(part3),cudaMemcpyDeviceToHost);

		if(count==0)save_vtk_bin_single_boundary(file_P1,file_P2,file_P3);
		if(count==0)save_vtk_bin_single_buildings(file_P1,file_P2,file_P3);
		//save_vtk_bin_single_fluid_CaseA(file_P1,file_P2,file_P3);     // Total fluid particles
		save_vtk_bin_single_fluid(file_P1,file_P2,file_P3);
		save_vtk_bin_single_fluid_condition(file_P1,file_P2,file_P3);    // Fluid particle with condition
		save_vtk_bin_single_vorticity(file_P1,file_P2,file_P3);
		
		// cudaMemcpy(file_P1,dev_SP1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
		// save_vtk_bin_single_test(file_P1,file_P2,file_P3);
		//cudaDeviceSynchronize();
	}

	if(count%100==0){
		cudaMemcpy(file_LP1,dev_LP1,num_part_LDM*sizeof(L_part1),cudaMemcpyDeviceToHost);
		save_vtk_bin_single_LDM(file_LP1);
	}

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
					printf("    [domain] Σk_turb=%.4e   Σ½|u|²=%.4e   avg_count=%u (AVG_START_TIME=%.2f)\n\n",
						sum_kt0, sum_uu0, h_ac, (double)AVG_START_TIME);
				}
			}
		}

		pthread_barrier_wait(&barrier);
		if(time>=time_end) break;

	}


	//-------------------------------------------------------------------------------------------------
	// ##. Save Restart File
	//-------------------------------------------------------------------------------------------------

	// if(ngpu==1) {
	//
	// 	save_restart(file_P1,file_P2,file_P3);
	// 	cudaMemcpy(file_P1,dev_SP1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
	// 	cudaMemcpy(file_P2,dev_SP2,num_part2*sizeof(part2),cudaMemcpyDeviceToHost);
	// 	cudaMemcpy(file_P3,dev_SP3,num_part2*sizeof(part3),cudaMemcpyDeviceToHost);
	//
	//
	// 	free(file_P2);
	// 	free(file_P3);
	// }

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
