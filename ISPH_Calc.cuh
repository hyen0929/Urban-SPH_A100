//-------------------------------------------------------------------------------------------------
// 모든 입자 계산을 수행하는 주함수
//-------------------------------------------------------------------------------------------------
void SOPHIA_single_ISPH(int_t*g_idx,int_t*p_idx,int_t*g_idx_in,int_t*p_idx_in,int_t*g_str,int_t*g_end,
	int_t*b_idx,int_t*b_idx_in,int_t*b_str,int_t*b_end,
	part1*dev_P1,part1*dev_SP1,part2*dev_P2,part2*dev_SP2,part3*dev_P3,
	int_t*p2p_af_in,int_t*p2p_idx_in,int_t*p2p_af,int_t*p2p_idx,
	void*dev_sort_storage,size_t*sort_storage_bytes,part1*file_P1,part2*file_P2,part3*file_P3,int tid)
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

		// b.x=(num_part2-1)/t.x+1;
		// if(dim==2) KERNEL_variable_smoothing_length2D<<<b,t>>>(g_str,g_end,dev_SP1);
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

//-------------------------------------------------------------------------------------------------
// 계산 준비: gradient correction, filter, reference density, p_type switch, penetration, density gradient
//-------------------------------------------------------------------------------------------------

	// 미분보정필터 계산 >> 진우 코드랑 비교 필요
	if(count==0){
		b.x=(num_part2-1)/t.x+1;
		//if(dim==2) KERNEL_clc_gradient_correction_2D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3);
		KERNEL_clc_gradient_correction_3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3);
		cudaDeviceSynchronize();
		KERNEL_clc_color_field3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3);
		cudaDeviceSynchronize();
		KERNEL_clc_IB_normal_vector3D<<<b,t>>>(g_str,g_end,dev_P1,dev_P3);
		cudaDeviceSynchronize();
	}

	// filter, reference density, p_type switch, penetration, normal gradient etc
	//if(dim==2) KERNEL_clc_prep2D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3,count);
	if(dim==3) KERNEL_clc_prep3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3,count,dt);
	cudaDeviceSynchronize();

	// Velocity condition for wall
	b.x=(num_part2-1)/t.x+1;
	if(noslip_bc==1){
		//if(dim==2) KERNEL_boundary2D<<<b,t>>>(g_str,g_end,dev_SP1);
		if(dim==3) KERNEL_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1);
		cudaDeviceSynchronize();
	}

	b.x=(num_part2-1)/t.x+1;
	if(dim==3) KERNEL_MOST_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3,count,dt);
	cudaDeviceSynchronize();

	// b.x=(num_part2-1)/t.x+1;
	// if(dim==3) KERNEL_turbulent_viscosity3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3,count,dt);
	// cudaDeviceSynchronize();

	if(count==0){
		//if(dim==2) KERNEL_Neumann_boundary2D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		if(dim==3) KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}

//-------------------------------------------------------------------------------------------------
// 압력힘을 제외한 힘 계산
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	if(dim==3) KERNEL_advection_force3D<<<b,t>>>(1,g_str,g_end,dev_SP1,dev_SP2,dev_P3,count);
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

	b.x=(num_part2-1)/t.x+1;
	if(dim==3) KERNEL_periodic_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	if(dim==3) KERNEL_freeslip_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// IBM_predictor<<<b,t>>>(dt_structure,time,dev_SP1,dev_SP2);
	// cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// 압력 계산
//-------------------------------------------------------------------------------------------------
	//For initial condition(Neumann boundary)

	for(int_t it=0;it<=0;it++){   // For iterative(implicit) pressure solver
		// PPE.cuh
		b.x=(num_part2-1)/t.x+1;
		if(dim==3) KERNEL_PPE3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
		// KERNEL_pressure_synchronize<<<b,t>>>(dev_SP1,dev_P1);
		// cudaDeviceSynchronize();

		// Pressure condition for wall(Neumann boundary)
		//if(dim==2) KERNEL_Neumann_boundary2D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		if(dim==3) KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();

		//if(dim==3) KERNEL_Neumann_boundary_buffer3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		if(dim==3) KERNEL_donothing_buffer3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		//if(dim==3) KERNEL_open_boundary_extrapolation_buffer3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
		cudaDeviceSynchronize();
	}

	// if(count%10==0){
	// 	// printf("save plot...........................\n");
	// 	cudaMemcpy(file_P1,dev_SP1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
	// 	cudaMemcpy(file_P2,dev_SP2,num_part2*sizeof(part2),cudaMemcpyDeviceToHost);
	// 	cudaMemcpy(file_P3,dev_P3,num_part2*sizeof(part3),cudaMemcpyDeviceToHost);
	// 	save_vtk_bin_single_flag(file_P1,file_P2,file_P3);
	//
	// 	// cudaMemcpy(file_P1,dev_SP1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
	// 	cudaDeviceSynchronize();
	//  }

	// b.x=(num_part2-1)/t.x+1;
	// if(dim==3) KERNEL_pressure_smoothing3D<<<b,t>>>(g_str,g_end,dev_SP1);
	// cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// 압력힘 계산
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	if(dim==3) KERNEL_pressureforce3D<<<b,t>>>(1,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	b.x=(num_part2-1)/t.x+1;
	// if(dim==3) KERNEL_MOST_boundary3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_P3,count,dt);
	cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// Immersed Boundary Method (IBM)
//-------------------------------------------------------------------------------------------------
	// IBM force interpolation (fb for IB)
	b.x=(num_part2-1)/t.x+1;
	if(dim==2) IBM_force_interpolation2D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	if(dim==3) IBM_force_interpolation3D<<<b,t>>>(dt,g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// IBM spreading interpolation (fb for fluid)
	if(dim==2) IBM_spreading_interpolation2D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	if(dim==3) IBM_spreading_interpolation3D<<<b,t>>>(g_str,g_end,dev_SP1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// 시간 적분 (Time Integration)
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	KERNEL_time_update_projection<<<b,t>>>(dt,dev_SP1,dev_P1,dev_SP2,dev_P2,dev_P3,count);
	cudaDeviceSynchronize();

	// KERNEL_time_update_buffer<<<b,t>>>(dt,dev_SP1,dev_P1,dev_SP2,dev_P2,dev_P3,count);
	// cudaDeviceSynchronize();

	IBM_corrector<<<b,t>>>(dt_structure,time,dev_SP1,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();

	// b.x=(num_part2-1)/t.x+1;
	// if(dim==3) KERNEL_xsph3D<<<b,t>>>(g_str,g_end,dev_P1);
	// cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// Open Boudnary && APR
//-------------------------------------------------------------------------------------------------

	b.x=(num_part2-1)/t.x+1;
	if(dim==3) KERNEL_Neumann_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
 	//if(dim==3) KERNEL_Neumann_boundary_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	if(dim==3) KERNEL_donothing_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	//if(dim==3) KERNEL_open_boundary_extrapolation_buffer3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	if(dim==3) KERNEL_periodic_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	if(dim==3) KERNEL_freeslip_boundary3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	cudaDeviceSynchronize();
	// if(dim==3) KERNEL_Recycling3D<<<b,t>>>(g_str,g_end,dev_P1,dev_SP2,dev_P3);
	// cudaDeviceSynchronize();

//-------------------------------------------------------------------------------------------------
// 출력
//-------------------------------------------------------------------------------------------------

	if((count%freq_output)==0){
		// printf("save plot...........................\n");
		cudaMemcpy(file_P1,dev_P1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P2,dev_SP2,num_part2*sizeof(part2),cudaMemcpyDeviceToHost);
		cudaMemcpy(file_P3,dev_P3,num_part2*sizeof(part3),cudaMemcpyDeviceToHost);
		save_vtk_bin_single_flag(file_P1,file_P2,file_P3);

		// cudaMemcpy(file_P1,dev_SP1,num_part2*sizeof(part1),cudaMemcpyDeviceToHost);
		// save_vtk_bin_single_test(file_P1,file_P2,file_P3);
		cudaDeviceSynchronize();
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

	//-------------------------------------------------------------------------------------------------
	// Device/GPU 변수 선언 및 메모리 할당
	//-------------------------------------------------------------------------------------------------

	// NNPS 관련 변수
	int_t*g_idx,*p_idx,*g_idx_in,*p_idx_in,*g_str,*g_end;
	int_t*b_idx,*b_idx_in,*b_str,*b_end;

	num_blocks=NBI*NBJ*NBK;

	// 주요 입자 변수
	part1*dev_P1,*dev_SP1;
	part2*dev_P2,*dev_SP2;
	part3*dev_SP3;

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

	// NNPS 메모리 초기화
	cudaMemset(g_idx_in,0,sizeof(int_t)*num_part2);
	cudaMemset(p_idx_in,0,sizeof(int_t)*num_part2);
	cudaMemset(g_idx,0,sizeof(int_t)*num_part2);
	cudaMemset(p_idx,0,sizeof(int_t)*num_part2);
	cudaMemset(g_str,cu_memset,sizeof(int_t)*num_cells);
	cudaMemset(g_end,0,sizeof(int_t)*num_cells);

	// Device 입자 메모리 초기화
	cudaMemset(dev_P1,0,sizeof(part1)*num_part2);
	cudaMemset(dev_SP1,0,sizeof(part1)*num_part2);
	cudaMemset(dev_P2,0,sizeof(part2)*num_part2);
	cudaMemset(dev_SP2,0,sizeof(part2)*num_part2);
	cudaMemset(dev_SP3,0,sizeof(part3)*num_part2);

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

	if(tid==0){
		printf("\n-----------------------------------------------------------\n");
		printf("GPU Domain Division Success\n");
		printf("-----------------------------------------------------------\n\n");
	}


	//-------------------------------------------------------------------------------------------------
	// 화면 출력용 기타 변수들 정의 및 메모리 할당 (최대속도, 최대힘 등)
	//-------------------------------------------------------------------------------------------------

	// host
	Real *max_umag0,*max_rho0,*max_ftotal0;
	max_umag0=(Real*)malloc(sizeof(Real));
	max_rho0=(Real*)malloc(sizeof(Real));
	max_ftotal0=(Real*)malloc(sizeof(Real));
	max_umag0[0]=max_ftotal0[0]=max_rho0[0]=0.0;

	// device
	Real*max_rho,*max_umag,*max_ft,*d_max_umag0,*d_max_rho0,*d_max_ftotal0;
	cudaMalloc((void**)&max_rho,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_umag,sizeof(Real)*num_part2);
	cudaMalloc((void**)&max_ft,sizeof(Real)*num_part2);
	cudaMalloc((void**)&d_max_umag0,sizeof(Real));
	cudaMalloc((void**)&d_max_rho0,sizeof(Real));
	cudaMalloc((void**)&d_max_ftotal0,sizeof(Real));
	cudaMemset(max_umag,0,sizeof(Real)*num_part2);
	cudaMemset(max_rho,0,sizeof(Real)*num_part2);
	cudaMemset(max_ft,0,sizeof(Real)*num_part2);
	cudaMemset(d_max_umag0,0,sizeof(Real));
	cudaMemset(d_max_rho0,0,sizeof(Real));
	cudaMemset(d_max_ftotal0,0,sizeof(Real));


	//-------------------------------------------------------------------------------------------------
	// 정렬(Sorting)을 위한 CUB 라이브러리 변수 준비
	//-------------------------------------------------------------------------------------------------

	// Sorting & Max variable to use CUB Library
	void*dev_sort_storage=NULL;
	void*dev_max_storage=NULL;
	size_t sort_storage_bytes=0;
	size_t max_storage_bytes=0;

	// Determine Sorting & Maximum Value Setting for Total Particle Data
	cub::DeviceRadixSort::SortPairs(dev_sort_storage,sort_storage_bytes,g_idx_in,g_idx,p_idx_in,p_idx,num_part2);
	cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_umag,d_max_umag0,num_part2);
	cudaDeviceSynchronize();
	cudaMalloc((void**)&dev_sort_storage,sort_storage_bytes);
	cudaMalloc((void**)&dev_max_storage,max_storage_bytes);
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
		SOPHIA_single_ISPH(g_idx,p_idx,g_idx_in,p_idx_in,g_str,g_end,b_idx,b_idx_in,b_str,b_end,dev_P1,dev_SP1,dev_P2,dev_SP2,dev_SP3,
					p2p_af_in,p2p_idx_in,p2p_af,p2p_idx,dev_sort_storage,&sort_storage_bytes,file_P1,file_P2,file_P3,tid);
		}

		//-------------------------------------------------------------------------------------------------
		// Time-step Control
		//-------------------------------------------------------------------------------------------------
		if(tid==0){
			time+=dt;
			count++;

			//timestep is updated every 10 steps ------------ estimate new timestep (Goswami & Pajarola(2011))
			if((count%freq_output)==0){
				dim3 t,b;
				t.x=128;
				b.x=(num_part2-1)/t.x+1;
				kernel_copy_max<<<b,t>>>(dev_P1,dev_SP2,dev_SP3,max_rho,max_ft,max_umag);
				cudaDeviceSynchronize();

				// Find Max Velocity & Force using CUB - TID=0
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_umag,d_max_umag0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_rho,d_max_rho0,num_part2);
				cub::DeviceReduce::Max(dev_max_storage,max_storage_bytes,max_ft,d_max_ftotal0,num_part2);
				cudaDeviceSynchronize();
				cudaMemcpy(max_umag0,d_max_umag0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_rho0,d_max_rho0,sizeof(Real),cudaMemcpyDeviceToHost);
				cudaMemcpy(max_ftotal0,d_max_ftotal0,sizeof(Real),cudaMemcpyDeviceToHost);

				end=clock();
				cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
				printf("Simulation time: %5.2f sec\n", time);
				printf("Computing time: %5.2f sec\n\n", cpu_time_used);
				printf("%d\t tu_max=%5.2f\tftotal_max=%5.5f\n\n",count,max_umag0[0],max_ftotal0[0]);
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
	cudaFree(g_idx);
	cudaFree(p_idx);
	cudaFree(g_idx_in);
	cudaFree(p_idx_in);
	cudaFree(g_str);
	cudaFree(g_end);
	cudaFree(dev_P1);
	cudaFree(dev_SP1);
	cudaFree(dev_P2);
	cudaFree(dev_SP2);
	cudaFree(dev_SP3);
	cudaFree(max_umag);
	cudaFree(max_rho);
	cudaFree(max_ft);
	cudaFree(d_max_umag0);
	cudaFree(d_max_rho0);
	cudaFree(d_max_ftotal0);
	cudaFree(dev_sort_storage);
	cudaFree(dev_max_storage);
	cudaStreamDestroy(str1[tid]);
	cudaStreamDestroy(str2[tid]);
	pthread_barrier_wait(&barrier);

	return 0;
}
