//---------------------------------------------------------------------------------------------------
// SOPHIA_gpu Version 2.0: Smoothed Particle Hydrodynamics code In Advanced nuclear safety
// Developed by Eung Soo Kim, Young Beom Jo, So Hyun Park, Hae Yoon Choi in 2017
// ENERGy SYSTEM LABORATORY, NUCLEAR ENGINEERING DEPARTIMENT, SEOUL NATIONAL UNIVERSITY, SOUTH KOREA
//---------------------------------------------------------------------------------------------------
// Optimized by Dong Hak Lee, Yong Woo Sim in 2018 (2018.01.08)
// Copyright 2018(C) CoColink Inc.
//---------------------------------------------------------------------------------------------------
// Multi-GPU Optimized by Dong Hak Lee, Yong Woo Sim in 2019 (2019.01.09)
// Copyright 2019(C) CoColink Inc.
//---------------------------------------------------------------------------------------------------
// Code Restructured by Eung Soo Kim, Hee Sang Yoo, Young Beom Jo, Hae Yoon Choi, Su-San Park, Jin Woo Kim, Yelyn Ahn, Tae Soo Choi in (2019.08.08)
// Copyright 2019(C) ESLAB, SEOUL NATIONAL UNIVERSITY, SOUTH KOREA.
//---------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <string>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <cub/cub.cuh>
#include <float.h>

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// #include "device_functions.h"
#include "Cuda_Error.cuh"

#include "Variable_Type.cuh"
#include "Parameters.cuh"
#include "class_Cuda_Particle_Array.cuh"


//---------------------------------------------------------------
// 전역변수 선언
//---------------------------------------------------------------

// Solver 파라미터 선언 (function_init.cuh)_______
// host
int_t vii[vii_size];
Real vif[vif_size];

// device - solver options
__constant__ int_t k_vii[vii_size];
__constant__ Real k_vif[vif_size];

// device - table
__constant__ Real k_Tab_T[table_size]; 	//table_size는 Parameter.cuh에 정의됨
__constant__ Real k_Tab_h[table_size];
__constant__ Real k_Tab_k[table_size];
__constant__ Real k_Tab_cp[table_size];
__constant__ Real k_Tab_vis[table_size];
__constant__ int k_table_index[10]; 		// table의 시작 주소
__constant__ int k_table_size[10]; 	 		// table별 데이터 수

// 입자 선언 _________
// for main domain
part1*HP1;							// host 전체입자
part1*DHP1[Max_GPU];		// gpu당 할당할 host 입자

// for data exchange
part1*send_P1[Max_GPU],*send_rSP1[Max_GPU],*send_lSP1[Max_GPU],*recv_P1[Max_GPU];
p2p_part3*send_P3[Max_GPU],*send_rSP3[Max_GPU],*send_lSP3[Max_GPU],*recv_P3[Max_GPU];

// table - host
Real host_Tab_T[table_size];
Real host_Tab_h[table_size];
Real host_Tab_k[table_size];
Real host_Tab_cp[table_size];
Real host_Tab_vis[table_size];

int host_table_index[10];
int host_table_size[10];

// 병렬 분기 _________
cudaStream_t str1[Max_GPU];
cudaStream_t str2[Max_GPU];
pthread_barrier_t barrier;

// Open Boundary _______
Real space;				  // open boundary inlet을 위한 가상 격자의 간격
int Nsx=0;					// x축 방향 격자 수
int Nsy=0;
int Nsz=0;					// y축 방향 격자 수
int buffer_size=0;  // 생성 입자룰 위한 여분의 메모리 크기

__device__ int num_buffer[1];	// buffer 입자의 총 수를 저장하는 전역 변수 (Open Boundary/APR 모델 사용시 필요)

// Plot Data ______
int num_plot_data;				// plot data 개수
char plot_data[20][20]; 	// plot 할 변수


//---------------------------------------------------------------


#include "function_init.cuh"
#include "functions_NNPS.cuh"
#include "functions_KNL.cuh"
#include "functions_PROP.cuh"
#include "functions_PREP.cuh"
#include "functions_IBM.cuh"
#include "functions_INTERACTION.cuh"
#include "functions_ENERGY.cuh"
#include "functions_TIME.cuh"
#include "functions_PPE.cuh"
#include "functions_OUTPUT.cuh"
#include "functions_BC.cuh"
#include "functions_MOST.cuh"
// #include "functions_APR.cuh"
//#include "functions_APR_etc.cuh"
#include "functions_TIME_ISPH.cuh"
#include "ISPH_Calc.cuh"
#include "ISPH.cuh"

////////////////////////////////////////////////////////////////////////
int main(int argc,char**argv)
{
	memset(vii,0,sizeof(int_t)*vii_size);
	memset(vif,0,sizeof(Real)*vif_size);

	ngpu=atoi(argv[1]);

	char fn[64],fn2[64];
	strcpy(fn,"./input/solv.txt");
	strcpy(fn2,"./input/data.txt");

	read_solv_input(vii,vif,fn);
	read_table(fn2);

	switch(solver_type){
		case Isph:
			ISPH(vii,vif);
			break;
		default:
			ISPH(vii,vif);
	}

	return 0;
}
