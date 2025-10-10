#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
__host__ __device__ int clc_idx_insert(Real xc, Real yc, Real zc, Real tdelx, int Nsx_, int Nsz_)
{
	int Ix,Iy, Iz;
	int idx_1D;

	Ix = round((xc-k_x_min)/tdelx);
	Iy = round((yc-k_y_min)/tdelx);
	Iz = round((zc-k_z_min)/tdelx);

	if (k_dim==2) idx_1D=(k_num_part2-1)-(Iy);
	if (k_dim==3) idx_1D=(k_num_part2-1)-(Iy*Nsz_+Iz);

	return idx_1D;
}
////////////////////////////////////////////////////////////////////////
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
__host__ __device__ int clc_idx_insert2(Real xc, Real yc, Real zc, Real tdelx, int Nsx_, int Nsz_)
{
	int Ix,Iy, Iz;
	int idx_1D;

	Ix = round((xc-k_x_min)/tdelx);
	Iy = round((yc-k_y_min)/tdelx);
	Iz = round((zc-k_z_min)/tdelx);

	if (k_dim==2) idx_1D=(k_num_part2-1)-(Iy);
	if (k_dim==3) idx_1D=(k_num_part2-1)-(Iy*Nsx_+Ix);

	return idx_1D;
}

////////////////////////////////////////////////////////////////////////
// open_boundary
// open_boundary
// __global__ void KERNEL_time_update_buffer(const Real tdt,part1*P1,part1*TP1,part2*P2,part3*P3,Real space_,int Nsx_, int Nsz_)
// {
// 	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part2) return;
// 	if(P1[i].i_type>2) return;
// 	if(P1[i].p_type<1) return;
// 	// if((P1[i].eli==1.0)&(P1[i].eli==1.0)) return;
//
// 	Real tx0,ty0,tz0,xc,yc,zc;									// position
// 	Real tux0,tuy0,tuz0,uxc,uyc,uzc;						// velocity
// 	Real t_dt=tdt;
// 	Real tmp_h=P1[i].h;
// 	int_t idx_insert;
//
// 	int_t p_type=P1[i].p_type;
// 	int_t i_type=P1[i].i_type;
// 	int_t buffer_type=P1[i].buffer_type;
//
// 	xc=TP1[i].x;												// correct x-directional position
// 	yc=TP1[i].y;												// correct Y-directional position
// 	zc=TP1[i].z;												// correct Z-directional position
//
// 	// inlet
// 	if ((buffer_type==1)&&(i_type==2)) {
// 		if (xc>=L2) {
// 			TP1[i].i_type=1;
// 			TP1[i].buffer_type=0;
//
// 			idx_insert=clc_idx_insert(xc,yc, zc,tmp_h/h_coeff,Nsx_,Nsz_);
//
// 			TP1[idx_insert]=P1[i];
//
// 			TP1[idx_insert].x=L1;
// 			TP1[idx_insert].y=yc;
// 			//TP1[idx_insert].y=yc-L2+L1-tmp_h/3.000;
// 			TP1[idx_insert].z=zc;
//
// 		}
// 	}
//
// 		// if (xc>=L3){
// 		// 		TP1[i].buffer_type=2;
// 		// 		TP1[i].i_type=2;
// 		// 	}
// 		//
// 		// if (yc>=L4 && TP1[i].buffer_type==Outlet){
// 			// 	TP1[i].i_type=3;
// 			// }
// 		if (xc>=L4){
// 			TP1[i].i_type=3;
// 		}
//
// 		if (TP1[i].buffer_type==0 && xc<L1){
// 			TP1[i].i_type=3;
// 		}
// 	}
////////////////////////////////////////////////////////////////////////
