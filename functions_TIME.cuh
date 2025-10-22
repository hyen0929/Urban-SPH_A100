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