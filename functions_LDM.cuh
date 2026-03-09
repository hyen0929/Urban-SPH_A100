//////////////////////// LDM 관련 커널 함수들 ////////////////////////
__device__ Real getRand(curandState *s){
    Real rand_int = curand_uniform(s);
    return rand_int;
}

__device__ Real GaussianRand(curandState *s, Real mu, Real stdv){
	Real u1 = getRand(s);
	Real u2 = getRand(s);
	
	Real mag = stdv*sqrt(-2.0*log(u1));

	return mag*cos(2*PI*u2)+mu;
}

__global__ void particle_releasing_Dispersion(Real tdt,int freq,Real ttime,Real tend,L_part1*LP1) {

	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==1) return;

	// if (ttime<=ttend){
	// 	for(j=k_num_part_LDM*(ttime)/(ttend); j<k_num_part_LDM*(ttime+freq*tdt)/(ttend); j++){
	// 		LP1[j].i_type=1;
	// 	}
	// }

	if(ttime>=LP1[i].t_release){
		LP1[i].i_type=1;
	}
}

__global__ void particle_releasing_Dispersion_Gaussian(Real tdt,int freq,Real ttime,Real tend,L_part1*LP1, Real src_radius) {

	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==1) return;

	if(ttime>=LP1[i].t_release){
		LP1[i].i_type=1;

		// RNG 초기화
		unsigned long long time_hash = (unsigned long long)(ttime * 1e6);
		if(LP1[i].seed < 1){
			LP1[i].seed = i * 321215 + time_hash + 36893211;
		}
		curandState ss;
		curand_init(LP1[i].seed, 0, 0, &ss);

		// 3σ = src_radius → 99.7% 입자가 반경 내 분포
		Real sigma = src_radius / 3.0;

		Real dx = GaussianRand(&ss, 0.0, sigma);
		Real dy = GaussianRand(&ss, 0.0, sigma);
		Real dz = GaussianRand(&ss, 0.0, sigma);

		// 반경 초과 시 clamp
		Real dist = sqrt(dx*dx + dy*dy + dz*dz) + 1e-20;
		if(dist > src_radius){
			Real scale = src_radius / dist;
			dx *= scale;
			dy *= scale;
			dz *= scale;
		}

		LP1[i].x += dx;
		LP1[i].y += dy;
		LP1[i].z += dz;

		LP1[i].seed = (unsigned long long)curand(&ss);
	}
}

__global__ void particle_releasing_Dispersion_Kernel(Real tdt,int freq,Real ttime,Real tend,L_part1*LP1) {

	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==1) return;

	if(ttime>=LP1[i].t_release){
		LP1[i].i_type=1;

		// RNG 초기화
		unsigned long long time_hash = (unsigned long long)(ttime * 1e6);
		if(LP1[i].seed < 1){
			LP1[i].seed = i * 321215 + time_hash + 36893211;
		}
		curandState ss;
		curand_init(LP1[i].seed, 0, 0, &ss);

		// SPH 커널 기반 rejection sampling
		Real h_src = 2.0*LP1[i].h;
		Real search_range = k_search_kappa * h_src;  // compact support 반경
		Real tmp_A = calc_tmpA(h_src);
		Real w_max = calc_kernel_wij(tmp_A, h_src, 1e-20);  // W(0) = peak

		Real dx, dy, dz, dist, wij;
		for(int iter = 0; iter < 100; iter++){
			dx = (2.0 * curand_uniform(&ss) - 1.0) * search_range;
			dy = (2.0 * curand_uniform(&ss) - 1.0) * search_range;
			dz = (2.0 * curand_uniform(&ss) - 1.0) * search_range;
			dist = sqrt(dx*dx + dy*dy + dz*dz);

			if(dist >= search_range) continue;

			wij = calc_kernel_wij(tmp_A, h_src, dist);
			if(curand_uniform(&ss) < wij / w_max) break;  // accept
		}

		LP1[i].x += dx;
		LP1[i].y += dy;
		LP1[i].z += dz;

		LP1[i].seed = (unsigned long long)curand(&ss);
	}
}

__global__ void KERNEL_clc_predictor_Dispersion(Real tdt,Real ttime,L_part1*LP1,L_part2*LP2)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;

	if(i>=k_num_part_LDM) return;
	//if(LP1[i].p_type==0) return;
	if(LP1[i].i_type==3) return;

	int_t p_typei;
	double tmp;

	p_typei=LP1[i].p_type;

	LP2[i].x0=LP1[i].x;															// update x-directional position
	LP2[i].y0=LP1[i].y;															// update y-directional position
	LP2[i].z0=LP1[i].z;															// update z-directional position

	LP2[i].ux0=LP1[i].ux;														// update x-directional velocity
	LP2[i].uy0=LP1[i].uy;														// update y-directional velocity
	LP2[i].uz0=LP1[i].uz;														// update z-directional velocity

}

__global__ void KERNEL_time_update_CFD_3D(const Real tdt, Real tt, Real tend, L_part1*LP1, L_part2*LP2, cfd_mesh* dev_cfd, Real xmin, Real xmax, Real xspan, Real ymin, Real ymax, Real yspan, Real zmin, Real zmax, Real zspan)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;


	unsigned long long seed;
	curandState ss;
	double tmp;

	// ========================================================================
	// Time-dependent RNG initialization (critical for stochastic dispersion)
	// ========================================================================
	// Add time-dependent hash to ensure different random sequences at each timestep
	// Without this, particles would follow identical stochastic paths every timestep
	unsigned long long time_hash = (unsigned long long)(tt * 1e6);

	if (LP1[i].seed < 1) {
		// First time: initialize with particle ID and current time
		LP1[i].seed = i * 321215 + time_hash + 36893211;
		seed = LP1[i].seed;
	}
	else {
		// Subsequent times: evolve seed with LCG + time mixing
		seed = LP1[i].seed;
		seed = (seed * 1103515245 + time_hash + 12345) % 2147483647;
		LP1[i].seed = seed;
	}

	// Initialize curand state with time-dependent offset
	curand_init(seed, i, (unsigned long long)(tt * 100), &ss);

	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==3) return;

	Real tx0,ty0,tz0,ts0,xc,yc,zc,sc;									// position
	Real tux0,tuy0,tuz0,uxc,uyc,uzc;						// velocity
	Real tux,tuy,tuz;														// velocity
	Real dux_dt,duy_dt,duz_dt;									// accleration (time derivative of velocity)

	Real k_turb, e_turb, yplus;
	Real xwind, ywind, zwind;

	tx0=LP2[i].x0;														// x-directional initial position
	ty0=LP2[i].y0;														// x-directional initial position
	tz0=LP2[i].z0;														// x-directional initial position
	ts0=LP2[i].sigma0;

	tux0=LP2[i].ux0;						// x-directional initial velocity
	tuy0=LP2[i].uy0;						// y-directional initial velocity
	tuz0=LP2[i].uz0;						// z-directional initial velocity

	int numx = round((xmax-xmin)/xspan);
	int numy = round((ymax-ymin)/yspan);
	int numz = round((zmax-zmin)/zspan);

	int xidx, yidx, zidx;

	//  tmpXYZ; (ex. tmLP101 = x-upside, y-downside, z-upside)

	// double tmp000, tmp001, tmp010, tmp011;
	// double tmLP100, tmLP101, tmLP110, tmLP111;

	double vec001, vec010, vec100;

	int xbasis, ybasis, zbasis;

	xidx = floor((tx0-xmin)/(xspan)+0.00001*SIZ);
	if (xidx==round((xmax)/xspan)) xidx = round((xmax)/xspan)-1;
	if (xidx<0) xidx=0;
	yidx = floor((ty0-ymin)/(yspan)+0.00001*SIZ);
	if (yidx==round((ymax)/yspan)) yidx = round((ymax)/yspan)-1;
	if (yidx<0) yidx=0;
	zidx = floor((tz0-zmin)/(zspan)+0.00001*SIZ);
	if (zidx==round((zmax)/zspan)) zidx = round((zmax)/zspan)-1;
	if (zidx<0) zidx=0;

	xbasis = (tx0>=(xidx+0.5)*xspan)-(tx0<(xidx+0.5)*xspan);
	ybasis = (ty0>=(yidx+0.5)*yspan)-(ty0<(yidx+0.5)*yspan);
	zbasis = (tz0>=(zidx+0.5)*zspan)-(tz0<(zidx+0.5)*zspan);


	if(zidx+1*zbasis>=numz||zidx+1*zbasis<0) vec001 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].u)<1e-9) vec001 = 0; // 0-value data handling
	else vec001 = (dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].u-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].u)/2;

	if(yidx+1*ybasis>=numy||yidx+1*ybasis<0) vec010 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].u)<1e-9) vec010 = 0; // 0-value data handling
	else vec010 = (dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].u-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].u)/2;

	if(xidx+1*xbasis>=numx||xidx+1*xbasis<0) vec100 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].u)<1e-9) vec100 = 0; // 0-value data handling
	else vec100 = (dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].u-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].u)/2;

	xwind = dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].u +
			vec001*(tz0-(zidx+0.5)*zspan)*zbasis/(0.5*zspan) +
			vec010*(ty0-(yidx+0.5)*yspan)*ybasis/(0.5*yspan) +
			vec100*(tx0-(xidx+0.5)*xspan)*xbasis/(0.5*xspan);


	if(zidx+1*zbasis>=numz||zidx+1*zbasis<0) vec001 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].v)<1e-9) vec001 = 0; // 0-value data handling
	else vec001 = (dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].v-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].v)/2;

	if(yidx+1*ybasis>=numy||yidx+1*ybasis<0) vec010 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].v)<1e-9) vec010 = 0; // 0-value data handling
	else vec010 = (dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].v-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].v)/2;

	if(xidx+1*xbasis>=numx||xidx+1*xbasis<0) vec100 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].v)<1e-9) vec100 = 0; // 0-value data handling
	else vec100 = (dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].v-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].v)/2;

	ywind = dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].v +
			vec001*(tz0-(zidx+0.5)*zspan)*zbasis/(0.5*zspan) +
			vec010*(ty0-(yidx+0.5)*yspan)*ybasis/(0.5*yspan) +
			vec100*(tx0-(xidx+0.5)*xspan)*xbasis/(0.5*xspan);


	if(zidx+1*zbasis>=numz||zidx+1*zbasis<0) vec001 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].w)<1e-9) vec001 = 0; // 0-value data handling
	else vec001 = (dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].w-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].w)/2;

	if(yidx+1*ybasis>=numy||yidx+1*ybasis<0) vec010 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].w)<1e-9) vec010 = 0; // 0-value data handling
	else vec010 = (dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].w-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].w)/2;

	if(xidx+1*xbasis>=numx||xidx+1*xbasis<0) vec100 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].w)<1e-9) vec100 = 0; // 0-value data handling
	else vec100 = (dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].w-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].w)/2;

	zwind = dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].w +
			vec001*(tz0-(zidx+0.5)*zspan)*zbasis/(0.5*zspan) +
			vec010*(ty0-(yidx+0.5)*yspan)*ybasis/(0.5*yspan) +
			vec100*(tx0-(xidx+0.5)*xspan)*xbasis/(0.5*xspan);


	if(zidx+1*zbasis>=numz||zidx+1*zbasis<0) vec001 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].k_turb)<1e-9) vec001 = 0; // 0-value data handling
	else vec001 = (dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].k_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].k_turb)/2;

	if(yidx+1*ybasis>=numy||yidx+1*ybasis<0) vec010 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].k_turb)<1e-9) vec010 = 0; // 0-value data handling
	else vec010 = (dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].k_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].k_turb)/2;

	if(xidx+1*xbasis>=numx||xidx+1*xbasis<0) vec100 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].k_turb)<1e-9) vec100 = 0; // 0-value data handling
	else vec100 = (dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].k_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].k_turb)/2;

	k_turb = dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].k_turb +
			vec001*(tz0-(zidx+0.5)*zspan)*zbasis/(0.5*zspan) +
			vec010*(ty0-(yidx+0.5)*yspan)*ybasis/(0.5*yspan) +
			vec100*(tx0-(xidx+0.5)*xspan)*xbasis/(0.5*xspan);


	if(zidx+1*zbasis>=numz||zidx+1*zbasis<0) vec001 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].e_turb)<1e-9) vec001 = 0; // 0-value data handling
	else vec001 = (dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx+1*zbasis)].e_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].e_turb)/2;

	if(yidx+1*ybasis>=numy||yidx+1*ybasis<0) vec010 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].e_turb)<1e-9) vec010 = 0; // 0-value data handling
	else vec010 = (dev_cfd[(xidx)*numy*numz+(yidx+1*ybasis)*numz+(zidx)].e_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].e_turb)/2;

	if(xidx+1*xbasis>=numx||xidx+1*xbasis<0) vec100 = 0; // Avoid Index error
	else if (abs(dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].e_turb)<1e-9) vec100 = 0; // 0-value data handling
	else vec100 = (dev_cfd[(xidx+1*xbasis)*numy*numz+(yidx)*numz+(zidx)].e_turb-dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].e_turb)/2;

	e_turb = dev_cfd[(xidx)*numy*numz+(yidx)*numz+(zidx)].e_turb +
			vec001*(tz0-(zidx+0.5)*zspan)*zbasis/(0.5*zspan) +
			vec010*(ty0-(yidx+0.5)*yspan)*ybasis/(0.5*yspan) +
			vec100*(tx0-(xidx+0.5)*xspan)*xbasis/(0.5*xspan);


#ifdef USE_LDM_ENHANCED_TURBULENCE
	// ========================================================================
	// Adaptive Lagrangian Timescales (Sawford 1991, Pope 2000)
	// ========================================================================
	Real C0 = 2.1;  // Sawford (1991) recommendation for atmospheric dispersion
#else
	Real C0 = 1.76;  // Original value
#endif

	// ========================================================================
	// Turbulent Velocity Scale and Lagrangian Timescales
	// ========================================================================
	Real Sigvel = sqrt(2*k_turb/3);

	// Real Tu = 1/(0.5+0.75*C0)*(k_turb/(e_turb+1e-10));
	// Real Tv = 1/(0.5+0.75*C0)*(k_turb/(e_turb+1e-10));
	// Real Tw = 1/(0.5+0.75*C0)*(k_turb/(e_turb+1e-10));

#ifdef USE_LDM_ENHANCED_TURBULENCE
	// Adaptive timescales with anisotropy (vertical vs horizontal)
	Real e_safe = fmax(e_turb, (Real)1e-10);
	Real k_safe = fmax(k_turb, (Real)1e-10);

	// Horizontal timescales (typically longer in stratified atmosphere)
	Real Tu_base = k_safe / (C0 * e_safe);
	Real Tv_base = k_safe / (C0 * e_safe);

	// Vertical timescale (faster dissipation due to stratification)
	Real Tw_base = k_safe / (C0 * e_safe) * (Real)0.5;

	// Apply limits for numerical stability
	Real Tu = fmax((Real)0.01, fmin((Real)1000.0, Tu_base));
	Real Tv = fmax((Real)0.01, fmin((Real)1000.0, Tv_base));
	Real Tw = fmax((Real)0.005, fmin((Real)100.0, Tw_base));
#else
	// Original isotropic timescales
	Real Tu = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);
	Real Tv = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);
	Real Tw = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);
#endif

	//printf("xidx = %d, yidx = %d, xwind = %f, ywind = %f\n", xidx, yidx, xwind, ywind);
	// printf("k_turb = %f, e_turb = %f\n", k_turb, e_turb);
	//printf("Time scale = %f\n", Tu);

	Real KH, KV;

	Real t_dt=tdt;

	Real Ru = exp(-t_dt/Tu);
	Real Rv = exp(-t_dt/Tv);
	Real Rw = exp(-t_dt/Tw);

	Real mu = 0.0;
	Real stdv = 1;

	int_t p_type_i=LP1[i].p_type;

	// Real Ru = 0;
	// Real Rv = 0;
	// Real Rw = 0;


	//printf("sig = %f\n", Sigvel);

	// ========================================================================
	// Langevin Equation Update with Thomson Drift Correction (Optional)
	// ========================================================================
#ifdef USE_LDM_ENHANCED_TURBULENCE
	// Thomson well-mixed condition: add drift term for non-homogeneous turbulence
	// This ensures the Fokker-Planck equation is satisfied in variable k-field

	// Compute gradients of k_turb using finite differences
	Real grad_k_x = (Real)0.0;
	Real grad_k_y = (Real)0.0;
	Real grad_k_z = (Real)0.0;

	if (xidx + 1 < numx && xidx > 0) {
		int idx_xp = (xidx + 1) * numy * numz + yidx * numz + zidx;
		int idx_xm = (xidx - 1) * numy * numz + yidx * numz + zidx;
		grad_k_x = (dev_cfd[idx_xp].k_turb - dev_cfd[idx_xm].k_turb) / ((Real)2.0 * xspan);
	}

	if (yidx + 1 < numy && yidx > 0) {
		int idx_yp = xidx * numy * numz + (yidx + 1) * numz + zidx;
		int idx_ym = xidx * numy * numz + (yidx - 1) * numz + zidx;
		grad_k_y = (dev_cfd[idx_yp].k_turb - dev_cfd[idx_ym].k_turb) / ((Real)2.0 * yspan);
	}

	if (zidx + 1 < numz && zidx > 0) {
		int idx_zp = xidx * numy * numz + yidx * numz + (zidx + 1);
		int idx_zm = xidx * numy * numz + yidx * numz + (zidx - 1);
		grad_k_z = (dev_cfd[idx_zp].k_turb - dev_cfd[idx_zm].k_turb) / ((Real)2.0 * zspan);
	}

	// Thomson drift correction: b = -T_L * ∂(σ²)/∂x = -T_L * (2/3) * ∂k/∂x
	Real drift_u = -Tu * ((Real)2.0 / (Real)3.0) * grad_k_x;
	Real drift_v = -Tv * ((Real)2.0 / (Real)3.0) * grad_k_y;
	Real drift_w = -Tw * ((Real)2.0 / (Real)3.0) * grad_k_z;
#else
	Real drift_u = (Real)0.0;
	Real drift_v = (Real)0.0;
	Real drift_w = (Real)0.0;
#endif

	// Stochastic diffusion term (standard Langevin)
	dux_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Ru*Ru);			// x-directional acceleration
	duy_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Rv*Rv);			// y-directional acceleration
	duz_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Rw*Rw);			// z-directional acceleration

	// Combine deterministic drift + Thomson correction + stochastic diffusion
	uxc=tux0*Ru+dux_dt + drift_u*t_dt;										// correct x-directional velocity
	uyc=tuy0*Rv+duy_dt + drift_v*t_dt;										// correct y-directional velocity
	uzc=tuz0*Rw+duz_dt + drift_w*t_dt;


	if(tx0>xmax || tx0<xmin || ty0>ymax || ty0<ymin || tz0>zmax || tz0<zmin){
		xwind = 0;
		ywind = 0;
		zwind = 0;
		uxc = 0;
		uyc = 0;
		uzc = 0;
	}

	// ========================================================================
	// Position Update with Ground Reflection Boundary Condition
	// ========================================================================
	xc=tx0+(xwind+uxc)*(t_dt);												// correct x-directional position
	yc=ty0+(ywind+uyc)*(t_dt);												// correct Y-directional position
	zc=tz0+(zwind+uzc)*(t_dt);

#ifdef USE_LDM_ENHANCED_TURBULENCE
	// Ground reflection with roughness-based dissipation
	Real z_ground = zmin;
	Real roughness_length = (Real)0.01;  // Typical urban roughness ~0.01-0.1 m

	// If particle hits ground, reflect with partial energy loss
	if (zc < z_ground + roughness_length && (zwind + uzc) < (Real)0.0) {
		// Perfect reflection with 90% coefficient (10% energy dissipation)
		Real reflection_coeff = (Real)0.9;
		uzc = -(zwind + uzc) * reflection_coeff;
		zc = z_ground + roughness_length + (z_ground + roughness_length - zc) * reflection_coeff;

		// Reset vertical position to just above ground
		zc = fmax(zc, z_ground + roughness_length);
	}
#endif
	//zc=tz0;											// correct Z-directional position

	LP1[i].x=xc;															// update x-directional position
	LP1[i].y=yc;															// update y-directional position
	LP1[i].z=zc;															// update z-directional position
	LP1[i].sigma=sc;

	LP1[i].ux=uxc;														// update x-directional velocity
	LP1[i].uy=uyc;														// update y-directional velocity
	LP1[i].uz=uzc;														// update z-directional velocity

	// LP1[i].uxr=uxc+xwind;
	// LP1[i].uyr=uyc+ywind;
	// LP1[i].uzr=uzc+zwind;

	// LP1[i].xwind = xwind;
	// LP1[i].ywind = ywind;
	// LP1[i].zwind = zwind;
	LP1[i].k_turb = k_turb;
	LP1[i].e_turb = e_turb;
	LP1[i].yplus = yplus;
	LP1[i].Tscale = Tu;
	LP1[i].Sigvel = Sigvel;
}

__global__ void KERNEL_time_update_LDM(const Real tdt, Real tt, Real tend, int_t*g_str, int_t*g_end, L_part1*LP1, L_part2*LP2, part1*P1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==3) return;

	unsigned long long seed;
	curandState ss;
	double tmp;

	// ========================================================================
	// Time-dependent RNG initialization (critical for stochastic dispersion)
	// ========================================================================
	// Add time-dependent hash to ensure different random sequences at each timestep
	// Without this, particles would follow identical stochastic paths every timestep
	unsigned long long time_hash = (unsigned long long)(tt * 1e6);

	if (LP1[i].seed < 1) {
		// First time: initialize with particle ID and current time
		LP1[i].seed = i * 321215 + time_hash + 36893211;
		seed = LP1[i].seed;
	}
	else {
		// Subsequent times: evolve seed with LCG + time mixing
		seed = LP1[i].seed;
		seed = (seed * 1103515245 + time_hash + 12345) % 2147483647;
		LP1[i].seed = seed;
	}

	// Initialize curand state with time-dependent offset
	curand_init(seed, i, (unsigned long long)(tt * 100), &ss);

	Real tx0,ty0,tz0,ts0,xc,yc,zc,sc;									// position
	Real tux0,tuy0,tuz0,uxc,uyc,uzc;						// velocity
	Real tux,tuy,tuz;														// velocity
	Real dux_dt,duy_dt,duz_dt;									// accleration (time derivative of velocity)

	Real k_turb, e_turb, yplus;
	Real xwind, ywind, zwind;
	xwind=ywind=zwind=0.0;
	k_turb=e_turb=0.0;

	Real tmp_flt=0.0;

	tx0=LP2[i].x0;														// x-directional initial position
	ty0=LP2[i].y0;														// x-directional initial position
	tz0=LP2[i].z0;														// x-directional initial position
	ts0=LP2[i].sigma0;

	tux0=LP2[i].ux0;						// x-directional initial velocity
	tuy0=LP2[i].uy0;						// y-directional initial velocity
	tuz0=LP2[i].uz0;						// z-directional initial velocity

	yplus=0.0;
	sc=ts0;

	Real xi,yi,zi,hi,search_range;
	Real tmp_A;

	xi=LP1[i].x;
	yi=LP1[i].y;
	zi=LP1[i].z;

	hi=LP1[i].h;
	tmp_A=calc_tmpA(hi);

	int_t icell,jcell,kcell;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t ncell=LP1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,hj,uxj,uyj,uzj,tdist;
						Real rhoj,mj,presj,tempj,ptypej,rho_ref_j;
						Real k_turbj, e_turbj;
						int itype;

						ptypej=P1[j].p_type;
						mj=P1[j].m;

						if(P1[j].p_type==1){     // only for SPH fluid particles
							xj=P1[j].x;
							yj=P1[j].y;
							zj=P1[j].z;

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;

							k_turbj=P1[j].k_turb;
							e_turbj=P1[j].e_turb;

							rhoj=P1[j].rho;
							hj=P1[j].h;

							search_range=k_search_kappa*fmax(hi,hj);   // search range for multi-resolution (KDH)

							tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;

							if(tdist<search_range){
								
								Real twij=calc_kernel_wij(tmp_A,hi,tdist);
								
								xwind+=uxj*mj/rhoj*twij;
								ywind+=uyj*mj/rhoj*twij;
								zwind+=uzj*mj/rhoj*twij;

								k_turb+=k_turbj*mj/rhoj*twij;
								e_turb+=e_turbj*mj/rhoj*twij;

								tmp_flt += (mj/rhoj)*twij;
							}
						}
					}
				}
			}
		}
	}
	if(tmp_flt<1e-6){
		xwind=0;
		ywind=0;
		zwind=0;
		k_turb=1e-10;
		e_turb=1e-10;
	}
	else{
		// Normalize by total weight
		xwind/=tmp_flt;
		ywind/=tmp_flt;
		zwind/=tmp_flt;
		k_turb/=tmp_flt;
		e_turb/=tmp_flt;
	}

	// Lagrangian Timescales
	Real C0 = 1.76;  // Original value

	// ========================================================================
	// Turbulent Velocity Scale and Lagrangian Timescales
	// ========================================================================
	Real Sigvel = sqrt(2*k_turb/3);

	// Original isotropic timescales
	Real Tu = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);
	Real Tv = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);
	Real Tw = 2*Sigvel*Sigvel/C0/(e_turb+1e-10);

	Real KH, KV;

	Real t_dt=tdt;

	Real Ru = exp(-t_dt/Tu);
	Real Rv = exp(-t_dt/Tv);
	Real Rw = exp(-t_dt/Tw);

	Real mu = 0.0;
	Real stdv = 1;

	int_t p_type_i=LP1[i].p_type;

	Real drift_u = (Real)0.0;
	Real drift_v = (Real)0.0;
	Real drift_w = (Real)0.0;

	// Stochastic diffusion term (standard Langevin)
	dux_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Ru*Ru);			// x-directional acceleration
	duy_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Rv*Rv);			// y-directional acceleration
	duz_dt=Sigvel*GaussianRand(&ss, mu, stdv)*sqrt(1-Rw*Rw);			// z-directional acceleration

	// Combine deterministic drift + Thomson correction + stochastic diffusion
	uxc=tux0*Ru+dux_dt + drift_u*t_dt;										// correct x-directional velocity
	uyc=tuy0*Rv+duy_dt + drift_v*t_dt;										// correct y-directional velocity
	uzc=tuz0*Rw+duz_dt + drift_w*t_dt;

	// ========================================================================
	// Position Update with Ground Reflection Boundary Condition
	// ========================================================================
	xc=tx0+(xwind+uxc)*(t_dt);												// correct x-directional position
	yc=ty0+(ywind+uyc)*(t_dt);												// correct Y-directional position
	zc=tz0+(zwind+uzc)*(t_dt);

	// Ground reflection (z <= L1)
	if(zc <= (Real)L1){
		zc = (Real)(2.0*L1) - zc;
		uzc = -uzc;
	}
	// Ceiling reflection (z >= L3)
	if(zc >= (Real)L3){
		zc = (Real)(2.0*L3) - zc;
		uzc = -uzc;
	}

	LP1[i].x=xc;															// update x-directional position
	LP1[i].y=yc;															// update y-directional position
	LP1[i].z=zc;															// update z-directional position
	LP1[i].sigma=sc;

	LP1[i].ux=uxc;														// update x-directional velocity
	LP1[i].uy=uyc;														// update y-directional velocity
	LP1[i].uz=uzc;														// update z-directional velocity

	// LP1[i].uxr=uxc+xwind;
	// LP1[i].uyr=uyc+ywind;
	// LP1[i].uzr=uzc+zwind;

	// LP1[i].xwind = xwind;
	// LP1[i].ywind = ywind;
	// LP1[i].zwind = zwind;
	LP1[i].k_turb = k_turb;
	LP1[i].e_turb = e_turb;
	LP1[i].yplus = yplus;
	LP1[i].Tscale = Tu;
	LP1[i].Sigvel = Sigvel;
}

__global__ void KERNEL_time_update_LDM_direct(const Real tdt, Real tt,int_t*g_str, int_t*g_end,L_part1*LP1, L_part2*LP2, part1*P1)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i >= k_num_part_LDM) return;
    if(LP1[i].i_type == 3) return;

    // ========================================================================
    // 직접 확산 파라미터 설정 (k, ε 대신)
    // ========================================================================
    const Real target_D = Scalar_Diffusion_coeff;  // m²/s - 목표 확산계수
    
    // 방법 1: sigma_u를 먼저 정하고 tau 계산 (추천)
    const Real sigma_u = 0.05;  // m/s - 속도 RMS
    const Real tau = target_D / (sigma_u * sigma_u);  // s - Lagrangian timescale
    
    // 방법 2: tau를 먼저 정하고 sigma_u 계산 (대안)
    // const Real tau = 1.0;  // s
    // const Real sigma_u = sqrt(target_D / tau);  // m/s
    
    // ========================================================================
    // RNG initialization (기존과 동일)
    // ========================================================================
    unsigned long long seed;
    curandState ss;
    unsigned long long time_hash = (unsigned long long)(tt * 1e6);

    if (LP1[i].seed < 1) {
        LP1[i].seed = i * 321215 + time_hash + 36893211;
        seed = LP1[i].seed;
    }
    else {
        seed = LP1[i].seed;
        seed = (seed * 1103515245 + time_hash + 12345) % 2147483647;
        LP1[i].seed = seed;
    }
    curand_init(seed, i, (unsigned long long)(tt * 100), &ss);

    // ========================================================================
    // 변수 초기화
    // ========================================================================
    Real tx0, ty0, tz0, ts0, xc, yc, zc, sc;
    Real tux0, tuy0, tuz0, uxc, uyc, uzc;
    Real dux_dt, duy_dt, duz_dt;

    Real xwind, ywind, zwind;
    xwind = ywind = zwind = 0.0;

    Real tmp_flt = 0.0;

    tx0 = LP2[i].x0;
    ty0 = LP2[i].y0;
    tz0 = LP2[i].z0;
    ts0 = LP2[i].sigma0;

    tux0 = LP2[i].ux0;
    tuy0 = LP2[i].uy0;
    tuz0 = LP2[i].uz0;

    sc = ts0;

    Real xi, yi, zi, hi, search_range;
    Real tmp_A;

    xi = LP1[i].x;
    yi = LP1[i].y;
    zi = LP1[i].z;

    hi = LP1[i].h;
    tmp_A = calc_tmpA(hi);

    // ========================================================================
    // SPH 입자로부터 평균 유동 보간
    // ========================================================================
    int_t icell, jcell, kcell;

    // Cell indexing
    if((k_x_max == k_x_min)) {icell = 0;}
    else {icell = min(floor((xi - k_x_min) / k_dcell), k_NI - 1);}
    if((k_y_max == k_y_min)) {jcell = 0;}
    else {jcell = min(floor((yi - k_y_min) / k_dcell), k_NJ - 1);}
    if((k_z_max == k_z_min)) {kcell = 0;}
    else {kcell = min(floor((zi - k_z_min) / k_dcell), k_NK - 1);}
    
    if(icell < 0) icell = 0;
    if(jcell < 0) jcell = 0;
    if(kcell < 0) kcell = 0;

    int_t ncell = LP1[i].ncell;
    
    // Neighbor loop - SPH 속도장 보간
    for(int_t z = -ncell; z <= ncell; z++) {
        for(int_t y = -ncell; y <= ncell; y++) {
            for(int_t x = -ncell; x <= ncell; x++) {
                int_t k = idx_cell(icell + x, jcell + y, kcell + z);

                if(((icell + x) < 0) || ((icell + x) > (k_NI - 1)) ||
                   ((jcell + y) < 0) || ((jcell + y) > (k_NJ - 1)) ||
                   ((kcell + z) < 0) || ((kcell + z) > (k_NK - 1))) continue;

                if(g_str[k] != cu_memset) {
                    int_t fend = g_end[k];
                    for(int_t j = g_str[k]; j < fend; j++) {
                        
                        Real ptypej = P1[j].p_type;
                        
                        if(P1[j].p_type == 1) {  // SPH fluid particles only
                            Real xj = P1[j].x;
                            Real yj = P1[j].y;
                            Real zj = P1[j].z;

                            Real uxj = P1[j].ux;
                            Real uyj = P1[j].uy;
                            Real uzj = P1[j].uz;

                            Real rhoj = P1[j].rho;
                            Real mj = P1[j].m;
                            Real hj = P1[j].h;

                            search_range = k_search_kappa * fmax(hi, hj);

                            Real tdist = sqrt((xi - xj) * (xi - xj) + 
                                            (yi - yj) * (yi - yj) + 
                                            (zi - zj) * (zi - zj)) + 1e-20;

                            if(tdist < search_range) {
                                Real twij = calc_kernel_wij(tmp_A, hi, tdist);

                                xwind += uxj * mj / rhoj * twij;
                                ywind += uyj * mj / rhoj * twij;
                                zwind += uzj * mj / rhoj * twij;

                                tmp_flt += (mj / rhoj) * twij;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Normalize interpolated velocity
    if(tmp_flt < 1e-6) {
        xwind = 0.0;
        ywind = 0.0;
        zwind = 0.0;
    }
    else {
        xwind /= tmp_flt;
        ywind /= tmp_flt;
        zwind /= tmp_flt;
    }

    // ========================================================================
    // Langevin 방정식 (직접 파라미터 사용)
    // ========================================================================
    
    Real t_dt = tdt;
    
    // Ornstein-Uhlenbeck process parameters
    Real R = exp(-t_dt / tau);
    Real sqrt_factor = sigma_u * sqrt(1.0 - R * R);
    
    // Stochastic velocity update (난류 변동)
    dux_dt = sqrt_factor * GaussianRand(&ss, 0.0, 1.0);
    duy_dt = sqrt_factor * GaussianRand(&ss, 0.0, 1.0);
    duz_dt = sqrt_factor * GaussianRand(&ss, 0.0, 1.0);

    // Velocity update (감쇠 + 확산)
    uxc = tux0 * R + dux_dt;
    uyc = tuy0 * R + duy_dt;
    uzc = tuz0 * R + duz_dt;

    // ========================================================================
    // Position update (평균 유동 + 난류 변동)
    // ========================================================================
    xc = tx0 + (xwind + uxc) * t_dt;
    yc = ty0 + (ywind + uyc) * t_dt;
    zc = tz0 + (zwind + uzc) * t_dt;

    // ========================================================================
    // Boundary conditions
    // ========================================================================
    if(zc >= 10.0 || zc <= 0.0) uzc = -uzc;  // 반사 경계

    // ========================================================================
    // Update particle data
    // ========================================================================
    LP1[i].x = xc;
    LP1[i].y = yc;
    LP1[i].z = zc;

    LP1[i].sigma = sc;

    LP1[i].ux = uxc;
    LP1[i].uy = uyc;
    LP1[i].uz = uzc;

    // 디버깅/출력용 저장
    LP1[i].Tscale = tau;
    LP1[i].Sigvel = sigma_u;
}

// ========================================================================
// LDM verification kernel: exact random walk in position space
// - Advection: x += U_mean * dt  (uniform flow, no SPH interpolation)
// - Diffusion: x += sqrt(2*D*dt) * N(0,1)  (exact for any dt)
// - U_mean: mean flow velocity (x-direction only)
// - D_coeff: molecular/turbulent diffusion coefficient [m^2/s]
// ========================================================================
__global__ void KERNEL_time_update_LDM_test(const Real tdt, Real tt, L_part1*LP1, L_part2*LP2, Real U_mean, Real D_coeff)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i >= k_num_part_LDM) return;
    if(LP1[i].i_type == 3) return;

    // RNG initialization
    unsigned long long seed;
    curandState ss;
    unsigned long long time_hash = (unsigned long long)(tt * 1e6);

    if (LP1[i].seed < 1) {
        LP1[i].seed = i * 321215 + time_hash + 36893211;
        seed = LP1[i].seed;
    }
    else {
        seed = LP1[i].seed;
        seed = (seed * 1103515245 + time_hash + 12345) % 2147483647;
        LP1[i].seed = seed;
    }
    curand_init(seed, i, (unsigned long long)(tt * 100), &ss);

    Real tx0 = LP2[i].x0;
    Real ty0 = LP2[i].y0;
    Real tz0 = LP2[i].z0;

    Real t_dt = tdt;

    // Random walk diffusion: sqrt(2*D*dt) * N(0,1)
    Real diff_scale = sqrt(2.0 * D_coeff * t_dt);
    Real dWx = diff_scale * GaussianRand(&ss, 0.0, 1.0);
    Real dWy = diff_scale * GaussianRand(&ss, 0.0, 1.0);
    Real dWz = diff_scale * GaussianRand(&ss, 0.0, 1.0);

    // Position update: advection + diffusion
    Real xc = tx0 + U_mean * t_dt + dWx;
    Real yc = ty0 + dWy;
    Real zc = tz0 + dWz;

    LP1[i].x = xc;
    LP1[i].y = yc;
    LP1[i].z = zc;

    LP1[i].ux = U_mean;
    LP1[i].uy = 0.0;
    LP1[i].uz = 0.0;
}

__global__ void KERNEL_clc_concentration(int_t*g_str, int_t*g_end, L_part1*LP1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part_LDM) return;
	if(LP1[i].i_type==3) return;

	int_t ptypei;
	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real pi;														// pressure
	//Real flt_si;
	Real tmpp;
	Real search_range,tmp_h,tmp_A;

	tmp_h=LP1[i].h;
	tmp_A=calc_tmpA(tmp_h);
	search_range=k_search_kappa*tmp_h;								// search range

	ptypei=LP1[i].p_type;

	xi=LP1[i].x;
	yi=LP1[i].y;
	zi=LP1[i].z;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpp=0.0;
	Real tmppp=0.0;
	int_t ncell=LP1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist;
						xj=LP1[j].x;
						yj=LP1[j].y;
						zj=LP1[j].z;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
						if(tdist<search_range){
							if(LP1[i].p_type==LP1[j].p_type){   // same species concentration
								Real twij,pj,mj,rhoj,ptype_j;
								twij=calc_kernel_wij(tmp_A,tmp_h,tdist);
								mj=LP1[j].m;
								rhoj=LP1[j].rho;
								tmpp+=mj/rhoj*twij;
								tmppp+=1;
							}
						}
					}
				}
			}
		}
	}
	LP1[i].concn=tmpp;
	LP1[i].concentration=tmppp;	
}

__global__ void KERNEL_compute_vorticity3D(int_t*g_str,int_t*g_end,part1*P1,part3*P3){
	uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= k_num_part2) return;
  if (P1[i].i_type > i_type_crt) return;
  if (P1[i].p_type >= 1000) return;  // IB 제외
  if (P1[i].p_type <= 0) return;     // 유체만

  // ---------------------
  // 로컬 변수
  // ---------------------
  int_t icell, jcell, kcell;
  Real xi = P1[i].x, yi = P1[i].y, zi = P1[i].z;
  Real uxi = P1[i].ux, uyi = P1[i].uy, uzi = P1[i].uz;
  Real rhoi = P1[i].rho;
  Real hi   = P1[i].h;

  // Changed: 다해상도/커널 준비
  Real tmp_A = calc_tmpA(hi);
  Real search_range = k_search_kappa * hi;

  // grid cell index
  if ((k_x_max==k_x_min)) icell=0; else icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);
  if ((k_y_max==k_y_min)) jcell=0; else jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);
  if ((k_z_max==k_z_min)) kcell=0; else kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);
  if (icell<0) icell=0; if (jcell<0) jcell=0; if (kcell<0) kcell=0;

  // 누적: omega = (omx,omy,omz)
  Real omx=0.0, omy=0.0, omz=0.0;

  int_t ncell = P1[i].ncell;

  for (int_t z=-ncell; z<=ncell; ++z){
    for (int_t y=-ncell; y<=ncell; ++y){
      for (int_t x=-ncell; x<=ncell; ++x){
        int_t k = idx_cell(icell+x, jcell+y, kcell+z);
        if (((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;

        if (g_str[k]!=cu_memset){
          int_t fend = g_end[k];
          for (int_t j=g_str[k]; j<fend; ++j){
            if (P1[j].p_type>=1000) continue;  // IB 제외

            // self-skip는 굳이 안해도 됨(gradW가 0)
            Real xj = P1[j].x, yj = P1[j].y, zj = P1[j].z;
            Real hj = P1[j].h;

            // Changed: 다해상도 탐색 반경
            search_range = k_search_kappa * fmax(hi, hj);

            Real rx = xi - xj;
            Real ry = yi - yj;
            Real rz = zi - zj;
            Real tdist = sqrt(rx*rx + ry*ry + rz*rz) + 1e-20;
            if (tdist >= search_range) continue;

            // 속도/밀도/질량
            Real uxj = P1[j].ux, uyj = P1[j].uy, uzj = P1[j].uz;
            Real mj  = P1[j].m;
            Real rhoj= P1[j].rho;

            // ---------------------
            // grad W_ij
            // ---------------------
            // 기본: i의 hi로 커널
            Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
            Real tdwx  = tdwij * rx / tdist;
            Real tdwy  = tdwij * ry / tdist;
            Real tdwz  = tdwij * rz / tdist;

            // Changed: gradient correction (기존 패턴과 동일)
            if (k_kgc_solve>0){
              Real twij_dummy = 0.0; // 함수 시그니처 유지용
              apply_gradient_correction_3D(P3[i].Cm, twij_dummy,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
            }

            // ---------------------
            // (u_j - u_i) × ∇W_ij
            // ---------------------
            Real dux = uxj - uxi;
            Real duy = uyj - uyi;
            Real duz = uzj - uzi;

            // cross = du × gradW
            Real cx =   duy*tdwz - duz*tdwy;
            Real cy =   duz*tdwx - dux*tdwz;
            Real cz =   dux*tdwy - duy*tdwx;

            // 누적
            Real w = mj / (rhoj + 1e-20);
            omx += w * cx;
            omy += w * cy;
            omz += w * cz;
          }
        }
      }
    }
  }

  // 출력
  P3[i].omx = omx;         // Changed
  P3[i].omy = omy;         // Changed
  P3[i].omz = omz;         // Changed
}