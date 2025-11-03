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

	int_t j;
	double ttend;
	ttend = tend;
	if (ttime<=ttend){
		for(j=k_num_part_LDM*(ttime)/(ttend); j<k_num_part_LDM*(ttime+freq*tdt)/(ttend); j++){
			LP1[j].i_type=1;
		}
	}

}

__global__ void KERNEL_clc_predictor_Dispersion(Real tdt,Real ttime,L_part1*LP1,L_part2*LP2,L_part3*LP3)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;

	if(i>=k_num_part_LDM) return;
	//if(LP1[i].p_type==0) return;
	if(LP1[i].i_type==3) return;

	int_t p_typei;
	double tmp;

	int_t buffer_type=LP1[i].buffer_type;

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

	int_t buffer_type=LP1[i].buffer_type;
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

	LP1[i].uxr=uxc+xwind;
	LP1[i].uyr=uyc+ywind;
	LP1[i].uzr=uzc+zwind;

	LP1[i].xwind = xwind;
	LP1[i].ywind = ywind;
	LP1[i].zwind = zwind;
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
								//if(i==135235) P1[j].PPE2=1.0;
								
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

	int_t buffer_type=LP1[i].buffer_type;
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

	LP1[i].x=xc;															// update x-directional position
	LP1[i].y=yc;															// update y-directional position
	LP1[i].z=zc;															// update z-directional position									
	if(zc<=0) uzc=-uzc;							// 바닥 반사 조건(kdh)								
	
	LP1[i].sigma=sc;

	LP1[i].ux=uxc;														// update x-directional velocity
	LP1[i].uy=uyc;														// update y-directional velocity
	LP1[i].uz=uzc;														// update z-directional velocity

	LP1[i].uxr=uxc+xwind;
	LP1[i].uyr=uyc+ywind;
	LP1[i].uzr=uzc+zwind;

	LP1[i].xwind = xwind;
	LP1[i].ywind = ywind;
	LP1[i].zwind = zwind;
	LP1[i].k_turb = k_turb;
	LP1[i].e_turb = e_turb;
	LP1[i].yplus = yplus;
	LP1[i].Tscale = Tu;
	LP1[i].Sigvel = Sigvel;
}