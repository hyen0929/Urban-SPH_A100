////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_clc_projection(int_t tcount, Real tdt,Real ttime,part1*P1,part2*P2,part3*P3)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>=3) return;
	if(P1[i].p_type>=1000) return;		// Immersed Boundary Method

	int_t p_typei;
	Real tx0,ty0,tz0;
	Real tux0,tuy0,tuz0,tuxp,tuyp,tuzp;
	Real tdux_dt0,tduy_dt0,tduz_dt0;
	Real t_dt;

	int_t buffer_type=P1[i].buffer_type;
	int_t MOST=P1[i].MOST_buffer;
	t_dt=tdt;
	p_typei=P1[i].p_type;

	if(p_typei==MOVING){
		P2[i].x0=P1[i].x;
		P2[i].y0=P1[i].y;
		P2[i].z0=P1[i].z;

		// P1[i].ux=-0.06*PI*sin(ttime*PI);
		P1[i].ux=Boundary_Velocity;
		P1[i].uy=0;
		P1[i].uz=0;
		P2[i].ux0=Boundary_Velocity;
		P2[i].uy0=0;
		P2[i].uz0=0;
	}else{
		tx0=P1[i].x;															// initial x-directional position
		ty0=P1[i].y;															// initial y-directional position
		tz0=P1[i].z;															// initial z-directional position
		tux0=P1[i].ux;															// initial x-directional position //YHS
		tuy0=P1[i].uy;															// initial y-directional position //YHS
		tuz0=P1[i].uz;															// initial z-directional position //YHS
		if(p_typei>0){
			tux0=P1[i].ux;													// initial x-directional velocity
			tuy0=P1[i].uy;													// initial y-directional velocity
			tuz0=P1[i].uz;													// initial z-directional velocity

			// tdux_dt0=P3[i].ftotalx*(buffer_type==0)+P1[i].fbx/P1[i].rho*(buffer_type==0);									// initial x-directional acceleration
			// tduy_dt0=P3[i].ftotaly*(buffer_type==0)+P1[i].fby/P1[i].rho*(buffer_type==0);									// initial y-directional acceleration
			// tduz_dt0=P3[i].ftotalz*(buffer_type==0)+P1[i].fbz/P1[i].rho*(buffer_type==0);									// initial z-directional acceleration

			// Geostrophic pressure gradient: balances MOST surface friction for steady ABL
#if CASE_ID >= 500 && CASE_ID < 600
			// Case CS: body force OFF (surface BC만 활성 — momentum source 없음)
			Real f_geo = 0.0;
#else
			Real f_geo = MOST_boundary * (us_ini * us_ini / ABL_h);
#endif
			tdux_dt0=P3[i].ftotalx*(buffer_type==0||MOST==1)+P1[i].fmx+f_geo*(buffer_type==0||MOST==1);									// initial x-directional acceleration

			//tdux_dt0=P3[i].ftotalx*(buffer_type==0||MOST==1)+P1[i].fmx;									// initial x-directional acceleration
			tduy_dt0=P3[i].ftotaly*(buffer_type==0||MOST==1)+P1[i].fmy;									// initial y-directional acceleration
			tduz_dt0=P3[i].ftotalz*(buffer_type==0||MOST==1);									// initial z-directional acceleration

			// Outlet buffer에서 속도 계산
			// tdux_dt0=P3[i].ftotalx*(buffer_type==0||buffer_type==2);									// initial x-directional acceleration
			// tduy_dt0=P3[i].ftotaly*(buffer_type==0||buffer_type==2);									// initial y-directional acceleration
			// tduz_dt0=P3[i].ftotalz*(buffer_type==0||buffer_type==2);									// initial z-directional acceleration

			tuxp=tux0+tdux_dt0*(t_dt);						// Predict x-directional velocity (dux_dt0 : acceleration of before time step)
			tuyp=tuy0+tduy_dt0*(t_dt);						// Predict y-directional velocity (duy_dt0 : acceleration of before time step)
			tuzp=tuz0+tduz_dt0*(t_dt);						// Predict z-directional velocity (duz_dt0 : acceleration of before time step)
																						// No need of Position update (Eulerian)

			//printf("ux0=%f, uy0=%f, uz0=%f, tdux=%f, tduy=%f, tduy=%f\n",tux0, tuy0, tuz0, tdux_dt0, tduy_dt0, tduz_dt0);

		}else{
			tuxp=tux0;tuyp=tuy0;tuzp=tuz0;
			tuxp=P1[i].ux;
			tuyp=P1[i].uy;
			tuzp=P1[i].uz;
		}

		P1[i].x=tx0;															// Update particle data by predicted x-directional position
		P1[i].y=ty0;															// Update particle data by predicted y-directional position
		P1[i].z=tz0;															// Update particle data by predicted z-directional position
		P1[i].ux=tuxp;														// Update particle data by predicted x-directional velocity
		P1[i].uy=tuyp;														// Update particle data by predicted y-directional velocity
		P1[i].uz=tuzp;														// Update particle data by predicted z-directional velocity

		P2[i].x0=tx0;															// update x-directional position
		P2[i].y0=ty0;															// update y-directional position
		P2[i].z0=tz0;															// update z-directional position
		P2[i].ux0=tux0;														// update x-directional velocity
		P2[i].uy0=tuy0;														// update y-directional velocity
		P2[i].uz0=tuz0;														// update z-directional velocity
	}

	if(ttime==0){
		P2[i].rho_ref=P1[i].rho;
	}
	if(k_con_solve==1){
		if(enthalpy_eqn){
			Real tenthalpyp=P1[i].enthalpy;
				P2[i].enthalpy0=tenthalpyp;
				P1[i].enthalpy=tenthalpyp;
				P1[i].temp=htoT(tenthalpyp,p_typei);
			}
			else{
				Real ttemp=P1[i].temp;
				P2[i].temp0=ttemp;
				P1[i].temp=ttemp;
			}
	}

	// if(k_concn_solve==1){
	// 	Real tconcn=P1[i].concn;
	// 	P2[i].concn0=tconcn;
	// 	P1[i].concn=tconcn+P3[i].dconcn*(t_dt);
	// }

	// if(isnan(P1[i].ux)) printf("ftotal=%f\n",P3[i].ftotalx);

	return;
}
////////////////////////////////////////////////////////////////////////
// corrector step for Predictor-Corrector time integration
__global__ void KERNEL_time_update_projection(const Real tdt,part1*P1,part1*TP1,part2*P2,part2*TP2,part3*P3,int_t tcount)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>=3) return;
	if(P1[i].p_type>=1000) return;		// Immersed Boundary Method
	if(P1[i].buffer_type==99) return;	// Interior fluid (inside buildings) — keep zeroed state from previous step
	// if(P1[i].buffer_type==1) return;

	Real tx0,ty0,tz0,xc,yc,zc;									// position
	Real tux0,tuy0,tuz0,uxc,uyc,uzc;						// velocity
	Real tux,tuy,tuz;														// velocity
	Real dux_dt,duy_dt,duz_dt;									// accleration (time derivative of velocity)
	Real t_dt=tdt;

	int_t buffer_type=P1[i].buffer_type;
	int_t p_type_i=P1[i].p_type;
	int_t MOST_buffer=P1[i].MOST_buffer;

	if(p_type_i==MOVING){
		tx0=P1[i].x;														// x-directional initial position
		ty0=P1[i].y;														// x-directional initial position
		tz0=P1[i].z;														// x-directional initial position
		tux=P1[i].ux;
		tuy=P1[i].uy;
		tuz=P1[i].uz;

		// TP1[i].x=xc;															// update x-directional position
		// TP1[i].y=yc;															// update y-directional position
		// TP1[i].z=zc;															// update z-directional position
		TP1[i].ux=tux;														// update x-directional velocity
		TP1[i].uy=tuy;														// update y-directional velocity
		TP1[i].uz=tuz;														// update z-directional velocity

	}else{
		tx0=P1[i].x;														// x-directional initial position
		ty0=P1[i].y;														// x-directional initial position
		tz0=P1[i].z;														// x-directional initial position
		if(p_type_i>0){
			tux0=P1[i].ux;						// x-directional initial velocity
			tuy0=P1[i].uy;						// y-directional initial velocity
			tuz0=P1[i].uz;						// z-directional initial velocity

			if(buffer_type==1||buffer_type==-1){        // Inlet buffer의 경우 advection force 배제
				tux0=P2[i].ux0;						// x-directional initial velocity
				tuy0=P2[i].uy0;						// y-directional initial velocity
				tuz0=P2[i].uz0;						// z-directional initial velocity
			}

			// // buffer를 제외한 속도 계산
			// dux_dt=P3[i].fpx*(buffer_type==0);			// x-directional acceleration
			// duy_dt=P3[i].fpy*(buffer_type==0);			// y-directional acceleration
			// duz_dt=P3[i].fpz*(buffer_type==0);			// z-directional acceleration

			// dux_dt=P3[i].fpx*(buffer_type==0)+P1[i].fbx/P1[i].rho+P1[i].fmx;			// x-directional acceleration
			// duy_dt=P3[i].fpy*(buffer_type==0)+P1[i].fby/P1[i].rho+P1[i].fmy;			// y-directional acceleration
			// duz_dt=P3[i].fpz*(buffer_type==0)+P1[i].fbz/P1[i].rho;			// z-directional acceleration

			dux_dt=P3[i].fpx*(buffer_type==0);			// x-directional acceleration (pressure only)
			duy_dt=P3[i].fpy*(buffer_type==0);			// y-directional acceleration (pressure only)
			duz_dt=P3[i].fpz*(buffer_type==0);			// z-directional acceleration
		}
		else{
			tux0=P2[i].ux0;
			tuy0=P2[i].uy0;
			tuz0=P2[i].uz0;
			dux_dt=duy_dt=duz_dt=0.0;
		}
		uxc=tux0+dux_dt*(t_dt);										// correct x-directional velocity
		uyc=tuy0+duy_dt*(t_dt);										// correct y-directional velocity
		uzc=tuz0+duz_dt*(t_dt);										// correct z-directional velocity

		if((uxc*uxc+uyc*uyc+uzc*uzc)>=k_u_limit*k_u_limit){
			uxc=tux0;
			uyc=tuy0;
			uzc=tuz0;
		}

		// Euler method
		// xc=tx0+uxc*(t_dt)*(p_type_i>0)*(P1[i].eli);												// correct x-directional position
		// yc=ty0+uyc*(t_dt)*(p_type_i>0)*(P1[i].eli);												// correct Y-directional position
		// zc=tz0+uzc*(t_dt)*(p_type_i>0)*(P1[i].eli);

		xc=tx0+(P2[i].ux0+uxc)/2.0*(t_dt)*(p_type_i>0)*(P1[i].eli);												// correct x-directional position
		yc=ty0+(P2[i].uy0+uyc)/2.0*(t_dt)*(p_type_i>0)*(P1[i].eli);												// correct Y-directional position
		zc=tz0+(P2[i].uz0+uzc)/2.0*(t_dt)*(p_type_i>0)*(P1[i].eli);
		
		Real umag=sqrt(uxc*uxc+uyc*uyc+uzc*uzc);
		if(umag>k_u_limit){
			TP1[i].ux=uxc*k_u_limit/umag;														// update x-directional velocity
		   	TP1[i].uy=uyc*k_u_limit/umag;														// update y-directional velocity
			TP1[i].uz=uzc*k_u_limit/umag;														// update z-directional velocity
		}
		else{
			TP1[i].ux=uxc;														// update x-directional velocity
		   	TP1[i].uy=uyc;														// update y-directional velocity
			TP1[i].uz=uzc;														// update z-directional velocity
		}
   }

	//update_properties_enthalpy-------------------------------
	if((k_con_solve==1)){
		if(enthalpy_eqn){
			//TP1[i].enthalpy=P2[i].enthalpy0+P3[i].denthalpy*t_dt*(P1[i].p_type!=-1);
			TP1[i].enthalpy=P2[i].enthalpy0+P3[i].denthalpy*t_dt;
			TP1[i].temp=P1[i].temp;
		}else{
			Real temp_new=P2[i].temp0+P3[i].dtemp*t_dt*(buffer_type==0);
			TP1[i].temp=temp_new;
		}
	}
	else{
		TP1[i].enthalpy=P1[i].enthalpy;
		TP1[i].temp=P1[i].temp;
	}

	//update_properties_concn----------------------------------
	// if(k_concn_solve==1) {
	// 	TP1[i].concn=P2[i].concn0+P3[i].dconcn*t_dt;
	// 	if(TP1[i].concn<0.0) TP1[i].concn=0.0;
	// }
	// else TP1[i].concn=P1[i].concn;

	TP1[i].pres=P1[i].pres;
	TP1[i].pres0=P1[i].pres;   // incremental PPE anchor: store previous-step pres
	TP1[i].flt_s=P1[i].flt_s;
	TP1[i].MOST_buffer=P1[i].MOST_buffer;
	TP1[i].m=P1[i].m;
	TP1[i].ncell=P1[i].ncell;
	TP1[i].h=P1[i].h;
	TP1[i].k_turb=P1[i].k_turb;
	TP1[i].e_turb=P1[i].e_turb;

	TP1[i].eli=P1[i].eli;

	TP1[i].vol=P1[i].vol;

	TP2[i].rho_ref=P2[i].rho_ref;
	TP1[i].vis_t=P1[i].vis_t;

	// TP1[i].ux_f=P1[i].ux_f;
	// TP1[i].uy_f=P1[i].uy_f;
	// TP1[i].uz_f=P1[i].uz_f;
	//
	// TP1[i].PPE1=P1[i].PPE1;
	// TP1[i].PPE2=P1[i].PPE2;

	TP1[i].fp_x=P1[i].fp_x;
	TP1[i].fp_y=P1[i].fp_y;
	TP1[i].fp_z=P1[i].fp_z;
	
	TP1[i].fe_x=P1[i].fe_x;
	TP1[i].fe_y=P1[i].fe_y;
	TP1[i].fe_z=P1[i].fe_z;
	
	TP1[i].fv_x=P1[i].fv_x;
	TP1[i].fv_y=P1[i].fv_y;
	TP1[i].fv_z=P1[i].fv_z;

	TP1[i].fbx=P1[i].fbx;
	TP1[i].fby=P1[i].fby;
	TP1[i].fbz=P1[i].fbz;

	TP1[i].fmx = P1[i].fmx;
	TP1[i].fmy = P1[i].fmy;
	TP1[i].hmz = P1[i].hmz;

	// Store results in particle

	TP1[i].concentration=P1[i].concentration;

	// if(TP1[i].source==1) TP1[i].concn+=0.5*t_dt;

	if(TP1[i].p_type==0){
		if(tdt*tcount>=T_start_t && tdt*tcount<T_con_t) TP1[i].temp=T_ini+tdt*tcount*(T_con-T_ini)/(T_con_t-T_start_t);
    else if(tdt*tcount<T_start_t) TP1[i].temp=T_ini;
    else TP1[i].temp=T_con;
	}
	//if(TP1[i].buffer_type==1) TP1[i].ux = us_update/k_vonKarman*log(TP1[i].z/MOST_z0);

	// if(tcount==0){
	// 	if(P1[i].p_type==1){
	// 		//TP1[i].ux=1.0;  // velocity perturbation on y direction
	// 		Real z=max(TP1[i].z,0.0);
	// 		Real z_ref = 300.0;
	//
	// 		TP1[i].ux=1.0*pow(z/z_ref,0.33);
	// 		TP1[i].ux=1.0;
	// 	}
	// }
}
//
// __global__ void KERNEL_time_update_buffer(const Real tdt,part1*P1,part1*TP1,part2*P2,part2*TP2,part3*P3,int_t tcount)
// {
// 	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part2) return;
// 	if(P1[i].buffer_type!=1) return;
//
// 	TP1[i].pres=P1[i].pres;
// }
////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_xsph3D(int_t*g_str,int_t*g_end,part1*P1)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>3) return;
	if(P1[i].p_type!=1) return;
	if(P1[i].buffer_type!=0) return;

	int_t ptypei;
	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real uxi,uyi,uzi;
	Real presi;
	Real tmpx,tmpy,tmpz;
	Real search_range,tmp_h,tmp_A;
	Real flt_s;

	tmp_h=XSPH_length*P1[i].h;
	tmp_A=calc_tmpA(tmp_h);
	search_range=k_search_kappa*tmp_h;								// search range

	ptypei=P1[i].p_type;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;
	presi=P1[i].pres;

	flt_s=0.0;

	Real c_xsph=0.003;      // XSPH coefficient
	Real c_RC=0.1;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpx=tmpy=tmpz=0.0;

	int_t ncell=P1[i].ncell;

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
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
						if(tdist<search_range){
							Real twij,tdwij,tdwx,tdwy,tdwz,uxj,uyj,uzj,mj,rhoj,presj;
							twij=calc_kernel_wij(tmp_A,tmp_h,tdist);
							tdwij=calc_kernel_dwij(tmp_A,tmp_h,tdist);

							tdwx=tdwij*(xi-xj)/tdist;
							tdwy=tdwij*(yi-yj)/tdist;
							tdwz=tdwij*(zi-zj)/tdist;

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;
							mj=P1[j].m;
							rhoj=P1[j].rho;
							presj=P1[j].pres;

							flt_s+=mj/rhoj*twij*(P1[j].p_type==1);

							// XSPH (KDH)
							tmpx+=c_xsph*mj/rhoj*(-uxi+uxj)*twij*(P1[j].p_type==1);
							tmpy+=c_xsph*mj/rhoj*(-uyi+uyj)*twij*(P1[j].p_type==1);
							tmpz+=c_xsph*mj/rhoj*(-uzi+uzj)*twij*(P1[j].p_type==1);

							// Rhie & Chow (KDH)
							// Real RC_ij=c_RC*(presi-presj)*mj/rhoj;
							// //if (isnan(RC_ij)==1) printf("RC=%f\n",RC_ij);
							// tmpx+=RC_ij*tdwx;
							// tmpy+=RC_ij*tdwy;
							// tmpz+=RC_ij*tdwz;

						}
					}
				}
			}
		}
	}
	if (flt_s>1e-6){
		// TP1[i].ux=uxi+tmpx/flt_s;				// update x-directional position
		// TP1[i].uy=uyi+tmpy/flt_s;				// update y-directional position
		// TP1[i].uz=uzi+tmpz/flt_s;				// update z-directional position

		P1[i].ux=uxi+tmpx/flt_s;				// update x-directional position
		P1[i].uy=uyi+tmpy/flt_s;				// update y-directional position
		P1[i].uz=uzi+tmpz/flt_s; 			// update z-directional position
		//printf("tmpx=%f\n",tmpx/flt_s);
	}
}

// =======================================================================
// Time-averaging accumulation (Eulerian SPH)
//   Adds current (ux,uy,uz,k_turb,e_turb) to running sums on each fluid
//   particle. Called from ISPH_Calc.cuh after the corrector step, gated by
//   simulation time >= AVG_START_TIME. The global counter avg_count is
//   incremented once per call by thread 0. Divide by avg_count at output.
// =======================================================================
__global__ void KERNEL_accumulate_avg(part1* P1)
{
	uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
	// Counter: single thread increments, independent of per-particle early-return
	if(i == 0) atomicAdd(&avg_count, 1u);

	if(i >= k_num_part2) return;
	if(P1[i].i_type > i_type_crt) return;   // skip dummy / outer-buffer
	if(P1[i].p_type != FLUID) return;       // fluid only

	P1[i].ux_avg     += P1[i].ux;
	P1[i].uy_avg     += P1[i].uy;
	P1[i].uz_avg     += P1[i].uz;
	P1[i].k_turb_avg += P1[i].k_turb;
	P1[i].e_turb_avg += P1[i].e_turb;
}

// =======================================================================
// Mean-field I/O for two-run transient-vs-mean LDM comparison.
//   Run A (LDM_MEAN_MODE=0): at time = AVG_END_TIME, save_mean_field dumps a
//     one-shot binary snapshot of <ux,uy,uz,k_turb,e_turb> averaged over
//     [AVG_START_TIME, AVG_END_TIME], keyed by particle index.
//   Run B (LDM_MEAN_MODE=1): at startup, load_mean_field reads the snapshot
//     and overwrites dev_P1.ux/uy/uz/k_turb/e_turb so LDM advects on the
//     frozen mean. SPH solver calls are skipped in ISPH_Calc.cuh.
// File layout: [mean_field_header][mean_field_record × n_part2]
// =======================================================================
typedef struct mean_field_record {
	float ux, uy, uz;
	float k_turb, e_turb;
} mean_field_record;

typedef struct mean_field_header {
	char         magic[16];           // "SOPHIA_MEAN_v1"
	int          n_part2;
	unsigned int avg_count_at_dump;
	float        us_update;
	float        OB_length_global;
	int          pad[4];              // reserved
} mean_field_header;

__global__ void KERNEL_apply_mean_field(part1* P1, const mean_field_record* mean)
{
	uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i >= k_num_part2) return;
	P1[i].ux     = mean[i].ux;
	P1[i].uy     = mean[i].uy;
	P1[i].uz     = mean[i].uz;
	P1[i].k_turb = mean[i].k_turb;
	P1[i].e_turb = mean[i].e_turb;
}

// file_P1 is a host buffer already copied from dev_P1 by the caller.
static void save_mean_field(const part1* file_P1, int n_part2_in, const char* path)
{
	unsigned int h_avg_count = 0;
	cudaMemcpyFromSymbol(&h_avg_count, avg_count, sizeof(unsigned int));
	if(h_avg_count == 0u){
		printf("[mean] save_mean_field skipped: avg_count == 0\n");
		return;
	}
	float h_us = 0.0f, h_OB = 0.0f;
	cudaMemcpyFromSymbol(&h_us, us_update,        sizeof(float));
	cudaMemcpyFromSymbol(&h_OB, OB_length_global, sizeof(float));

	FILE* fp = fopen(path, "wb");
	if(!fp){
		printf("[mean] save_mean_field: cannot open %s for write\n", path);
		return;
	}

	mean_field_header hdr;
	memset(&hdr, 0, sizeof(hdr));
	memcpy(hdr.magic, "SOPHIA_MEAN_v1", 14);
	hdr.n_part2           = n_part2_in;
	hdr.avg_count_at_dump = h_avg_count;
	hdr.us_update         = h_us;
	hdr.OB_length_global  = h_OB;
	fwrite(&hdr, sizeof(hdr), 1, fp);

	Real inv_n = (Real)1.0 / (Real)h_avg_count;
	for(int i = 0; i < n_part2_in; i++){
		mean_field_record rec;
		rec.ux     = (float)(file_P1[i].ux_avg     * inv_n);
		rec.uy     = (float)(file_P1[i].uy_avg     * inv_n);
		rec.uz     = (float)(file_P1[i].uz_avg     * inv_n);
		rec.k_turb = (float)(file_P1[i].k_turb_avg * inv_n);
		rec.e_turb = (float)(file_P1[i].e_turb_avg * inv_n);
		fwrite(&rec, sizeof(rec), 1, fp);
	}
	fclose(fp);
	printf("[mean] saved %d records, avg_count=%u, us=%.4f, OB=%.4e -> %s\n",
	       n_part2_in, h_avg_count, h_us, h_OB, path);
}

static int load_mean_field(part1* dev_P1, int n_part2_in, const char* path)
{
	FILE* fp = fopen(path, "rb");
	if(!fp){
		printf("[mean] load_mean_field: cannot open %s\n", path);
		return 0;
	}
	mean_field_header hdr;
	if(fread(&hdr, sizeof(hdr), 1, fp) != 1){
		printf("[mean] load_mean_field: header read failed\n");
		fclose(fp);
		return 0;
	}
	if(memcmp(hdr.magic, "SOPHIA_MEAN_v1", 14) != 0){
		printf("[mean] load_mean_field: bad magic in %s\n", path);
		fclose(fp);
		return 0;
	}
	if(hdr.n_part2 != n_part2_in){
		printf("[mean] load_mean_field: n_part2 mismatch (file=%d, runtime=%d)\n",
		       hdr.n_part2, n_part2_in);
		fclose(fp);
		return 0;
	}

	mean_field_record* host_rec =
		(mean_field_record*)malloc(sizeof(mean_field_record) * (size_t)n_part2_in);
	size_t nread = fread(host_rec, sizeof(mean_field_record), (size_t)n_part2_in, fp);
	fclose(fp);
	if((int)nread != n_part2_in){
		printf("[mean] load_mean_field: record read short (%zu / %d)\n",
		       nread, n_part2_in);
		free(host_rec);
		return 0;
	}

	mean_field_record* dev_rec = NULL;
	cudaMalloc((void**)&dev_rec, sizeof(mean_field_record) * (size_t)n_part2_in);
	cudaMemcpy(dev_rec, host_rec,
	           sizeof(mean_field_record) * (size_t)n_part2_in,
	           cudaMemcpyHostToDevice);

	dim3 t(128), b((n_part2_in - 1) / 128 + 1);
	KERNEL_apply_mean_field<<<b, t>>>(dev_P1, dev_rec);
	cudaDeviceSynchronize();

	cudaFree(dev_rec);
	free(host_rec);

	// Restore reference scalars updated by MOST at freeze time
	cudaMemcpyToSymbol(us_update,        &hdr.us_update,        sizeof(float));
	cudaMemcpyToSymbol(OB_length_global, &hdr.OB_length_global, sizeof(float));

	printf("[mean] loaded %d records, avg_count_at_dump=%u, us=%.4f, OB=%.4e <- %s\n",
	       n_part2_in, hdr.avg_count_at_dump, hdr.us_update, hdr.OB_length_global, path);
	return 1;
}
