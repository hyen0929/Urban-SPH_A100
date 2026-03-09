////////////////////////////////////////////////////////////////////////
__global__ void IBM_predictor(Real t_dt,Real ttime,part1*P1,part2*P2)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type<1000) return;

	int_t p_typei;
	Real tx0,ty0,tz0,txp,typ,tzp;
	Real tux0,tuy0,tuz0,tuxp,tuyp,tuzp;
	Real tdux_dt0,tduy_dt0,tduz_dt0;

	Real freq = 0.8*0.196*0.1/0.025;
	Real A = 0.55*0.02;
	P1[i].ux=0.0;
	// P1[i].uy=-0.02*sin(2.0*PI*freq*ttime)*(P1[i].p_type==1000);
	// P1[i].uy=-2.0*PI*freq*A*sin(2.0*PI*freq*ttime);
	P1[i].uy=0.0;
	P1[i].uz=0.0;


	tx0 = P1[i].x;
	ty0 = P1[i].y;
	tz0 = P1[i].z;

	tux0= P1[i].ux;
	tuy0= P1[i].uy;
	tuz0= P1[i].uz;

	// tdux_dt0=P1[i].ax;
	// tduy_dt0=P1[i].ay;
	// tduz_dt0=P1[i].az;

	tdux_dt0=0.0;
	tduy_dt0=0.0;
	tduz_dt0=0.0;    // 김도현 수정

	//update position

	txp = tx0 + (tux0*0.5)*t_dt;
	typ = ty0 + (tuy0*0.5)*t_dt;
	tzp = tz0 + (tuz0*0.5)*t_dt;

	// Static point

	// if(P1[i].ccc==3){
	// 	tux0=0.0;
	// 	tuy0=0.0;
	// 	tuz0=0.0;
	// 	tdux_dt0=tduy_dt0=tduz_dt0=0.0;
	// }

	//update velocity

	tuxp = tux0 + (tdux_dt0*0.5)*t_dt;
	tuyp = tuy0 + (tduy_dt0*0.5)*t_dt;
	tuzp = tuz0 + (tduz_dt0*0.5)*t_dt;


	P1[i].x = txp;
	P1[i].y = typ;
	P1[i].z = tzp;

	P1[i].ux = tuxp;
	P1[i].uy = tuyp;
	P1[i].uz = tuzp;

	P2[i].x0 = tx0;
	P2[i].y0 = ty0;
	P2[i].z0 = tz0;

	P2[i].ux0 = tux0;
	P2[i].uy0 = tuy0;
	P2[i].uz0 = tuz0;

	P2[i].rho_ref=1.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void IBM_force_interpolation3D(Real t_dt,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type<1000)	return;		// Immersed Boundary Method

	// Stationary building: enforce zero velocity (no-slip)

	int_t icell,jcell,kcell;
	Real xi,yi,zi,uxi,uyi,uzi;
	Real rhoi;
	Real search_range,hi,tmp_A;
	// Real tmpx,tmpy,tmppx,tmppy,tmp_R;
	Real tmpux, tmpuy, tmpuz, filt;
	int p_type_i;

	p_type_i=P1[i].p_type;
	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
  	uzi=P1[i].uz;
	rhoi=P1[i].rho;

	hi=IBM_length*P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	// tmpx=tmpy=tmppx=tmppy=tmp_R=0.0;
	tmpux=tmpuy=tmpuz=filt=0.0;

	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist;
						int itype;
						int p_type_j, buffer_typej;

						xj=P1[j].x;
						yj=P1[j].y;
          				zj=P1[j].z;
						itype=P1[j].i_type;
						//if(abs(xj-0.76)<1e-3 && abs(zj-0.156)<1e-3) continue;    // for AIJ Case H
						// if(abs(xj-95)<1e-3 && abs(zj-19.5)<1e-3) continue;    // for AIJ Case HS

						if(itype!=4){
							if(P1[j].p_type==1 && P1[j].buffer_type!=99){

            				tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
								if(tdist<search_range){
									Real mj,uxj,uyj,uzj,rhoj,tmprho;
									Real hj;
									Real twij=calc_kernel_wij(tmp_A,hi,tdist);
									hj=P1[i].h;  // IBM 입자 기준으로

									// apply_MLS_filter_2D(P3[i].A,twij,&twij,xi,xj,yi,yj,hj);

									p_type_j=P1[j].p_type;

									mj=P1[j].m;
									rhoj=P1[j].rho;
									uxj=P1[j].ux;
									uyj=P1[j].uy;
            						uzj=P1[j].uz;

									tmprho=mj/rhoj;
									tmpux+=(uxj)*tmprho*twij;
									tmpuy+=(uyj)*tmprho*twij;
            						tmpuz+=(uzj)*tmprho*twij;
									filt+=tmprho*twij;
								}
							}
						}
          			}
				}
			}
		}
	}

	// preliminary velocity
	P1[i].concentration = filt;
	// if(filt<1e-3){
	// 	P1[i].fbx = 0.0;
	// 	P1[i].fby = 0.0;
	// 	P1[i].fbz = 0.0;
	// }
	// else{
	// 	P1[i].fbx = 1.225 * (uxi-tmpux/(filt+1e-10))/t_dt;
	// 	P1[i].fby = 1.225 * (uyi-tmpuy/(filt+1e-10))/t_dt;
	// 	P1[i].fbz = 1.225 * (uzi-tmpuz/(filt+1e-10))/t_dt;
	// }
	// Normal/tangential IBM forcing
	if(filt<1e-3){
		P1[i].fbx = 0.0;
		P1[i].fby = 0.0;
		P1[i].fbz = 0.0;
	}
	else{
		// Interpolated fluid velocity at IBM particle
		Real ufx = tmpux/(filt+1e-10);
		Real ufy = tmpuy/(filt+1e-10);
		Real ufz = tmpuz/(filt+1e-10);

		// Surface normal vector (precomputed, outward from building)
		Real nx = P3[i].nx;
		Real ny = P3[i].ny;
		Real nz = P3[i].nz;
		Real nmag = sqrt(nx*nx + ny*ny + nz*nz);

		if(nmag > 0.5){
			nx /= nmag;  ny /= nmag;  nz /= nmag;

			// Normal/tangential decomposition
			Real un  = ufx*nx + ufy*ny + ufz*nz;   // normal component (scalar)
			Real utx = ufx - un*nx;                  // tangential velocity
			Real uty = ufy - un*ny;
			Real utz = ufz - un*nz;

			// Normal: full correction (impenetrability)
			// Tangential: relaxed correction toward zero (partial no-slip)
			Real beta_t = 0.3;
			P1[i].fbx = (-un)*nx + beta_t*(uxi - utx);
			P1[i].fby = (-un)*ny + beta_t*(uyi - uty);
			P1[i].fbz = (-un)*nz + beta_t*(uzi - utz);
		}
		else{
			// Fallback: isotropic direct forcing (no valid normal)
			P1[i].fbx = uxi - ufx;
			P1[i].fby = uyi - ufy;
			P1[i].fbz = uzi - ufz;
		}
	}
}

__global__ void IBM_spreading_interpolation3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type!=1)	return;		// Immersed Boundary Method
	if(P1[i].buffer_type==Outlet) return;
	if(P1[i].buffer_type==99) return;	// Skip interior fluid (inside buildings)

	int_t icell,jcell,kcell;
	Real xi,yi,zi,uxi,uyi,uzi;
	Real rhoi;
	Real search_range,hi,tmp_A;
	Real tmpux, tmpuy,tmpuz,filt;
	int p_type_i;
	Real fb_x, fb_y, fb_z;
	Real tmp_twx, fb_so;

	p_type_i=P1[i].p_type;
	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
  	uzi=P1[i].uz;
	rhoi=P1[i].rho;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);

	//if(abs(xi-0.76)<1e-3 && abs(zi-0.156)<1e-3) return;      // for AIJ Case H
	//if(abs(xi-95)<1e-3 && abs(zi-19.5)<1e-3) return;         // for AIJ Case HS

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpux=tmpuy=tmpuz=filt=0.0;
	fb_x = fb_y = fb_z = 0.0;
	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist;
						int itype;
						int p_type_j, buffer_typej;

						xj=P1[j].x;
						yj=P1[j].y;
          				zj=P1[j].z;
						itype=P1[j].i_type;

						if(itype!=4){
							if(P1[j].p_type==1000){
								Real hj=IBM_length*P1[j].h;
								Real Aj=calc_tmpA(hj);
								search_range=k_search_kappa*hj;    // IBM 커널 기준 (force_interpolation과 대칭)
								
            					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
								if(tdist<search_range){
									Real mj,tdwx,tdwy,uxj,uyj,uzj,rhoj,tmprho;
									Real fb_xm, fb_ym, fb_zm;
									Real twij=calc_kernel_wij(Aj,hj,tdist);    // IBM 커널 사용 (force_interpolation과 대칭)
 
									p_type_j=P1[j].p_type;

									mj=P1[j].m;
									rhoj=P1[j].rho;
									uxj=P1[j].ux;
									uyj=P1[j].uy;
              						uzj=P1[j].uz;

									fb_xm = P1[j].fbx;
									fb_ym = P1[j].fby;
									fb_zm = P1[j].fbz;

									// if(k_kgc_solve>0){
									// 	Real twij=calc_kernel_wij(tmp_A,hi,tdist);
									// 	apply_gradient_correction_2D(P3[i].Cm,twij,tdwx,tdwy,&tdwx,&tdwy);
									// }

									tmprho=mj/rhoj;
									fb_x+=fb_xm*twij*tmprho;
									fb_y+=fb_ym*twij*tmprho;
    								fb_z+=fb_zm*twij*tmprho;
									filt+=tmprho*twij;

									// Kernel-weighted average (no volume element — IBM is 2D surface)
									// fb_x+=fb_xm*twij;
									// fb_y+=fb_ym*twij;
    								// fb_z+=fb_zm*twij;
									// filt+=twij;
								}
							}
						}
					}
				}
			}
		}
	}

	// const Real filt_min = 1e-3;   // Changed: align with force kernel
	// const Real filt_max = 2.0;    // Changed: optional upper clamp (tune)

	// Real filt_eff = fmin(fmax(filt, filt_min), filt_max); // Changed

	if(filt < 1e-12){
        P1[i].fbx = 0.0f;
        P1[i].fby = 0.0f;
        P1[i].fbz = 0.0f;
       return;
    }
	// Direct sum (not normalized): use with single-layer surface IBM only
	P1[i].fbx = fb_x;
	P1[i].fby = fb_y;
  	P1[i].fbz = fb_z;
}
////////////////////////////////////////////////////////////////////////
__global__ void IBM_corrector(Real t_dt,Real ttime,part1*P1,part1*TP1,part2*P2,part3*P3)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type<1000) return;

	int_t p_typei;
	Real tx0,ty0,tz0,txc,tyc,tzc;
	Real tux0,tuy0,tuz0,tuxc,tuyc,tuzc;
	Real tdux_dt,tduy_dt,tduz_dt;

	// Get previous step
	tx0 = P2[i].x0;
	ty0 = P2[i].y0;
	tz0 = P2[i].z0;

	tux0=P2[i].ux0;
	tuy0=P2[i].uy0;
	tuz0=P2[i].uz0;

	//////// force to solid for no-slip boundary condition ////////
	// tdux_dt=P3[i].ftotalx-P1[i].fbx/P1[i].rho;
	// tduy_dt=P3[i].ftotaly-P1[i].fby/P1[i].rho;
	// tduz_dt=P3[i].ftotalz-P1[i].fbz/P1[i].rho;
	// two way coupling
	// tdux_dt=P1[i].ax-P1[i].fbx/P1[i].rho;
	// tduy_dt=P1[i].ay-P1[i].fby/P1[i].rho;
	// tduz_dt=P1[i].az-P1[i].fbz/P1[i].rho;
	// one way coupling


	// tdux_dt=P1[i].ax;
	// tduy_dt=P1[i].ay;
	// tduz_dt=P1[i].az;

	tdux_dt=0.0;
	tduy_dt=0.0;
	tduz_dt=0.0;


	// if(P1[i].ccc==3){
	// 	tux0=0.0;
	// 	tuy0=0.0;
	// 	tuz0=0.0;
	// 	tdux_dt=tduy_dt=tduz_dt=0.0;
	// }

	//update velocity

	tuxc = tux0 + tdux_dt*t_dt;
	tuyc = tuy0 + tduy_dt*t_dt;
	tuzc = tuz0 + tduz_dt*t_dt;

	// update position

	txc = tx0 + tuxc*(t_dt);
	tyc = ty0 + tuyc*(t_dt);
	tzc = tz0 + tuzc*(t_dt);

	// Stationary building: copy position from sorted array, velocity = 0
	TP1[i].x = P1[i].x;
	TP1[i].y = P1[i].y;
	TP1[i].z = P1[i].z;

	TP1[i].ux = 0.0;
	TP1[i].uy = 0.0;
	TP1[i].uz = 0.0;

	// if(P1[i].p_type>=2000)	TP1[i].rho=P2[i].rho_ref/(P1[i].jacob-1e-20);

	TP1[i].fbx = P1[i].fbx;
	TP1[i].fby = P1[i].fby;
	TP1[i].fbz = P1[i].fbz;
	//
	// TP1[i].fb_x += P1[i].fbx;
	// TP1[i].fb_y += P1[i].fby;
	// TP1[i].fb_z += P1[i].fbz;

	// TP1[i].ux_i = P1[i].ux_i;
	// TP1[i].uy_i = P1[i].uy_i;
	// TP1[i].uz_i = P1[i].uz_i;


	// TP1[i].ax = P1[i].ax;
	// TP1[i].ay = P1[i].ay;
	// TP1[i].az = P1[i].az;

	TP1[i].concentration=P1[i].concentration;

	TP1[i].vol=P1[i].vol;
	TP1[i].pres=P1[i].pres;

	//Structure modeling

	// TP1[i].F11=P1[i].F11;
	// TP1[i].F12=P1[i].F12;
	// TP1[i].F13=P1[i].F13;
	// TP1[i].F21=P1[i].F21;
	// TP1[i].F22=P1[i].F22;
	// TP1[i].F23=P1[i].F23;
	// TP1[i].F31=P1[i].F31;
	// TP1[i].F32=P1[i].F32;
	// TP1[i].F33=P1[i].F33;

	// TP1[i].px=P1[i].px;
	// TP1[i].py=P1[i].py;
	// TP1[i].pz=P1[i].pz;

	// TP1[i].jacob=P1[i].jacob;


}

__global__ void IBM_corrector_fluid(Real t_dt,Real ttime,part1*P1,part2*P2,part3*P3)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type!=1) return;
	if(P1[i].buffer_type==99) return;	// Skip interior fluid

	// Direct forcing: velocity correction [m/s]
	P1[i].ux += P1[i].fbx;
	P1[i].uy += P1[i].fby;
	P1[i].uz += P1[i].fbz;

}

// Zero velocity for fluid particles flagged as interior (inside buildings)
__global__ void KERNEL_zero_interior_velocity(part1* P1)
{
	int_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= k_num_part2) return;
	if (P1[i].buffer_type != 99) return;

	P1[i].ux = 0.0;
	P1[i].uy = 0.0;
	P1[i].uz = 0.0;
	P1[i].pres = 0.0;
}

// Monitor a specific particle each timestep
__global__ void KERNEL_monitor_particle(part1* P1, part3* P3, int step, int phase)
{
	int_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= k_num_part2) return;

	// Target: (377.5, 32.5, 2.5)
	Real tx = 377.5, ty = 32.5, tz = 2.5;
	Real dist = (P1[i].x-tx)*(P1[i].x-tx) + (P1[i].y-ty)*(P1[i].y-ty) + (P1[i].z-tz)*(P1[i].z-tz);
	if (dist > 1.0) return;

	// phase: 0=after_prep, 1=after_predictor, 2=after_IBM, 3=after_PPE, 4=after_corrector
	const char* ph[] = {"PREP", "PRED", "IBM ", "PPE ", "CORR"};
	printf("[MON] stp=%d %s  pos=(%.1f,%.1f,%.1f)  v=(%.4f,%.4f,%.4f)  p=%.2f  "
	       "fb=(%.4f,%.4f,%.4f)  fv=(%.4f,%.4f,%.4f)  fp=(%.4f,%.4f,%.4f)  bt=%d\n",
	       step, ph[phase],
	       P1[i].x, P1[i].y, P1[i].z,
	       P1[i].ux, P1[i].uy, P1[i].uz,
	       P1[i].pres,
	       P1[i].fbx, P1[i].fby, P1[i].fbz,
	       P1[i].fv_x, P1[i].fv_y, P1[i].fv_z,
	       P1[i].fp_x, P1[i].fp_y, P1[i].fp_z,
	       P1[i].buffer_type);
}

__global__ void KERNEL_change_fluid_to_IBM(part1* P1){
	int_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= k_num_part2) return;
	if (P1[i].p_type == 1000) P1[i].i_type=3;
	if (P1[i].buffer_type != 99) return;

	P1[i].p_type = 1000;
}
