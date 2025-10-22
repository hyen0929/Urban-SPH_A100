// ////////////////////////////////////////////////////////////////////////
#define DEBUG_INDEX 11111

__global__ void KERNEL_clc_prep3D(int_t*g_str,int_t*g_end,part1*P1, part2*P2, part3*P3, int_t tcount, Real tdt)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	//if(k_open_boundary>0 && P1[i].buffer_type>0) return;
	if(P1[i].p_type>=1000) return;
	if(P1[i].p_type<=0)	return;

	int_t icell,jcell,kcell;
	int_t ptypei=P1[i].p_type;

	Real xi,yi,zi;
	Real uxi,uyi,uzi;
	Real mi=P1[i].m;
	Real rhoi=P1[i].rho;
	Real cci=P3[i].cc;
	Real hi,tmp_A,search_range;
	Real tmp_flt,tmp_SR;
	Real tmp_rhox,tmp_rhoy,tmp_rhoz;
	Real tvis_t=0.0,th;
	Real tempi=P1[i].temp;

	// for LES turbulence model (CTS)
	Real visi;
	Real tmp_sigma;
	Real gc11, gc12, gc22, gc23, gc31, gc33;
	Real st11, st12, st22, st23, st31, st33;

  	Real ki = fmax(P1[i].k_turb, 0.0);
  	Real nui = P1[i].vis_t;
  	Real akti = nui/Prt_sgs;
  	Real dTdx, dTdy, dTdz;
  	Real Tk;

  	dTdx=dTdy=dTdz=Tk=0.0;

	gc11=P3[i].inv_cm_xx, gc12=P3[i].inv_cm_xy;
	gc22=P3[i].inv_cm_yy, gc23=P3[i].inv_cm_yz;
	gc31=P3[i].inv_cm_zx, gc33=P3[i].inv_cm_zz;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	th=hi*L_SPS;
	search_range=k_search_kappa*hi;	// search range
	visi=viscosity2(ptypei,tempi);

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	// 초기화
	tmp_flt=tmp_SR=tmp_rhox=tmp_rhoy=tmp_rhoz=0.0;

	tmp_sigma=0.0;
	st11=st12=st22=st23=st31=st33=0.0;

  	int_t ncell=P1[i].ncell;

	// 계산
	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
        		if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,uxj,uyj,uzj,uij2,mj,rhoj,ccj,tdwx,tdwy,tdwz,tmp_wij,tmp_dwij,tdist,tmp_val;
						int_t ptypej;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						mj=P1[j].m;
						rhoj=P1[j].rho;
						ptypej=P1[j].p_type;
						ccj=P3[j].cc;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;

						if(tdist<search_range){
							tmp_wij=calc_kernel_wij(tmp_A,hi,tdist);
							tmp_dwij=calc_kernel_dwij(tmp_A,hi,tdist);

							tdwx=tmp_dwij*(xi-xj)/tdist;
							tdwy=tmp_dwij*(yi-yj)/tdist;
							tdwz=tmp_dwij*(zi-zj)/tdist;

							// filter
							if((tcount%k_freq_filt)==0) tmp_flt+=mj/rhoj*tmp_wij;

							// strain rate
							if((k_fv_solve!=0)&&(k_turbulence_model!=Laminar))
							{
								// 벽면 Strain 계산 부정확;
								// uij2=(uxi-uxj)*(uxi-uxj)+(uyi-uyj)*(uyi-uyj)+(uzi-uzj)*(uzi-uzj);
								// tmp_val=-0.5*mj*(rhoi+rhoj)*uij2;
								// tmp_val/=(rhoi*rhoj*tdist*tdist);
								// tmp_SR+=tmp_val*(xi-xj)*tdwx+tmp_val*(yi-yj)*tdwy+tmp_val*(zi-zj)*tdwz;

								tmp_sigma+=tdist*tdist*tmp_wij*mj/rhoj;

								if(ptypej==1){

									Real tdwxc, tdwyc, tdwzc;
                  					Real kj, nuj, akk, tmprr;

									tdwxc=gc11*tdwx+gc12*tdwy+gc31*tdwz;
									tdwyc=gc12*tdwx+gc22*tdwy+gc23*tdwz;
									tdwzc=gc31*tdwx+gc23*tdwy+gc33*tdwz;

									st11+=mj/rhoj*(uxj-uxi)*tdwxc;
									st12+=0.5*mj/rhoj*((uxj-uxi)*tdwyc+(uyj-uyi)*tdwxc);
									st22+=mj/rhoj*(uyj-uyi)*tdwyc;
									st23+=0.5*mj/rhoj*((uyj-uyi)*tdwzc+(uzj-uzi)*tdwyc);
									st31+=0.5*mj/rhoj*((uxj-uxi)*tdwzc+(uzj-uzi)*tdwxc);
									st33+=mj/rhoj*(uzj-uzi)*tdwzc;

                  					dTdx += mj/rhoj*(P1[j].temp - tempi)*tdwxc;                       // Changed
                  					dTdy += mj/rhoj*(P1[j].temp - tempi)*tdwyc;                       // Changed
                  					dTdz += mj/rhoj*(P1[j].temp - tempi)*tdwzc;

                  					kj = fmax(P1[j].k_turb,0.0);
                  					nuj = P1[j].vis_t;
                  					akk = 0.5*(nui+nuj)/Sigk_sgs;
                  					tmprr = ((xi-xj)*tdwxc + (yi-yj)*tdwyc + (zi-zj)*tdwzc)/(tdist*tdist+ 1e-12);

                  					//Tk += akk*(ki-kj)*tmprr*mj/rhoj;
								}
							}
						}
					}
				}
			}
		}
	}

	// strain_rate
	if((k_fv_solve!=0)&&(k_turbulence_model!=Laminar)){
		// tmp_SR=max(1e-20,tmp_SR);
		// tmp_SR=sqrt(tmp_SR);
		tmp_SR=2*(st11*st11+st22*st22+st33*st33+2*st12*st12+2*st31*st31+2*st23*st23);
		tmp_SR=sqrt(tmp_SR);
		P2[i].SR=tmp_SR;

   		tmp_sigma=0.3333*tmp_sigma;
    	tmp_sigma=sqrt(tmp_sigma);
    	Real delta_m = cbrt(mi/rhoi);

    	Real tmp_cl=max(tmp_sigma, hi);
    	//printf("tmp_sigma=%f  hi=%f  tmp_cl=%f\n",tmp_sigma,hi,tmp_cl);

		if(k_turbulence_model==SSM){
			tvis_t=Cs_sgs*Cs_sgs*tmp_cl*tmp_cl*tmp_SR;
			if(P1[i].buffer_type!=1)P1[i].vis_t=tvis_t;
		}
    	else if(k_turbulence_model==DDF){
      		Real p = nui*tmp_SR*tmp_SR;
      		Real b = 0.0;
      		// Real b = -(Gravitational_CONST/Theta0_ref)*(akti)*dTdz;
      		Real eps = Ce_sgs*ki*sqrt(ki)/tmp_cl;
      		Real dk = p + b - eps + Tk;
      		ki += dk*tdt;
      		ki=max(ki,1e-4);    // Clipping for negative TKE
      		tvis_t = Ck_sgs*tmp_cl*sqrt(ki);
      		if(isnan(tvis_t)) printf("p=%f eps=%f dk=%f ki=%f tmp_cl=%f\n",p,eps,dk,ki,tmp_cl);

      		if(P1[i].buffer_type!=1){
        		P1[i].vis_t=tvis_t;
        		P1[i].k_turb=ki;
        		P1[i].e_turb=eps;
      		}
		}
	}

	// reference density
	if(tcount==0){
		P2[i].rho_ref=P1[i].rho;
	}

	// filter
	if((tcount%k_freq_filt)==0) P1[i].flt_s=tmp_flt;
}

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_clc_correction_KGC_3D(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;

	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real search_range,hi,tmp_A;;
	Real tmpxx,tmpyy,tmpzz,tmpxy,tmpyz,tmpzx;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpxx=tmpyy=tmpzz=0;
	tmpxy=tmpyz=tmpzx=0;
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
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
						if(tdist>0&&tdist<search_range){
							Real tdwij,mj,rhoj,txx,txy,tyy,tzx,tyz,tzz,rtd,mtd;
							tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
							mj=P1[j].m;
							rhoj=P1[j].rho;

							mtd=mj*tdwij;
							rtd=1.0/(rhoj*tdist);

							txx=-mtd*(xi-xj)*(xi-xj);
							txx*=rtd;
							txy=-mtd*(yi-yj)*(xi-xj);
							txy*=rtd;
							tyy=-mtd*(yi-yj)*(yi-yj);
							tyy*=rtd;
							tzx=-mtd*(xi-xj)*(zi-zj);
							tzx*=rtd;
							tyz=-mtd*(yi-yj)*(zi-zj);
							tyz*=rtd;
							tzz=-mtd*(zi-zj)*(zi-zj);
							tzz*=rtd;
							tmpxx+=txx;
							tmpxy+=txy;
							tmpyy+=tyy;
							tmpzx+=tzx;
							tmpyz+=tyz;
							tmpzz+=tzz;
						}
					}
				}
			}
		}
	}
	// save values to particle array
	Real tmpcmd;
	tmpcmd=tmpxx*(tmpyy*tmpzz-tmpyz*tmpyz);
	tmpcmd-=tmpxy*(tmpxy*tmpzz-tmpyz*tmpzx);
	tmpcmd+=tmpzx*(tmpxy*tmpyz-tmpyy*tmpzx);

	if(abs(tmpcmd)>Min_det){
		Real rtcmd=1.0/tmpcmd;
		P3[i].Cm[0][0]=(tmpyy*tmpzz-tmpyz*tmpyz)*rtcmd;
		P3[i].Cm[0][1]=(tmpzx*tmpyz-tmpxy*tmpzz)*rtcmd;
		P3[i].Cm[0][2]=(tmpxy*tmpyz-tmpzx*tmpyy)*rtcmd;
		P3[i].Cm[1][0]=(tmpzx*tmpyz-tmpxy*tmpzz)*rtcmd;
		P3[i].Cm[1][1]=(tmpxx*tmpzz-tmpzx*tmpzx)*rtcmd;
		P3[i].Cm[1][2]=(tmpzx*tmpxy-tmpxx*tmpyz)*rtcmd;
		P3[i].Cm[2][0]=(tmpxy*tmpyz-tmpzx*tmpyy)*rtcmd;
		P3[i].Cm[2][1]=(tmpzx*tmpxy-tmpxx*tmpyz)*rtcmd;
		P3[i].Cm[2][2]=(tmpxx*tmpyy-tmpxy*tmpxy)*rtcmd;
	}
	else{
		P3[i].Cm[0][0]=1;
		P3[i].Cm[0][1]=0;
		P3[i].Cm[0][2]=0;
		P3[i].Cm[1][0]=0;
		P3[i].Cm[1][1]=1;
		P3[i].Cm[1][2]=0;
		P3[i].Cm[2][0]=0;
		P3[i].Cm[2][1]=0;
		P3[i].Cm[2][2]=1;
	}
}

void gradient_correction(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
{
	dim3 b,t;
	t.x=128;
	b.x=(num_part2-1)/t.x+1;
	switch (kgc_solve){
		case KGC:
		if(dim==2) KERNEL_clc_correction_KGC_2D<<<b,t>>>(g_str,g_end,P1,P3);
		if(dim==3) KERNEL_clc_correction_KGC_3D<<<b,t>>>(g_str,g_end,P1,P3);
		cudaDeviceSynchronize();
		break;
		case FPM:
		// if(dim==2) KERNEL_clc_correction_FPM_2D<<<b,t>>>(g_str,g_end,P1,P3);
		// if(dim==3) KERNEL_clc_correction_FPM_3D<<<b,t>>>(inout,g_str,g_end,P1,P3);
		// cudaDeviceSynchronize();
		break;
		case DFPM:
		// if(dim==2) KERNEL_clc_correction_DFPM_2D<<<b,t>>>(g_str,g_end,P1,P3);
		// if(dim==3) KERNEL_clc_correction_DFPM_3D<<<b,t>>>(inout,g_str,g_end,P1,P3);
		// cudaDeviceSynchronize();
		break;
		case KGF:
		// if(dim==2) KERNEL_clc_correction_KGF_2D<<<b,t>>>(g_str,g_end,P1,P3);
		// if(dim==3) KERNEL_clc_correction_KGF_3D<<<b,t>>>(inout,g_str,g_end,P1,P3);
		// cudaDeviceSynchronize();
		break;
		default:
		if(dim==2) KERNEL_clc_correction_KGC_2D<<<b,t>>>(g_str,g_end,P1,P3);
		if(dim==3) KERNEL_clc_correction_KGC_3D<<<b,t>>>(g_str,g_end,P1,P3);
		cudaDeviceSynchronize();
		break;
	}
}

////////////////////////////////////////////////////////////////////////
// calcuate color field for two-phase flow surface tension model (2017.05.08 jyb)
__global__ void KERNEL_clc_color_field3D(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;


	int_t ptypei;
	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real tmpn,tmpd;
	Real search_range,hi,tmp_A;

	ptypei=P1[i].p_type;
	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpn=tmpd=0.0;

	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				//int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist; //rr
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));

						if(tdist<search_range){
							Real twij,mj,rhoj;
							twij=calc_kernel_wij(tmp_A,hi,tdist);

							mj=P1[j].m;
							rhoj=P1[j].rho;

							tmpn+=mj*twij*(ptypei==P1[j].p_type)/rhoj;
							tmpd+=mj*twij/rhoj;
						}
					}
				}
			}
		}
	}
	P3[i].cc=tmpn/tmpd;
}

__global__ void KERNEL_clc_IB_normal_vector3D(int_t* g_str, int_t* g_end, part1* P1, part3* P3)
{
    int_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;

    // Process only IBM marker particles
    if (P1[i].p_type < 1000) return;

    int_t icell, jcell, kcell;
    Real xi, yi, zi;
    Real hi, tmp_A, search_range;
    Real nx, ny, nz; // Normal vector components

    xi = P1[i].x;
    yi = P1[i].y;
    zi = P1[i].z;
    hi = P1[i].h;
    tmp_A = calc_tmpA(hi);
    search_range = k_search_kappa * hi;

    // Find cell index for particle i
    if ((k_x_max == k_x_min)) { icell = 0; }
    else { icell = min(floor((xi - k_x_min) / k_dcell), k_NI - 1); }
    if ((k_y_max == k_y_min)) { jcell = 0; }
    else { jcell = min(floor((yi - k_y_min) / k_dcell), k_NJ - 1); }
    if ((k_z_max == k_z_min)) { kcell = 0; }
    else { kcell = min(floor((zi - k_z_min) / k_dcell), k_NK - 1); }

    if (icell < 0) icell = 0; if (jcell < 0) jcell = 0; if (kcell < 0) kcell = 0;

    // Initialize normal vector
    nx = ny = nz = 0.0;

    int_t ncell = P1[i].ncell;

	if(abs(xi-0.764)<1e-5) nx=-1.0;
	if(abs(xi-0.836)<1e-5) nx=1.0;
	if(abs(yi+0.036)<1e-5) ny=-1.0;
	if(abs(yi-0.036)<1e-5) ny=1.0;
	if(abs(zi-0.0)<1e-5) nz=-1.0;
	if(abs(zi-0.152)<1e-5) nz=1.0;

    // // Loop over neighboring cells
    // for (int_t z = -ncell; z <= ncell; z++) {
    //     for (int_t y = -ncell; y <= ncell; y++) {
    //         for (int_t x = -ncell; x <= ncell; x++) {
    //             int_t k = idx_cell(icell + x, jcell + y, kcell + z);

    //             if (((icell + x) < 0) || ((icell + x) > (k_NI - 1)) || ((jcell + y) < 0) || ((jcell + y) > (k_NJ - 1)) || ((kcell + z) < 0) || ((kcell + z) > (k_NK - 1))) continue;
    //             if (g_str[k] != cu_memset) {
    //                 int_t fend = g_end[k];
    //                 for (int_t j = g_str[k]; j < fend; j++) {
    //                     // Interact only with fluid particles
    //                     if (P1[j].p_type < 1000) {
    //                         Real xj, yj, zj, tdist;
    //                         xj = P1[j].x;
    //                         yj = P1[j].y;
    //                         zj = P1[j].z;

    //                         tdist = sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj));

    //                         if (tdist > 1e-9 && tdist < search_range) {
    //                             Real tmp_dwij, mj, rhoj;
    //                             tmp_dwij = calc_kernel_dwij(tmp_A, hi, tdist);
    //                             mj = P1[j].m;
    //                             rhoj = P1[j].rho;

    //                             // Accumulate the gradient of the color field
    //                             // This is equivalent to summing the kernel gradients from fluid particles
    //                             Real vol_j = mj / rhoj;
    //                             nx += vol_j * tmp_dwij * (xi - xj) / tdist;
    //                             ny += vol_j * tmp_dwij * (yi - yj) / tdist;
    //                             nz += vol_j * tmp_dwij * (zi - zj) / tdist;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // Normalize the vector
    Real mag = sqrt(nx * nx + ny * ny + nz * nz);
    if (mag > 1e-9) {
        P3[i].nx = nx / mag;
        P3[i].ny = ny / mag;
        P3[i].nz = nz / mag;
    } else {
        P3[i].nx = 0.0;
        P3[i].ny = 0.0;
        P3[i].nz = 0.0;
    }
}