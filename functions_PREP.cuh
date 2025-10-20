// ////////////////////////////////////////////////////////////////////////
#define DEBUG_INDEX 11111

// __global__ void KERNEL_clc_prep2D(int_t*g_str,int_t*g_end,part1*P1, part2*P2, part3*P3, int_t tcount)
// {
// 	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part2) return;
// 	if(P1[i].i_type>i_type_crt) return;
// 	// if(P1[i].eli>1) return;
// 	// if((P1[i].eli>1)&(P1[i].eli>1)) return;
// 	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
//
// 	int_t icell,jcell;
// 	int_t ptypei=P1[i].p_type;
//
// 	Real xi,yi;
// 	Real uxi,uyi;
// 	Real mi=P1[i].m;
// 	Real rhoi=P1[i].rho;
// 	Real cci=P3[i].cc;
// 	Real hi,tmp_Ai,search_range;
// 	Real tmp_flt,tmp_SR;
// 	Real tmp_rhox,tmp_rhoy;
// 	Real tvis_t=0.0,th;
// 	Real tmp_uxx, tmp_uxy, tmp_uyx, tmp_uyy;
//
// 	Real st11, st12, st22;
//
// 	xi=P1[i].x;
//   yi=P1[i].y;
//
// 	uxi=P1[i].ux;
// 	uyi=P1[i].uy;
//
// 	hi=P1[i].h;
// 	tmp_Ai=calc_tmpA(hi);
// 	th=hi*L_SPS;
// 	search_range=k_search_kappa*hi;	// search range
//
// 	// calculate I,J,K in cell
// if((k_x_max==k_x_min)){icell=0;}
// else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
// if((k_y_max==k_y_min)){jcell=0;}
// else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
//
// 	// out-of-range handling
// if(icell<0) icell=0;	if(jcell<0) jcell=0;
//
// 	// 변수 초기화
// tmp_flt=tmp_SR=tmp_rhox=tmp_rhoy=0.0;
// tmp_uxx=tmp_uxy=tmp_uyx=tmp_uyy=0.0;
//
// 	// 계산
// for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
// 	for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
// 			// int_t k=(icell+x)+k_NI*(jcell+y);
// 		int_t k=idx_cell(icell+x,jcell+y,0);
//
// 		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
//
// 		if(g_str[k]!=cu_memset){
// 			int_t fend=g_end[k];
// 				for(int_t j=g_str[k];j<fend;j++){
// 					Real xj,yj,uxj,uyj,uij2,mj,rhoj,ccj,tdwx,tdwy,tmp_wij,tmp_dwij,tdist,tmp_val,hj;
// 					int_t ptypej;
// 					int itype;
// 					Real volj;
//
// 					xj=P1[j].x;
// 					yj=P1[j].y;
// 					uxj=P1[j].ux;
// 					uyj=P1[j].uy;
// 					mj=P1[j].m;
// 					rhoj=P1[j].rho;
// 					ptypej=P1[j].p_type;
// 					ccj=P3[j].cc;
// 					itype=P1[j].i_type;
// 					volj=P1[j].vol;
// 					hj=P1[j].h;
//
// 					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;
// 					if(itype!=4){
// 						if(tdist<search_range){
//
// 						// tmp_wij=calc_kernel_wij(tmp_A,hi,tdist);
// 						// tmp_dwij=calc_kernel_dwij(tmp_A,hi,tdist);
//
// 						// for multi resolution (KDH)
// 							Real tmp_Aj=calc_tmpA(hj);
// 							tmp_wij=0.5*(calc_kernel_wij(tmp_Ai,hi,tdist)+calc_kernel_wij(tmp_Aj,hj,tdist));
// 							tmp_dwij=0.5*(calc_kernel_dwij(tmp_Ai,hi,tdist)+calc_kernel_dwij(tmp_Aj,hj,tdist));
//
// 						// dwij
// 							tdwx=tmp_dwij*(xi-xj)/tdist;
// 							tdwy=tmp_dwij*(yi-yj)/tdist;
//
// 						// filter
// 							if((tcount%k_freq_filt)==0) tmp_flt+=mj/rhoj*tmp_wij;
//
// 						// strain rate
// 							if((k_fv_solve!=0)&&(k_turbulence_model!=Laminar)){
// 								uij2=(uxi-uxj)*(uxi-uxj)+(uyi-uyj)*(uyi-uyj);
// 								tmp_val=-0.5*mj*(rhoi+rhoj)*uij2;
// 								tmp_val/=(rhoi*rhoj*tdist*tdist);
// 								tmp_SR+=tmp_val*(xi-xj)*tdwx+tmp_val*(yi-yj)*tdwy;
// 							}
// 							if(k_fv_solve==Cleary){
//
// 								tmp_uxx+=-volj*(uxi-uxj)*tdwx;
// 								tmp_uxy+=-volj*(uxi-uxj)*tdwy;
// 								tmp_uyx+=-volj*(uyi-uyj)*tdwx;
// 								tmp_uyy+=-volj*(uyi-uyj)*tdwy;
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
//
// 	// reference density
// 	P2[i].rho_ref=1000.0;
//
// 	// filter
// 	if((tcount%k_freq_filt)==0) P1[i].flt_s=tmp_flt;
// 	if(k_fv_solve==Cleary){
//
// 		P3[i].Gxx=tmp_uxx;
// 		P3[i].Gxy=tmp_uxy;
// 		P3[i].Gyx=tmp_uyx;
// 		P3[i].Gyy=tmp_uyy;
// 	}
//}

__global__ void KERNEL_clc_prep2D(int_t* g_str, int_t* g_end, part1* P1, part2* P2, part3* P3, int_t tcount)
{
    int_t i=threadIdx.x+blockIdx.x*blockDim.x;
    if (i>=k_num_part2) return;
    if (P1[i].i_type>i_type_crt) return;
    if (k_open_boundary>0 && P1[i].buffer_type>0) return;
    if (P1[i].p_type>=1000) return;
    if (P1[i].p_type<=0) return;

    int_t icell, jcell;
    int_t ptypei=P1[i].p_type;

    Real xi, yi;
    Real uxi, uyi;
    Real mi=P1[i].m;
    Real rhoi=P1[i].rho;
    Real cci=P3[i].cc;
    Real hi, tmp_A, search_range;
    Real tmp_flt, tmp_SR;
    Real tmp_rhox, tmp_rhoy;
    Real tvis_t=0.0, th;
    Real tempi=P1[i].temp;

    // For LES turbulence model
    Real visi;
    Real tmp_sigma;
    Real gc11, gc12, gc22;
    Real st11, st12, st22;

    gc11=P3[i].inv_cm_xx;
    gc12=P3[i].inv_cm_xy;
    gc22=P3[i].inv_cm_yy;

    xi=P1[i].x;
    yi=P1[i].y;

    uxi=P1[i].ux;
    uyi=P1[i].uy;

    hi=P1[i].h;
    tmp_A=calc_tmpA(hi);
    th=hi*L_SPS;
    search_range=k_search_kappa*hi; // search range
    visi=viscosity2(ptypei, tempi);

    // Calculate I, J in cell
    if ((k_x_max==k_x_min)){icell=0;}
    else {icell=min(floor((xi-k_x_min)/k_dcell), k_NI-1);}
    if ((k_y_max==k_y_min)){jcell=0;}
    else {jcell=min(floor((yi-k_y_min)/k_dcell), k_NJ-1);}
    // Out-of-range handling
    if (icell<0) icell=0;
    if (jcell<0) jcell=0;

    // Initialization
    tmp_flt=tmp_SR=tmp_rhox=tmp_rhoy=0.0;

    tmp_sigma=0.0;
    st11=st12=st22=0.0;

    // Computation
    for (int_t y=-P1[i].ncell; y<=P1[i].ncell; y++){
        for (int_t x=-P1[i].ncell; x<=P1[i].ncell; x++){
            int_t k=idx_cell(icell+x, jcell+y,0);

            if (((icell+x)<0) || ((icell+x)>(k_NI-1)) || ((jcell+y)<0) || ((jcell+y)>(k_NJ-1))) continue;
            if (g_str[k]!=cu_memset){
                int_t fend=g_end[k];
                for (int_t j=g_str[k]; j<fend; j++){
                    Real xj, yj, uxj, uyj, uij2, mj, rhoj, ccj, tdwx, tdwy, tmp_wij, tmp_dwij, tdist, tmp_val;
                    int_t ptypej;

                    xj=P1[j].x;
                    yj=P1[j].y;
                    uxj=P1[j].ux;
                    uyj=P1[j].uy;
                    mj=P1[j].m;
                    rhoj=P1[j].rho;
                    ptypej=P1[j].p_type;
                    ccj=P3[j].cc;

                    if(fabs(rhoj) < 1e-8) continue;

                    tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;

                    if (tdist<search_range){
                        tmp_wij=calc_kernel_wij(tmp_A, hi, tdist);
                        tmp_dwij=calc_kernel_dwij(tmp_A, hi, tdist);

                        tdwx=tmp_dwij*(xi-xj)/tdist;
                        tdwy=tmp_dwij*(yi-yj)/tdist;

                        // Filter
                        if ((tcount%k_freq_filt)==0) tmp_flt+=mj/rhoj*tmp_wij;

                        // Strain rate
                        if ((k_fv_solve!=0) && (k_turbulence_model!=Laminar)){
                            tmp_sigma+=tdist*tdist*tmp_wij*mj/rhoj;

                            if (ptypej==1){
                                Real tmp_xc, tmp_yc;
                                Real tdwxc, tdwyc;

                                tmp_xc=tmp_dwij*(xi-xj)/tdist;
                                tmp_yc=tmp_dwij*(yi-yj)/tdist;

                                tdwxc=gc11*tmp_xc+gc12*tmp_yc;
                                tdwyc=gc12*tmp_xc+gc22*tmp_yc;

                                st11+=mj/rhoj*(uxj-uxi)*tdwxc;
                                st12+=0.5*mj/rhoj*((uxj-uxi)*tdwyc+(uyj-uyi)*tdwxc);
                                st22+=mj/rhoj*(uyj-uyi)*tdwyc;
                           }
                       }
                   }
               }
           }
       }
   }

    // Strain rate magnitude
    if ((k_fv_solve!=0) && (k_turbulence_model!=Laminar)){
        tmp_SR=2*(st11*st11+st22*st22+2*st12*st12);
        tmp_SR=sqrt(tmp_SR);
        P2[i].SR=tmp_SR;

        if (k_turbulence_model==SSM){
            tmp_sigma=0.5*tmp_sigma; // Adjusted for 2D
            tmp_sigma=sqrt(tmp_sigma);

            Real tmp_cl, tmp_length;

            tmp_cl=0.0;
            tmp_length=100.0;

            tmp_cl=min(Cs_SPS*tmp_sigma, 0.41*tmp_length);
            tvis_t=0.2*tmp_cl*tmp_cl*tmp_SR;   // Cs=0.1

            P1[i].vis_t=tvis_t;
       }
   }

    // Reference density
    if (tcount==0){
        P2[i].rho_ref=P1[i].rho;
   }

    // Filter
    if ((tcount%k_freq_filt)==0) P1[i].flt_s=tmp_flt;
}

////////////////////////////////////////////////////////////////////////
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
__global__ void KERNEL_clc_correction_KGC_2D(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;

	int_t icell,jcell;
	Real search_range,hi,tmp_A;
	Real xi,yi;
	Real tmpxx,tmpyy,tmpxy;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;

	tmpxx=tmpyy=tmpxy=0;

	for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
			//int_t k=(icell+x)+k_NI*(jcell+y);
			int_t k=idx_cell(icell+x,jcell+y,0);
			if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
			if(g_str[k]!=cu_memset){
				int_t fend=g_end[k];
				for(int_t j=g_str[k];j<fend;j++){
					if(P1[j].p_type<1000){

						Real xj,yj,tdist;
						xj=P1[j].x;
						yj=P1[j].y;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
						if(tdist>0&&tdist<search_range){
							Real tdwij,mj,rhoj,txx,txy,tyy,rtd,mtd;
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

							tmpxx+=txx;
							tmpxy+=txy;
							tmpyy+=tyy;
						}
					}
				}
			}
		}
	}
	// save values to particle array

	Real tmpcmd=tmpxx*tmpyy-tmpxy*tmpxy;
	if(abs(tmpcmd)>Min_det){
		Real rtcmd=1.0/tmpcmd;
		P3[i].Cm[0][0]=tmpyy*rtcmd;
		P3[i].Cm[0][1]=-tmpxy*rtcmd;
		P3[i].Cm[1][0]=-tmpxy*rtcmd;
		P3[i].Cm[1][1]=tmpxx*rtcmd;
	}
	else{
		P3[i].Cm[0][0]=1;
		P3[i].Cm[0][1]=0;
		P3[i].Cm[1][0]=0;
		P3[i].Cm[1][1]=1;
	}
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
__global__ void KERNEL_clc_color_field2D(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;

	int_t ptypei;
	int_t icell,jcell;
	Real xi,yi;
	Real tmpn,tmpd;
	Real search_range,hi,tmp_A;

	ptypei=P1[i].p_type;
	xi=P1[i].x;
	yi=P1[i].y;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;

	// calculate I,J,K in cell
if((k_x_max==k_x_min)){icell=0;}
else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
if((k_y_max==k_y_min)){jcell=0;}
else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;

	tmpn=tmpd=0.0;

	for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
			// int_t k=(icell+x)+k_NI*(jcell+y);
			int_t k=idx_cell(icell+x,jcell+y,0);

			if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
			if(g_str[k]!=cu_memset){
				int_t fend=g_end[k];
				for(int_t j=g_str[k];j<fend;j++){
					Real xj,yj,tdist; //rr
					xj=P1[j].x;
					yj=P1[j].y;

					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
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
	P3[i].cc=tmpn/tmpd;
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

__global__ void KERNEL_variable_smoothing_length2D(int_t*g_str,int_t*g_end,part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if(i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  if(P1[i].p_type<1) return;

  int_t icell,jcell;
  Real xi,yi,hi,h_ref_i;
  Real search_range,tmp_A,tmp_R,tmp_flt,tmp_RR;

  tmp_R=tmp_flt=0.0;

  Real scale = 1.0;
  hi=scale*P1[i].h;
  h_ref_i=P1[i].h_ref;
  tmp_A=calc_tmpA(hi);
  search_range=k_search_kappa*hi;	// search range

  xi=P1[i].x;
  yi=P1[i].y;

  // calculate I,J,K in cell
  if((k_x_max==k_x_min)){icell=0;}
  else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  if((k_y_max==k_y_min)){jcell=0;}
  else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  // out-of-range handling
  if(icell<0) icell=0;	if(jcell<0) jcell=0;

  tmp_RR=0.0;
  for(int_t y=-scale*P1[i].ncell;y<=scale*P1[i].ncell;y++){
		for(int_t x=-scale*P1[i].ncell;x<=scale*P1[i].ncell;x++){
      // int_t k=(icell+x)+k_NI*(jcell+y);
      int_t k=idx_cell(icell+x,jcell+y,0);
      if (k<0) continue;
      if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
      if(g_str[k]!=cu_memset){
        int_t fend=g_end[k];
        for(int_t j=g_str[k];j<fend;j++){
          Real xj,yj,mj,hj,rhoj,tdist,ptypej;

          xj=P1[j].x;
          yj=P1[j].y;
          mj=P1[j].m;
          hj=P1[j].h;
          rhoj=P1[j].rho;
          ptypej=P1[j].p_type;
          tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;

          //search_range=k_search_kappa*fmax(hi,hj);

          if(tdist<search_range){
            Real twij;

            twij=calc_kernel_wij(tmp_A,hi,tdist);
            tmp_R+=twij;
            tmp_RR+=hj*mj/rhoj*twij;
            tmp_flt+=mj/rhoj*twij;
          }
        }
      }
    }
  }
  //if (tcount%k_freq_h_reset==0)P1[i].h=1.5*(sqrt(1/tmp_R));
  //else P1[i].h=tmp_RR/tmp_flt;
  P1[i].h=tmp_RR/tmp_flt;
}

////////////////////////////////////////////////////////////////////////////////////////////
// __global__ void KERNEL_clc_correction_preLaplacian(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
// {
// 	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part2) return;
// 	if(P1[i].i_type>i_type_crt) return;
//
// 	int_t icell,jcell;
// 	Real xi, yi;
// 	Real search_range,hi,tmp_Ai;
// 	Real tA111, tA112, tA122, tA211, tA212, tA222;
// 	Real tflt_s;
//
// 	hi=P1[i].h;
// 	tmp_Ai=calc_tmpA(hi);
// 	search_range=k_search_kappa*hi;	// search range
//
// 	xi=P1[i].x;
// 	yi=P1[i].y;
//
// 	tA111=tA112=tA122=tA211=tA212=tA222=0.0;
//
// 	// calculate I,J,K in cell
// if((k_x_max==k_x_min)){icell=0;}
// else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
// if((k_y_max==k_y_min)){jcell=0;}
// else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
// 	// out-of-range handling
// 	if(icell<0) icell=0;	if(jcell<0) jcell=0;
//
// 	for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
// 		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
// 			//int_t k=(icell+x)+k_NI*(jcell+y);
// 			int_t k=idx_cell(icell+x,jcell+y,0);
//
// 			if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
// 			if(g_str[k]!=cu_memset){
// 				int_t fend=g_end[k];
// 				for(int_t j=g_str[k];j<fend;j++){
// 					if(P1[j].p_type<1000){
//
// 					Real xj,yj,tdist;
//
// 					xj=P1[j].x;
// 					yj=P1[j].y;
//
// 					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
//
// 					if(tdist>0&&tdist<search_range){
// 						Real twij,mj,rhoj,hj,rt,mtd;
// 						Real tdwij, tdwx, tdwy;
// 						Real tmprho;
//
// 						hj=P1[j].h;
//
// 						// twij=calc_kernel_wij(tmp_A,hi,tdist);
// 						// tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
//
// 						// for multi resolution (KDH)
// 						Real tmp_Aj=calc_tmpA(hj);
// 						twij=0.5*(calc_kernel_wij(tmp_Ai,hi,tdist)+calc_kernel_wij(tmp_Aj,hj,tdist));
// 						tdwij=0.5*(calc_kernel_dwij(tmp_Ai,hi,tdist)+calc_kernel_dwij(tmp_Aj,hj,tdist));
//
// 						apply_gradient_correction_2D(P3[i].Cm,twij,tdwx,tdwy,&tdwx,&tdwy);
//
// 						tdwx=tdwij*(xi-xj)/tdist;
// 						tdwy=tdwij*(yi-yj)/tdist;
//
// 						mj=P1[j].m;
// 						rhoj=P1[j].rho;
// 						tmprho=mj/rhoj;
//
// 						tA111+=(xi-xj)*(xi-xj)*tdwx*tmprho;
// 						tA112+=(xi-xj)*(yi-yj)*tdwx*tmprho;
// 						tA122+=(yi-yj)*(yi-yj)*tdwx*tmprho;
// 						tA211+=(xi-xj)*(xi-xj)*tdwy*tmprho;
// 						tA212+=(xi-xj)*(yi-yj)*tdwy*tmprho;
// 						tA222+=(yi-yj)*(yi-yj)*tdwy*tmprho;
//
// 					}
// 				}
// 			}
// 			}
// 		}
// 	}
// 	P3[i].A111=tA111;
// 	P3[i].A112=tA112;
// 	P3[i].A122=tA122;
// 	P3[i].A211=tA211;
// 	P3[i].A212=tA212;
// 	P3[i].A222=tA222;
//
// }

////////////////////////////////////////////////////////////////////////////////////////////
// __global__ void KERNEL_clc_correction_Laplacian(int_t*g_str,int_t*g_end,part1*P1,part3*P3)
// {
// 	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part2) return;
// 	if(P1[i].i_type>i_type_crt) return;
//
// 	int_t icell,jcell;
// 	Real xi, yi;
// 	Real search_range,hi,tmp_Ai;
// 	Real tB11, tB12, tB13, tB21, tB22, tB23, tB31, tB32, tB33;
// 	Real tA111, tA112, tA122, tA211, tA212, tA222;
// 	Real tflt_s;
//
// 	hi=P1[i].h;
// 	tmp_Ai=calc_tmpA(hi);
// 	search_range=k_search_kappa*hi;	// search range
//
// 	xi=P1[i].x;
// 	yi=P1[i].y;
//
// 	tB11=tB12=tB13=tB21=tB22=tB23=tB31=tB32=tB33=0.0;
//
// 	tA111=P3[i].A111;
// 	tA112=P3[i].A112;
// 	tA122=P3[i].A122;
// 	tA211=P3[i].A211;
// 	tA212=P3[i].A212;
// 	tA222=P3[i].A222;
//
// 	// calculate I,J,K in cell
// if((k_x_max==k_x_min)){icell=0;}
// else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
// if((k_y_max==k_y_min)){jcell=0;}
// else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
// 	// out-of-range handling
// 	if(icell<0) icell=0;	if(jcell<0) jcell=0;
//
// 	for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
// 		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
// 			//int_t k=(icell+x)+k_NI*(jcell+y);
// 			int_t k=idx_cell(icell+x,jcell+y,0);
//
// 			if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
// 			if(g_str[k]!=cu_memset){
// 				int_t fend=g_end[k];
// 				for(int_t j=g_str[k];j<fend;j++){
// 					if(P1[j].p_type<1000){
//
// 					Real xj,yj,tdist;
// 					xj=P1[j].x;
// 					yj=P1[j].y;
//
// 					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
// 					if(tdist>0&&tdist<search_range){
// 						Real twij,mj,rhoj,rt,mtd;
// 						Real tdwij, tdwx, tdwy;
// 						Real tmprho;
// 						Real e1,e2,r1,r2;
//
// 						// twij=calc_kernel_wij(tmp_A,hi,tdist);
// 						// tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
//
// 						// for multi resolution
// 						Real hj=P1[j].h;
// 						Real tmp_Aj=calc_tmpA(hj);
// 						twij=0.5*(calc_kernel_wij(tmp_Ai,hi,tdist)+calc_kernel_wij(tmp_Aj,hj,tdist));
// 						tdwij=0.5*(calc_kernel_dwij(tmp_Ai,hi,tdist)+calc_kernel_dwij(tmp_Aj,hj,tdist));
//
// 						tdwx=tdwij*(xi-xj)/tdist;
// 						tdwy=tdwij*(yi-yj)/tdist;
//
// 						e1=(xi-xj)/tdist;
// 						e2=(yi-yj)/tdist;
//
// 						r1=(xi-xj);
// 						r2=(yi-yj);
//
// 						mj=P1[j].m;
// 						rhoj=P1[j].rho;
// 						tmprho=mj/rhoj;
//
// 						tB11+=tmprho*(tA111*e1+tA211*e2+r1*e1)*(e1*tdwx);
// 						tB12+=tmprho*(tA111*e1+tA211*e2+r1*e1)*(e1*tdwy+e2*tdwx);
// 						tB13+=tmprho*(tA111*e1+tA211*e2+r1*e1)*(e2*tdwy);
// 						tB21+=tmprho*(tA112*e1+tA212*e2+r1*e2)*(e1*tdwx);
// 						tB22+=tmprho*(tA112*e1+tA212*e2+r1*e2)*(e1*tdwy+e2*tdwx);
// 						tB23+=tmprho*(tA112*e1+tA212*e2+r1*e2)*(e2*tdwy);
// 						tB31+=tmprho*(tA122*e1+tA222*e2+r2*e2)*(e1*tdwx);
// 						tB32+=tmprho*(tA122*e1+tA222*e2+r2*e2)*(e1*tdwy+e2*tdwx);
// 						tB33+=tmprho*(tA122*e1+tA222*e2+r2*e2)*(e2*tdwy);
// 					}
// 				}
// 			}
// 			}
// 		}
// 	}
//
// 	Real det_B;
// 	Real tL11, tL12, tL22;
// 	tL11=tL12=tL22=0.0;
//
// 	det_B=tB11*(tB22*tB33-tB23*tB32)-tB12*(tB21*tB33-tB23*tB31)+tB13*(tB21*tB32-tB22*tB31);
//
// 	if(det_B>1e-10){
// 		tL11=-1.0/det_B*(tB22*tB33-tB23*tB32+tB12*tB23-tB13*tB22);
// 		tL12=-1.0/det_B*(-tB21*tB33+tB23*tB31-tB11*tB23+tB13*tB21);
// 		tL22=-1.0/det_B*(tB21*tB32-tB22*tB31+tB11*tB22-tB12*tB21);
// 	}
// 	else{
// 	 	tL11=1.0;
// 	 	tL12=0.0;
// 	 	tL22=1.0;
// 	}
//
//
// 	P3[i].L11=tL11;
// 	P3[i].L12=tL12;
// 	P3[i].L21=tL12;
// 	P3[i].L22=tL22;
//
// }
