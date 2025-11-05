////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_PPE3D(Real tdt, int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;

	int_t icell,jcell,kcell;
	Real xi,yi,zi,uxi,uyi,uzi,rhoi,tempi;
	Real drhostar;
	Real search_range,hi,tmp_A;
	Real flt;
	Real bi,bix,biy,biz;
	Real Aij;
	Real AijPij;

	hi=Pressure_length*P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;
	rhoi=P1[i].rho;
	tempi=P1[i].temp;

	Real gc11, gc12, gc22, gc23, gc31, gc33;
	gc11=P3[i].inv_cm_xx, gc12=P3[i].inv_cm_xy;
	gc22=P3[i].inv_cm_yy, gc23=P3[i].inv_cm_yz;
	gc31=P3[i].inv_cm_zx, gc33=P3[i].inv_cm_zz;

	Real betai,Fb_i, divFb;
	betai=thermal_expansion(tempi,1);
 	Fb_i= -Gravitational_CONST*betai*(tempi-tempb);
	divFb=0.0;

	drhostar=0.0;
	bi = Aij = AijPij = 0.0;
	bix = biy = biz = 0.0;
	flt = 0.0;

	Real unity=0.0;
	Real uijsum=0.0;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t ncell=P1[i].ncell;
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
						int itype;

						itype=P1[j].i_type;
						ptypej=P1[j].p_type;
						mj=P1[j].m;
						tempj=P1[j].temp;

						if(P1[j].p_type<1000){
							if(itype!=4){
								xj=P1[j].x;
								yj=P1[j].y;
								zj=P1[j].z;

								uxj=P1[j].ux;
								uyj=P1[j].uy;
								uzj=P1[j].uz;

								rhoj=P1[j].rho;
								rho_ref_j=P2[j].rho_ref;
								presj=P1[j].pres;
								hj=P1[j].h;

								search_range=k_search_kappa*fmax(hi,hj);   // search range for multi-resolution (KDH)

								tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;

								if(tdist<search_range){
									//if(i==135235) P1[j].PPE2=1.0;
									Real twij=calc_kernel_wij(tmp_A,hi,tdist);
									Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);    // 원본
									Real tdwx=tdwij*(xi-xj)/tdist;
									Real tdwy=tdwij*(yi-yj)/tdist;
									Real tdwz=tdwij*(zi-zj)/tdist;

									Real tdwxc, tdwyc, tdwzc;      // gradient correction
									tdwxc=gc11*tdwx+gc12*tdwy+gc31*tdwz;
									tdwyc=gc12*tdwx+gc22*tdwy+gc23*tdwz;
									tdwzc=gc31*tdwx+gc23*tdwy+gc33*tdwz;

									// Real hij=0.5*(hi+hj);
									// Real tmp_Aij=calc_tmpA(hij);
									// Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
									// Real tdwij=calc_kernel_dwij(tmp_Aij,hij,tdist);          // For multi-resolution (KDH)

									// apply_gradient_correction_3D(P3[i].Cm,twij,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);

									// bix += mj/rhoj*(uxi-uxj)*tdwx*(2.0*rhoi*rhoj/(rhoi+rhoj));
									// biy += mj/rhoj*(uyi-uyj)*tdwy*(2.0*rhoi*rhoj/(rhoi+rhoj));
									// biz += mj/rhoj*(uzi-uzj)*tdwz*(2.0*rhoi*rhoj/(rhoi+rhoj));

									// bix += mj/rhoj*(uxj-uxi)*tdwx*rhoi;
									// biy += mj/rhoj*(uyj-uyi)*tdwy*rhoi;
									// biz += mj/rhoj*(uzj-uzi)*tdwz*rhoi;

									// bix += mj/rhoj*(uxj-uxi)*tdwx*(2.0*rhoi*rhoj/(rhoi+rhoj));
									// biy += mj/rhoj*(uyj-uyi)*tdwy*(2.0*rhoi*rhoj/(rhoi+rhoj));
									// biz += mj/rhoj*(uzj-uzi)*tdwz*(2.0*rhoi*rhoj/(rhoi+rhoj));

									bix += mj/rhoj*(uxj-uxi)*tdwxc*(2.0*rhoi*rhoj/(rhoi+rhoj));       // gradient correction
									biy += mj/rhoj*(uyj-uyi)*tdwyc*(2.0*rhoi*rhoj/(rhoi+rhoj));
									biz += mj/rhoj*(uzj-uzi)*tdwzc*(2.0*rhoi*rhoj/(rhoi+rhoj));

									// Real betaj=thermal_expansion(tempj,ptypej);
									// Real Fb_j= -Gravitational_CONST*betaj*(tempj-tempb);
									// divFb += mj/rhoj*(Fb_i-Fb_j)*tdwz;

									flt += (mj/rhoj)*twij;
									unity += twij;

									// Real c = 2.0*mj/rhoj*((xi-xj)*tdwx+(yi-yj)*tdwy+(zi-zj)*tdwz)/tdist/tdist;
									Real c = 2.0*mj/rhoj*((xi-xj)*tdwxc+(yi-yj)*tdwyc+(zi-zj)*tdwzc)/tdist/tdist;    // gradient correction

									Aij += c;
									AijPij += c*(presj-P1[j].pres0);
								}
							}
						}
					}
				}
			}
		}
	}
	Real weight=1.0;
	//if(P1[i].buffer_type==-1||P1[i].buffer_type==1) weight=0.05*xi;

	bi = weight*(bix+biy+biz+divFb)/tdt;

	// bi = -rhoi*(1.0-drhostar/flt)/tdt/tdt;
	//P1[i].PPE1 = rhoi*(bix+biy+biz)/tdt;
	//P1[i].PPE2 = -rhoi*(1.0-drhostar/flt)/tdt/tdt;  // YHS

	//P1[i].PPE1 = unity;
	P1[i].pres=P1[i].pres0+(bi+AijPij)/Aij;
	//if(xi>4.5 && xi<15.5 && yi>4 && yi<6 && zi>7 && zi<8) printf("xi=%f bix=%f biy=%f biz=%f Aij=%f AijPij=%f\n",xi,bix,biy,biz,Aij,AijPij);

	//if(isnan(P1[i].pres)==1) printf("xi=%f bi=%f bix=%f biy=%f biz=%f Aij=%f AijPij=%f\n",xi,bi,bix,biy,biz,Aij,AijPij);

	// P1[i].fv_x=bi;
	// P1[i].fv_y=Aij;
	// P1[i].fv_z=AijPij;

}

__global__ void KERNEL_pressure_synchronize(part1*P1,part1*TP1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;

	P1[i].pres=TP1[i].pres;
	printf("P1[i].pres=%f %f\n",P1[i].pres, TP1[i].pres);

}

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_advection_force3D(int_t inout,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3,int_t tcount)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type==1000) return;		// Immersed Boundary Method
	if(P1[i].p_type<=0) return;		// Immersed Boundary Method
	// if(P1[i].buffer_type>1) return;      // inlet

	int_t ptypei;
	int_t buffer_type_i;
	int_t icell,jcell,kcell;
	Real xi,yi,zi,uxi,uyi,uzi,kci,eta;
	Real hi,mi,rhoi,tempi,visi,betai;
	Real diffi,concni;
	Real search_range,tmp_A,tmp_Rc,tmp_Rd,tmp_Rd2;
	Real tmpx,tmpy,tmpz;
	Real eulerx, eulery, eulerz, eulert;
	Real cpi;
	Real tmp_temp,tmp_flt;

	Real fvx,fvy,fvz,fex,fey,fez,fbs;

	ptypei=P1[i].p_type;
	buffer_type_i=P1[i].buffer_type;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;
	hi=P1[i].h;
	tempi=P1[i].temp;
	mi=P1[i].m;
	rhoi=P1[i].rho;

	tmp_A=calc_tmpA(hi);

	if(k_con_solve){
		// kci=conductivity(tempi,ptypei);
		// cpi=specific_heat(tempi,ptypei);
		kci=Solid_k*(ptypei==0||ptypei==1000)+Air_k*(ptypei==1||ptypei==-1);
		cpi=Solid_cp*(ptypei==0||ptypei==1000)+Air_cp*(ptypei==1||ptypei==-1);
		eta=0.001*hi;
	}

	//printf("kci, cpi=%f %f\n",kci,cpi);
	if(k_concn_solve){
		concni=P1[i].concn;
		diffi=diffusion_coefficient(tempi,ptypei);
	}

	visi=viscosity2(tempi,ptypei)+P1[i].vis_t*rhoi;
	betai=thermal_expansion(tempi,ptypei);

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpx=tmpy=tmpz=0.0;
	tmp_Rc=0.0;
	tmp_Rd=0.0;
	tmp_Rd2=0.0;
	eulerx=eulery=eulerz=eulert=0.0;

	fvx=fvy=fvz=0.0;
	fex=fey=fez=0.0;
	tmp_temp=tmp_flt=0.0;

	// For artificial viscosity
	int_t fva_solve=1;
	Real Alpha=0.05;
	Real Beta=0.05;
	Real soundspeed=30.0;
	Real Turb_dT=0.0;

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
						Real xj,yj,zj,hj,ptypej,tdist;
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						hj=P1[j].h;
						ptypej=P1[j].p_type;

						if(P1[j].p_type<1000){

							tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
							search_range=k_search_kappa*fmax(hi,hj);	// search range for multi-resoltuion(KDH)

							if(tdist<search_range){
								int_t ptypej;
								Real tdwx,tdwy,tdwz,uxj,uyj,uzj,mj,tempj,rhoj,kcj,sum_con_H,diffj,concnj,tmprd,tmprd1,tmprd2;

								Real twij=calc_kernel_wij(tmp_A,hi,tdist);
								Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);

								// for multi resolution2 (KDH)
								// Real hij=0.5*(hi+hj);
								// Real tmp_Aij=calc_tmpA(hij);
								// Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
								// Real tdwij=calc_kernel_dwij(tmp_Aij,hij,tdist);

								tdwx=tdwij*(xi-xj)/tdist;
								tdwy=tdwij*(yi-yj)/tdist;
								tdwz=tdwij*(zi-zj)/tdist;

								apply_gradient_correction_3D(P3[i].Cm,twij,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);

								ptypej=P1[j].p_type;
								uxj=P1[j].ux;
								uyj=P1[j].uy;
								uzj=P1[j].uz;
								mj=P1[j].m;
								tempj=P1[j].temp;
								rhoj=P1[j].rho;
								hj=P1[j].h;

								if(k_fv_solve==Cleary){
									Real visj,C_v;
									visj=viscosity2(tempj,ptypej)+P1[j].vis_t*rhoj;        // mu (kg/(ms))
									//C_v=4*(mj/(rhoi*rhoj))*((visi*visj)/(visi+visj))*((xi-xj)*tdwx+(yi-yj)*tdwy+(zi-zj)*tdwz)/tdist/tdist;
									C_v=(xi-xj)*tdwx+(yi-yj)*tdwy+(zi-zj)*tdwz;
									C_v*=(visi*visj)/(visi+visj);
									C_v*=4*(mj/(rhoi*rhoj));
									C_v/=tdist;
									C_v/=tdist;

									tmpx+=C_v*(uxi-uxj);
									tmpy+=C_v*(uyi-uyj);
									tmpz+=C_v*(uzi-uzj);

									fvx+=C_v*(uxi-uxj);
									fvy+=C_v*(uyi-uyj);
									fvz+=C_v*(uzi-uzj);

								}else if(k_fv_solve==Monaghan){

									Real uij_xij=(uxi-uxj)*(xi-xj)+(uyi-uyj)*(yi-yj)+(uzi-uzj)*(zi-zj);
									Real phi_ij,P_ij;
									Real visj=viscosity2(tempj,ptypej)+P1[j].vis_t*rhoj;

									// Monaghan type
									phi_ij=uij_xij;
									phi_ij/=(tdist*tdist+0.01*hi*hi);
									// P_ij=-2.0*(k_dim+2.0)*phi_ij*mj/rhoj;
									P_ij=-2.0*(k_dim+2.0)*phi_ij*mj/rhoj;
									Real nui,C_v;
									C_v=2.0*(visi*visj)/(visi+visj-1e-20);
									nui = C_v/rhoi;
									P_ij*= nui;

									if(P1[i].p_type<=0 || P1[i].p_type>=1000)	P_ij=0.0;
									if(P1[j].p_type<=0 || P1[i].p_type>=1000)	P_ij=0.0;

									tmpx+=-(P_ij)*tdwx;
									tmpy+=-(P_ij)*tdwy;
									tmpz+=-(P_ij)*tdwz;

									fvx+=-(P_ij)*tdwx;
									fvy+=-(P_ij)*tdwy;
									fvz+=-(P_ij)*tdwz;

								}
								if(P1[i].eli<1.0){

									eulerx += (uxi*(uxj-uxi)*tdwx+uyi*(uxj-uxi)*tdwy+uzi*(uxj-uxi)*tdwz)*mj/rhoj*(1-P1[i].eli);
									eulery += (uxi*(uyj-uyi)*tdwx+uyi*(uyj-uyi)*tdwy+uzi*(uyj-uyi)*tdwz)*mj/rhoj*(1-P1[i].eli);
									eulerz += (uxi*(uzj-uzi)*tdwx+uyi*(uzj-uzi)*tdwy+uzi*(uzj-uzi)*tdwz)*mj/rhoj*(1-P1[i].eli);

									if(k_con_solve){
										//eulert -= (uxi*(tempj-tempi)*tdwx+uyi*(tempj-tempi)*tdwy+uzi*(tempj-tempi)*tdwz)*mj/rhoj*(1-P1[i].eli);
										eulert -= 2.0*mj/(rhoi+rhoj)*(tempj-tempi)*((uxi-uxj)*tdwx+(uyi-uyj)*tdwy+(uzi-uzj)*tdwz)*(1-P1[i].eli);
									}
									if(k_concn_solve){
										concnj=P1[j].concn;
										tmp_Rd-= (uxi*(concnj-concni)*tdwx+uyi*(concnj-concni)*tdwy+uzi*(concnj-concni)*tdwz)*mj/rhoj*(1-P1[i].eli);
									}
								}
								if(k_con_solve){
									//kcj=conductivity(tempj, ptypej);
									kcj=Solid_k*(ptypej==0||ptypej==1000)+Air_k*(ptypej==1||ptypej==-1);
									Real rhoj_temp, mj_temp;
									if(ptypej==1||ptypej==-1) rhoj_temp=rhoj;
									else rhoj_temp=Solid_rho;
									mj_temp=mj*rhoj_temp/rhoj;
									sum_con_H=4.0*mj_temp*kcj*kci*(tempi-tempj)*tdwij;
									//sum_con_H/=(tdist+1e-10)*rhoi*rhoj*(kci+kcj);
									sum_con_H/=(tdist+1e-10)*rhoi*rhoj_temp*(kci+kcj);
									tmp_Rc+=sum_con_H;
								}
								//if(k_turbulence_model==SPS){
									// tmp_Rc+=0.5*mj*(P1[i].vis_t+P1[j].vis_t)/(rhoi*rhoj)*((uxi-uxj)*tdwx+(uyi-uyj)*tdwy+(uzi-uzj)*tdwz);
									// Turb_dT+=0.5*mj*(P1[i].vis_t+P1[j].vis_t)/(rhoi*rhoj)*((uxi-uxj)*tdwx+(uyi-uyj)*tdwy+(uzi-uzj)*tdwz);
									//tmp_Rc += mj/(rhoi*rhoj)*((P1[i].vis_t+P1[j].vis_t)/Pr_t)*(tempi-tempj)*tdwij/(tdist+1e-20);
								//}
								if(k_concn_solve){
									concnj=P1[j].concn;
									diffj=diffusion_coefficient(tempj,ptypej);
									tmprd=mi*(diffi*rhoi+diffj*rhoj);
									tmprd*=((xi-xj)*tdwx+(yi-yj)*tdwy+(zi-zj)*tdwz);
									tmprd/=(rhoi*rhoj*(tdist*tdist+0.01*hi*hi));
									tmprd1=tmprd*(concni-concnj)*(ptypei==ptypej);

									tmp_Rd+=tmprd1;
								}
								if(k_boussinesq_solve){
									tmp_temp+=(mj/rhoj)*tempj*twij;
									tmp_flt+=(mj/rhoj)*twij;
								}
								if(fva_solve){    // Artificial viscosity solve
									Real uij_xij=(uxi-uxj)*(xi-xj)+(uyi-uyj)*(yi-yj)+(uzi-uzj)*(zi-zj);
									if(uij_xij<0){
										Real h_ij,phi_ij,P_ij;
										//
										h_ij=(hi+hj)*0.5;
										phi_ij=h_ij*uij_xij;
										phi_ij/=(tdist*tdist+0.01*h_ij*h_ij);
										P_ij=mi*phi_ij*(-Alpha*soundspeed+Beta*phi_ij);
										P_ij/=(rhoi+rhoj);
										P_ij*=0.5;
										//P_ij=mi*(-Alpha*k_soundspeed*phi_ij+Beta*phi_ij*phi_ij)/(rhoi+rhoj)*0.5;
										tmpx+=-(P_ij)*tdwx;
										tmpy+=-(P_ij)*tdwy;
										tmpz+=-(P_ij)*tdwz;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	// z-directional gravitational force
	Real fg=0.0;
	if(k_fg_solve) tmpz+=-Gravitational_CONST;
	if((k_boussinesq_solve)&&(ptypei>0)&&(ptypei<1000)) tmpz=betai*Gravitational_CONST*(tmp_temp/tmp_flt-tempb);

	//P1[i].PPE1=betai*Gravitational_CONST*(tmp_temp/tmp_flt-tempb);
	P3[i].ftotalx=tmpx-eulerx;
	P3[i].ftotaly=tmpy-eulery;
	P3[i].ftotalz=tmpz-eulerz;

	// P1[i].fe_x=-eulerx;
	// P1[i].fe_y=-eulery;
	// P1[i].fe_z=-eulerz;
	//
	// P1[i].fv_x=fvx;
	// P1[i].fv_y=fvy;
	// P1[i].fv_z=fvz;

	if(k_con_solve){
		//if (enthalpy_eqn) P3[i].denthalpy=tmp_Rc*(ptypei!=-1);
		if(enthalpy_eqn) P3[i].denthalpy=tmp_Rc*(ptypei!=-1);
		else{
			P3[i].dtemp+=tmp_Rc/cpi-eulert;
			//P1[i].PPE1=Turb_dT;
			//P1[i].PPE2=-eulert;
			//if (tmp_Rc>1e-6) printf("dtemp=%f\n",P3[i].dtemp);
		}
	}
	// P3[i].fsx=tmpsx;
	// P3[i].fsy=tmpsy;
	// P3[i].fsz=tmpsz;

	P3[i].ftotal=sqrt(P3[i].ftotalx*P3[i].ftotalx+P3[i].ftotaly*P3[i].ftotaly+P3[i].ftotalz*P3[i].ftotalz);
	if(k_concn_solve) P3[i].dconcn=tmp_Rd;

	// if(isnan(P3[i].ftotalx)) printf("ftotalx:%f eulerx=%f\n",P3[i].ftotalx, eulerx);

	return;
}
////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_pressureforce3D(int_t inout,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3){
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	//if(k_open_boundary>0 && P1[i].buffer_type>0) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;		// Immersed Boundary Method
	if(P1[i].buffer_type!=0)	return;

	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real pi,hi,rhoi;
	Real search_range,tmp_A;
	Real tmpx,tmpy,tmpz;
	int_t ptypei;

	Real fpx,fpy,fpz;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	hi=P1[i].h;
	pi=P1[i].pres-P1[i].pres0;
	rhoi=P1[i].rho;
	ptypei=P1[i].p_type;

	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpx=tmpy=tmpz=0.0;
	fpx=fpy=fpz=0.0;

	// Real pj_max = -FLT_MAX;  // 초기값 설정
  // Real pj_min = FLT_MAX;
	// Real pi_sum = 0.0;
	// Real flt_pi = 0.0;

  // **1단계: pj_max와 pj_min 계산**
  // for(int_t z=-P1[i].ncell; z<=P1[i].ncell; z++) {
  //   for(int_t y=-P1[i].ncell; y<=P1[i].ncell; y++) {
  //     for(int_t x=-P1[i].ncell; x<=P1[i].ncell; x++) {
  //       int_t k=idx_cell(icell+x, jcell+y, kcell+z);
  //       // 셀 범위를 벗어나는 경우 무시
  //       if(((icell+x)<0) || ((icell+x)>=k_NI) || ((jcell+y)<0) || ((jcell+y)>=k_NJ) || ((kcell+z)<0) || ((kcell+z)>=k_NK)) continue;
  //       if(g_str[k]!=cu_memset){
  //         int_t fend=g_end[k];
  //         for(int_t j=g_str[k]; j<fend; j++) {
  //           // 입자 j에 대한 물리량 가져오기
  //           Real pj=P1[j].pres;
	// 					Real xj,yj,zj,mj,hj,rhoj,tdist;
	// 					xj=P1[j].x;
	// 					yj=P1[j].y;
	// 					zj=P1[j].z;
	// 					hj=P1[j].h;
	// 					rhoj=P1[j].rho;
	// 					mj=P1[j].m;
	//
	// 					if(P1[j].p_type<1000){
	// 					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
	//
	// 						if(tdist<search_range){
	//
	// 							Real twij=calc_kernel_wij(tmp_A,hi,tdist);
	// 							Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
	//
	// 							if(pj>pj_max) pj_max=pj;
  //           		if(pj<pj_min) pj_min=pj;
	// 							pi_sum+=pj*(mj/rhoj)*twij;
	// 							flt_pi+=(mj/rhoj)*twij;
	// 						}
	// 					}
	// 				}
  //     	}
  //   	}
  // 	}
  // }

	// **2단계: 압력 힘 계산 및 적용**
	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,hj,tdist;
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						hj=P1[j].h;
						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
						search_range=k_search_kappa*fmax(hi,hj);

						if(P1[j].p_type<1000){
							if(tdist<search_range){
								Real tdwx,tdwy,tdwz,mj,rhoj,pj;
								int_t ptypej;

								ptypej=P1[j].p_type;

								// Real twij=calc_kernel_wij(tmp_A,hi,tdist);
								// Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);

								// for multi resolution2 (KDH)
								Real hij=0.5*(hi+hj);
								Real tmp_Aij=calc_tmpA(hij);
								Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
								Real tdwij=calc_kernel_dwij(tmp_Aij,hij,tdist);

								tdwx=tdwij*(xi-xj)/tdist;
								tdwy=tdwij*(yi-yj)/tdist;
								tdwz=tdwij*(zi-zj)/tdist;

								if(k_kgc_solve>0){
									apply_gradient_correction_3D(P3[i].Cm,twij,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
								}

								mj=P1[j].m;
								rhoj=P1[j].rho;
								pj=P1[j].pres-P1[j].pres0;

								if(k_fp_solve){

									// p_stab 계산  (for pressure stabilization)
	               	// Real p_stab=min((pj_max-pj_min)/2,pi)-0.0;
									//Real p_stab=min((pj_max-pj_min)/2,pi_sum/flt_pi);

                	// **C_p 계산 및 힘 적용**
									// Real C_p=-mj/(rhoi*rhoj)*(pi*rhoj+pj*rhoi)/(rhoi+rhoj);   // two phase
									// Real C_p=-mj*(pi+pj)/(rhoi*rhoj);      // Standard Monaghan
									Real C_p=-2.0*mj/(rhoi*rhoj)*((pi+Pb)*rhoj+(pj+Pb)*rhoi)/(rhoi+rhoj);      // Monaghan variation
                	// Real C_p=-mj*(pj-pi)/(rhoi*rhoj)-mj*(2.0*p_stab)/(rhoi*rhoj)*(ptypei==ptypej);     // Pressure stabilization(KDH)
									tmpx+=C_p*tdwx;
									tmpy+=C_p*tdwy;
									tmpz+=C_p*tdwz;
								}
							}
						}
					}
				}
			}
		}
	}
	// P1[i].fp_x=tmpx;
	// P1[i].fp_y=tmpy;
	// P1[i].fp_z=tmpz;

	P3[i].fpx=tmpx;
	P3[i].fpy=tmpy;
	P3[i].fpz=tmpz;

	// P1[i].uy_f=tmpx+tmpy+tmpz;
}

__global__ void KERNEL_pressure_smoothing3D(int_t*g_str,int_t*g_end,part1*P1)
{
	int_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type!=1) return;
	if(P1[i].buffer_type!=0) return;

	int_t ptypei;
	int_t icell,jcell,kcell;
	Real xi,yi,zi;
	Real pi;														// pressure
	Real flt_si;
	Real tmpp;
	Real search_range,tmp_h,tmp_A;
	Real pc=0.0;

	tmp_h=P1[i].h;
	tmp_A=calc_tmpA(tmp_h);
	search_range=k_search_kappa*tmp_h;								// search range

	ptypei=P1[i].p_type;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	pi=P1[i].pres;

	flt_si=P1[i].flt_s;

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

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj));
						if(tdist<search_range){
							Real twij,pj,mj,rhoj,ptype_j;
							twij=calc_kernel_wij(tmp_A,tmp_h,tdist);
							mj=P1[j].m;
							rhoj=P1[j].rho;
							pj=P1[j].pres;
							ptype_j=P1[j].p_type;

							tmpp+=C_pres*mj/rhoj*(-pi+pj)*twij*(ptype_j==1);
						}
					}
				}
			}
		}
	}
	//printf("flt_si=%f \n", flt_si);
	pc=pi+tmpp*(ptypei>0);			// correct pressure

	P1[i].pres=pc;				// update pressure
}

__global__ void KERNEL_variable_smoothing_length3D(int_t*g_str,int_t*g_end,part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if (i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  if(P1[i].p_type<1) return;

  int_t icell,jcell,kcell;
  Real xi,yi,zi,hi;
  Real search_range,tmp_A,tmp_R,tmp_RR,tmp_flt;

  tmp_R=tmp_RR=tmp_flt=0.0;
  Real smooth_range=3.0;
  hi=smooth_range*P1[i].h;
  tmp_A=calc_tmpA(hi);
  search_range=k_search_kappa*hi;	// search range

  xi=P1[i].x;
  yi=P1[i].y;
  zi=P1[i].z;

  // calculate I,J,K in cell
  if((k_x_max==k_x_min)){icell=0;}
  else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  if((k_y_max==k_y_min)){jcell=0;}
  else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  if((k_z_max==k_z_min)){kcell=0;}
  else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
  // out-of-range handling
  if(icell<0) icell=0;	if(jcell<0) jcell=0;  if(kcell<0) kcell=0;

  int_t ncell=smooth_range*P1[i].ncell;

  for(int_t z=-ncell;z<=ncell;z++){
    for(int_t y=-ncell;y<=ncell;y++){
      for(int_t x=-ncell;x<=ncell;x++){
        // int_t k=(icell+x)+k_NI*(jcell+y);
        int_t k=idx_cell(icell+x,jcell+y,kcell+z);
        if (k<0) continue;
        if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
        if(g_str[k]!=cu_memset){
          int_t fend=g_end[k];
          for(int_t j=g_str[k];j<fend;j++){
            Real xj,yj,zj,hj,mj,rhoj,tdist;

            xj=P1[j].x;
            yj=P1[j].y;
            zj=P1[j].z;
            mj=P1[j].m;
            hj=P1[j].h;
            rhoj=P1[j].rho;

            tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;

            if(k_apr_solv) search_range=k_search_kappa*fmax(hi,hj);

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
  }
  //if(tcount%k_freq_h_reset==0)P1[i].h=1.5*(cbrt(1/tmp_R));
  // if(tcount%k_freq_h_reset==0)P1[i].h=P1[i].h_ref*cbrt(P1[i].m/P1[i].m_ref);
  // else P1[i].h=tmp_RR/tmp_flt;
  P1[i].h=tmp_RR/tmp_flt;
}
