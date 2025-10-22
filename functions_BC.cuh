////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_boundary3D(int_t*g_str,int_t*g_end,part1*P1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;

	int_t icell,jcell,kcell;
	uint_t p_type_i;
	p_type_i=P1[i].p_type;
	if(p_type_i>0) return;

	Real xi,yi,zi,uxi,uyi,uzi;
	Real search_range,hi,tmp_A;
	Real tmpx,tmpy,tmpz,flt;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;

	// calculate I,J,K in cell
if((k_x_max==k_x_min)){icell=0;}
else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
if((k_y_max==k_y_min)){jcell=0;}
else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
if((k_z_max==k_z_min)){kcell=0;}
else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tmpx=0.0;
	tmpy=0.0;
	tmpz=0.0;
	flt=1.0e-10;
	int_t ncell=P1[i].ncell;

	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist;
						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
						if(tdist<search_range){
							uint_t p_type_j;
							Real twij,uxj,uyj,uzj;
							twij=calc_kernel_wij(tmp_A,hi,tdist);
							p_type_j=P1[j].p_type;
							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;

							// noslip boundary condition
							if(k_noslip_bc==1){
								if((p_type_j>0)&(p_type_j!=MOVING)){//YHS
									tmpx+=uxj*twij;
									tmpy+=uyj*twij;
									tmpz+=uzj*twij;
									flt+=twij;
								}
							}
						}
					}
				}
			}
		}
	}

	// noslip boundary condition
	if(k_noslip_bc==1){
		if ((p_type_i == BOUNDARY) || (p_type_i==MOVING)){
			P1[i].ux=2*uxi*(p_type_i==MOVING)-tmpx/flt;
			P1[i].uy=2*uyi*(p_type_i==MOVING)-tmpy/flt;
			P1[i].uz=2*uzi*(p_type_i==MOVING)-tmpz/flt;
		}
	}
}

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_Neumann_boundary3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type==1) return;
	//if(P1[i].p_type==1 && P1[i].buffer_type!=-1 && P1[i].buffer_type!=1) return;
	//if(P1[i].p_type==1 && P1[i].buffer_type!=-1) return;
	//if(P1[i].p_type>0 && P1[i].buffer_type!=Outlet) return;

	int_t icell,jcell,kcell;
	Real xi, yi, zi;
	Real uxi, uyi, uzi;
	Real rhoi, presi;
	Real rho_ref_i;
	Real search_range, hi, tmp_A;
	Real tflt_s;
	Real tux, tuy, tuz;
	Real tpres, thpres, pres_g, grad_px, grad_py, grad_pz, aix, aiy, aiz;
	Real search_coeff=2.0;      // Search increment factor (KDH)
	hi=search_coeff*P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;
	rhoi=P1[i].rho;
	presi=P1[i].pres;

	grad_px=grad_py=grad_pz=0.0;
	tpres=thpres=tflt_s=pres_g=0.0;
	tux=tuy=tuz=0.0;

	Real xgnode_i, ygnode_i, zgnode_i;

	if(P1[i].buffer_type==1 || P1[i].buffer_type==-1) xgnode_i=2.0*L2-xi;
	else xgnode_i=xi;
	ygnode_i=yi;
	zgnode_i=zi;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t flag=0;

	for(int_t z=-1;z<=1;z++){
		for(int_t y=-1;y<=1;y++){
			for(int_t x=-1;x<=1;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,presj;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						presj=P1[j].pres;

						Real dist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						if(dist<0.1){
							pres_g=presj;      // ghost node pressure
							flag=1;
						}
					}
				}
			}
		}
	}

	Real ncell=P1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,uxj,uyj,uzj,tdist,mj,rhoj,presj;
						Real ajx,ajy,ajz;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						ajx=P3[j].ftotalx;
						ajy=P3[j].ftotaly;
						ajz=P3[j].ftotalz;
						mj=P1[j].m;
						rhoj=P1[j].rho;
						presj=P1[j].pres;

						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						if(P1[j].p_type<1000){

							if(tdist<search_range){
								// if(tdist<search_range){
								Real twij, tdwij, tdwx, tdwy, tdwz;
								int_t ptype_j=P1[j].p_type;
								int_t buffer_type_j=P1[j].buffer_type;

								twij=calc_kernel_wij(tmp_A,hi,tdist);
								tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
								tdwx=(xgnode_i-xj)/tdist * tdwij;
								tdwy=(ygnode_i-yj)/tdist * tdwij;
								tdwz=(zgnode_i-zj)/tdist * tdwij;

								apply_gradient_correction_3D(P3[i].Cm,twij,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);

								Real neighbor = (ptype_j==1)*(buffer_type_j==0);

								tpres+=mj/rhoj*presj*twij*neighbor;
								thpres+=((P3[j].ftotalx)*(xj-xi)+(P3[j].ftotaly)*(yj-yi)+(P3[j].ftotalz)*(zj-zi))*mj/rhoj*twij*(ptype_j>=1);
								tux+=mj/rhoj*uxj*twij*neighbor;
								tuy+=mj/rhoj*uyj*twij*neighbor;
								tuz+=mj/rhoj*uzj*twij*neighbor;
								tflt_s+=mj/rhoj*twij*neighbor;

								grad_px += mj/rhoj*(presj-pres_g)*tdwx*neighbor;
								grad_py += mj/rhoj*(presj-pres_g)*tdwy*neighbor;
								grad_pz += mj/rhoj*(presj-pres_g)*tdwz*neighbor;
							}
						}
					}
				}
			}
		}
	}
	if(tflt_s<1e-6){
		P1[i].pres=0.0;
	}
	else{
		P1[i].pres = (tpres)/tflt_s;
		if(P1[i].p_type==-1){
			P1[i].ux=tux/tflt_s;
			P1[i].uy=tuy/tflt_s;
			P1[i].uz=tuz/tflt_s;
		}
		//if(P1[i].p_type==0) P1[i].pres=0.0;
		if(P1[i].buffer_type==1 || P1[i].buffer_type==-1) P1[i].pres = pres_g+(xi-xgnode_i)*grad_px+(yi-ygnode_i)*grad_py+(zi-zgnode_i)*grad_pz;
	}
}

// __global__ void KERNEL_Neumann_boundary3D_YHS(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
// {
// 	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(i>=k_num_part3) return;
// 	if((P1[i].p_type>0)&&(P1[i].p_type!=MOVING))	return;		// Immersed Boundary Method
//
//
// 	int_t icell,jcell,kcell;
// 	Real xi, yi, zi;
// 	Real uxi, uyi, uzi;
// 	Real rhoi, presi;
// 	Real rho_ref_i;
// 	Real search_range, tmp_h, tmp_A;
// 	Real tflt_s;
// 	Real tux, tuy, tuz;
// 	Real tpres, thpresx, thpresy, thpresz, ttemp, aix, aiy, aiz;
// 	Real t1, t2;
// 	tmp_h=P1[i].h;
// 	tmp_A=calc_tmpA(tmp_h);
// 	search_range=k_search_kappa*tmp_h;	// search range
//
// 	xi=P1[i].x;
// 	yi=P1[i].y;
// 	zi=P1[i].z;
// 	uxi=P1[i].ux;
// 	uyi=P1[i].uy;
// 	uzi=P1[i].uz;
// 	rhoi=P1[i].rho;
// 	presi=P1[i].pres;
//
// 	tpres=thpresx=thpresy=thpresz=ttemp=tflt_s=0.0;
// 	tux=tuy=tuz=0.0;
// 	t1=t2=0.0;
//
// 	// calculate I,J,K in cell
// 	if((k_x_max==k_x_min)){icell=0;}
// 	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
// 	if((k_y_max==k_y_min)){jcell=0;}
// 	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
// 	if((k_z_max==k_z_min)){kcell=0;}
// 	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
// 	// out-of-range handling
// 	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;
//
// 	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
// 		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
// 			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
// 				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
// 				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
// 				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
// 				if(g_str[k]!=cu_memset){
// 				int_t fend=g_end[k];
// 				for(int_t j=g_str[k];j<fend;j++){
// 					Real xj,yj,zj,uxj,uyj,uzj,tdist,mj,rhoj,presj,tempj;
// 					Real ajx,ajy,ajz;
//
// 					xj=P1[j].x;
// 					yj=P1[j].y;
// 					zj=P1[j].z;
// 					uxj=P1[j].ux;
// 					uyj=P1[j].uy;
// 					uzj=P1[j].uz;
// 					ajx=P3[j].ftotalx;
// 					ajy=P3[j].ftotaly;
// 					ajz=P3[j].ftotalz;
// 					mj=P1[j].m;
// 					rhoj=P1[j].rho;
// 					presj=P1[j].pres;
// 					tempj=P1[j].temp;
//
// 					tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
// 					if(P1[j].p_type<1000){
//
// 					if(tdist<search_range){
// 					// if(tdist<search_range){
// 						Real twij, tdwij, tdwx, tdwy, tdwz;
// 						int_t ptype_j=P1[j].p_type;
// 						int_t buffer_type_j=P1[j].buffer_type;
// 						twij=calc_kernel_wij(tmp_A,tmp_h,tdist);
// 						tdwij=calc_kernel_dwij(tmp_A,tmp_h,tdist);
//
// 						tpres+=presj*twij*(ptype_j>0);
// 						thpresx+=rhoj*(xi-xj)*twij*(ptype_j>0);
// 						thpresy+=rhoj*(yi-yj)*twij*(ptype_j>0);
// 						thpresz+=rhoj*(zi-zj)*twij*(ptype_j>0);
// 						tflt_s+=twij*(ptype_j>0);
// 						t1+=twij*(ptype_j==1);
// 						t2+=twij*(ptype_j==2);
//
// 						// ttemp+=mj/rhoj*tempj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
// 						// thpres+=((P3[j].ftotalx)*(xj-xi)+(P3[j].ftotaly)*(yj-yi)+(P3[j].ftotalz)*(zj-zi))*mj/rhoj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
// 						// tux+=mj/rhoj*uxj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
// 						// tuy+=mj/rhoj*uyj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
// 						// tuz+=mj/rhoj*uzj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
// 						// tflt_s+=mj/rhoj*twij*((ptype_j>0)&&(ptype_j!=MOVING));
//
// 					}
// 				}
// 			}
// 				}
// 			}
// 		}
// 	}
//
//
// 	if(tflt_s<1e-6){
// 		P1[i].rho=P2[i].rho_ref;
// 		P1[i].pres=0.0;
// 	}else{
// 		P1[i].pres = (tpres + (-0.0*thpresx -0.0*thpresy + (-Gravitational_CONST-0.0)*thpresz))/tflt_s;
// 		if(t2>1e-6&&t1<1e-6)		P1[i].pres = (0.0 + (-0.0*thpresx -0.0*thpresy + (-Gravitational_CONST-0.0)*thpresz))/tflt_s;
//
// 	if(P1[i].p_type==0)	P1[i].temp = ttemp/tflt_s;
// 	if(k_solver_type==Wcsph){
// 		Real rho0=fmax(1e-6,k_rho0_eos);
// 		double B = rho0*k_soundspeed*k_soundspeed/k_gamma;
// 		double K = P1[i].pres/B + 1.0;
// 		P1[i].rho = P2[i].rho_ref*pow(K, 1.0/k_gamma);
// 		}
// 	}
// 	if(P1[i].pres<0.0){
// 		P1[i].pres=0.0;
// }
// }

__global__ void KERNEL_Neumann_boundary_buffer3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type!=2 && P1[i].buffer_type!=-1) return;		// Outlet boundary condition
	if(P1[i].MOST_buffer==1 && P1[i].buffer_type!=-1) return;    // MOST boundary

	int_t icell,jcell,kcell;
	Real xi, yi, zi;
	Real rhoi, presi;
	Real rho_ref_i;
	Real search_range, hi, tmp_A;
	Real tflt_s;
	Real tux, tuy, tuz, tpres;
	hi=2.0*P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	rhoi=P1[i].rho;
	presi=P1[i].pres;

	tux=tuy=tuz=tpres=tflt_s=0.0;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t ncell=2.0*P1[i].ncell;

	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,uxj,uyj,uzj,tdist,mj,rhoj,presj;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						mj=P1[j].m;
						rhoj=P1[j].rho;
						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						presj=P1[j].pres;

						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
						if(P1[j].p_type<1000){

							if(tdist<search_range){
								// if(tdist<search_range){
								Real twij;
								int_t ptype_j=P1[j].p_type;
								int_t buffer_type_j=P1[j].buffer_type;
								twij=calc_kernel_wij(tmp_A,hi,tdist);

								tux+=mj/rhoj*uxj*twij*(ptype_j==1)*(buffer_type_j==0)*(P1[j].MOST_buffer!=1);
								tuy+=mj/rhoj*uyj*twij*(ptype_j==1)*(buffer_type_j==0)*(P1[j].MOST_buffer!=1);
								tuz+=mj/rhoj*uzj*twij*(ptype_j==1)*(buffer_type_j==0)*(P1[j].MOST_buffer!=1);
								tpres+=mj/rhoj*presj*twij*(ptype_j==1)*(buffer_type_j==0)*(P1[j].MOST_buffer!=1);
								tflt_s+=mj/rhoj*twij*(ptype_j==1)*(buffer_type_j==0)*(P1[j].MOST_buffer!=1);
							}
						}
					}
				}
			}
		}
	}
	if(tflt_s<1e-6){
		P1[i].pres=0.0;
		P1[i].ux=0.0;
		P1[i].uy=0.0;
		P1[i].uz=0.0;
	}
	else{
		if(P1[i].buffer_type==2){
			P1[i].pres = tpres/tflt_s;
			P1[i].ux=tux/tflt_s;
			P1[i].uy=tuy/tflt_s;
			P1[i].uz=tuz/tflt_s;
		}
		else{
			P1[i].pres = tpres/tflt_s;
		}
	}
}

__global__ void KERNEL_donothing_buffer3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type!=2) return;		// Outlet boundary condition
	//if(P1[i].MOST_buffer==1) return;    // MOST boundary

	int_t icell,jcell,kcell;
	Real xi, yi, zi;
	Real rhoi, presi;
	Real rho_ref_i;
	Real search_range, hi, tmp_A;
	Real tflt_s;
	Real tux, tuy, tuz, tpres, ttemp;
	hi=2.0*P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	rhoi=P1[i].rho;
	presi=P1[i].pres;

	tux=tuy=tuz=tpres=ttemp=tflt_s=0.0;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t ncell=2*P1[i].ncell;

	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,uxj,uyj,uzj,tdist,mj,rhoj,presj,tempj;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						mj=P1[j].m;
						rhoj=P1[j].rho;
						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						presj=P1[j].pres;
						tempj=P1[j].temp;
						Real xg=L3-0.008;

						tdist=sqrt((xg-xj)*(xg-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
						if(P1[j].p_type<1000){

							if(tdist<0.001){
								// if(tdist<search_range){
								tux=uxj;
								tuy=uyj;
								tuz=uzj;
								tpres=presj;
								ttemp=tempj;
							}
						}
					}
				}
			}
		}
	}
	if(P1[i].buffer_type==2){
			P1[i].pres = tpres;
			P1[i].ux=tux;
			P1[i].uy=tuy;
			P1[i].uz=tuz;
			P1[i].temp=ttemp;
		}
	else{
		P1[i].pres = tpres;
	}
}

__global__ void KERNEL_open_boundary_extrapolation_buffer3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type!=Outlet)	return;		// Outlet boundary condition
	if(P1[i].MOST_buffer==1) return;

	int_t icell,jcell,kcell;
	Real xi, yi, zi, xgnode_i, ygnode_i, zgnode_i;
	Real uxi, uyi, uzi;
	Real rhoi;
	Real search_range, hi, tmp_A;
	Real tflt_s;
	Real tux, tuy, tuz;
	Real tpres;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;
	rhoi=P1[i].rho;

	tpres=tflt_s=0.0;
	tux=tuy=tuz=0.0;

	if(P1[i].buffer_type==Inlet){
		xgnode_i=2*L2-xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}else if(P1[i].buffer_type==Outlet){
		xgnode_i=2*L3-xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}else if(P1[i].buffer_type==Left){
		xgnode_i=xi;
		ygnode_i=2*(L_left)-yi;
		zgnode_i=zi;
	}else if(P1[i].buffer_type==Right){
		xgnode_i=xi;
		ygnode_i=2*(L_right)-yi;
		zgnode_i=zi;
	}

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xgnode_i-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((ygnode_i-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zgnode_i-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y)+k_NI*k_NJ*(kcell+z);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,uxj,uyj,uzj,tdist,mj,rhoj,presj;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						mj=P1[j].m;
						rhoj=P1[j].rho;
						presj=P1[j].pres;

						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						if(P1[j].p_type<1000){

							if(tdist<search_range){
								// if(tdist<search_range){
								Real twij;
								int_t ptype_j=P1[j].p_type;
								int_t buffer_type_j=P1[j].buffer_type;
								twij=calc_kernel_wij(tmp_A,hi,tdist);

								// tux+=mj/rhoj*uxj*twij*(ptype_j==1&&buffer_type_j!=Outlet);
								// tuy+=mj/rhoj*uyj*twij*(ptype_j==1&&buffer_type_j!=Outlet);
								// tuz+=mj/rhoj*uzj*twij*(ptype_j==1&&buffer_type_j!=Outlet);
								// tflt_s+=mj/rhoj*twij*(ptype_j==1&&buffer_type_j!=Outlet);
								tpres+=mj/rhoj*presj*twij*(ptype_j==1)*(P1[j].MOST_buffer!=1);
								tux+=mj/rhoj*uxj*twij*(ptype_j==1)*(P1[j].MOST_buffer!=1);
								tuy+=mj/rhoj*uyj*twij*(ptype_j==1)*(P1[j].MOST_buffer!=1);
								tuz+=mj/rhoj*uzj*twij*(ptype_j==1)*(P1[j].MOST_buffer!=1);
								tflt_s+=mj/rhoj*twij*(ptype_j==1)*(P1[j].MOST_buffer!=1);
							}
						}
					}
				}
			}
		}
	}
	if(tflt_s<1e-6){
		return;
	}
	else{
		P1[i].pres = 0.0;
		P1[i].ux = tux/tflt_s;
		P1[i].uy = tuy/tflt_s;
		P1[i].uz = tuz/tflt_s;
	}
}

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_periodic_boundary3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type==0) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method

	int_t icell,jcell,kcell;
	Real xi, xgnode_i, yi,ygnode_i, zi,zgnode_i;
	Real rho_ref_i;
	Real search_range, hi, tmp_A, trho, tdrho_x, tdrho_y, tux, tuy, tuz, ttemp;
	Real drho_x, drho_y;
	Real ux_if, uy_if, uz_if, pres_if;
	Real tflt_s;
	Real tpres, tdpres_x, tdpres_y, tdpres_z, dpres_x, dpres_y;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	Real pi_ref, ux_ref, uy_ref, rhoi_ref;
	Real J1, J2, J3, J4;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	// calculate ghost node
	if(P1[i].buffer_type==Left){       // periodic boundary
		xgnode_i=xi;
		ygnode_i=yi+L_right-L_left;
		zgnode_i=zi;
	}
	else if(P1[i].buffer_type==Right){      // periodic boundary
		xgnode_i=xi;
		ygnode_i=yi-L_right+L_left;
		zgnode_i=zi;
	}
	else return;

	trho=tdrho_x=tdrho_y=0.0;
	tux=tuy=tuz=0.0;
	tflt_s=ttemp=0.0;
	tpres=tdpres_x=tdpres_y=tdpres_z=0.0;
	ux_if=uy_if=uz_if=pres_if=0.0;
	pi_ref=uy_ref=ux_ref=0.0;
	rhoi_ref=Inlet_Density;

	int_t ncell=P1[i].ncell;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}

	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	int_t flag=0;

	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(flag==1) continue;
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
						Real volj;
						int_t buffer_type_j, ptype_j;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						mj=P1[j].m;
						hj=P1[j].h;

						rhoj=P1[j].rho;
						volj=P1[j].m/P1[j].rho;

						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						presj=P1[j].pres;
						tempj=P1[j].temp;
						buffer_type_j=P1[j].buffer_type;
						ptype_j=P1[j].p_type;
						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						// if(tdist<search_range){       // for SPH interpolation
						//
						// 	Real hij=0.5*(hi+hj);
						// 	Real tmp_Aij=calc_tmpA(hij);
						// 	Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
						// 	Real neighbor=(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
						//
						// 	tux += uxj*(mj/rhoj)*twij*neighbor;
						// 	tuy += uyj*(mj/rhoj)*twij*neighbor;
						// 	tuz += uzj*(mj/rhoj)*twij*neighbor;
						// 	tpres += presj*(mj/rhoj)*twij*neighbor;
						// 	ttemp += tempj*(mj/rhoj)*twij*neighbor;
						//
						// 	tflt_s += (mj/rhoj)*twij*neighbor;
						// }
						if(tdist<0.001){      // For exact position (KDH)
							tux=uxj;
							tuy=uyj;
							tuz=uzj;
							tpres=presj;
							ttemp=tempj;
							flag=1;
						}
					}
				}
			}
		}
		//P1[i].ux=tux/tflt_s;
		//P1[i].uy=tuy/tflt_s;
		//P1[i].uz=tuz/tflt_s;
	}
	// For SPH interpolation (KDH)
	// P1[i].ux=tux/tflt_s;
	// P1[i].uy=tuy/tflt_s;
	// P1[i].uz=tuz/tflt_s;
	// P1[i].pres=tpres/tflt_s;
	// P1[i].temp=ttemp/tflt_s;
	// P1[i].flt_s = P1[i].flt_s;

	// For exact value (KDH)
	P1[i].ux=tux;
	P1[i].uy=tuy;
	P1[i].uz=tuz;
	P1[i].pres=tpres;
	P1[i].temp=ttemp;
	P1[i].flt_s = P1[i].flt_s;
	return;
}

__global__ void KERNEL_freeslip_boundary3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type!=-1) return;

	int_t icell,jcell,kcell;
	Real xi, xgnode_i, yi,ygnode_i, zi,zgnode_i;
	Real rho_ref_i;
	Real search_range, hi, tmp_A, trho, tdrho_x, tdrho_y, tux, tuy, tuz, ttemp;
	Real drho_x, drho_y;
	Real ux_if, uy_if, uz_if, pres_if;
	Real tflt_s;
	Real tpres, tdpres_x, tdpres_y, tdpres_z, dpres_x, dpres_y;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	Real pi_ref, ux_ref, uy_ref, rhoi_ref;
	Real J1, J2, J3, J4;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	// calculate ghost node
	xgnode_i=xi;
	ygnode_i=yi;
	//zgnode_i=2*(L_ceiling)-zi;
	zgnode_i=zi;

	trho=tdrho_x=tdrho_y=0.0;
	tux=tuy=tuz=0.0;
	tflt_s=ttemp=0.0;
	tpres=tdpres_x=tdpres_y=tdpres_z=0.0;
	ux_if=uy_if=uz_if=pres_if=0.0;
	pi_ref=uy_ref=ux_ref=0.0;
	rhoi_ref=Inlet_Density;

	int_t ncell=P1[i].ncell;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}

	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	for(int_t z=-ncell;z<=ncell;z++){
		for(int_t y=-ncell;y<=ncell;y++){
			for(int_t x=-ncell;x<=ncell;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
						Real volj;
						int_t buffer_type_j, ptype_j;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						mj=P1[j].m;
						hj=P1[j].h;

						rhoj=P1[j].rho;
						volj=P1[j].m/P1[j].rho;

						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						presj=P1[j].pres;
						tempj=P1[j].temp;
						buffer_type_j=P1[j].buffer_type;
						ptype_j=P1[j].p_type;
						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						if(tdist<search_range){       // for SPH interpolation

							// Real hij=0.5*(hi+hj);
							// Real tmp_Aij=calc_tmpA(hij);
							// Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
							Real neighbor=(ptype_j==1);

							Real twij=calc_kernel_wij(tmp_A,hi,tdist);
							Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);

							tux += uxj*(mj/rhoj)*twij*neighbor;
							tuy += uyj*(mj/rhoj)*twij*neighbor;
							tuz += uzj*(mj/rhoj)*twij*neighbor;
							tpres += presj*(mj/rhoj)*twij*neighbor;
							ttemp += tempj*(mj/rhoj)*twij*neighbor;

							tflt_s += (mj/rhoj)*twij*neighbor;
						}

					}
				}
			}
		}
	}
	// For SPH interpolation (KDH)
	if(tflt_s>1e-5){
		P1[i].ux=tux/tflt_s;
		P1[i].uy=tuy/tflt_s;
		P1[i].uz=tuz/tflt_s;
		P1[i].pres=tpres/tflt_s;
		P1[i].temp=ttemp/tflt_s;
		P1[i].flt_s = P1[i].flt_s;
	}
		return;
}

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_inout_boundary3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type==0) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method

	int_t icell,jcell,kcell;
	Real xi, xgnode_i, yi,ygnode_i, zi,zgnode_i;
	Real rho_ref_i;
	Real search_range, hi, tmp_A, trho, tdrho_x, tdrho_y, tux, tuy, tuz, ttemp;
	Real drho_x, drho_y;
	Real ux_if, uy_if, uz_if, pres_if;
	Real tflt_s;
	Real tpres, tdpres_x, tdpres_y, tdpres_z, dpres_x, dpres_y;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	Real pi_ref, ux_ref, uy_ref, rhoi_ref;
	Real J1, J2, J3, J4;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	// calculate ghost node
	if(P1[i].buffer_type==Inlet){
		xgnode_i=xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}
	else if(P1[i].buffer_type==Outlet){	 // To ensure that gradient=0
		xgnode_i=2*L3-xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}
	else if(P1[i].buffer_type==Left){       // periodic boundary
		xgnode_i=xi;
		ygnode_i=yi+L_right-L_left;
		zgnode_i=zi;
	}
	else if(P1[i].buffer_type==Right){      // periodic boundary
		xgnode_i=xi;
		ygnode_i=yi-L_right+L_left;
		zgnode_i=zi;
	}
	else return;

	// else if(P1[i].buffer_type==Left){        // Mirroring boundary
	// 	ygnode_i=2*(L_left)-yi;
	// 	xgnode_i=xi;
	// 	//ygnode_i=L3;
	// }else if(P1[i].buffer_type==Right){
	// 	ygnode_i=2*(L_right)-yi;
	// 	xgnode_i=xi;
	// }

	trho=tdrho_x=tdrho_y=0.0;
	tux=tuy=tuz=0.0;
	tflt_s=ttemp=0.0;
	tpres=tdpres_x=tdpres_y=tdpres_z=0.0;
	ux_if=uy_if=uz_if=pres_if=0.0;
	pi_ref=uy_ref=ux_ref=0.0;
	rhoi_ref=Inlet_Density;

	int_t ncell=P1[i].ncell;

	// calculate I,J,K in cell
	if(P1[i].buffer_type==Left || P1[i].buffer_type==Right){
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}
		// out-of-range handling
		if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

		int_t flag=0;

		for(int_t z=-ncell;z<=ncell;z++){
			for(int_t y=-ncell;y<=ncell;y++){
				for(int_t x=-ncell;x<=ncell;x++){
					int_t k=idx_cell(icell+x,jcell+y,kcell+z);
					if(flag==1) continue;
					if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
					if(g_str[k]!=cu_memset){
						int_t fend=g_end[k];
						for(int_t j=g_str[k];j<fend;j++){
							Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
							Real volj;
							int_t buffer_type_j, ptype_j;

							xj=P1[j].x;
							yj=P1[j].y;
							zj=P1[j].z;
							mj=P1[j].m;
							hj=P1[j].h;

							rhoj=P1[j].rho;
							volj=P1[j].m/P1[j].rho;

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;
							presj=P1[j].pres;
							tempj=P1[j].temp;

							buffer_type_j=P1[j].buffer_type;
							ptype_j=P1[j].p_type;

							tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
							// if(tdist<search_range){
							//
							// 	Real hij=0.5*(hi+hj);
							// 	Real tmp_Aij=calc_tmpA(hij);
							// 	Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
							// 	Real neighbor=(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
							//
							// 	tux += uxj*(mj/rhoj)*twij*neighbor;
							// 	tuy += uyj*(mj/rhoj)*twij*neighbor;
							// 	tuz += uzj*(mj/rhoj)*twij*neighbor;
							// 	tpres += presj*(mj/rhoj)*twij*neighbor;
							// 	ttemp += tempj*(mj/rhoj)*twij*neighbor;
							//
							// 	tflt_s += (mj/rhoj)*twij*neighbor;
							// }
							if(tdist<0.001){      // For exact position (KDH)
								tux=uxj;
								tuy=uyj;
								tuz=uzj;
								tpres=presj;
								ttemp=tempj;
								flag=1;
							}
						}
					}
				}
			}
			//P1[i].ux=tux/tflt_s;
			//P1[i].uy=tuy/tflt_s;
			//P1[i].uz=tuz/tflt_s;
		}
		// For SPH interpolation (KDH)
		// P1[i].ux=tux/tflt_s;
		// P1[i].uy=tuy/tflt_s;
		// P1[i].uz=tuz/tflt_s;
		// P1[i].pres=tpres/tflt_s;
		// P1[i].temp=ttemp/tflt_s;
		// P1[i].flt_s = P1[i].flt_s;

		// For exact value (KDH)
		P1[i].ux=tux;
		P1[i].uy=tuy;
		P1[i].uz=tuz;
		P1[i].pres=tpres;
		P1[i].temp=ttemp;
		P1[i].flt_s = P1[i].flt_s;
	}
	else if(P1[i].buffer_type==Outlet){         // Outlet일 경우 interface와의 gradient로 교정
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((L3-k_x_min)/k_dcell),k_NI-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}     // 계면 위치에서 값 추정

		// out-of-range handling
		if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

		tux=P1[i].ux;
		tuy=P1[i].uy;
		tuz=P1[i].uz;
		tpres=P1[i].pres;

		for(int_t z=-ncell;z<=ncell;z++){
			for(int_t y=-ncell;y<=ncell;y++){
				for(int_t x=-ncell;x<=ncell;x++){
					int_t k=idx_cell(icell+x,jcell+y,kcell+z);

					if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
					if(g_str[k]!=cu_memset){
						int_t fend=g_end[k];
						for(int_t j=g_str[k];j<fend;j++){
							Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
							Real volj;
							int_t buffer_type_j, ptype_j;

							xj=P1[j].x;
							yj=P1[j].y;
							zj=P1[j].z;

							rhoj=P1[j].rho;
							mj=P1[j].m;
							volj=P1[j].m/P1[j].rho;

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;
							hj=P1[j].h;
							presj=P1[j].pres;
							tempj=P1[j].temp;

							buffer_type_j=P1[j].buffer_type;
							ptype_j=P1[j].p_type;

							tdist=sqrt((L3-xj)*(L3-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;

							if(tdist<search_range){

								Real hij=0.5*(hi+hj);
								Real tmp_Aij=calc_tmpA(hij);
								Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);

								ux_if += uxj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								uy_if += uyj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								uz_if += uzj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								pres_if += presj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);

								tflt_s += (mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
							}
						}
					}
				}
			}
		}
		P1[i].ux = -tux+2.0*ux_if/tflt_s;
		P1[i].uy = -tuy+2.0*uy_if/tflt_s;
		P1[i].uz = -tuz+2.0*uz_if/tflt_s;
		P1[i].pres = -tpres+2.0*pres_if/tflt_s;
	}
	// else if(P1[i].buffer_type==Inlet){       // Inlet에서는 압력만 조절
	// 	if((k_x_max==k_x_min)){icell=0;}
	// 	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
	// 	if((k_y_max==k_y_min)){jcell=0;}
	// 	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
	// 	if((k_z_max==k_z_min)){kcell=0;}
	// 	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}     // 계면 위치에서 값 추정
	//
	// 	// out-of-range handling
	// 	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;
	//
	// 	tpres=P1[i].pres;
	//
	// 	for(int_t z=-ncell;z<=ncell;z++){
	// 		for(int_t y=-ncell;y<=ncell;y++){
	// 			for(int_t x=-ncell;x<=ncell;x++){
	// 				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
	//
	// 				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
	// 				if(g_str[k]!=cu_memset){
	// 					int_t fend=g_end[k];
	// 					for(int_t j=g_str[k];j<fend;j++){
	// 						Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
	// 						Real volj;
	// 						int_t buffer_type_j, ptype_j;
	//
	// 						xj=P1[j].x;
	// 						yj=P1[j].y;
	// 						zj=P1[j].z;
	//
	// 						rhoj=P1[j].rho;
	// 						mj=P1[j].m;
	//
	// 						hj=P1[j].h;
	// 						presj=P1[j].pres;
	// 						tempj=P1[j].temp;
	//
	// 						buffer_type_j=P1[j].buffer_type;
	// 						ptype_j=P1[j].p_type;
	//
	// 						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
	//
	// 						if(tdist<search_range){
	//
	// 							Real hij=0.5*(hi+hj);
	// 							Real tmp_Aij=calc_tmpA(hij);
	// 							Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);
	//
	// 							tpres += presj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
	// 							tflt_s += (mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// 	if(tflt_s>1e-5) P1[i].pres=tpres/tflt_s;
	// 	else P1[i].pres=0.0;
	// 	P1[i].flt_s = tflt_s;
	// }
	return;
}

__global__ void KERNEL_periodic_boundary3D_ori(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type==0) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method

	int_t icell,jcell,kcell;
	Real xi, xgnode_i, yi,ygnode_i, zi,zgnode_i;
	Real rho_ref_i;
	Real search_range, hi, tmp_A, trho, tdrho_x, tdrho_y, tux, tuy, tuz, ttemp;
	Real drho_x, drho_y;
	Real ux_if, uy_if, uz_if, pres_if;
	Real tflt_s;
	Real tpres, tdpres_x, tdpres_y, tdpres_z, dpres_x, dpres_y;
	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	Real pi_ref, ux_ref, uy_ref, rho_ref;
	Real J1, J2, J3, J4;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;

	// calculate ghost node

	if(P1[i].buffer_type==Inlet){
		xgnode_i=2*L2-xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}
	if(P1[i].buffer_type==Outlet){	 // To ensure that gradient=0
		xgnode_i=2*L3-xi;
		ygnode_i=yi;
		zgnode_i=zi;
	}
	if(P1[i].buffer_type==Left){       // periodic boundary
		xgnode_i=xi;
		//ygnode_i=yi-46.2;
		ygnode_i=yi+L_right-L_left;
		zgnode_i=zi;
	}
	if(P1[i].buffer_type==Right){
		xgnode_i=xi;
		//ygnode_i=yi+46.2;
		ygnode_i=yi-L_right+L_left;
		zgnode_i=zi;
	}
	// else if(P1[i].buffer_type==Left){        // Mirroring boundary
	// 	ygnode_i=2*(L_left)-yi;
	// 	xgnode_i=xi;
	// 	//ygnode_i=L3;
	// }else if(P1[i].buffer_type==Right){
	// 	ygnode_i=2*(L_right)-yi;
	// 	xgnode_i=xi;
	// }

	trho=tdrho_x=tdrho_y=tux=tuy=tuz=tflt_s=tpres=tdpres_x=tdpres_y=tdpres_z=ttemp=0.0;
	ux_if=uy_if=uz_if=pres_if=0.0;
	pi_ref=0.0; uy_ref=1.0;  ux_ref=0.0; rho_ref=Inlet_Density;

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
		for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
			for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
						Real volj;
						int_t buffer_type_j, ptype_j;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;
						mj=P1[j].m;
						hj=P1[j].h;

						rhoj=P1[j].rho;
						volj=P1[j].m/P1[j].rho;

						uxj=P1[j].ux;
						uyj=P1[j].uy;
						uzj=P1[j].uz;
						presj=P1[j].pres;
						tempj=P1[j].temp;

						buffer_type_j=P1[j].buffer_type;
						ptype_j=P1[j].p_type;

						tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
						if(tdist<search_range){

							Real hij=0.5*(hi+hj);
							Real tmp_Aij=calc_tmpA(hij);
							Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);

							tux += uxj*(mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
							tuy += uyj*(mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
							tuz += uzj*(mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
							tpres += presj*(mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
							ttemp += tempj*(mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);

							tflt_s += (mj/rhoj)*twij*(buffer_type_j!=Left)*(buffer_type_j!=Right)*(ptype_j==1);
						}
					}
				}
			}
		}

		//P1[i].ux=tux/tflt_s;
		//P1[i].uy=tuy/tflt_s;
		//P1[i].uz=tuz/tflt_s;
		P1[i].ux=tux/tflt_s;
		P1[i].uy=tuy/tflt_s;
		P1[i].uz=tuz/tflt_s;
		P1[i].pres=tpres/tflt_s;
		P1[i].temp=ttemp/tflt_s;

		//printf("ux=%f pres=%f flt=%f\n",P1[i].ux, P1[i].pres, tflt_s);
		//if (P1[i].buffer_type==Outlet) P1[i].pres=0.0;   // Pressure boundary condition for outlet
	}

	tux=tuy=tuz=tpres=ttemp=tflt_s=0.0;

	if(P1[i].buffer_type==Outlet){         // Outlet일 경우 interface와의 gradient로 교정
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((L3-k_x_min)/k_dcell),k_NI-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}     // 계면 위치에서 값 추정

		tux=P1[i].ux;
		tuy=P1[i].uy;
		tuz=P1[i].uz;
		tpres=P1[i].pres;

		for(int_t z=-P1[i].ncell;z<=P1[i].ncell;z++){
			for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
				for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
					int_t k=idx_cell(icell+x,jcell+y,kcell+z);

					if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;

					if(g_str[k]!=cu_memset){
						int_t fend=g_end[k];
						for(int_t j=g_str[k];j<fend;j++){
							Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,hj,presj,tempj,hij,tmp_Aij;
							Real volj;
							int_t buffer_type_j, ptype_j;

							xj=P1[j].x;
							yj=P1[j].y;
							zj=P1[j].z;

							rhoj=P1[j].rho;
							mj=P1[j].m;
							volj=P1[j].m/P1[j].rho;

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;
							hj=P1[j].h;
							presj=P1[j].pres;
							tempj=P1[j].temp;

							buffer_type_j=P1[j].buffer_type;
							ptype_j=P1[j].p_type;

							tdist=sqrt((L3-xj)*(L3-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;

							if(tdist<search_range){

								Real hij=0.5*(hi+hj);
								Real tmp_Aij=calc_tmpA(hij);
								Real twij=calc_kernel_wij(tmp_Aij,hij,tdist);

								ux_if += uxj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								uy_if += uyj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								uz_if += uzj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
								pres_if += presj*(mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);

								tflt_s += (mj/rhoj)*twij*(buffer_type_j==0)*(ptype_j==1);
							}
						}
					}
				}
			}
		}

		P1[i].ux = -tux+2.0*ux_if/tflt_s;
		P1[i].uy = -tuy+2.0*uy_if/tflt_s;
		P1[i].uz = -tuz+2.0*uz_if/tflt_s;
		P1[i].pres = -tpres+2.0*pres_if/tflt_s;
	}
	return;
}

////////////////////////////////////////////////////////////////////////
__global__ void set_inlet_buffer(part1*P1, Real ttime)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type!=1) return;

	Real xi = 0;
	Real H = 0.2;
	Real d = 0.3;
	Real C = sqrt(Gravitational_CONST*(H+d));
	Real T = 2*d/C*sqrt(4*d/3/H)*(3.8+H/d);
	// Real T=0;
	Real X = sqrt(3*H/4/d/d/d)*(xi-C*(ttime-T));

	// Real limit = H*(1.0/(exp(2*X)+exp(-2*X)+2.0))+d;
	Real k, w;
	k = PI/1.0;
	w = PI/1.0;

	Real limit = H/cos(PI/2.0-w*ttime)+d;


	if(P1[i].y>limit){
		P1[i].i_type=3;
	}else{
		P1[i].i_type=2;
	}
}

__global__ void KERNEL_inlet_boundary(part1*P1){
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	Real a=0.0;
	Real z_ref=0.0;

	switch(terrain_type){
		case 0:
			a=0.14;
			z_ref=300.0;
			printf("terrain type: flat\n");
			break;
		case 1:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: Other\n");
			break;
		case 2:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		case 3:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		case 4:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		default:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: flat\n");
			break;
	}
	Real u_ref=1.0;
	if(P1[i].buffer_type==1) P1[i].ux=u_ref*pow(P1[i].y/z_ref,a);
	else P1[i].ux=0.0001;

	return;
}

__global__ void KERNEL_Recycling3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].buffer_type!=Inlet&&P1[i].buffer_type!=-1) return;

	int_t icell,jcell,kcell;
	Real xi,yi,zi,uxi,uyi,uzi;
	Real tux, tuy, tuz;
	Real search_range,hi,tmp_A;

	hi=P1[i].h;
	tmp_A=calc_tmpA(hi);
	search_range=k_search_kappa*hi;	// search range

	xi=P1[i].x+980.0;    // Length of the recycling station
	yi=P1[i].y;
	zi=P1[i].z;
	uxi=P1[i].ux;
	uyi=P1[i].uy;
	uzi=P1[i].uz;

	Real z_ref=300.0;
	Real z_i = max(0.0, zi);

	// calculate I,J,K in cell
	if((k_x_max==k_x_min)){icell=0;}
	else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
	if((k_y_max==k_y_min)){jcell=0;}
	else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
	if((k_z_max==k_z_min)){kcell=0;}
	else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
	// out-of-range handling
	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

	tux=tuy=tuz=0.0;

	for(int_t z=-1;z<=1;z++){
		for(int_t y=-1;y<=1;y++){
			for(int_t x=-1;x<=1;x++){
				// int_t k=(icell+x)+k_NI*(jcell+y);
				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
				if(g_str[k]!=cu_memset){
					int_t fend=g_end[k];
					for(int_t j=g_str[k];j<fend;j++){
						Real xj,yj,zj,tdist,uxj,uyj,uzj,mj,rhoj,presj;
						Real volj;

						xj=P1[j].x;
						yj=P1[j].y;
						zj=P1[j].z;


						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-1e-20;
						if(tdist<1e-2){

							uxj=P1[j].ux;
							uyj=P1[j].uy;
							uzj=P1[j].uz;

							//printf("uxj=%f\n", uxj);
							tux=uxj;
							tuy=uyj;
							tuz=uzj;
							break;
						}
					}
				}
			}
		}
	}

	if(P1[i].buffer_type==Inlet||P1[i].buffer_type==-1){
		// P1[i].ux=tux;
		// P1[i].uy=tuy;
		// P1[i].uz=tuz;
		P1[i].ux=(us_update/k_vonKarman)*log(P1[i].z/MOST_z0);

		// P1[i].ux_f=0.1*tux;
		// P1[i].uy_f=tuy;
		// P1[i].uz_f=tuz;

		return;
	}
}
