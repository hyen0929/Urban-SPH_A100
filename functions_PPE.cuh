// PPE_WALL_NB(P,j): may neighbour j enter the pressure sums?
// Production (IBM_PPE_WALL==0): fluid + ground only; building markers and
// interior fluid are both dropped, which leaves the operator truncated and
// unclosed at building surfaces. With IBM_PPE_WALL==1 the markers are kept,
// so the wall contributes div(u) with u=0 and a Neumann pressure term, the
// same way the ground (p_type==0) already does. The i_type guard keeps
// DEACTIVATE_IBM'd and out-of-subdomain markers out. See Parameters.cuh.
#if IBM_PPE_WALL
#define PPE_WALL_NB(P,j) (((P)[j].p_type<1000 || ((P)[j].p_type==1000 && (P)[j].i_type<=i_type_crt)) && (P)[j].buffer_type!=99)
#else
#define PPE_WALL_NB(P,j) ((P)[j].p_type<1000 && (P)[j].buffer_type!=99)
#endif

// g_bi (optional, may be NULL): per-particle |bi| = |rho*div(u*)/dt|, the PPE
// source magnitude. Written only for fluid particles; caller must zero it first
// (skipped particles early-return). Used by the PPE_DIAG residual diagnostic.
__global__ void KERNEL_PPE_prep3D(Real tdt, int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3, Real* g_bi)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;
	if(P1[i].buffer_type==99)	return;	// Interior fluid (inside building) — skip PPE

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

						if(PPE_WALL_NB(P1,j)){
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
									AijPij += c*presj;
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

	bi = weight*(bix+biy+biz)/tdt;

	if(g_bi!=NULL) g_bi[i]=fabsf(bi);   // PPE source magnitude for PPE_DIAG (|rho*div(u*)/dt|)

	// bi = -rhoi*(1.0-drhostar/flt)/tdt/tdt;
	//P1[i].PPE1 = rhoi*(bix+biy+biz)/tdt;
	//P1[i].PPE2 = -rhoi*(1.0-drhostar/flt)/tdt/tdt;  // YHS

	P1[i].PPE1 = bi+AijPij;
	P1[i].PPE2 = Aij;
 	//P1[i].pres=(bi+AijPij)/Aij;
	//if(xi>4.5 && xi<15.5 && yi>4 && yi<6 && zi>7 && zi<8) printf("xi=%f bix=%f biy=%f biz=%f Aij=%f AijPij=%f\n",xi,bix,biy,biz,Aij,AijPij);

	//if(isnan(P1[i].pres)==1) printf("xi=%f bi=%f bix=%f biy=%f biz=%f Aij=%f AijPij=%f\n",xi,bi,bix,biy,biz,Aij,AijPij);

	// P1[i].fv_x=bi;
	// P1[i].fv_y=Aij;
	// P1[i].fv_z=AijPij;

}

__global__ void KERNEL_PPE_calc3D(part1*P1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;

	P1[i].pres=(P1[i].PPE1)/(P1[i].PPE2+1e-20);
}

// True PPE residual per particle: |A*p - b| = |PPE2*pres - PPE1|  (~ |div(u^{n+1})|).
// MUST be called after a fresh KERNEL_PPE_prep3D evaluated with the post-solve
// pressure; otherwise pres = PPE1/PPE2 straight from PPE_calc makes this identically 0.
// resid must be zeroed by the caller (non-fluid particles early-return).
__global__ void KERNEL_compute_ppe_residual(part1*P1, Real* resid)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000)	return;		// Immersed Boundary Method
	if(P1[i].p_type<=0)	return;
	if(P1[i].buffer_type==99)	return;	// Interior fluid (inside building)

	resid[i]=fabsf(P1[i].PPE2*P1[i].pres - P1[i].PPE1);
}

// Snapshot current pressure into pres_old (so the first Jacobi delta is well-defined).
__global__ void KERNEL_snapshot_pres(part1*P1, Real* pres_old)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	pres_old[i]=P1[i].pres;
}

// Per-particle Jacobi delta |p_new - p_old| (early-exit gate), then snapshot p_new -> pres_old
// for the next sweep. Non-fluid particles get delta=0 so they don't pollute the Max reduce.
__global__ void KERNEL_compute_pres_delta(part1*P1, Real* pres_old, Real* pres_delta)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;

	Real pnew=P1[i].pres;
	bool is_fluid = !(P1[i].i_type>i_type_crt || P1[i].p_type>=1000 || P1[i].p_type<=0 || P1[i].buffer_type==99);
	pres_delta[i] = is_fluid ? fabsf(pnew - pres_old[i]) : 0.0f;
	pres_old[i]=pnew;   // snapshot for next iteration
}

//=================================================================================================
// IISPH consistent-operator PPE (Parameters.cuh CONSISTENT_PPE). Plain gradW (no KGC) so the
// matrix A == divergence of the applied symmetric pressure force a_p=-sum m_j(p_i/rho_i^2+p_j/rho_j^2)gradW.
//   prep  : G_i=sum m_j gradW, H_i=sum m_j|gradW|^2, D_i=sum m_j (v_i-v_j).gradW
//           -> d_ii=-dt^2/rho_i^2 G_i,  a_ii=-dt^2/rho_i^2(|G_i|^2+m_i H_i),  s_i=rho0-(rho_i+dt D_i)
//   pass1 : Delta x_i = d_ii p_i + sum_j d_ij p_j,  d_ij=-dt^2 m_j/rho_j^2 gradW_ij
//   pass2 : Ap_i = sum_j m_j (dx_i-dx_j).gradW_ij ;  p_i += (omega/a_ii)(s_i-Ap_i), clamp >=0
//=================================================================================================
#define IISPH_NEI_BOILERPLATE \
	int_t icell,jcell,kcell; \
	if((k_x_max==k_x_min)){icell=0;}else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);} \
	if((k_y_max==k_y_min)){jcell=0;}else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);} \
	if((k_z_max==k_z_min)){kcell=0;}else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);} \
	if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

__global__ void KERNEL_iisph_prep(Real tdt,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000) return;
	if(P1[i].p_type<=0) return;
	if(P1[i].buffer_type==99) return;

	Real xi=P1[i].x, yi=P1[i].y, zi=P1[i].z;
	Real uxi=P1[i].ux, uyi=P1[i].uy, uzi=P1[i].uz;
	Real rhoi=P1[i].rho, mi=P1[i].m;
	Real hi=Pressure_length*P1[i].h, tmp_A=calc_tmpA(hi);
	Real search_range=k_search_kappa*hi;

	Real Gx=0,Gy=0,Gz=0,H=0,D=0;
	IISPH_NEI_BOILERPLATE
	int_t ncell=P1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++) for(int_t y=-ncell;y<=ncell;y++) for(int_t x=-ncell;x<=ncell;x++){
		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
		int_t k=idx_cell(icell+x,jcell+y,kcell+z);
		if(g_str[k]==cu_memset) continue;
		int_t fend=g_end[k];
		for(int_t j=g_str[k];j<fend;j++){
			if(P1[j].i_type==4) continue;
			if(P1[j].p_type>=1000||P1[j].buffer_type==99) continue;
			Real xj=P1[j].x,yj=P1[j].y,zj=P1[j].z,hj=P1[j].h;
			Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
			Real sr=k_search_kappa*fmax(hi,hj);
			if(tdist>=sr) continue;
			Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
			Real tdwx=tdwij*(xi-xj)/tdist, tdwy=tdwij*(yi-yj)/tdist, tdwz=tdwij*(zi-zj)/tdist;
#if IISPH_KGC
			apply_gradient_correction_3D(P3[i].Cm,0.0f,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
#endif
			Real mj=P1[j].m;
			Gx+=mj*tdwx; Gy+=mj*tdwy; Gz+=mj*tdwz;
			H +=mj*(tdwx*tdwx+tdwy*tdwy+tdwz*tdwz);
			D +=mj*((uxi-P1[j].ux)*tdwx+(uyi-P1[j].uy)*tdwy+(uzi-P1[j].uz)*tdwz);
		}
	}
	Real dt2=tdt*tdt, inv_r2=1.0f/(rhoi*rhoi);
	P3[i].dii_x=-dt2*inv_r2*Gx; P3[i].dii_y=-dt2*inv_r2*Gy; P3[i].dii_z=-dt2*inv_r2*Gz;
	P3[i].iisph_aii=-dt2*inv_r2*(Gx*Gx+Gy*Gy+Gz*Gz + mi*H);
	Real rho0=P2[i].rho_ref;
	P3[i].iisph_src=(rho0-rhoi)-tdt*D;     // s_i = rho0 - rho*_i
}

__global__ void KERNEL_iisph_pass1(Real tdt,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000) return;
	if(P1[i].p_type<=0) return;
	if(P1[i].buffer_type==99) return;

	Real xi=P1[i].x, yi=P1[i].y, zi=P1[i].z;
	Real pi=P1[i].pres;
	Real hi=Pressure_length*P1[i].h, tmp_A=calc_tmpA(hi);
	Real search_range=k_search_kappa*hi;
	Real dt2=tdt*tdt;

	Real axx=0,axy=0,axz=0;   // sum_j d_ij p_j
	IISPH_NEI_BOILERPLATE
	int_t ncell=P1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++) for(int_t y=-ncell;y<=ncell;y++) for(int_t x=-ncell;x<=ncell;x++){
		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
		int_t k=idx_cell(icell+x,jcell+y,kcell+z);
		if(g_str[k]==cu_memset) continue;
		int_t fend=g_end[k];
		for(int_t j=g_str[k];j<fend;j++){
			if(P1[j].i_type==4) continue;
			if(P1[j].p_type>=1000||P1[j].buffer_type==99) continue;
			Real xj=P1[j].x,yj=P1[j].y,zj=P1[j].z,hj=P1[j].h;
			Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
			Real sr=k_search_kappa*fmax(hi,hj);
			if(tdist>=sr) continue;
			Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
			Real tdwx=tdwij*(xi-xj)/tdist, tdwy=tdwij*(yi-yj)/tdist, tdwz=tdwij*(zi-zj)/tdist;
#if IISPH_KGC
			apply_gradient_correction_3D(P3[i].Cm,0.0f,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
#endif
			Real rhoj=P1[j].rho, mj=P1[j].m, pj=P1[j].pres;
			Real dij=-dt2*mj/(rhoj*rhoj)*pj;   // d_ij . (unit) * pj
			axx+=dij*tdwx; axy+=dij*tdwy; axz+=dij*tdwz;
		}
	}
	P3[i].dx_x=P3[i].dii_x*pi+axx;
	P3[i].dx_y=P3[i].dii_y*pi+axy;
	P3[i].dx_z=P3[i].dii_z*pi+axz;
}

// resid (optional): per-particle density residual |s_i - Ap_i| for the early-exit
// convergence check (avg over fluid = post-projection density error). Caller zeroes it.
__global__ void KERNEL_iisph_pass2(Real tdt,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3, Real* resid)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type>i_type_crt) return;
	if(P1[i].p_type>=1000) return;
	if(P1[i].p_type<=0) return;
	if(P1[i].buffer_type==99) return;

	Real xi=P1[i].x, yi=P1[i].y, zi=P1[i].z;
	Real dxi_x=P3[i].dx_x, dxi_y=P3[i].dx_y, dxi_z=P3[i].dx_z;
	Real hi=Pressure_length*P1[i].h, tmp_A=calc_tmpA(hi);
	Real search_range=k_search_kappa*hi;

	Real Ap=0;   // sum_j m_j (dx_i - dx_j).gradW_ij
	IISPH_NEI_BOILERPLATE
	int_t ncell=P1[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++) for(int_t y=-ncell;y<=ncell;y++) for(int_t x=-ncell;x<=ncell;x++){
		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
		int_t k=idx_cell(icell+x,jcell+y,kcell+z);
		if(g_str[k]==cu_memset) continue;
		int_t fend=g_end[k];
		for(int_t j=g_str[k];j<fend;j++){
			if(P1[j].i_type==4) continue;
			if(P1[j].p_type>=1000||P1[j].buffer_type==99) continue;
			Real xj=P1[j].x,yj=P1[j].y,zj=P1[j].z,hj=P1[j].h;
			Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
			Real sr=k_search_kappa*fmax(hi,hj);
			if(tdist>=sr) continue;
			Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
			Real tdwx=tdwij*(xi-xj)/tdist, tdwy=tdwij*(yi-yj)/tdist, tdwz=tdwij*(zi-zj)/tdist;
#if IISPH_KGC
			apply_gradient_correction_3D(P3[i].Cm,0.0f,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
#endif
			Real mj=P1[j].m;
			Ap+=mj*((dxi_x-P3[j].dx_x)*tdwx+(dxi_y-P3[j].dx_y)*tdwy+(dxi_z-P3[j].dx_z)*tdwz);
		}
	}
	Real src=P3[i].iisph_src;
	if(resid!=NULL) resid[i]=fabsf(src-Ap);   // density residual of the current iterate
	Real aii=P3[i].iisph_aii;
	if(fabsf(aii)<1e-20f) return;          // isolated particle: leave pressure unchanged
	Real pnew=P1[i].pres+IISPH_OMEGA/aii*(src-Ap);
	P1[i].pres=fmaxf(pnew,0.0f);           // clamp negative pressure (IISPH)
}

// Divergence of the CORRECTED velocity u^{n+1}: out[i] = |rho*div(u^{n+1})/dt|.
// Same operator/units as bi in KERNEL_PPE_prep3D (KGC-corrected gradient) so it is
// directly comparable to bi_max. Geometry/rho/m read from Pg (=dev_SP1, fully valid),
// velocity from Pv (=dev_P1, holds u^{n+1} after the corrector). Caller zeroes out[].
// use_kgc: 1 = KGC-corrected gradient (matches explicit's PPE operator), 0 = plain gradW
// (matches the IISPH operator). Lets each run be graded on its own divergence operator.
__global__ void KERNEL_compute_div_un1(Real tdt,int_t*g_str,int_t*g_end,part1*Pg,part1*Pv,part2*P2,part3*P3,Real* out,int use_kgc)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(Pg[i].i_type>i_type_crt) return;
	if(Pg[i].p_type>=1000) return;
	if(Pg[i].p_type<=0) return;
	if(Pg[i].buffer_type==99) return;

	Real xi=Pg[i].x, yi=Pg[i].y, zi=Pg[i].z;
	Real uxi=Pv[i].ux, uyi=Pv[i].uy, uzi=Pv[i].uz;     // corrected velocity u^{n+1}
	Real rhoi=Pg[i].rho;
	Real hi=Pressure_length*Pg[i].h, tmp_A=calc_tmpA(hi);
	Real search_range=k_search_kappa*hi;
	Real gc11=P3[i].inv_cm_xx, gc12=P3[i].inv_cm_xy, gc22=P3[i].inv_cm_yy,
	     gc23=P3[i].inv_cm_yz, gc31=P3[i].inv_cm_zx, gc33=P3[i].inv_cm_zz;

	Real bix=0,biy=0,biz=0;
	IISPH_NEI_BOILERPLATE
	int_t ncell=Pg[i].ncell;
	for(int_t z=-ncell;z<=ncell;z++) for(int_t y=-ncell;y<=ncell;y++) for(int_t x=-ncell;x<=ncell;x++){
		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
		int_t k=idx_cell(icell+x,jcell+y,kcell+z);
		if(g_str[k]==cu_memset) continue;
		int_t fend=g_end[k];
		for(int_t j=g_str[k];j<fend;j++){
			if(Pg[j].i_type==4) continue;
			if(Pg[j].p_type>=1000||Pg[j].buffer_type==99) continue;
			Real xj=Pg[j].x,yj=Pg[j].y,zj=Pg[j].z,hj=Pg[j].h;
			Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
			Real sr=k_search_kappa*fmax(hi,hj);
			if(tdist>=sr) continue;
			Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
			Real tdwx=tdwij*(xi-xj)/tdist, tdwy=tdwij*(yi-yj)/tdist, tdwz=tdwij*(zi-zj)/tdist;
			Real gx=tdwx, gy=tdwy, gz=tdwz;            // plain gradW (IISPH operator)
			if(use_kgc){                               // KGC-corrected (explicit operator)
				gx=gc11*tdwx+gc12*tdwy+gc31*tdwz;
				gy=gc12*tdwx+gc22*tdwy+gc23*tdwz;
				gz=gc31*tdwx+gc23*tdwy+gc33*tdwz;
			}
			Real mj=Pg[j].m, rhoj=Pg[j].rho;
			Real w=2.0*rhoi*rhoj/(rhoi+rhoj);
			bix += mj/rhoj*(Pv[j].ux-uxi)*gx*w;
			biy += mj/rhoj*(Pv[j].uy-uyi)*gy*w;
			biz += mj/rhoj*(Pv[j].uz-uzi)*gz*w;
		}
	}
	out[i]=fabsf((bix+biy+biz)/tdt);
}

// Fluid-particle indicator (1.0 fluid / 0.0 else) for averaging counts. Caller zeroes mask[].
__global__ void KERNEL_fluid_mask(part1*P1, Real* mask)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	bool fl = !(P1[i].i_type>i_type_crt || P1[i].p_type>=1000 || P1[i].p_type<=0 || P1[i].buffer_type==99);
	mask[i]= fl?1.0f:0.0f;
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
// minmod: 부호가 다르면 0, 같으면 절댓값 작은 쪽
__device__ inline Real minmod2(Real a, Real b){ return (a*b<=0.0) ? 0.0 : (fabsf(a)<fabsf(b) ? a : b); }
#if RIEMANN_EULERV==3
#define RECON(a,b) (a)               // 선형 upwind: 리미터 없음
#else
#define RECON(a,b) minmod2(a,b)
#endif

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
	Real search_range,tmp_A,tmp_Rc,tmp_Rc_diff,tmp_Rd,tmp_Rd2;
	Real tmpx,tmpy,tmpz;
	Real eulerx, eulery, eulerz, eulert;
	Real cpi;
	Real tmp_temp,tmp_flt;

	Real fvx,fvy,fvz,fex,fey,fez,fbs;

	ptypei=P1[i].p_type;
	buffer_type_i=P1[i].buffer_type;

	Real gc11, gc12, gc22, gc23, gc31, gc33;

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
	tmp_Rc_diff=0.0;
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

	Real wx=0.0, wy=0.0, wz=0.0;

	Real a_conv_x = 0.0;
	Real a_conv_y = 0.0;
	Real a_conv_z = 0.0;

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

						// IBM walls treated as ADIABATIC (no fluid-IBM conduction).
						// Reason: previous Dirichlet T=T_con on all IBM particles (incl. roof, top
						// edges, side walls) injected unphysical heat at corners where the SPH
						// operator has severe boundary truncation. Combined with Boussinesq
						// buoyancy, this drove a ±60K checkerboard at the windward roof corner.
						// Ground heat is already handled by MOST (KERNEL_MOST_boundary3D).
						// If wall heat balance is needed later, implement PALM-style USM
						// (radiation + sensible + wall conduction → θ* feedback).

						if(P1[j].p_type<1000 && P1[j].buffer_type!=99){

							tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
							search_range=k_search_kappa*fmax(hi,hj);	// search range for multi-resoltuion(KDH)

							if(tdist<search_range){
								int_t ptypej;
								Real tmp_tdwx,tmp_tdwy,tmp_tdwz,tdwx,tdwy,tdwz,uxj,uyj,uzj,mj,tempj,rhoj,kcj,sum_con_H,diffj,concnj,tmprd,tmprd1,tmprd2;
								
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

								// 위 함수 형태를 밖으로 꺼내면
								// Real tdwxc=tdwx*P3[i].Cm[0][0]+tdwy*P3[i].Cm[0][1]+tdwz*P3[i].Cm[0][2];
								// Real tdwyc=tdwx*P3[i].Cm[1][0]+tdwy*P3[i].Cm[1][1]+tdwz*P3[i].Cm[1][2];
								// Real tdwzc=tdwx*P3[i].Cm[2][0]+tdwy*P3[i].Cm[2][1]+tdwz*P3[i].Cm[2][2];
								// KGC - CTS (KDH)
								// tdwx=gc11*tmp_tdwx+gc12*tmp_tdwy+gc31*tmp_tdwz;
								// tdwy=gc12*tmp_tdwx+gc22*tmp_tdwy+gc23*tmp_tdwz;
								// tdwz=gc31*tmp_tdwx+gc23*tmp_tdwy+gc33*tmp_tdwz;

								ptypej=P1[j].p_type;
								uxj=P1[j].ux;
								uyj=P1[j].uy;
								uzj=P1[j].uz;
								mj=P1[j].m;
								tempj=P1[j].temp;
								rhoj=P1[j].rho;
								hj=P1[j].h;

								Real dvx = uxj-uxi;
								Real dvy = uyj-uyi;
								Real dvz = uzj-uzi;

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

									wx += (mj/rhoj) * (dvy*tdwz - dvz*tdwy);
    								wy += (mj/rhoj) * (dvz*tdwx - dvx*tdwz);
    								wz += (mj/rhoj) * (dvx*tdwy - dvy*tdwx);

#if RIEMANN_EULERV
									{
									Real vb_dw=0.5*((uxi+uxj)*tdwx+(uyi+uyj)*tdwy+(uzi+uzj)*tdwz);
									Real cv=4.0*mj/(rhoi+rhoj)*vb_dw*(1-P1[i].eli);
									Real dux=uxj-uxi, duy=uyj-uyi, duz=uzj-uzi, sx, sy, sz;
									if(vb_dw<=0.0){ sx=dux; sy=duy; sz=duz; }      // j 가 upwind: 계면값 = j 쪽
									else{ sx=sy=sz=0.0; }                           // i 가 upwind: 계면값 = i 쪽
#if RIEMANN_EULERV>=2
									Real rx=xj-xi, ry=yj-yi, rz=zj-zi;
									if(vb_dw<=0.0){ if(ptypej==1){                  // j 쪽 재구성 u_j − ½∇u_j·r (경계입자는 기울기 없음 → 1차)
										sx-=0.5*RECON(P3[j].gux[0]*rx+P3[j].gux[1]*ry+P3[j].gux[2]*rz, dux);
										sy-=0.5*RECON(P3[j].guy[0]*rx+P3[j].guy[1]*ry+P3[j].guy[2]*rz, duy);
										sz-=0.5*RECON(P3[j].guz[0]*rx+P3[j].guz[1]*ry+P3[j].guz[2]*rz, duz); } }
									else{                                           // i 쪽 재구성 u_i + ½∇u_i·r
										sx=0.5*RECON(P3[i].gux[0]*rx+P3[i].gux[1]*ry+P3[i].gux[2]*rz, dux);
										sy=0.5*RECON(P3[i].guy[0]*rx+P3[i].guy[1]*ry+P3[i].guy[2]*rz, duy);
										sz=0.5*RECON(P3[i].guz[0]*rx+P3[i].guz[1]*ry+P3[i].guz[2]*rz, duz); }
#endif
#if RIEMANN_EULERV==3
									Real vc=mj/rhoj*(1-P1[i].eli), uidw=uxi*tdwx+uyi*tdwy+uzi*tdwz;   // 중앙차분 항 (RIEMANN_EULERV==0 과 같은 식)
									eulerx+=(1-UPWIND_BLEND)*uidw*dux*vc+UPWIND_BLEND*cv*sx;
									eulery+=(1-UPWIND_BLEND)*uidw*duy*vc+UPWIND_BLEND*cv*sy;
									eulerz+=(1-UPWIND_BLEND)*uidw*duz*vc+UPWIND_BLEND*cv*sz;
#else
									eulerx+=cv*sx; eulery+=cv*sy; eulerz+=cv*sz;
#endif
									}
#else
									eulerx += (uxi*(uxj-uxi)*tdwx+uyi*(uxj-uxi)*tdwy+uzi*(uxj-uxi)*tdwz)*mj/rhoj*(1-P1[i].eli);
									eulery += (uxi*(uyj-uyi)*tdwx+uyi*(uyj-uyi)*tdwy+uzi*(uyj-uyi)*tdwz)*mj/rhoj*(1-P1[i].eli);
									eulerz += (uxi*(uzj-uzi)*tdwx+uyi*(uzj-uzi)*tdwy+uzi*(uzj-uzi)*tdwz)*mj/rhoj*(1-P1[i].eli);
#endif

									if(k_con_solve && ptypej<1000 && P1[j].buffer_type!=99){
#if RIEMANN_EULERT
										// Vila 1999 conservative Eulerian advection, 1st-order upwind (Godunov).
										// eulert = u.gradT = sum 2mj/(rhoi+rhoj)(T*-Ti)(u*.dW)(1-eli),
										//   u* = interface velocity = 0.5(ui+uj),
										//   T* = upwind state: Tj if neighbour j is upwind (u*.dW<=0), else Ti.
										// Monotone by construction -> suppresses corner temperature overshoot.
										Real us_dw = 0.5*((uxi+uxj)*tdwx+(uyi+uyj)*tdwy+(uzi+uzj)*tdwz);
										Real tstar = (us_dw<=0.0) ? tempj : tempi;
										eulert += 2.0*mj/(rhoi+rhoj)*(tstar-tempi)*us_dw*(1-P1[i].eli);
#else
										eulert -= 2.0*mj/(rhoi+rhoj)*(tempj-tempi)*((uxi-uxj)*tdwx+(uyi-uyj)*tdwy+(uzi-uzj)*tdwz)*(1-P1[i].eli);
#endif
									}
									if(k_concn_solve){
										concnj=P1[j].concn;
										tmp_Rd-= (uxi*(concnj-concni)*tdwx+uyi*(concnj-concni)*tdwy+uzi*(concnj-concni)*tdwz)*mj/rhoj*(1-P1[i].eli);
									}
								}
								if(k_con_solve && P1[j].buffer_type!=99){
									//kcj=conductivity(tempj, ptypej);
									kcj=Solid_k*(ptypej==0)+Air_k*(ptypej==1||ptypej==-1);  // IBM walls adiabatic (ptypej==1000 excluded)
									Real rhoj_temp, mj_temp;
									if(ptypej==1||ptypej==-1) rhoj_temp=rhoj;
									else rhoj_temp=Solid_rho;
									mj_temp=mj*rhoj_temp/rhoj;
									sum_con_H=4.0*mj_temp*kcj*kci*(tempi-tempj)*tdwij;
									//sum_con_H/=(tdist+1e-10)*rhoi*rhoj*(kci+kcj);
									sum_con_H/=(tdist+1e-10)*rhoi*rhoj_temp*(kci+kcj);
									tmp_Rc+=sum_con_H;
								}
								if(k_turbulence_model==DDF && P1[j].buffer_type!=99){
									// SGS thermal diffusion (Brookshaw form, K/s units — does not divide by cp)
									tmp_Rc_diff += mj/rhoj*((P1[i].vis_t+P1[j].vis_t)/Pr_t)*(tempi-tempj)*tdwij/(tdist+1e-20);
								}
								if(k_con_solve){
									// Artificial thermal diffusion: damp temperature instability near building corners.
									// For converging flow (uij.xij < 0): always apply (standard Monaghan AV).
									// For diverging flow (uij.xij >= 0): also apply near IBM walls (flt_s < 0.7),
									//   where separation-point positive feedback causes eulert runaway.
									Real uij_xij_t=(uxi-uxj)*(xi-xj)+(uyi-uyj)*(yi-yj)+(uzi-uzj)*(zi-zj);
									if(uij_xij_t<0 || P1[i].flt_s < 0.7f){
										Real h_ij=(hi+hj)*0.5f;
										Real u_sig=fabsf(uij_xij_t)/(tdist+1e-20f);
										Real alpha_art=0.5f*h_ij*u_sig;  // m2/s, artificial thermal diffusivity
										// Brookshaw form with factor 2 (single-particle alpha): K/s units
										tmp_Rc_diff += 2.0f*mj/rhoj*alpha_art*(tempi-tempj)*tdwij/(tdist+1e-20f);
									}
								}
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
	if((k_boussinesq_solve)&&(ptypei>0)&&(ptypei<1000)&&(tmp_flt>1e-6f)) tmpz+=betai*Gravitational_CONST*(tmp_temp/tmp_flt-tempb);

	//P1[i].PPE1=betai*Gravitational_CONST*(tmp_temp/tmp_flt-tempb);
	P3[i].ftotalx=tmpx-eulerx;
	P3[i].ftotaly=tmpy-eulery;
	P3[i].ftotalz=tmpz-eulerz;

	// a_conv_x = uyi*wz - uzi*wy;
	// a_conv_y = uzi*wx - uxi*wz;
	// a_conv_z = uxi*wy - uyi*wx;

	// P3[i].ftotalx=tmpx + a_conv_x;
	// P3[i].ftotaly=tmpy + a_conv_y;
	// P3[i].ftotalz=tmpz + a_conv_z;

	P1[i].fe_x=-eulerx;
	P1[i].fe_y=-eulery;
	P1[i].fe_z=-eulerz;
	
	P1[i].fv_x=fvx;
	P1[i].fv_y=fvy;
	P1[i].fv_z=fvz;

	if(k_con_solve){
		//if (enthalpy_eqn) P3[i].denthalpy=tmp_Rc*(ptypei!=-1);
		if(enthalpy_eqn) P3[i].denthalpy=tmp_Rc*(ptypei!=-1);
		else{
			// tmp_Rc: conduction in W/kg form (Cleary harmonic, divided by cp -> K/s)
			// tmp_Rc_diff: SGS + artificial diffusion in K/s (Brookshaw form, no cp division)
			P3[i].dtemp+=tmp_Rc/cpi+tmp_Rc_diff-eulert;
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

	Real gc11, gc12, gc22, gc23, gc31, gc33;

	gc11=P3[i].inv_cm_xx, gc12=P3[i].inv_cm_xy;
	gc22=P3[i].inv_cm_yy, gc23=P3[i].inv_cm_yz;
	gc31=P3[i].inv_cm_zx, gc33=P3[i].inv_cm_zz;

	xi=P1[i].x;
	yi=P1[i].y;
	zi=P1[i].z;
	hi=P1[i].h;
	pi=P1[i].pres;
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
						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
						search_range=k_search_kappa*fmax(hi,hj);

						if(PPE_WALL_NB(P1,j)){
							if(tdist<search_range){
								Real tmp_tdwx,tmp_tdwy,tmp_tdwz,tdwx,tdwy,tdwz,mj,rhoj,pj;
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

								// Gradient correction 적용 (CONSISTENT_PPE: skip KGC -> plain gradW so force==IISPH operator).
								// IISPH_KGC overrides: re-apply KGC so force==KGC-IISPH operator (targets accurate div).
#if !CONSISTENT_PPE || IISPH_KGC
								apply_gradient_correction_3D(P3[i].Cm,twij,tdwx,tdwy,tdwz,&tdwx,&tdwy,&tdwz);
#endif
								// 위 함수 형태를 밖으로 꺼내면
								// Real tdwxc=tdwx*P3[i].Cm[0][0]+tdwy*P3[i].Cm[0][1]+tdwz*P3[i].Cm[0][2];
								// Real tdwyc=tdwx*P3[i].Cm[1][0]+tdwy*P3[i].Cm[1][1]+tdwz*P3[i].Cm[1][2];
								// Real tdwzc=tdwx*P3[i].Cm[2][0]+tdwy*P3[i].Cm[2][1]+tdwz*P3[i].Cm[2][2];
								// KGC - CTS (KDH)
								// tdwx=gc11*tmp_tdwx+gc12*tmp_tdwy+gc31*tmp_tdwz;
								// tdwy=gc12*tmp_tdwx+gc22*tmp_tdwy+gc23*tmp_tdwz;
								// tdwz=gc31*tmp_tdwx+gc23*tmp_tdwy+gc33*tmp_tdwz;

								mj=P1[j].m;
								rhoj=P1[j].rho;
								pj=P1[j].pres;

								if(k_fp_solve){

									// p_stab 계산  (for pressure stabilization)
	               					// Real p_stab=min((pj_max-pj_min)/2,pi)-0.0;
									//Real p_stab=min((pj_max-pj_min)/2,pi_sum/flt_pi);

                					// **C_p 계산 및 힘 적용**
									// Real C_p=-mj/(rhoi*rhoj)*(pi*rhoj+pj*rhoi)/(rhoi+rhoj);   // two phase
									// Real C_p=-mj*(pi+pj)/(rhoi*rhoj);      // Standard Monaghan
									// Real C_p=-2.0*mj/(rhoi*rhoj)*((pi+Pb)*rhoj+(pj+Pb)*rhoi)/(rhoi+rhoj);      // Monaghan variation
                					// Real C_p=-mj*(pj-pi)/(rhoi*rhoj)-mj*(2.0*p_stab)/(rhoi*rhoj)*(ptypei==ptypej);     // Pressure stabilization(KDH)
#if CONSISTENT_PPE
									// IISPH-matching symmetric form: a_p=-sum m_j(p_i/rho_i^2+p_j/rho_j^2)gradW.
									// This is exactly the force whose divergence is the IISPH operator A.
									Real C_p=-mj*(pi/(rhoi*rhoi)+pj/(rhoj*rhoj));
#else
									Real C_p=-2.0*mj/(rhoi*rhoj)*(pi*rhoj+pj*rhoi)/(rhoi+rhoj);      // Origianl (KDH, 251226)
									// Real C_p = -mj / rhoj * (pj - pi);
#endif
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
	P1[i].fp_x=tmpx;
	P1[i].fp_y=tmpy;
	P1[i].fp_z=tmpz;

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
