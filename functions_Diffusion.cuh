__global__ void KERNEL_clc_concn_transport3D(int_t*g_str,int_t*g_end,part1*P1,part2*P2,part3*P3, Real Diff_Coeff)
{
    uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
    if(i>=k_num_part2) return;
    if(P1[i].i_type>i_type_crt) return; 
    if(P1[i].p_type!=1)  return; 

    // 1. 기본 변수 로드
    Real xi=P1[i].x;  Real yi=P1[i].y;  Real zi=P1[i].z;
    Real uxi=P1[i].ux; Real uyi=P1[i].uy; Real uzi=P1[i].uz;
    Real concn_i = P1[i].concn;
    Real eli = P1[i].eli;
    Real hi=P1[i].h; 
    Real rhoi=P1[i].rho; // 밀도 필요

    Real tmp_A=calc_tmpA(hi);
    Real search_range=k_search_kappa*hi; 
    Real concn_advection = 0.0; 
    Real concn_diffusion = 0.0; 

    // Gradient Correction Matrix 로드
    Real gc11=P3[i].inv_cm_xx; Real gc12=P3[i].inv_cm_xy;
    Real gc22=P3[i].inv_cm_yy; Real gc23=P3[i].inv_cm_yz;
    Real gc31=P3[i].inv_cm_zx; Real gc33=P3[i].inv_cm_zz;

    // --- Neighbor Loop ---
    int_t icell,jcell,kcell;
    
    // Cell Indexing
    if((k_x_max==k_x_min)){icell=0;}
    else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
    if((k_y_max==k_y_min)){jcell=0;}
    else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
    if((k_z_max==k_z_min)){kcell=0;}
    else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
    
    if(icell<0) icell=0;    if(jcell<0) jcell=0;    if(kcell<0) kcell=0;

    int_t ncell=P1[i].ncell;
    for(int_t z=-ncell;z<=ncell;z++){
        for(int_t y=-ncell;y<=ncell;y++){
            for(int_t x=-ncell;x<=ncell;x++){
                int_t k=idx_cell(icell+x,jcell+y,kcell+z);
                
                if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
                
                if(g_str[k]!=cu_memset){
                    int_t fend=g_end[k];
                    for(int_t j=g_str[k];j<fend;j++){
                        if(P1[j].p_type<1000){
                            Real xj=P1[j].x; Real yj=P1[j].y; Real zj=P1[j].z;
                            Real rhoj=P1[j].rho; Real mj=P1[j].m;
                            Real concn_j = P1[j].concn;

                            // 거리 계산
                            Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
                            
                            if(tdist<search_range){
                                Real twij=calc_kernel_wij(tmp_A,hi,tdist);
                                Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
                                
                                Real tdwx=tdwij*(xi-xj)/tdist;
                                Real tdwy=tdwij*(yi-yj)/tdist;
                                Real tdwz=tdwij*(zi-zj)/tdist;

                                // Gradient Correction 적용 (GC)
                                Real tdwxc=gc11*tdwx+gc12*tdwy+gc31*tdwz;
                                Real tdwyc=gc12*tdwx+gc22*tdwy+gc23*tdwz;
                                Real tdwzc=gc31*tdwx+gc23*tdwy+gc33*tdwz;

                                // 1. Advection (이송) 계산
                                Real dc = concn_j - concn_i;
                                concn_advection += (uxi*tdwxc + uyi*tdwyc + uzi*tdwzc) * dc * (mj/rhoj);
                                
                                // 2. Diffusion (확산) 계산
                                // [수정됨] KERNEL_advection_force3D 참고하여 분모 안전장치 추가
                                Real dot_r_gradW = ((xi-xj)*tdwxc + (yi-yj)*tdwyc + (zi-zj)*tdwzc);
                                
                                Real denominator = tdist*tdist + 0.01*hi*hi;

                                // Laplacian Term (Cleary & Monaghan Type approximation)
                                Real geom_term = 2.0 * (mj/rhoj) * dot_r_gradW / denominator;

                                Real concn_j_eff = concn_j;
                                if(P1[j].p_type!=1) concn_j_eff = concn_i; // 경계면 입자 처리
                                
                                concn_diffusion += (concn_i - concn_j_eff) * geom_term;
                                Real diff_term = (concn_i - concn_j_eff) * geom_term;

                                if(isnan(diff_term)) printf("concn_j=%f, concn_i=%f, geom_term=%f, dot_r_gradW=%f, mj=%f, rhoj=%f, tdist=%f \n",concn_j,concn_i,geom_term,dot_r_gradW,mj,rhoj,tdist);
                            }
                        }
                    }
                }
            }
        }
    }

    // 최종 결과 저장 (확산 - 이송)
    P3[i].dconcn = (Diff_Coeff * concn_diffusion) - concn_advection * (1.0 - eli);
}


__global__ void KERNEL_apply_source_concn(part1* P1, Real tdt, Real cx, Real cy, Real cz, Real radius, Real intensity)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;

    // 1. 유효한 유체 입자인지 확인
    if (P1[i].i_type>i_type_crt) return; 
    if (P1[i].p_type>=1000) return; // 벽면 입자 제외

    // 2. 입자와 소스 중심 사이의 거리 계산 (Squared distance)
    Real dx = P1[i].x-cx;
    Real dy = P1[i].y-cy;
    Real dz = P1[i].z-cz;
    Real dist2 = dx*dx+dy*dy+dz*dz;

    // 3. 반경 안에 있는 경우 농도 추가
    if (dist2 <= radius*radius){
        P1[i].concn += intensity;
        P1[i].source = 1;
        printf("Source applied to particle %d: new concn = %f\n", i, P1[i].concn);
        // [옵션] 만약 농도가 특정 값 이상으로 커지는 것을 막고 싶다면(Saturation):
        // if(P1[i].concn > 1.0) P1[i].concn = 1.0; 
    }
}

__global__ void KERNEL_update_concn(part1* P1, part3* P3, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return; // 필요한 경우 조건 추가

    // Euler 적분: C_new = C_old + (dC/dt) * dt
    P1[i].concn += P3[i].dconcn * tdt;

    // (옵션) 농도가 음수가 되는 것을 방지 (수치적 오차 보정)
    if (P1[i].concn < 0.0) P1[i].concn = 0.0;
}

__global__ void KERNEL_clc_concn_LDM(int_t*g_str_LDM, int_t*g_end_LDM, part1* P1, L_part1* LP1)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return; // 필요한 경우 조건 추가
    if (P1[i].p_type != 1) return; // 오직 오염물질 입자에 대해서만 계산

    // 1. 기본 변수 로드
    Real xi=P1[i].x;
    Real yi=P1[i].y;
    Real zi=P1[i].z;
    Real uxi=P1[i].ux;
    Real uyi=P1[i].uy;
    Real uzi=P1[i].uz;
    Real hi=2.0*P1[i].h; 
    Real rhoi=P1[i].rho; // 밀도 필요

    Real tmp_A=calc_tmpA(hi);
    Real search_range=k_search_kappa*hi; 

    Real tmp_concn = 0.0;

    // --- Neighbor Loop ---
    int_t icell,jcell,kcell;
    
    // Cell Indexing
    if((k_x_max==k_x_min)){icell=0;}
    else{icell=min(floor((xi-k_x_min)/k_dcell),k_NI-1);}
    if((k_y_max==k_y_min)){jcell=0;}
    else{jcell=min(floor((yi-k_y_min)/k_dcell),k_NJ-1);}
    if((k_z_max==k_z_min)){kcell=0;}
    else{kcell=min(floor((zi-k_z_min)/k_dcell),k_NK-1);}
    
    if(icell<0) icell=0;    if(jcell<0) jcell=0;    if(kcell<0) kcell=0;

    int_t ncell=P1[i].ncell;
    for(int_t z=-ncell;z<=ncell;z++){
        for(int_t y=-ncell;y<=ncell;y++){
            for(int_t x=-ncell;x<=ncell;x++){
                int_t k=idx_cell(icell+x,jcell+y,kcell+z);
                
                if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
                
                if(g_str_LDM[k]!=cu_memset){
                    int_t fend=g_end_LDM[k];
                    for(int_t j=g_str_LDM[k];j<fend;j++){
                        if(P1[j].p_type<1000){
                            Real xj=LP1[j].x;
                            Real yj=LP1[j].y;
                            Real zj=LP1[j].z;
                            Real rhoj=LP1[j].rho;
                            Real mj=LP1[j].m;
                            
                            // 거리 계산
                            Real tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
                            
                            if(tdist<search_range){
                                Real twij=calc_kernel_wij(tmp_A,hi,tdist);

                                tmp_concn += twij * (mj/rhoj);

                            }
                        }
                    }
                }
            }
        }
    }
    P1[i].concn_LDM = tmp_concn;
}