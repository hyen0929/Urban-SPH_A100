// ------------------------------
// Trilinear interpolation at query point (FDM-style, ~8 particles)
// Uses tent function weights: w = max(0, 1-|dx|/sp) * max(0, 1-|dy|/sp) * max(0, 1-|dz|/sp)
// where sp = particle spacing = hi / h_coeff
// Numerical diffusion: O(sp^2) per interpolation
// ------------------------------
__device__ __forceinline__ Real interp_C_at_trilinear(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, Real hi, int ncell)
{
    Real sp = hi / h_coeff;   // particle spacing
    Real search_range = 1.01 * sp;  // just over 1*sp to catch the 8 nearest

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;

    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real dx = fabs(xq - P1[j].x);
                    Real dy = fabs(yq - P1[j].y);
                    Real dz = fabs(zq - P1[j].z);

                    if(dx < search_range && dy < search_range && dz < search_range){
                        Real wx = 1.0 - dx / sp;
                        Real wy = 1.0 - dy / sp;
                        Real wz = 1.0 - dz / sp;
                        Real w = wx * wy * wz;

                        num += w * P1[j].concn;
                        den += w;
                    }
                }
            }
        }
    }

    return (den > 0.0) ? (num/den) : 0.0;
}

// ------------------------------
// Tensor-product cubic B-spline interpolation (~64 particles on regular grid)
// 1D basis M4(q): q=|dx|/sp
//   q<1: (4 - 6q^2 + 3q^3) / 6
//   q<2: (2-q)^3 / 6
// 3D weight: M4(qx) * M4(qy) * M4(qz)
// Numerical diffusion: O(sp^4) per interpolation — 100x less than trilinear
// ------------------------------
__device__ __forceinline__ Real interp_C_at_cubic(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, Real hi, int ncell)
{
    Real sp = hi / h_coeff;   // particle spacing

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;

    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real qx = fabs(xq - P1[j].x) / sp;
                    Real qy = fabs(yq - P1[j].y) / sp;
                    Real qz = fabs(zq - P1[j].z) / sp;

                    if(qx < 2.0 && qy < 2.0 && qz < 2.0){
                        // 1D cubic B-spline basis (tensor product)
                        Real wx, wy, wz;

                        if(qx < 1.0) wx = (4.0 - 6.0*qx*qx + 3.0*qx*qx*qx) / 6.0;
                        else          wx = (2.0 - qx)*(2.0 - qx)*(2.0 - qx) / 6.0;

                        if(qy < 1.0) wy = (4.0 - 6.0*qy*qy + 3.0*qy*qy*qy) / 6.0;
                        else          wy = (2.0 - qy)*(2.0 - qy)*(2.0 - qy) / 6.0;

                        if(qz < 1.0) wz = (4.0 - 6.0*qz*qz + 3.0*qz*qz*qz) / 6.0;
                        else          wz = (2.0 - qz)*(2.0 - qz)*(2.0 - qz) / 6.0;

                        Real w = wx * wy * wz;
                        num += w * P1[j].concn;
                        den += w;
                    }
                }
            }
        }
    }

    return (den > 0.0) ? (num/den) : 0.0;
}

// ------------------------------
// Interp scalar at query point (normalized)
// ------------------------------
__device__ __forceinline__ Real interp_C_at(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, Real hi, int ncell)
{
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;

    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;

                    if(r < search_range){
                        Real w = calc_kernel_wij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;  // volume weight
                        num += V * P1[j].concn * w;
                        den += V * w;
                    }
                }
            }
        }
    }

    return (den > 0.0) ? (num/den) : 0.0;
}

// ------------------------------
// Interp scalar at query point with 1st order gradient reconstruction
// ------------------------------
__device__ __forceinline__ Real interp_C_at_1st(Real xq, Real yq, Real zq, int_t*g_str, int_t*g_end, part1*P1, Real hi, int ncell)
{
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;
    Real grad_x = 0.0, grad_y = 0.0, grad_z = 0.0;
    Real xc = 0.0, yc = 0.0, zc = 0.0;  // Weighted center

    // Pass 1: 0th order interpolation + weighted center
    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;

                    if(r < search_range){
                        Real w = calc_kernel_wij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        Real Vw = V * w;

                        num += Vw * P1[j].concn;
                        den += Vw;
                        xc += Vw * P1[j].x;
                        yc += Vw * P1[j].y;
                        zc += Vw * P1[j].z;
                    }
                }
            }
        }
    }

    if(den <= 0.0) return 0.0;

    Real C0 = num / den;
    xc /= den;  yc /= den;  zc /= den;

    // Pass 2: Gradient calculation using SPH gradient
    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;

                    if(r < search_range){
                        Real dwij = calc_kernel_dwij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        Real dC = P1[j].concn - C0;

                        grad_x += V * dC * dwij * dx / r;
                        grad_y += V * dC * dwij * dy / r;
                        grad_z += V * dC * dwij * dz / r;
                    }
                }
            }
        }
    }

    // Apply 1st order correction: C = C0 + grad·(xq - xc)
    Real dx_q = xq - xc;
    Real dy_q = yq - yc;
    Real dz_q = zq - zc;

    return C0 + grad_x * dx_q + grad_y * dy_q + grad_z * dz_q;
}

// ------------------------------
// Interp velocity at query point (normalized)  // Changed
// ------------------------------
__device__ __forceinline__ void interp_u_at(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, Real hi, Real& uqx, Real& uqy, Real& uqz, int ncell)
{
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real numx=0.0, numy=0.0, numz=0.0, den=0.0;

    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;

                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;

                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;

                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;

                    if(r < search_range){
                        Real w = calc_kernel_wij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        numx += V * P1[j].ux * w;
                        numy += V * P1[j].uy * w;
                        numz += V * P1[j].uz * w;
                        den  += V * w;
                    }
                }
            }
        }
    }

    if(den > 0.0){
        uqx = numx/den; uqy = numy/den; uqz = numz/den;
    } else {
        uqx = 0.0; uqy = 0.0; uqz = 0.0;
    }
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

__global__ void KERNEL_apply_source_concn_Gaussian(part1* P1, Real tdt, Real cx, Real cy, Real cz, Real radius, Real intensity)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;

    if (P1[i].i_type>i_type_crt) return;
    if (P1[i].p_type>=1000) return;

    Real dx = P1[i].x-cx;
    Real dy = P1[i].y-cy;
    Real dz = P1[i].z-cz;
    Real dist = sqrt(dx*dx+dy*dy+dz*dz)+1e-20;

    // radius = 커널 compact support 반경, h_src = smoothing length
    Real h_src = radius/k_search_kappa;
    h_src = 2.0*P1[i].h;
    Real tmp_A = calc_tmpA(h_src);
    Real search_range = k_search_kappa * h_src;

    if (dist < search_range){
        Real wij = calc_kernel_wij(tmp_A, h_src, dist);
        // SPH 정규화: Σ_j V_j·W ≈ 1 → 총 스칼라 질량 ≈ intensity
        P1[i].concn += intensity * wij;
        P1[i].source = 1;
    }
}

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
    Real dCdx=0.0, dCdy=0.0, dCdz=0.0; // Changed (루프 밖 초기화)

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

                                // // Gradient Correction 적용 (GC)
                                // Real tdwxc=gc11*tdwx+gc12*tdwy+gc31*tdwz;
                                // Real tdwyc=gc12*tdwx+gc22*tdwy+gc23*tdwz;
                                // Real tdwzc=gc31*tdwx+gc23*tdwy+gc33*tdwz;

                                Real tdwxc=tdwx;
                                Real tdwyc=tdwy;
                                Real tdwzc=tdwz;

                                // 1. Advection (이송) 계산
                                Real uxj = P1[j].ux;  Real uyj = P1[j].uy;  Real uzj = P1[j].uz; // Changed

                                // (c u) 벡터의 차이
                                Real fx = concn_j*uxj - concn_i*uxi; // Changed
                                Real fy = concn_j*uyj - concn_i*uyi; // Changed
                                Real fz = concn_j*uzj - concn_i*uzi; // Changed

                                Real dC = concn_j - concn_i;
                                Real w = (mj/rhoj);
                                dCdx += dC * tdwxc * w;
                                dCdy += dC * tdwyc * w;
                                dCdz += dC * tdwzc * w;

                                
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

    // 루프 종료 후 gradient -> advection 계산
    concn_advection = uxi*dCdx + uyi*dCdy + uzi*dCdz;

    // 최종 결과 저장 (확산 - 이송)
    P3[i].dconcn = (Diff_Coeff * concn_diffusion) - concn_advection;
}

// ------------------------------
// Semi-Lagrangian advection: C_adv = C(x - u dt)
// ------------------------------
__global__ void KERNEL_advect_concn_SL(int_t* g_str,int_t* g_end,part1* P1, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i>=k_num_part2) return;
    if(P1[i].i_type>i_type_crt) return;
    if(P1[i].p_type!=1) return;

    Real xi=P1[i].x, yi=P1[i].y, zi=P1[i].z;
    Real uxi=P1[i].ux, uyi=P1[i].uy, uzi=P1[i].uz;
    Real hi=P1[i].h;
    Real concn_n = P1[i].concn;
    int_t ncell = P1[i].ncell;

    // Forward SL backtrace + 1st order interpolation
    Real xb = xi - uxi*tdt;
    Real yb = yi - uyi*tdt;
    Real zb = zi - uzi*tdt;

    Real C_advected = interp_C_at_1st(xb, yb, zb, g_str, g_end, P1, hi, ncell);
    if(C_advected < 0.0) C_advected = 0.0;

    // Diffusion (Cleary-Monaghan Laplacian)
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;
    Real diffusion_sum = 0.0;

    int_t icell,jcell,kcell;
    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xi-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yi-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zi-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii=icell+x, jj=jcell+y, kk=kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t k=idx_cell(ii,jj,kk);
                if(g_str[k]==cu_memset) continue;
                for(int_t j=g_str[k]; j<g_end[k]; j++){
                    if(P1[j].p_type>=1000) continue;
                    Real rx=xi-P1[j].x, ry=yi-P1[j].y, rz=zi-P1[j].z;
                    Real tdist=sqrt(rx*rx+ry*ry+rz*rz)+1e-20;
                    if(tdist<search_range){
                        Real tdwij=calc_kernel_dwij(tmp_A,hi,tdist);
                        Real dwx=tdwij*rx/tdist, dwy=tdwij*ry/tdist, dwz=tdwij*rz/tdist;
                        Real concn_j=(P1[j].p_type!=1) ? concn_n : P1[j].concn;
                        Real dot_r_gradW=rx*dwx+ry*dwy+rz*dwz;
                        Real denom=tdist*tdist+0.01*hi*hi;
                        diffusion_sum+=(concn_n-concn_j)*2.0*(P1[j].m/P1[j].rho)*dot_r_gradW/denom;
                    }
                }
            }
        }
    }

    Real C_final = C_advected + tdt * Diff_Coeff * diffusion_sum;
    if(C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
}

// MacCormack 방식을 적용한 고품질 이송 함수
// (기존 interp_C_at 함수가 정의되어 있어야 합니다)
__global__ void KERNEL_advect_concn_MacCormack(int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return;
    if (P1[i].p_type != 1) return;

    Real xi = P1[i].x;   Real yi = P1[i].y;   Real zi = P1[i].z;
    Real uxi = P1[i].ux; Real uyi = P1[i].uy; Real uzi = P1[i].uz;
    Real hi = P1[i].h;
    Real concn_n = P1[i].concn; // 현재 농도

    // ---------------------------------------------------------
    // 1단계: Forward Advection (Predictor)
    // "과거 위치(xb)에서 값을 가져옴" -> 일반적인 Semi-Lagrangian
    // ---------------------------------------------------------
    Real xb = xi - uxi * tdt;
    Real yb = yi - uyi * tdt;
    Real zb = zi - uzi * tdt;
    
    // C_star: 1차 예측값 (약간 뭉개진 값)
    Real C_star = interp_C_at(xb, yb, zb, g_str, g_end, P1, 0.65*hi, P1[i].ncell);

    // ---------------------------------------------------------
    // 2단계: Backward Advection (Corrector)
    // "예측된 위치에서, 다시 미래(현재)로 오면 원래 값이 나와야 함"
    // 원래는 미래 필드(C_star field)에서 가져와야 하지만, 
    // 근사적으로 현재 필드(P1)에서 미래 위치(xf)를 조회하여 오차를 추정합니다.
    // ---------------------------------------------------------
    Real xf = xi + uxi * tdt;
    Real yf = yi + uyi * tdt;
    Real zf = zi + uzi * tdt;

    // C_reverse: 1차 예측값을 다시 되돌렸을 때의 값 (2번 뭉개진 값)
    Real C_reverse = interp_C_at(xf, yf, zf, g_str, g_end, P1, 0.65*hi, P1[i].ncell);

    // ---------------------------------------------------------
    // 3단계: Error Correction & Clamping
    // 오차 = (원래값 - 되돌려진값) / 2
    // ---------------------------------------------------------
    Real correction = (concn_n - C_reverse) * 0.5;
    Real C_final = C_star + correction;

    // [중요] 오버슈팅 방지 (Clamp)
    // 보정 과정에서 값이 튀어(음수 or 최대값 초과) 불안정해지는 것을 막습니다.
    // 여기서는 간단히 Min/Max Limit을 주변 값 기준으로 잡아야 하지만,
    // 약식으로 음수 방지 및 기존 값 범위 내 유지로 처리합니다.
    
    if (C_final < 0.0) C_final = 0.0;
    
    // (선택) 보정이 너무 과하면 원래 1차 예측값 사용 (안정성 확보)
    // if (fabs(correction) > fabs(C_star) * 0.5) C_final = C_star;

    P3[i].dconcn = C_final;
}

// =====================================================================
// MacCormack advection + Cleary-Monaghan diffusion (combined)
// dC/dt = -u·∇C + D∇²C
// Operator splitting: advection first, then diffusion
// =====================================================================
__global__ void KERNEL_Scalar_transport_MacCormack(int_t* g_str, int_t* g_end, part1* P1, part2* P2, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return;
    if (P1[i].p_type != 1) return;

    Real xi = P1[i].x;   Real yi = P1[i].y;   Real zi = P1[i].z;
    Real uxi = P1[i].ux; Real uyi = P1[i].uy; Real uzi = P1[i].uz;
    Real hi = P1[i].h;
    Real concn_n = P1[i].concn; // Current concentration

    // =========================================================
    // STEP 1: MacCormack Advection
    // =========================================================

    // 1-1. Forward Advection (Predictor): go back in time
    Real xb = xi - uxi * tdt;
    Real yb = yi - uyi * tdt;
    Real zb = zi - uzi * tdt;

    Real C_star = interp_C_at(xb, yb, zb, g_str, g_end, P1, hi, P1[i].ncell);

    // 1-2. Backward Advection (Corrector): estimate error
    Real xf = xi + uxi * tdt;
    Real yf = yi + uyi * tdt;
    Real zf = zi + uzi * tdt;

    Real C_reverse = interp_C_at(xf, yf, zf, g_str, g_end, P1, hi, P1[i].ncell);

    // 1-3. Error Correction
    Real correction = (concn_n - C_reverse) * 0.5;
    Real C_advected = C_star + correction;

    // Clamp to prevent negative
    if (C_advected < 0.0) C_advected = 0.0;

    // =========================================================
    // STEP 2: Diffusion (Cleary & Monaghan Laplacian)
    // Apply diffusion to the current concentration field
    // =========================================================

    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;
    Real diffusion_sum = 0.0;

    int_t icell, jcell, kcell;
    if (k_x_max == k_x_min) icell = 0; else icell = min((int_t)floor((xi - k_x_min) / k_dcell), k_NI - 1);
    if (k_y_max == k_y_min) jcell = 0; else jcell = min((int_t)floor((yi - k_y_min) / k_dcell), k_NJ - 1);
    if (k_z_max == k_z_min) kcell = 0; else kcell = min((int_t)floor((zi - k_z_min) / k_dcell), k_NK - 1);
    if (icell < 0) icell = 0; if (jcell < 0) jcell = 0; if (kcell < 0) kcell = 0;

    int_t ncell = P1[i].ncell;
    for (int_t z = -ncell; z <= ncell; z++) {
        for (int_t y = -ncell; y <= ncell; y++) {
            for (int_t x = -ncell; x <= ncell; x++) {
                int_t k = idx_cell(icell + x, jcell + y, kcell + z);
                if ((icell + x) < 0 || (icell + x) > (k_NI - 1) || (jcell + y) < 0 || (jcell + y) > (k_NJ - 1) || (kcell + z) < 0 || (kcell + z) > (k_NK - 1)) continue;
                if (g_str[k] != cu_memset) {
                    int_t fend = g_end[k];
                    for (int_t j = g_str[k]; j < fend; j++) {
                        if (P1[j].p_type < 1000) {
                            Real xj = P1[j].x, yj = P1[j].y, zj = P1[j].z;
                            Real mj = P1[j].m, rhoj = P1[j].rho;
                            Real concn_j = P1[j].concn;
                            Real vol_j = mj / rhoj;

                            Real rx = xi - xj, ry = yi - yj, rz = zi - zj;
                            Real tdist = sqrt(rx * rx + ry * ry + rz * rz) + 1e-20;

                            if (tdist < search_range) {
                                Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
                                Real tdwx = tdwij * rx / tdist;
                                Real tdwy = tdwij * ry / tdist;
                                Real tdwz = tdwij * rz / tdist;

                                // Diffusion term
                                Real concn_j_eff = (P1[j].p_type != 1) ? concn_n : concn_j;
                                Real dot_r_gradW = rx * tdwx + ry * tdwy + rz * tdwz;
                                Real denom = tdist * tdist + 0.01 * hi * hi;
                                diffusion_sum += (concn_n - concn_j_eff) * 2.0 * vol_j * dot_r_gradW / denom;
                            }
                        }
                    }
                }
            }
        }
    }

    // =========================================================
    // STEP 3: Combine advection + diffusion
    // =========================================================
    Real C_final = C_advected + tdt * Diff_Coeff * diffusion_sum;

    // Final clamp
    if (C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
}

// =====================================================================
// MacCormack advection with RK2 back-tracing + Cleary-Monaghan diffusion
// RK2 provides more accurate departure point tracking
// =====================================================================
__global__ void KERNEL_Scalar_transport_MacCormack_RK2(int_t* g_str, int_t* g_end, part1* P1, part2* P2, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return;
    if (P1[i].p_type != 1) return;

    Real xi = P1[i].x;   Real yi = P1[i].y;   Real zi = P1[i].z;
    Real uxi = P1[i].ux; Real uyi = P1[i].uy; Real uzi = P1[i].uz;
    Real hi = P1[i].h;
    Real concn_n = P1[i].concn;

    // CFL-based adaptive kernel radius
    Real u_mag = sqrt(uxi*uxi + uyi*uyi + uzi*uzi);
    Real CFL = u_mag * tdt / hi;
    Real kernel_factor = 0.7 + min(CFL * 0.6, 0.4);  // Range: 0.7~1.1
    Real interp_radius = kernel_factor * hi;

    // =========================================================
    // STEP 1: MacCormack Advection with RK2 Back-tracing
    // =========================================================

    // ----------------------
    // 1-1. Forward Advection (Predictor) with RK2
    // ----------------------
    
    // RK2 Stage 1: Half-step backward to get midpoint
    Real x_mid_b = xi - uxi * (tdt * 0.5);
    Real y_mid_b = yi - uyi * (tdt * 0.5);
    Real z_mid_b = zi - uzi * (tdt * 0.5);

    // RK2 Stage 2: Interpolate velocity at midpoint
    Real u_mid_bx, u_mid_by, u_mid_bz;
    interp_u_at(x_mid_b, y_mid_b, z_mid_b, g_str, g_end, P1, hi, u_mid_bx, u_mid_by, u_mid_bz, P1[i].ncell);

    // RK2 Stage 3: Full-step backward using midpoint velocity
    Real xb = xi - u_mid_bx * tdt;
    Real yb = yi - u_mid_by * tdt;
    Real zb = zi - u_mid_bz * tdt;

    // Interpolate concentration at departure point
    Real C_star = interp_C_at(xb, yb, zb, g_str, g_end, P1, interp_radius, P1[i].ncell);

    // ----------------------
    // 1-2. Backward Advection (Corrector) with RK2
    // ----------------------
    
    // RK2 Stage 1: Half-step forward to get midpoint
    Real x_mid_f = xi + uxi * (tdt * 0.5);
    Real y_mid_f = yi + uyi * (tdt * 0.5);
    Real z_mid_f = zi + uzi * (tdt * 0.5);

    // RK2 Stage 2: Interpolate velocity at midpoint
    Real u_mid_fx, u_mid_fy, u_mid_fz;
    interp_u_at(x_mid_f, y_mid_f, z_mid_f, g_str, g_end, P1, hi, u_mid_fx, u_mid_fy, u_mid_fz, P1[i].ncell);

    // RK2 Stage 3: Full-step forward using midpoint velocity
    Real xf = xi + u_mid_fx * tdt;
    Real yf = yi + u_mid_fy * tdt;
    Real zf = zi + u_mid_fz * tdt;

    // Interpolate concentration at forward point
    Real C_reverse = interp_C_at(xf, yf, zf, g_str, g_end, P1, interp_radius, P1[i].ncell);

    // ----------------------
    // 1-3. Error Correction
    // ----------------------
    Real correction = (concn_n - C_reverse) * 0.75;
    Real C_advected = C_star + correction;

    // Clamp to prevent negative
    if (C_advected < 0.0) C_advected = 0.0;

    // =========================================================
    // STEP 2: Diffusion (Cleary & Monaghan Laplacian)
    // =========================================================

    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;
    Real diffusion_sum = 0.0;

    int_t icell, jcell, kcell;
    if (k_x_max == k_x_min) icell = 0; else icell = min((int_t)floor((xi - k_x_min) / k_dcell), k_NI - 1);
    if (k_y_max == k_y_min) jcell = 0; else jcell = min((int_t)floor((yi - k_y_min) / k_dcell), k_NJ - 1);
    if (k_z_max == k_z_min) kcell = 0; else kcell = min((int_t)floor((zi - k_z_min) / k_dcell), k_NK - 1);
    if (icell < 0) icell = 0; if (jcell < 0) jcell = 0; if (kcell < 0) kcell = 0;

    int_t ncell = P1[i].ncell;
    for (int_t z = -ncell; z <= ncell; z++) {
        for (int_t y = -ncell; y <= ncell; y++) {
            for (int_t x = -ncell; x <= ncell; x++) {
                int_t k = idx_cell(icell + x, jcell + y, kcell + z);
                if ((icell + x) < 0 || (icell + x) > (k_NI - 1) || (jcell + y) < 0 || (jcell + y) > (k_NJ - 1) || (kcell + z) < 0 || (kcell + z) > (k_NK - 1)) continue;
                if (g_str[k] != cu_memset) {
                    int_t fend = g_end[k];
                    for (int_t j = g_str[k]; j < fend; j++) {
                        if (P1[j].p_type < 1000) {
                            Real xj = P1[j].x, yj = P1[j].y, zj = P1[j].z;
                            Real mj = P1[j].m, rhoj = P1[j].rho;
                            Real concn_j = P1[j].concn;
                            Real vol_j = mj / rhoj;

                            Real rx = xi - xj, ry = yi - yj, rz = zi - zj;
                            Real tdist = sqrt(rx * rx + ry * ry + rz * rz) + 1e-20;

                            if (tdist < search_range) {
                                Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
                                Real tdwx = tdwij * rx / tdist;
                                Real tdwy = tdwij * ry / tdist;
                                Real tdwz = tdwij * rz / tdist;

                                // Diffusion term
                                Real concn_j_eff = (P1[j].p_type != 1) ? concn_n : concn_j;
                                Real dot_r_gradW = rx * tdwx + ry * tdwy + rz * tdwz;
                                Real denom = tdist * tdist + 0.01 * hi * hi;
                                diffusion_sum += (concn_n - concn_j_eff) * 2.0 * vol_j * dot_r_gradW / denom;
                            }
                        }
                    }
                }
            }
        }
    }

    // =========================================================
    // STEP 3: Combine advection + diffusion
    // =========================================================
    Real C_final = C_advected + tdt * Diff_Coeff * diffusion_sum;

    // Final clamp
    if (C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
}

// =====================================================================
// MacCormack advection with RK2 back-tracing + 1st order interpolation + Cleary-Monaghan diffusion
// Uses gradient reconstruction for reduced numerical diffusion
// =====================================================================
__global__ void KERNEL_Scalar_transport_MacCormack_RK2_1st(int_t* g_str, int_t* g_end, part1* P1, part2* P2, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part2) return;
    if (P1[i].i_type > i_type_crt) return;
    if (P1[i].p_type != 1) return;

    Real xi = P1[i].x;   Real yi = P1[i].y;   Real zi = P1[i].z;
    Real uxi = P1[i].ux; Real uyi = P1[i].uy; Real uzi = P1[i].uz;
    Real hi = P1[i].h;
    Real concn_n = P1[i].concn;

    // CFL-based adaptive kernel radius
    Real u_mag = sqrt(uxi*uxi + uyi*uyi + uzi*uzi);
    Real CFL = u_mag * tdt / hi;
    Real kernel_factor = 0.7 + min(CFL * 0.6, 0.4);  // Range: 0.7~1.1
    kernel_factor=0.6;
    Real interp_radius = kernel_factor * hi;

    // =========================================================
    // STEP 1: MacCormack Advection with RK2 Back-tracing
    // =========================================================

    // ----------------------
    // 1-1. Forward Advection (Predictor) with RK2
    // ----------------------

    // RK2 Stage 1: Half-step backward to get midpoint
    Real x_mid_b = xi - uxi * (tdt * 0.5);
    Real y_mid_b = yi - uyi * (tdt * 0.5);
    Real z_mid_b = zi - uzi * (tdt * 0.5);

    // RK2 Stage 2: Interpolate velocity at midpoint
    Real u_mid_bx, u_mid_by, u_mid_bz;
    interp_u_at(x_mid_b, y_mid_b, z_mid_b, g_str, g_end, P1, hi, u_mid_bx, u_mid_by, u_mid_bz, P1[i].ncell);

    // RK2 Stage 3: Full-step backward using midpoint velocity
    Real xb = xi - u_mid_bx * tdt;
    Real yb = yi - u_mid_by * tdt;
    Real zb = zi - u_mid_bz * tdt;

    // Interpolate concentration at departure point with 1st order
    Real C_star = interp_C_at_1st(xb, yb, zb, g_str, g_end, P1, interp_radius, P1[i].ncell);

    // ----------------------
    // 1-2. Backward Advection (Corrector) with RK2
    // ----------------------

    // RK2 Stage 1: Half-step forward to get midpoint
    Real x_mid_f = xi + uxi * (tdt * 0.5);
    Real y_mid_f = yi + uyi * (tdt * 0.5);
    Real z_mid_f = zi + uzi * (tdt * 0.5);

    // RK2 Stage 2: Interpolate velocity at midpoint
    Real u_mid_fx, u_mid_fy, u_mid_fz;
    interp_u_at(x_mid_f, y_mid_f, z_mid_f, g_str, g_end, P1, hi, u_mid_fx, u_mid_fy, u_mid_fz, P1[i].ncell);

    // RK2 Stage 3: Full-step forward using midpoint velocity
    Real xf = xi + u_mid_fx * tdt;
    Real yf = yi + u_mid_fy * tdt;
    Real zf = zi + u_mid_fz * tdt;

    // Interpolate concentration at forward point with 1st order
    Real C_reverse = interp_C_at_1st(xf, yf, zf, g_str, g_end, P1, interp_radius, P1[i].ncell);

    // ----------------------
    // 1-3. Error Correction
    // ----------------------
    Real correction = (concn_n - C_reverse) * 0.0;
    Real C_advected = C_star + correction;

    // Clamp to prevent negative
    if (C_advected < 0.0) C_advected = 0.0;

    // =========================================================
    // STEP 2: Diffusion (Cleary & Monaghan Laplacian)
    // =========================================================

    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;
    Real diffusion_sum = 0.0;

    int_t icell, jcell, kcell;
    if (k_x_max == k_x_min) icell = 0; else icell = min((int_t)floor((xi - k_x_min) / k_dcell), k_NI - 1);
    if (k_y_max == k_y_min) jcell = 0; else jcell = min((int_t)floor((yi - k_y_min) / k_dcell), k_NJ - 1);
    if (k_z_max == k_z_min) kcell = 0; else kcell = min((int_t)floor((zi - k_z_min) / k_dcell), k_NK - 1);
    if (icell < 0) icell = 0; if (jcell < 0) jcell = 0; if (kcell < 0) kcell = 0;

    int_t ncell = P1[i].ncell;
    for (int_t z = -ncell; z <= ncell; z++) {
        for (int_t y = -ncell; y <= ncell; y++) {
            for (int_t x = -ncell; x <= ncell; x++) {
                int_t k = idx_cell(icell + x, jcell + y, kcell + z);
                if ((icell + x) < 0 || (icell + x) > (k_NI - 1) || (jcell + y) < 0 || (jcell + y) > (k_NJ - 1) || (kcell + z) < 0 || (kcell + z) > (k_NK - 1)) continue;
                if (g_str[k] != cu_memset) {
                    int_t fend = g_end[k];
                    for (int_t j = g_str[k]; j < fend; j++) {
                        if (P1[j].p_type < 1000) {
                            Real xj = P1[j].x, yj = P1[j].y, zj = P1[j].z;
                            Real mj = P1[j].m, rhoj = P1[j].rho;
                            Real concn_j = P1[j].concn;
                            Real vol_j = mj / rhoj;

                            Real rx = xi - xj, ry = yi - yj, rz = zi - zj;
                            Real tdist = sqrt(rx * rx + ry * ry + rz * rz) + 1e-20;

                            if (tdist < search_range) {
                                Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
                                Real tdwx = tdwij * rx / tdist;
                                Real tdwy = tdwij * ry / tdist;
                                Real tdwz = tdwij * rz / tdist;

                                // Diffusion term
                                Real concn_j_eff = (P1[j].p_type != 1) ? concn_n : concn_j;
                                Real dot_r_gradW = rx * tdwx + ry * tdwy + rz * tdwz;
                                Real denom = tdist * tdist + 0.01 * hi * hi;
                                diffusion_sum += (concn_n - concn_j_eff) * 2.0 * vol_j * dot_r_gradW / denom;
                            }
                        }
                    }
                }
            }
        }
    }

    // =========================================================
    // STEP 3: Combine advection + diffusion
    // =========================================================
    Real C_final = C_advected + tdt * Diff_Coeff * diffusion_sum;

    // Final clamp
    if (C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
}

__global__ void KERNEL_copy_concn_adv(part1* P1, part3* P3)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i>=k_num_part2) return;
    if(P1[i].i_type>i_type_crt) return;
    if(P1[i].p_type!=1) return;
    P1[i].concn = P3[i].dconcn;
}

// =====================================================================
// Proper BiMacCormack: two-pass split kernels
//
//   Pass 1  KERNEL_BiMC_forward   : C_star = SL(C^n, x - u*dt)
//                                   result → P3->concn_adv
//                                   (cudaDeviceSynchronize between passes)
//   Pass 2  KERNEL_BiMC_correct   : reads C_star field from P3->concn_adv
//                                   C_reverse = SL(C_star, x + u*dt)
//                                   C_mc = C_star + 0.5*(C^n - C_reverse)
//                                   monotone limiter (Selle 2008 style)
//                                   + Cleary-Monaghan diffusion
//                                   result → P3->dconcn
//
//   Calling sequence (ISPH_Calc):
//     KERNEL_BiMC_forward<<<b,t>>>(g_str,g_end,dev_P1,dev_P3,dt);
//     cudaDeviceSynchronize();
//     KERNEL_BiMC_correct<<<b,t>>>(g_str,g_end,dev_P1,dev_P3,diff_coeff,dt);
//     cudaDeviceSynchronize();
//     KERNEL_copy_concn_adv<<<b,t>>>(dev_P1,dev_P3);
//     cudaDeviceSynchronize();
// =====================================================================

// Helper: SPH-weighted interpolation from P3->concn_adv (the C_star field)
__device__ __forceinline__ Real interp_Cstar_at(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real hi, int ncell)
{
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    int_t icell,jcell,kcell;

    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;

    for(int_t z = -ncell; z <= ncell; z++){
        for(int_t y = -ncell; y <= ncell; y++){
            for(int_t x = -ncell; x <= ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t cell = idx_cell(ii, jj, kk);
                if(g_str[cell] == cu_memset) continue;
                for(int_t j = g_str[cell]; j < g_end[cell]; j++){
                    if(P1[j].p_type != 1) continue;
                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx + dy*dy + dz*dz) + 1e-20;
                    if(r < search_range){
                        Real w = calc_kernel_wij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        num += V * P3[j].concn_adv * w;  // reads C_star field
                        den += V * w;
                    }
                }
            }
        }
    }
    return (den > 0.0) ? (num / den) : 0.0;
}

// Helper: 1st order interpolation from P3->concn_adv (the C_star field)
__device__ __forceinline__ Real interp_Cstar_at_1st(Real xq, Real yq, Real zq, int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real hi, int ncell)
{
    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    int_t icell,jcell,kcell;
    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xq-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yq-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zq-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real num = 0.0, den = 0.0;
    Real grad_x = 0.0, grad_y = 0.0, grad_z = 0.0;
    Real xc = 0.0, yc = 0.0, zc = 0.0;

    // Pass 1: 0th order + weighted center
    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;
                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;
                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;
                    if(r < search_range){
                        Real w = calc_kernel_wij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        Real Vw = V * w;
                        num += Vw * P3[j].concn_adv;
                        den += Vw;
                        xc += Vw * P1[j].x;
                        yc += Vw * P1[j].y;
                        zc += Vw * P1[j].z;
                    }
                }
            }
        }
    }

    if(den <= 0.0) return 0.0;

    Real C0 = num / den;
    xc /= den;  yc /= den;  zc /= den;

    // Pass 2: Gradient calculation
    for(int_t z=-ncell; z<=ncell; z++){
        for(int_t y=-ncell; y<=ncell; y++){
            for(int_t x=-ncell; x<=ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t cell = idx_cell(ii,jj,kk);
                if(g_str[cell]==cu_memset) continue;
                for(int_t j=g_str[cell]; j<g_end[cell]; j++){
                    if(P1[j].p_type!=1) continue;
                    Real dx = xq - P1[j].x;
                    Real dy = yq - P1[j].y;
                    Real dz = zq - P1[j].z;
                    Real r  = sqrt(dx*dx+dy*dy+dz*dz) + 1e-20;
                    if(r < search_range){
                        Real dwij = calc_kernel_dwij(tmp_A, hi, r);
                        Real V = P1[j].m / P1[j].rho;
                        Real dC = P3[j].concn_adv - C0;
                        grad_x += V * dC * dwij * dx / r;
                        grad_y += V * dC * dwij * dy / r;
                        grad_z += V * dC * dwij * dz / r;
                    }
                }
            }
        }
    }

    return C0 + grad_x * (xq - xc) + grad_y * (yq - yc) + grad_z * (zq - zc);
}

// -------------------------------------------------------
// Pass 1: Forward SL step  (C^n → C_star stored in P3->concn_adv)
// -------------------------------------------------------
__global__ void KERNEL_BiMC_forward(int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i >= k_num_part2) return;
    if(P1[i].i_type > i_type_crt) return;
    if(P1[i].p_type != 1) return;

    Real xi = P1[i].x, yi = P1[i].y, zi = P1[i].z;
    Real uxi = P1[i].ux, uyi = P1[i].uy, uzi = P1[i].uz;
    Real hi  = P1[i].h;

    Real xb = xi - uxi * tdt;
    Real yb = yi - uyi * tdt;
    Real zb = zi - uzi * tdt;

    P3[i].concn_adv = interp_C_at_1st(xb, yb, zb, g_str, g_end, P1, hi, P1[i].ncell);
}

// -------------------------------------------------------
// Pass 2: Backward advection + BiMacCormack correction
//         + monotone limiter + Cleary-Monaghan diffusion
//         result → P3->dconcn
// -------------------------------------------------------
__global__ void KERNEL_BiMC_correct(int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i >= k_num_part2) return;
    if(P1[i].i_type > i_type_crt) return;
    if(P1[i].p_type != 1) return;

    Real xi = P1[i].x, yi = P1[i].y, zi = P1[i].z;
    Real uxi = P1[i].ux, uyi = P1[i].uy, uzi = P1[i].uz;
    Real hi  = P1[i].h;
    Real concn_n = P1[i].concn;
    Real C_star  = P3[i].concn_adv;
    int_t ncell  = P1[i].ncell;

    // ----------------------------------------------------------
    // STEP 1: BiMacCormack correction
    // Backward step: from C_star field, advect forward to (x + u*dt)
    // ----------------------------------------------------------
    Real xf = xi + uxi * tdt;
    Real yf = yi + uyi * tdt;
    Real zf = zi + uzi * tdt;

    Real C_reverse = interp_Cstar_at_1st(xf, yf, zf, g_str, g_end, P1, P3, hi, ncell);

    Real C_mc = C_star + 0.5 * (concn_n - C_reverse);

    // ----------------------------------------------------------
    // STEP 2: Monotone limiter (Selle et al. 2008)
    // Find min/max of C^n in the neighborhood; if C_mc out of
    // range, fall back to the first-order C_star.
    // ----------------------------------------------------------
    int_t icell,jcell,kcell;
    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xi-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yi-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zi-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    Real C_min = concn_n, C_max = concn_n;
    Real search_range = k_search_kappa * hi;

    for(int_t z = -ncell; z <= ncell; z++){
        for(int_t y = -ncell; y <= ncell; y++){
            for(int_t x = -ncell; x <= ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t cell = idx_cell(ii, jj, kk);
                if(g_str[cell] == cu_memset) continue;
                for(int_t j = g_str[cell]; j < g_end[cell]; j++){
                    if(P1[j].p_type != 1) continue;
                    Real dx = xi - P1[j].x, dy = yi - P1[j].y, dz = zi - P1[j].z;
                    if(sqrt(dx*dx + dy*dy + dz*dz) < search_range){
                        Real cj = P1[j].concn;
                        if(cj < C_min) C_min = cj;
                        if(cj > C_max) C_max = cj;
                    }
                }
            }
        }
    }

    Real C_advected = (C_mc < C_min || C_mc > C_max) ? C_star : C_mc;
    if(C_advected < 0.0) C_advected = 0.0;

    // ----------------------------------------------------------
    // STEP 3: Diffusion (Cleary & Monaghan Laplacian)
    // ----------------------------------------------------------
    Real tmp_A = calc_tmpA(hi);
    Real diffusion_sum = 0.0;

    for(int_t z = -ncell; z <= ncell; z++){
        for(int_t y = -ncell; y <= ncell; y++){
            for(int_t x = -ncell; x <= ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t k = idx_cell(ii, jj, kk);
                if(g_str[k] == cu_memset) continue;
                for(int_t j = g_str[k]; j < g_end[k]; j++){
                    if(P1[j].p_type >= 1000) continue;
                    Real xj = P1[j].x, yj = P1[j].y, zj = P1[j].z;
                    Real mj = P1[j].m, rhoj = P1[j].rho;
                    Real vol_j = mj / rhoj;
                    Real rx = xi-xj, ry = yi-yj, rz = zi-zj;
                    Real tdist = sqrt(rx*rx + ry*ry + rz*rz) + 1e-20;
                    if(tdist < search_range){
                        Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
                        Real tdwx = tdwij * rx / tdist;
                        Real tdwy = tdwij * ry / tdist;
                        Real tdwz = tdwij * rz / tdist;
                        Real concn_j_eff = (P1[j].p_type != 1) ? concn_n : P1[j].concn;
                        Real dot_r_gradW = rx*tdwx + ry*tdwy + rz*tdwz;
                        Real denom = tdist*tdist + 0.01*hi*hi;
                        diffusion_sum += (concn_n - concn_j_eff) * 2.0 * vol_j * dot_r_gradW / denom;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------
    // STEP 4: Combine
    // ----------------------------------------------------------
    Real C_final = C_advected + tdt * Diff_Coeff * diffusion_sum;
    if(C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
}

// -------------------------------------------------------
// Direct SPH advection + diffusion (no interpolation)
// Advection: -u_i · ∇C_i  (anti-symmetric SPH gradient)
// Diffusion: Cleary-Monaghan Laplacian
// Result: P3->dconcn = C^{n+1} (final value, not derivative)
// -------------------------------------------------------
__global__ void KERNEL_Scalar_transport_direct(int_t* g_str, int_t* g_end, part1* P1, part3* P3, Real Diff_Coeff, Real tdt)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i >= k_num_part2) return;
    if(P1[i].i_type > i_type_crt) return;
    if(P1[i].p_type != 1) return;

    Real xi = P1[i].x, yi = P1[i].y, zi = P1[i].z;
    Real uxi = P1[i].ux, uyi = P1[i].uy, uzi = P1[i].uz;
    Real hi  = P1[i].h;
    Real concn_i = P1[i].concn;
    int_t ncell = P1[i].ncell;

    Real tmp_A = calc_tmpA(hi);
    Real search_range = k_search_kappa * hi;

    // Gradient Correction Matrix
    Real gc11=P3[i].inv_cm_xx; Real gc12=P3[i].inv_cm_xy;
    Real gc22=P3[i].inv_cm_yy; Real gc23=P3[i].inv_cm_yz;
    Real gc31=P3[i].inv_cm_zx; Real gc33=P3[i].inv_cm_zz;

    int_t icell,jcell,kcell;
    if(k_x_max==k_x_min) icell=0; else icell=min((int_t)floor((xi-k_x_min)/k_dcell), k_NI-1);
    if(k_y_max==k_y_min) jcell=0; else jcell=min((int_t)floor((yi-k_y_min)/k_dcell), k_NJ-1);
    if(k_z_max==k_z_min) kcell=0; else kcell=min((int_t)floor((zi-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0; if(jcell<0) jcell=0; if(kcell<0) kcell=0;

    // ∇C (SPH gradient with GC) + Laplacian (Cleary-Monaghan) in single neighbor loop
    Real grad_C_x = 0.0, grad_C_y = 0.0, grad_C_z = 0.0;
    Real diffusion_sum = 0.0;

    for(int_t z = -ncell; z <= ncell; z++){
        for(int_t y = -ncell; y <= ncell; y++){
            for(int_t x = -ncell; x <= ncell; x++){
                int_t ii = icell+x, jj = jcell+y, kk = kcell+z;
                if(ii<0||ii>k_NI-1||jj<0||jj>k_NJ-1||kk<0||kk>k_NK-1) continue;
                int_t k = idx_cell(ii, jj, kk);
                if(g_str[k] == cu_memset) continue;
                for(int_t j = g_str[k]; j < g_end[k]; j++){
                    if(P1[j].p_type >= 1000) continue;

                    Real rx = xi - P1[j].x;
                    Real ry = yi - P1[j].y;
                    Real rz = zi - P1[j].z;
                    Real tdist = sqrt(rx*rx + ry*ry + rz*rz) + 1e-20;

                    if(tdist < search_range){
                        Real V_j = P1[j].m / P1[j].rho;
                        Real tdwij = calc_kernel_dwij(tmp_A, hi, tdist);
                        Real dwx = tdwij * rx / tdist;
                        Real dwy = tdwij * ry / tdist;
                        Real dwz = tdwij * rz / tdist;

                        // Gradient Correction 적용
                        Real dwxc = gc11*dwx + gc12*dwy + gc31*dwz;
                        Real dwyc = gc12*dwx + gc22*dwy + gc23*dwz;
                        Real dwzc = gc31*dwx + gc23*dwy + gc33*dwz;

                        Real concn_j = (P1[j].p_type != 1) ? concn_i : P1[j].concn;
                        Real dC = concn_j - concn_i;

                        // SPH gradient with GC: ∇C_i = Σ V_j (C_j - C_i) ∇W_ij^c
                        grad_C_x += V_j * dC * dwxc;
                        grad_C_y += V_j * dC * dwyc;
                        grad_C_z += V_j * dC * dwzc;

                        // Cleary-Monaghan diffusion (GC 적용)
                        Real dot_r_gradW = rx*dwxc + ry*dwyc + rz*dwzc;
                        Real denom = tdist*tdist + 0.01*hi*hi;
                        diffusion_sum += (concn_i - concn_j) * 2.0 * V_j * dot_r_gradW / denom;
                    }
                }
            }
        }
    }

    // Advection: -u · ∇C
    Real advection = -(uxi * grad_C_x + uyi * grad_C_y + uzi * grad_C_z);

    // Explicit Euler: C^{n+1} = C^n + dt * (-u·∇C + D·∇²C)
    Real C_final = concn_i + tdt * (advection + Diff_Coeff * diffusion_sum);
    if(C_final < 0.0) C_final = 0.0;

    P3[i].dconcn = C_final;
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
    int Scale=1;
    Real xi=P1[i].x;
    Real yi=P1[i].y;
    Real zi=P1[i].z;
    Real uxi=P1[i].ux;
    Real uyi=P1[i].uy;
    Real uzi=P1[i].uz;
    Real hi=Scale*P1[i].h; 
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

    int_t ncell=Scale*P1[i].ncell;
    for(int_t z=-ncell;z<=ncell;z++){
        for(int_t y=-ncell;y<=ncell;y++){
            for(int_t x=-ncell;x<=ncell;x++){
                int_t k=idx_cell(icell+x,jcell+y,kcell+z);
                
                if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
                
                if(g_str_LDM[k]!=cu_memset){
                    int_t fend=g_end_LDM[k];
                    for(int_t j=g_str_LDM[k];j<fend;j++){
                        if(LP1[j].p_type<1000){
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
