// =============================================================================
// functions_DD.cuh  —  Dry Deposition & Resuspension for SOPHIA_gpu LDM
// =============================================================================
//
// Deposition velocity model (Wesely 1989, Slinn 1982, Zhang et al. 2001):
//   Vd = 1 / (Ra + Rb) + Vs
//   Ra = [ln(z/z0) - psi_h(z/L) + psi_h(z0/L)] / (kappa * u*)   [aerodynamic]
//   Rb = 1 / (u* * eps0 * (E_B + E_IM))                          [sublayer]
//   Vs = rho_p * dp^2 * g * Cc / (18*mu)                         [settling]
//
// Stochastic deposition (Stohl et al. 2005, FLEXPART):
//   P_dep = 1 - exp(-Vd * dt / dz)
//
// Stochastic resuspension (Biasi et al. 2001, Loosmore 2003):
//   lambda_res = C_res * max(0, u*^2 - u*_t^2)   [s^-1]
//   P_res      = 1 - exp(-lambda_res * dt)
//
// Particle state (L_part1.buffer_type):
//   0  —  active    : subject to Lagrangian advection + deposition check
//   1  —  deposited : immobile, subject to resuspension check only
//
// Requires (included before this file):
//   functions_MOST.cuh  → psi_h(), Phi_h()
//   Parameters.cuh      → us_update, OB_length_global, MOST_z0, L1, k_vonKarman
// =============================================================================

#pragma once

// -----------------------------------------------------------------------------
// Physical constants — Cs-137 aerosol (dp ≈ 1 μm, accumulation mode)
// References: Seinfeld & Pandis (2006), Zhang et al. (2001)
// -----------------------------------------------------------------------------
#define DD_dp           1.0e-6f     // particle diameter [m]
#define DD_rho_p        2500.0f     // particle density  [kg/m³]
#define DD_Cc           1.17f       // Cunningham slip correction (dp = 1 μm)
#define DD_DB           2.8e-11f    // Brownian diffusivity [m²/s]
#define DD_tau_p        9.0e-6f     // particle relaxation time [s]
#define DD_Vs_const     8.9e-5f     // gravitational settling velocity [m/s]
#define DD_nu_air       1.5e-5f     // air kinematic viscosity [m²/s] at ~20°C
#define DD_eps0         3.0f        // empirical sublayer constant (Slinn 1982)

// -----------------------------------------------------------------------------
// Resuspension parameters
// Biasi et al. (2001) "Use of simple models in the first two STORM campaigns"
// Loosmore (2003) "Temporal and spatial trends for resuspension"
//   lambda = C_res * max(0, u*^2 - u*_t^2)   [s^-1]
//   At u* = 0.3 m/s, u*_t = 0.2 m/s : lambda ≈ 5e-5 s^-1  (half-life ~4 h)
//   At u* = 0.5 m/s                  : lambda ≈ 2e-4 s^-1  (half-life ~1 h)
//   At u* = 1.0 m/s                  : lambda ≈ 1e-3 s^-1  (half-life ~12 min)
// -----------------------------------------------------------------------------
#define DD_ustar_thresh  0.2f       // threshold friction velocity [m/s]
#define DD_C_res         1.0e-3f    // resuspension rate coefficient [s^-1 / (m/s)^2]

// -----------------------------------------------------------------------------
// Deposition capture height above ground surface [m]
// Particles within this distance from ground are candidates for deposition.
// Should be O(h_particle) — a few metres for typical urban SPH resolution.
// -----------------------------------------------------------------------------
#define DD_h_dep         2.0f

// 1 = 침적·재부유 끔 (비침적 기체 추적자 검증용, 예: CEDVAL 에탄)
#ifndef DD_DISABLE
#define DD_DISABLE       0
#endif

// --- 벽면·지붕 침적 (건물 IBM 표면) ---
// 0 = 종전대로 지면(z=L1)에만 침적, 1 = IBM 표면 근처 입자도 침적 검사
#ifndef DD_WALL_DEP
#define DD_WALL_DEP      0
#endif
#ifndef DD_wall_dist
#define DD_wall_dist     2.0f     // 표면으로부터 이 거리 안이면 침적 후보 [m]
#endif


// =============================================================================
// DD_calc_Vd: deposition velocity Vd [m/s]
//   z_particle : absolute z of LDM particle [m]
//   z0         : aerodynamic roughness length [m]
//   u_star, OB_len : 이 위치의 국소 값 (DD_local_most 가 준다)
// =============================================================================
__device__ Real DD_calc_Vd(Real z_particle, Real z0, Real u_star, Real OB_len)
{

    // Ra: aerodynamic resistance  Ra = Phi_h(L, z, z0) / (kappa * u*)
    // Phi_h = ln(z/z0) - psi_h(z/L) + psi_h(z0/L)  defined in functions_MOST.cuh
    Real z_eff = fmaxf(z_particle, z0 * 2.0f);
    Real Ra = Phi_h(OB_len, z_eff, z0) / (k_vonKarman * fmaxf(u_star, 0.01f));
    Ra = fmaxf(Ra, 1.0f);  // physical lower bound ~1 s/m

    // Rb: viscous sublayer resistance (Zhang et al. 2001)
    Real Sc   = (Real)DD_nu_air / (Real)DD_DB;                      // Schmidt number ≈ 5.4e5
    Real E_B  = powf(Sc, -2.0f / 3.0f);                             // Brownian efficiency ≈ 1.4e-3
    Real St   = (Real)DD_tau_p * u_star * u_star / (Real)DD_nu_air; // Stokes number
    Real E_IM = (St / (St + 0.5f)) * (St / (St + 0.5f));           // impaction efficiency
    Real Rb   = 1.0f / (fmaxf(u_star, 0.01f) * (Real)DD_eps0 * (E_B + E_IM) + 1.0e-20f);

    return 1.0f / (Ra + Rb) + (Real)DD_Vs_const;
}


// =============================================================================
// DD_calc_Vd_wall: 벽면·지붕용 침적속도 [m/s]
//   dist_to_wall : 가장 가까운 IBM 표면까지의 거리 [m]
//   중력침강 항 없음(수직면), 저항은 표면까지 거리로 평가
// =============================================================================
__device__ Real DD_calc_Vd_wall(Real dist_to_wall, Real z0, Real u_star, Real OB_len)
{
    Real d_eff  = fmaxf(dist_to_wall, z0 * 2.0f);
    Real Ra = Phi_h(OB_len, d_eff, z0) / (k_vonKarman * fmaxf(u_star, 0.01f));
    Ra = fmaxf(Ra, 1.0f);
    Real Sc   = (Real)DD_nu_air / (Real)DD_DB;
    Real E_B  = powf(Sc, -2.0f / 3.0f);
    Real St   = (Real)DD_tau_p * u_star * u_star / (Real)DD_nu_air;
    Real E_IM = (St / (St + 0.5f)) * (St / (St + 0.5f));
    Real Rb   = 1.0f / (fmaxf(u_star, 0.01f) * (Real)DD_eps0 * (E_B + E_IM) + 1.0e-20f);
    return 1.0f / (Ra + Rb);
}


// =============================================================================
// DD_wall_distance: 가장 가까운 IBM 표면 입자까지의 거리 [m]
//   KERNEL_time_update_LDM 의 셀 검색과 같은 방식 (27 셀 스캔)
// =============================================================================
__device__ Real DD_wall_distance(Real xi, Real yi, Real zi, int_t ncell,
                                 int_t* g_str, int_t* g_end, part1* P1)
{
    int_t icell = (k_x_max==k_x_min) ? 0 : min((int_t)floor((xi-k_x_min)/k_dcell), k_NI-1);
    int_t jcell = (k_y_max==k_y_min) ? 0 : min((int_t)floor((yi-k_y_min)/k_dcell), k_NJ-1);
    int_t kcell = (k_z_max==k_z_min) ? 0 : min((int_t)floor((zi-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0;  if(jcell<0) jcell=0;  if(kcell<0) kcell=0;

    Real near_d = (Real)1e20;
    for(int_t z=-ncell; z<=ncell; z++)
    for(int_t y=-ncell; y<=ncell; y++)
    for(int_t x=-ncell; x<=ncell; x++){
        if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||
           ((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
        int_t k = idx_cell(icell+x, jcell+y, kcell+z);
        if(g_str[k]==cu_memset) continue;
        int_t fend = g_end[k];
        for(int_t j=g_str[k]; j<fend; j++){
            if(P1[j].p_type < 1000) continue;                 // IBM 표면만
            Real dx=xi-P1[j].x, dy=yi-P1[j].y, dz=zi-P1[j].z;
            Real d = sqrtf(dx*dx + dy*dy + dz*dz);
            if(d < near_d) near_d = d;
        }
    }
    return near_d;
}


// =============================================================================
// DD_calc_P_dep: stochastic deposition probability (Stohl et al. 2005)
//   P = 1 - exp(-Vd * tstep / dz)
// =============================================================================
__device__ Real DD_calc_P_dep(Real Vd, Real tstep, Real dz)
{
    return 1.0f - expf(-Vd * tstep / fmaxf(dz, 1.0e-4f));
}


// =============================================================================
// DD_calc_lambda_res: resuspension rate [s^-1]  (Biasi 2001 / Loosmore 2003)
//   lambda = C_res * max(0, u*^2 - u*_t^2)
// =============================================================================
__device__ Real DD_calc_lambda_res(Real u_star)
{
    Real excess = u_star * u_star - (Real)(DD_ustar_thresh * DD_ustar_thresh);
    return (excess > 0.0f) ? (Real)DD_C_res * excess : 0.0f;
}


// =============================================================================
// DD_calc_P_res: stochastic resuspension probability
//   P = 1 - exp(-lambda_res * tstep)
//   u_star : 이 위치의 국소 마찰속도
// =============================================================================
__device__ Real DD_calc_P_res(Real tstep, Real u_star)
{
    return 1.0f - expf(-DD_calc_lambda_res(u_star) * tstep);
}


// =============================================================================
// DD_local_most: 이 위치의 표면층 u* 와 Obukhov 길이
//   27 셀(ncell) 안의 MOST 밴드 유체 입자에서 수평풍속·온도·높이를 먼저 평균내고,
//   그 평균 상태로 Newton 을 한 번만 돌린다. KERNEL_MOST_boundary3D 와 같은 식이다.
//   이웃이 없으면 전역 밴드 평균(us_update / OB_length_global)으로 폴백한다.
//   지면 침적·재부유는 z<DD_h_dep 에서만 부르므로, 셀 크기(dcell)가 밴드 두께보다
//   커서 입자 자기 셀이 이미 밴드를 포함한다.
//   theta_surface 는 T_con_t=0 (즉시 적용) 기준이라 MOST 커널의 램프와 동일하다.
// =============================================================================
__device__ void DD_local_most(Real xi, Real yi, Real zi, int_t ncell,
                              int_t* g_str, int_t* g_end, part1* P1,
                              Real* u_star_out, Real* OB_out)
{
    int_t icell = (k_x_max==k_x_min) ? 0 : min((int_t)floor((xi-k_x_min)/k_dcell), k_NI-1);
    int_t jcell = (k_y_max==k_y_min) ? 0 : min((int_t)floor((yi-k_y_min)/k_dcell), k_NJ-1);
    int_t kcell = (k_z_max==k_z_min) ? 0 : min((int_t)floor((zi-k_z_min)/k_dcell), k_NK-1);
    if(icell<0) icell=0;  if(jcell<0) jcell=0;  if(kcell<0) kcell=0;

    Real sU=0.0f, sT=0.0f, sZ=0.0f;
    int  n=0;
    for(int_t z=-ncell; z<=ncell; z++)
    for(int_t y=-ncell; y<=ncell; y++)
    for(int_t x=-ncell; x<=ncell; x++){
        if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||
           ((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
        int_t k = idx_cell(icell+x, jcell+y, kcell+z);
        if(g_str[k]==cu_memset) continue;
        int_t fend = g_end[k];
        for(int_t j=g_str[k]; j<fend; j++){
            if(P1[j].MOST_buffer != 1) continue;
            sU += sqrtf(P1[j].ux*P1[j].ux + P1[j].uy*P1[j].uy);
            sT += P1[j].temp;
            sZ += P1[j].z;
            n++;
        }
    }

    if(n == 0){                       // 이웃 없음 → 전역 밴드 평균으로 폴백
        *u_star_out = (Real)us_update;
        *OB_out     = (Real)OB_length_global;
        return;
    }

    Real U  = fmaxf(sU/(Real)n, 0.001f);
    Real T  = sT/(Real)n;
    Real zm = fmaxf(sZ/(Real)n, (Real)MOST_z0*2.0f);
    Real z0 = (Real)MOST_z0, z0h = (Real)MOST_z0h;

    Real Ri_b = (Gravitational_CONST/T)*(T-(Real)T_con)*zm/(U*U+1e-10f);
    Real L, ustar;
    if(fabsf(Ri_b) > 1e-3f){
        L = clamp(zm/Ri_b, -500.0f, 500.0f);
        for(int it=0; it<MAX_ITER; it++){
            Real lm=L, ll=0.999f*lm, lu=1.001f*lm;
            Real fL  = Ri_b-(zm/lm)*Phi_h(lm,zm,z0h)/powf(Phi_m(lm,zm,z0),2);
            Real dfL = ((-(zm/lu)*Phi_h(lu,zm,z0h)/(Phi_m(lu,zm,z0)*Phi_m(lu,zm,z0)))
                       -(-(zm/ll)*Phi_h(ll,zm,z0h)/(Phi_m(ll,zm,z0)*Phi_m(ll,zm,z0))))/(lu-ll);
            if(fabsf(dfL) < 1e-10f) break;
            L = clamp(L - fL/dfL, -500.0f, 500.0f);
            if(fabsf(L) > 1e-6f && fabsf((L-lm)/L) < TOLERANCE) break;
        }
        ustar = k_vonKarman*U/Phi_m(L,zm,z0);
    }else{                            // 중립
        L = 1e6f;
        ustar = k_vonKarman*U/logf(zm/z0);
    }
    *u_star_out = fminf(ustar, 2.0f);
    *OB_out     = L;
}


// =============================================================================
// KERNEL_dry_deposition
//   — Active particles (buffer_type == 0) within DD_h_dep of ground:
//       stochastic deposition check → buffer_type = 1 on hit.
//   — Deposited particles (buffer_type == 1):
//       stochastic resuspension check → buffer_type = 0 on lift.
//
//   Called once per freq_LDM steps after KERNEL_time_update_LDM.
//   Uses LP1[i].seed with curand offset +17 to stay independent
//   from the Langevin step (offset 0) within the same timestep.
//
//   Parameters:
//     tdt_in : effective LDM timestep (= freq_LDM * dt)  [s]
//     tt     : current simulation time                    [s]
//     LP1    : LDM particle array
//     z0_in  : aerodynamic roughness length               [m]
// =============================================================================
__global__ void KERNEL_dry_deposition(Real tdt_in, Real tt, L_part1* LP1, Real z0_in,
                                      int_t* g_str, int_t* g_end, part1* P1)
{
    uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= k_num_part_LDM) return;
    if (LP1[i].i_type == 3) return;   // slot never used
    if (LP1[i].i_type == 0) return;   // not yet released

    // -------------------------------------------------------------------------
    // RNG: advance seed from value left by KERNEL_time_update_LDM,
    // then curand_init with offset +17 for independence from Langevin step.
    // -------------------------------------------------------------------------
    curandState ss;
    unsigned long long time_hash = (unsigned long long)(tt * 1.0e6);
    unsigned long long seed_val = LP1[i].seed;
    seed_val = (seed_val * 1103515245ULL + time_hash + 54321ULL) % 2147483647ULL;
    LP1[i].seed = seed_val;
    curand_init(seed_val, i, (unsigned long long)(tt * 100.0) + 17ULL, &ss);

    Real r = curand_uniform(&ss);
    LP1[i].seed = (unsigned long long)curand(&ss);  // advance for next call

    Real z_ground = (Real)L1;
    Real zi       = LP1[i].z;

    // =========================================================================
    // Case 1: active particle — ground deposition check
    // =========================================================================
    if (LP1[i].buffer_type == 0) {
        Real dz = zi - z_ground;
        if (dz < (Real)DD_h_dep) {
            Real us_l, OB_l;
            DD_local_most(LP1[i].x, LP1[i].y, zi, LP1[i].ncell, g_str, g_end, P1, &us_l, &OB_l);
            Real Vd    = DD_calc_Vd(zi, z0_in, us_l, OB_l);
            Real P_dep = DD_calc_P_dep(Vd, tdt_in, fmaxf(dz, 1.0e-4f));
            if (r < P_dep) {
                LP1[i].buffer_type = 1;
                LP1[i].z  = z_ground + 1.0e-4f;
                LP1[i].ux = 0.0f;
                LP1[i].uy = 0.0f;
                LP1[i].uz = 0.0f;
            }
        }
#if DD_WALL_DEP
        // 지면에 붙지 않은 입자만 벽면 검사. 난수 r 을 그대로 재사용하므로
        // 침적하지 않는 입자의 궤적은 DD_WALL_DEP=0 런과 비트 단위로 같다.
        if (LP1[i].buffer_type == 0) {
            Real near_d = DD_wall_distance(LP1[i].x, LP1[i].y, zi, LP1[i].ncell, g_str, g_end, P1);
            if (near_d < (Real)DD_wall_dist) {
                Real us_w, OB_w;
                DD_local_most(LP1[i].x, LP1[i].y, zi, LP1[i].ncell, g_str, g_end, P1, &us_w, &OB_w);
                Real Vdw = DD_calc_Vd_wall(near_d, z0_in, us_w, OB_w);
                if (r < DD_calc_P_dep(Vdw, tdt_in, near_d)) {
                    LP1[i].buffer_type = 2;   // 2 = 벽면·지붕 침적 (위치 보존)
                    LP1[i].ux = 0.0f;
                    LP1[i].uy = 0.0f;
                    LP1[i].uz = 0.0f;
                }
            }
        }
#endif
    }
    // =========================================================================
    // Case 2: deposited particle — resuspension check
    // =========================================================================
    else if (LP1[i].buffer_type == 1) {
        Real us_l, OB_l;
        DD_local_most(LP1[i].x, LP1[i].y, zi, LP1[i].ncell, g_str, g_end, P1, &us_l, &OB_l);
        (void)OB_l;
        Real P_res = DD_calc_P_res(tdt_in, us_l);
        if (r < P_res) {
            LP1[i].buffer_type = 0;
            LP1[i].uz = us_l * 0.3f;   // upward kick ~ 0.3 u* (국소값)
        }
    }
}
