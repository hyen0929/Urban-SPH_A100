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


// =============================================================================
// DD_calc_Vd: deposition velocity Vd [m/s]
//   z_particle : absolute z of LDM particle [m]
//   z0         : aerodynamic roughness length [m]
//   reads us_update and OB_length_global from global device vars (Parameters.cuh)
// =============================================================================
__device__ Real DD_calc_Vd(Real z_particle, Real z0)
{
    Real u_star = (Real)us_update;
    Real OB_len  = (Real)OB_length_global;

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
//   reads us_update from global device var (Parameters.cuh)
// =============================================================================
__device__ Real DD_calc_P_res(Real tstep)
{
    Real u_star = (Real)us_update;
    return 1.0f - expf(-DD_calc_lambda_res(u_star) * tstep);
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
__global__ void KERNEL_dry_deposition(Real tdt_in, Real tt, L_part1* LP1, Real z0_in)
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
            Real Vd    = DD_calc_Vd(zi, z0_in);
            Real P_dep = DD_calc_P_dep(Vd, tdt_in, fmaxf(dz, 1.0e-4f));
            if (r < P_dep) {
                LP1[i].buffer_type = 1;
                LP1[i].z  = z_ground + 1.0e-4f;
                LP1[i].ux = 0.0f;
                LP1[i].uy = 0.0f;
                LP1[i].uz = 0.0f;
            }
        }
    }
    // =========================================================================
    // Case 2: deposited particle — resuspension check
    // =========================================================================
    else if (LP1[i].buffer_type == 1) {
        Real P_res = DD_calc_P_res(tdt_in);
        if (r < P_res) {
            LP1[i].buffer_type = 0;
            LP1[i].uz = (Real)us_update * 0.3f;  // upward kick ~ 0.3 u*
        }
    }
}
