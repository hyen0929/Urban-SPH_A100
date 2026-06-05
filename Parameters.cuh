#define		cu_memset				0xffffffff
#define 	vii_size				64
#define 	vif_size				32
#define		Max_GPU					10

#define 	Wcsph 					0
#define 	Isph 					1

#define 	FLUID				    1
#define     AIR                     2
#define 	BOUNDARY				0
#define 	MOVING					9

#define 	Liquid 					0
#define 	Gas 					1
#define 	Solid 					2

#define 	WRONG_INDEX		    	1e8			// Limitation of WRONG_INDEX: 2.2e9 (signed integer)

#define 	Gaussian				0
#define 	Quintic					1
#define 	Quartic					2
#define 	Wendland2				3
#define 	Wendland4				4
#define 	Wendland6				5

#define 	Euler					0				// Euler Explicit Time Stepping
#define 	Pre_Cor					1				// Predictor-Corrector Time Stepping

#define 	PI						3.141592653

#define 	Gravitational_CONST     9.810	// gravitational constant [m/s2]

#define 	Cs_SPS					0.12			// check please (by esk)
#define 	CI_SPS					0.00066
#define 	L_SPS					0.01			// scale of length scale (for test)
#define     C_pres                  0.2       // Pressure_smoothing_coefficient

// viscous model
#define     Cleary                  1
#define     Monaghan                2
#define     Violeau                  3

#define 	Lm					0.01
#define 	Laminar				0
#define 	K_LM				1
#define 	K_E					2
#define 	HB      			3
#define     SSM                 4
#define     DDF                 5

#define 	DIFF_DENSITY		    0
#define     Correction_Matrix_Size  3

#define 	KGC		    1
#define 	FPM		    2
#define 	DFPM		  3
#define 	KGF		    4

#define   Min_det   1e-7      // minimum of determinent

#define   Pb    0
#define   i_type_crt  2
#define   C_p2p   1.5

#define table_size 50

#define h_coeff 1.5

#define Eulerian  0
#define Lagrangian  1
#define ALE 2

#define Boundary_Velocity 0.0 // No-slip condition for moving boundary

/// APR ///
#define APR_Direct 1
#define APR_Activate 2
#define APR_Block 3

//// open boudary /////
#define Inlet 1
#define Outlet 2
#define Left 3
#define Right 4
#define Dummy 5
#define Ceiling 6
#define Inlet_Density 1.0
#define Inlet_Velocity 1.0
#define Wall_Velocity 0.5
#define MIRROR_LENGTH 0.5
#define Extrapolation_Length 0.4
////////////////////////

//// energy ////////////
#define enthalpy_eqn 0          // energy equation described as enthalpy
#define h_CONV 15.0
#define T_SUR 300.0
#define EMISSIVITY 0.65
#define sigma_SB 5.670374e-8
////////////////////////

//// temperature advection scheme ////
// eulert (SPH temperature advection) discretization, switched at compile time.
//   0 = central-diff (current production: 2mj/(rhoi+rhoj)(Tj-Ti)((ui-uj).dW)).
//       Non-monotone; drives the windward roof-corner temperature runaway.
//   1 = Vila 1999 conservative Eulerian advection, 1st-order upwind (Godunov).
//       eulert = sum 2mj/(rhoi+rhoj)(T*-Ti)(u*.dW), u*=0.5(ui+uj),
//       T* = upwind state (Tj if neighbour upwind, Ti otherwise). Monotone by construction.
// See Temperature_advection.md and functions_PPE.cuh:KERNEL_advection_force3D.
#define RIEMANN_EULERT 1
//////////////////////////////////////

//// PPE solver mode ////
// 0 = EISPH (Explicit ISPH, 1-sweep Jacobi: current production).
//     The PPE loop in ISPH_Calc.cuh runs exactly once per step.
// 1 = Implicit PPE (Jacobi N-sweep with early-exit on residual).
//     Iterates PPE_prep + PPE_calc + Neumann up to maxIteration (solv vii[21]),
//     stopping when max |p_new - p_old| < PPE_TOLERANCE. Iter count and final
//     residual are logged per step to ./<CASE_NAME>/PPE_iter_log.txt and
//     printed to stdout every 50 steps. minIteration (vii[20]) is honoured
//     as a lower bound (no early exit before that). Costs N* the per-step
//     PPE time (typically 100-1000* slower per step on large grids).
#define IMPLICIT_PPE     0
#define PPE_TOLERANCE    1e-5f
/////////////////////////

//// building (IBM) deactivation ////
// 0 = IBM building active (normal run).
// 1 = building-free run: KERNEL_deactivate_IBM sets i_type=3 on all
//     p_type==1000 particles at init, so the building is fully removed
//     (skipped as target and as neighbour). Interior fluid stays active
//     -> the former building volume becomes ordinary fluid.
//     Used for clean SGS-model (DDF vs SSM) buoyancy comparisons.
#define DEACTIVATE_IBM 1
/////////////////////////////////////

#define ncell_init 2
#define Cell_division_factor 4
#define IBM_length 1.0
#define Pressure_length 1.0
#define XSPH_length 0.66667

/// MOST ///
#define MAX_ITER 200
#define TOLERANCE 1e-5
#define k_vonKarman 0.4
#define Businger_const 16.0   // For MOST psi function
#define MOST_z0 0.1            // Roughness length [m] (0.1=grassland, 1.0=urban)
#define MOST_z0h 0.01          // Roughness length for heat [m]
#define us_ini 0.5             // Initial friction velocity [m/s]
#define ABL_h 125.0            // ABL height [m]
#define Pr_t 0.71             // Turbulent pr number

#define Solid_cp 900.0
#define Air_cp 1005         // J/kgK
#define Solid_k 1.0
#define Air_k 0.025
#define Solid_rho 2300.0

// ---- Temperature settings (change T_con per case) ----
// 211 Unstable: T_con=315.0  (ΔT=+15K, L<0)
// 212 Neutral:  T_con=300.0  (ΔT=  0K, L→∞)
// 213 Stable:   T_con=285.0  (ΔT=-15K, L>0)
#define T_ini 300.0       // Initial fluid temperature [K]
#define T_con 305       // Surface temperature [K]  ← DDF/SSM 부력 비교: 지표 330K, T_ini=300K 성층
#define T_start_t 0       // Start time [s]
#define T_con_t 0.0       // 0=즉시 T_con 적용
#define Pb 0.0
#define tempb 300.0       // Reference potential temperature [K]

#define Scalar_Diffusion_coeff 0.01  // Scalar diffusion coefficient (m^2/s) — matched to LDM D=0.01

// Deardorff SGS model constants
#define Ck_sgs 0.1
#define Ce_sgs 0.93
#define Cs_sgs 0.15
#define Prt_sgs 0.7
#define Sigk_sgs 1.0

// Standard k-epsilon (Launder-Spalding) constants, used by functions_TURB.cuh
#define C_mu     0.09
#define C_e1     1.44
#define C_e2     1.92
#define sigma_k  1.0
#define sigma_e  1.3
#define kappa_t  0.41

#define k_ib_cx 0.0
#define k_ib_cy 0.0
#define k_ib_cz 0.0

#define MOST_boundary 0   // 1: MOST boundary condition 적용, 0: 적용 안함 (DDF/SSM 부력 비교: ON)
#define Rho_air 1.225

// Channel-flow validation (CASE_ID 1500): streamwise body force [m/s2]
// f_x = u_tau^2 / delta = 0.5^2 / 80 = 0.003125  (drives periodic channel)
#define CHANNEL_FX 0.003125
#define CHANNEL_U0 10.0    // initial uniform streamwise velocity [m/s]

#define SIZ 0.01

#define INPUT_CASE_NAME "CaseCS1"
#define OUTPUT_CASE_NAME "CaseCS1"

#define INPUT_SPH_PATH  "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME ".txt"
#define INPUT_LDM_PATH  "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME "_LDM.txt"
#define SOLV_SPH_PATH   "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME "_solv.txt"
#define OUTPUT_DIR_NAME OUTPUT_CASE_NAME

__device__ float us_update = us_ini;
__device__ float OB_length_global = 1e6f;  // Obukhov length for LDM dry deposition [m], updated by MOST kernel

// Time-averaging (Eulerian SPH)
// Per-case override in case_config.cuh. 0.0 = average from start of run.
#ifndef AVG_START_TIME
#define AVG_START_TIME 30.0
#endif
__device__ unsigned int avg_count = 0;  // # of accumulated samples (incremented by KERNEL_accumulate_avg)