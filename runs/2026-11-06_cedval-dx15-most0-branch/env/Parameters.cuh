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

// LDM 입자가 IBM 표면에서 거울반사되는 임계거리 [m] (functions_LDM.cuh no-penetration).
// 기본 <=0 이면 종전대로 0.5*dx_sp (= 0.5*h/h_coeff, Leipzig 에서 1.25 m).
// 하한은 IBM 표면 입자 간격의 절반(Leipzig 측정 중앙 1.06 m → 0.53 m)과 LDM 스텝 이동거리
// (0.1 s × 3~5 m/s = 0.3~0.5 m) — 그보다 작으면 점군 틈으로 새거나 한 스텝에 벽을 건너뛴다.
#ifndef LDM_WALL_REFLECT_D
#define LDM_WALL_REFLECT_D 0.0
#endif

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
// See docs/Temperature_advection.md and functions_PPE.cuh:KERNEL_advection_force3D.
#ifndef RIEMANN_EULERT
#define RIEMANN_EULERT 1
#endif
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
#ifndef IMPLICIT_PPE
#define IMPLICIT_PPE     0
#endif
// Early-exit gate: stop sweeping once max|p^{k+1}-p^k| < PPE_TOLERANCE [Pa].
// Capped at maxIteration (solv-file maximum-iteration = 200). Early in the run
// it will hit the 200 cap; once the flow develops the warm-started per-step
// pressure increment shrinks and this gate fires earlier. Retune from
// ./<CASE>/PPE_iter_log.txt (iter_used, delta_p columns).
#ifndef PPE_TOLERANCE
#define PPE_TOLERANCE    1e-3f
#endif
/////////////////////////

//// PPE residual diagnostic (explicit & implicit) ////
// 1 = every PPE_DIAG_FREQ steps, recompute the PPE source bi and the TRUE
//     residual |A*p - b| using the post-solve pressure, max-reduce both, and
//     append a CSV row to ./<OUTPUT_CASE_NAME>/PPE_diag_{explicit,implicit}.csv
//     (file name picked by IMPLICIT_PPE so the two runs don't clobber).
//     resid_max ~ |div(u^{n+1})| is the real incompressibility metric: lets the
//     explicit(1-sweep) and implicit(N-sweep) runs be compared on equal footing.
//     Cost: one extra PPE_prep + 2 CUB reduces per logged step. tid==0 only
//     (single-GPU; multi-GPU logs the GPU-0 subdomain max).
#ifndef PPE_DIAG
#define PPE_DIAG       1
#endif
#ifndef PPE_DIAG_FREQ
#define PPE_DIAG_FREQ  10
#endif
///////////////////////////////////////////////////////

//// Consistent-operator PPE (IISPH) ////
// 0 = production Brookshaw PPE operator (A != div of applied force; inconsistent).
// 1 = IISPH consistent operator: the matrix A is built as the divergence of the
//     ACTUALLY-applied pressure force, so converging it is a true projection.
//     Replaces PPE_prep/calc with KERNEL_iisph_prep + (pass1,pass2) relaxed-Jacobi
//     sweeps, source = density error (rho_0 - rho*), and switches pressureforce3D
//     to the matching symmetric form a_p = -sum_j m_j(p_i/rho_i^2 + p_j/rho_j^2)gradW
//     (plain gradW, no KGC, so operator==force and A is symmetric). Iterated via the
//     IMPLICIT_PPE loop; validate with PPE_DIAG resid_max (should drop & stay low).
#ifndef CONSISTENT_PPE
#define CONSISTENT_PPE   0
#endif
#ifndef IISPH_OMEGA
#define IISPH_OMEGA      0.5f   // relaxed-Jacobi factor (Ihmsen 2014)
#endif
// Early-exit gate in DIVERGENCE units: stop once mean |rho*div(u^{n+1})/dt|
// (= |residual|/dt^2, same quantity as div_un1) drops below IISPH_ETA. This is the
// actual incompressibility, NOT density error — important because this solver holds
// rho=rho0 (density never evolves), so a density-error gate is misleadingly loose.
// Capped at maxIteration (solv). Logged to ./<CASE>/IISPH_log.txt.
// 0.09 sits just above the relaxed-Jacobi floor (~0.08) so it early-exits in ~tens
// of sweeps (vs 200 cap at 0.05), while still well below explicit (~0.16).
#ifndef IISPH_ETA
#define IISPH_ETA        0.09f
#endif
// Stagnation gate: also stop when the per-sweep RELATIVE residual drop falls below this,
// i.e. the solver has hit THIS step's floor and further sweeps won't help. Without it, when
// the developed-flow floor rises above IISPH_ETA the loop spins to maxIteration (200-cap)
// even though it converged in ~30-50 sweeps. Convergence = (avg_err<IISPH_ETA) OR (rel_drop<this).
#ifndef IISPH_STAG_TOL
#define IISPH_STAG_TOL   1e-3f
#endif
// Every IISPH_CONV_FREQ steps, dump that step's per-iteration div curve (no early-exit
// on those steps, run full maxIteration) to ./<CASE>/IISPH_conv.txt  (step, iter, avg_div).
// 0 = disabled: every step uses the 0.09 early-exit (no full-200 conv steps -> faster).
#ifndef IISPH_CONV_FREQ
#define IISPH_CONV_FREQ  0
#endif
// KGC-consistent IISPH (smoke test): 0 = plain gradW (symmetric, stable, but solves the
//   LESS-accurate plain divergence -> KGC-div stays ~0.27 > explicit). 1 = apply the SAME
//   KGC correction (Cm_i) to gradW in EVERY IISPH term (prep G/H/D, pass1 d_ij, pass2 Ap)
//   AND the applied pressure force, so the solver actually targets the accurate KGC
//   divergence. RISK (see investigation summary section 2): KGC breaks gradW_ij=-gradW_ji
//   -> operator loses symmetry/diagonal-dominance -> relaxed-Jacobi may blow up exactly
//   like the Brookshaw attempt. This flag exists to TEST that: run a few hundred steps,
//   watch bi_max. Requires CONSISTENT_PPE=1 and solv kgc=KGC.
#ifndef IISPH_KGC
#define IISPH_KGC        1
#endif
/////////////////////////////////////////

//// building (IBM) deactivation ////
// 0 = IBM building active (normal run).
// 1 = building-free run: KERNEL_deactivate_IBM sets i_type=3 on all
//     p_type==1000 particles at init, so the building is fully removed
//     (skipped as target and as neighbour). Interior fluid stays active
//     -> the former building volume becomes ordinary fluid.
//     Used for clean SGS-model (DDF vs SSM) buoyancy comparisons.
#ifndef DEACTIVATE_IBM
#define DEACTIVATE_IBM 1
#endif
/////////////////////////////////////

#define ncell_init 2
#define Cell_division_factor 4
#ifndef IBM_length
#define IBM_length 1.0
#endif
#define Pressure_length 1.0
#define XSPH_length 0.66667

// IBM_PPE_WALL: treat building markers (p_type==1000) as Neumann pressure
//     neighbours in the PPE, exactly like ground particles (p_type==0).
//     0 = production (markers excluded from every SPH operator)
//     1 = markers enter KERNEL_PPE_prep3D and KERNEL_pressureforce3D
//   Why: the IBM velocity forcing is applied to u** BEFORE the projection
//   (ISPH_Calc.cuh IBM_corrector_fluid -> PPE -> pressureforce -> corrector),
//   but the PPE excludes p_type>=1000, so the pressure solve has no boundary
//   condition at a building surface. The projection then undoes the forcing
//   and fluid is pushed through the wall. Markers already carry an
//   extrapolated pressure from KERNEL_Neumann_boundary3D (called every PPE
//   sweep) and have u=0; this flag just lets the PPE read them.
//   Caveat: a single marker layer only fills a ~2.5 m shell of a 15 m kernel
//   support (Kori: dx_IBM=2.5 m, kappa*h=15 m), so this is a partial fix.
//   Full closure needs the interior particles (buffer_type==99) as wall
//   particles too. EISPH path only; CONSISTENT_PPE/IISPH kernels untouched.
#ifndef IBM_PPE_WALL
#define IBM_PPE_WALL 0
#endif

/// MOST ///
#define MAX_ITER 200
#define TOLERANCE 1e-5
#define k_vonKarman 0.4
#define Businger_const 16.0   // For MOST psi function
#ifndef MOST_z0
#define MOST_z0 0.1            // Roughness length [m] (0.1=grassland, 1.0=urban)
#endif
#ifndef MOST_z0h
#define MOST_z0h 0.01          // Roughness length for heat [m]
#endif
#ifndef us_ini
#define us_ini 0.5             // Initial friction velocity [m/s]
#endif
#ifndef ABL_h
#define ABL_h 125.0            // ABL height [m]
#endif
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
#ifndef T_ini
#define T_ini 300.0       // Initial fluid temperature [K]
#endif
#ifndef T_con
#define T_con 300.0       // Surface temperature [K]  ← DDF/SSM 부력 비교: 지표 330K, T_ini=300K 성층
#endif
#define T_start_t 0       // Start time [s]
#define T_con_t 0.0       // 0=즉시 T_con 적용
#define Pb 0.0
#define tempb 300.0       // Reference potential temperature [K]

#ifndef Scalar_Diffusion_coeff
#define Scalar_Diffusion_coeff 0.01  // Scalar diffusion coefficient (m^2/s) — matched to LDM D=0.01
#endif

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

#ifndef MOST_boundary
#define MOST_boundary 0   // 1: MOST boundary condition 적용, 0: 적용 안함 (Kori: OFF)
#endif
#define Rho_air 1.225

// Convective outflow (Orlanski) phase speed [m/s] for the outlet buffer.
// <=0 : use the local kernel-interpolated streamwise velocity (default, z-dependent)
// >0  : pin to this constant (e.g. -D OUTFLOW_UC=10.0 reproduces the old behaviour)
#ifndef OUTFLOW_UC
#define OUTFLOW_UC 0.0
#endif

// Channel-flow validation (CASE_ID 1500): streamwise body force [m/s2]
// f_x = u_tau^2 / delta = 0.5^2 / 80 = 0.003125  (drives periodic channel)
#ifndef CHANNEL_FX
#define CHANNEL_FX 0.003125
#endif
#ifndef CHANNEL_U0
#define CHANNEL_U0 10.0    // initial uniform streamwise velocity [m/s]
#endif

// Rayleigh damping (sponge) layer near the ceiling.
//   γ(z) = SPONGE_GAMMA_MAX × ramp(z)²,  ramp = (z - Z_BOTTOM)/(z_top - Z_BOTTOM)
//   uz   *= (1 - γ·dt)   [uz always; uy too when SPONGE_DAMP_UY=1; ux untouched to protect mean flow]
//   z_top = k_z_max (domain ceiling)
// Set SPONGE_GAMMA_MAX <= 0 to disable.
#ifndef SPONGE_Z_BOTTOM
#define SPONGE_Z_BOTTOM 90.0
#endif
#ifndef SPONGE_GAMMA_MAX
#define SPONGE_GAMMA_MAX 0.2f    // [1/s], e-folding time = 10 s at the ceiling
#endif
// 1 → sponge also damps uy with the same γ(z). Only for cases with zero cross-stream mean wind.
// Default 0 = previous behaviour (uz only).
#ifndef SPONGE_DAMP_UY
#define SPONGE_DAMP_UY 0
#endif

#define SIZ 0.01

#ifndef INPUT_CASE_NAME
#define INPUT_CASE_NAME "Kori"
#endif
#ifndef OUTPUT_CASE_NAME
#define OUTPUT_CASE_NAME "Kori_vy1"
#endif

#define INPUT_SPH_PATH  "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME ".txt"
#define INPUT_LDM_PATH  "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME "_LDM.txt"
#define SOLV_SPH_PATH   "./input/" INPUT_CASE_NAME "/" INPUT_CASE_NAME "_solv.txt"
#define OUTPUT_DIR_NAME OUTPUT_CASE_NAME

__device__ float us_update = us_ini;
__device__ float OB_length_global = 1e6f;  // Obukhov length for LDM dry deposition [m], updated by MOST kernel

// MOST 밴드 누적기 — 매 스텝 KERNEL_MOST_boundary3D 가 더하고
// KERNEL_MOST_band_reduce 가 평균내어 us_update/OB_length_global 에 넣은 뒤 0 으로 되돌린다.
// 이 값은 이웃을 못 찾은 LDM 입자의 폴백으로만 쓰인다 (DD_local_most 참고).
// 1/L 을 누적하는 이유: L 은 중립에서 발산하지만 1/L 은 0 을 지나며 연속이다.
__device__ float us_band_sum  = 0.0f;
__device__ float invL_band_sum= 0.0f;
__device__ int   us_band_cnt  = 0;

// Time-averaging (Eulerian SPH)
//
// Two paths to set the averaging-window start time:
//   AVG_START_TIME == 0 (default) → quasi-steady auto-detect from Σ½|u|²
//                                    time series (see QS_DETECT_* below)
//   AVG_START_TIME >  0           → fixed time, bypass detection
#ifndef AVG_START_TIME
#define AVG_START_TIME 0.0
#endif
// End of the averaging window.
//   AVG_END_TIME == 0 → auto = (detected/fixed start) + AVG_WINDOW
//   AVG_END_TIME >  0 → fixed absolute time
//   AVG_WINDOW   <=0 → no end (accumulate to time_end, no freeze/dump)
#ifndef AVG_END_TIME
#define AVG_END_TIME 0.0
#endif
// Averaging window length added to the (detected/fixed) start to derive
// the end time when AVG_END_TIME == 0.
// ------------------------------------------------------------
// 초기 횡류(y) 속도 [m/s] — Leipzig 계열(CASE_ID 1000~1100) 전역 초기조건.
// 운동량 방정식에 y 방향 강제항이 없으므로(f_geo 는 x 전용) 솔버가 이 값을
// 유지해주지 않는다. 순수한 초기조건일 뿐이며 기본값 0 = 횡류 없음.
// 2026-09-10 이전 기본값은 -1.0 이었다.
// ------------------------------------------------------------
#ifndef INIT_CROSSFLOW_UY
#define INIT_CROSSFLOW_UY 0.0
#endif

// ------------------------------------------------------------
// VTK 출력 케이던스 [step]. 기본값은 종전 하드코딩과 동일.
// 검증런처럼 전체장을 자주 봐야 할 때만 -D FULL_DUMP_EVERY=1000 식으로 낮춘다
// (전체장 1개 = 43M 입자 x 6 필드 ~ 2.2 GB).
// ------------------------------------------------------------
#ifndef SLICE_DUMP_EVERY
#define SLICE_DUMP_EVERY 100
#endif
#ifndef FULL_DUMP_EVERY
#define FULL_DUMP_EVERY 10000
#endif
#ifndef VORT_DUMP_EVERY
#define VORT_DUMP_EVERY 1000
#endif
#ifndef LDM_DUMP_EVERY
#define LDM_DUMP_EVERY 25      // LDM VTK 출력 주기. 25의 배수여야 한다 (출력 블록이 25 step 게이트 안에 있음)
#endif
// 가벼운 순간 속도장 (functions_OUTPUT.cuh save_transient_velocity). 0 = 끔, 25의 배수.
//   z < TRANSIENT_DUMP_ZMAX 유체 입자의 (ux,uy,uz) float32 만 쓰고 좌표는 처음 한 번만 쓴다.
//   time >= TRANSIENT_DUMP_START 부터. 읽기: python/read_transient.py
#ifndef TRANSIENT_DUMP_EVERY
#define TRANSIENT_DUMP_EVERY 0
#endif
#ifndef TRANSIENT_DUMP_ZMAX
#define TRANSIENT_DUMP_ZMAX 50.0
#endif
#ifndef TRANSIENT_DUMP_START
#define TRANSIENT_DUMP_START 0.0
#endif

#ifndef AVG_WINDOW
#define AVG_WINDOW 30.0
#endif

// Quasi-steady detection (Σ½|u|² sliding-window relative std test).
// Effective window = QS_DETECT_N × QS_SAMPLE_STRIDE × 50 simulation steps.
//   QS_DETECT_N: ring-buffer size (# of samples kept)
//   QS_SAMPLE_STRIDE: how many 50-step diagnostic blocks per sample push
//   (default N=60, STRIDE=4 → 60 × 4 × 50 = 12000 steps = 60s at dt=0.005)
#ifndef QS_DETECT_N
#define QS_DETECT_N 60
#endif
#ifndef QS_SAMPLE_STRIDE
#define QS_SAMPLE_STRIDE 4
#endif
#ifndef QS_DETECT_EPS
#define QS_DETECT_EPS 0.01f       // relative std threshold (1 %)
#endif
// Minimum time below which QS auto-detect cannot latch. Prevents the
// detector from firing on the initial ring-buffer fill window when the
// flow has not yet developed. Set 0 to disable.
#ifndef QS_MIN_TIME
#define QS_MIN_TIME 49.0
#endif
// Fallback: if QS is not detected by this absolute time, force-latch
// h_avg_start = time. Set <=0 to disable fallback.
#ifndef AVG_FORCE_START_TIME
#define AVG_FORCE_START_TIME 50.0
#endif

// Wall-clock-style stop: end the simulation this many seconds after the
// averaging window starts (h_avg_start latched). Overrides time_end.
// Set <=0 to disable (then time_end governs).
#ifndef LDM_DURATION_AFTER_AVG
#define LDM_DURATION_AFTER_AVG 0.0   // disabled: time_end governs (averaging off)
#endif

__device__ unsigned int avg_count = 0;  // # of accumulated samples (incremented by KERNEL_accumulate_avg)

// Dynamic averaging window bounds, updated at runtime.
//   avg_start_dyn = ∞ → not yet started (LDM dormant, no accumulation)
//   avg_end_dyn   = ∞ → no freeze yet (or never, when AVG_WINDOW<=0)
__device__ Real avg_start_dyn = (Real)1e30;
__device__ Real avg_end_dyn   = (Real)1e30;

// Mean-flow LDM mode.
//   0 → 평소대로 SPH transient + LDM (Run A)
//   1 → SPH 솔버 skip, MEAN_FIELD_PATH 로드해서 P1.ux/k_turb/e_turb를 frozen
//        mean field로 채우고 LDM만 풀이 (Run B). AVG_START_TIME gate 우회.
#ifndef LDM_MEAN_MODE
#define LDM_MEAN_MODE 0
#endif
// Run B 가 다른 실험 디렉토리의 mean_field.bin 을 읽을 수 있도록 -D 오버라이드 허용.
#ifndef MEAN_FIELD_PATH
#define MEAN_FIELD_PATH "./" OUTPUT_CASE_NAME "/mean_field.bin"
#endif

// Checkpoint / restart (ISPH_Calc.cuh save_checkpoint / load_checkpoint).
//   CKPT_EVERY > 0 → count%CKPT_EVERY==0 인 스텝 끝마다 CKPT_PATH 에 전체 상태 저장
//                    (이전 파일은 CKPT_PATH.prev 로 한 벌 보관). 단일 GPU 전용.
//   재개: 같은 바이너리를 환경변수 SOPHIA_RESTART=<체크포인트 경로> 로 실행.
#ifndef CKPT_EVERY
#define CKPT_EVERY 0
#endif
#ifndef CKPT_PATH
#define CKPT_PATH "./" OUTPUT_CASE_NAME "/restart.bin"
#endif