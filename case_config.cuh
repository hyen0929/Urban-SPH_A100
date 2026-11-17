// case_config.h
#pragma once

// ============================================================
// Case selector
// ============================================================
#ifndef CASE_ID
#define CASE_ID 1010
#endif

// ============================================================
// Case definitions
//  규칙:
//   - Case A 계열: 1, 2, 3 ...
//   - Case H 계열: 101, 102, 103 ...
//   - Case C 계열: 201 ...
//   - Multi block 계열: 301 ...
//   - Case CS (Single building block) 계열: 501 ...
//   - Leipzig 계열: 1001 ...
//   - Transport test 계열: 2001 ...
//   - MOST Validation 계열: 3001 ...
// ============================================================

#if CASE_ID == 101
// ------------------------------------------------------------
// Case A, n=10
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.028
  #define L3  13.304
  #define L4  2210.0

  #define L_left   -0.412
  #define L_right   0.412
  #define L_ceiling 1.910

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 102
// ------------------------------------------------------------
// Case A, n=12
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.02333
  #define L3  13.304
  #define L4  2210.0

  #define L_left   -0.416667
  #define L_right   0.416667
  #define L_ceiling 1.910

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 103
// ------------------------------------------------------------
// Case A_n16 (dx=0.005m)
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.0125
  #define L3  10.0
  #define L4  100.0

  #define L_left   -0.4275
  #define L_right   0.4275
  #define L_ceiling 0.965

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 201
// ------------------------------------------------------------
// Case AS, n=10, Scale up (dx=1.0m, MOST at 3rd layer z=2.5m)
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  3.5
  #define L3  2013.304
  #define L4  2210.0

  #define L_left   -62.5
  #define L_right   62.5
  #define L_ceiling 121.0

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0
  #define L_MOST_z_lo 0.0
  #define L_MOST_z_hi 1.0

// ------------------------------------------------------------
// Single Building - Stability comparison (domain same as 201)
//   211: Unstable  (T_con=315K, ΔT=+15K)
//   212: Neutral   (T_con=300K, ΔT=  0K)
//   213: Stable    (T_con=285K, ΔT=-15K)
// dx=1.0m, MOST at 3rd fluid layer (z=2.5m)
// ------------------------------------------------------------
#elif CASE_ID == 211 || CASE_ID == 212 || CASE_ID == 213
  #define L1  0.0
  #define L2  3.5
  #define L3  2013.304
  #define L4  2210.0

  #define L_left   -76.5
  #define L_right   76.5
  #define L_ceiling 121.0

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0
  #define L_MOST_z_lo 0.0
  #define L_MOST_z_hi 1.0

#elif CASE_ID == 301
// ------------------------------------------------------------
// Case C (dx=0.01m)
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.035
  #define L3  10.0
  #define L4  100.0

  #define L_left   -1.46
  #define L_right   1.46
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 302
// ------------------------------------------------------------
// Case C (dx=0.0133333m)
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.033333333
  #define L3  10.0
  #define L4  100.0

  #define L_left   -1.46
  #define L_right   1.46
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 401
// ------------------------------------------------------------
// Multi block (Case RS, n=10, Scale up)  // 나머지 계열: 301로 임의 배정
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  3.5
  #define L3  2013.304
  #define L4  2210.0

  #define L_left   -132.5
  #define L_right   132.5
  #define L_ceiling 148.0

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 501
// ------------------------------------------------------------
// Case CS (Single building block, dx=1.0m)
//   - MOST off; inlet sustained via KERNEL_Recycling3D (log profile w/ us_ini)
//   - Outlet buffer 미사용 (L3/L4 도메인 밖) → 도메인 충분히 길게 빼서 wake가 끝에 못 닿게
//   - Ceiling은 input의 p_type=-1 layer가 free-slip 처리
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  3.5
  #define L3  246.5
  #define L4  251.0

  #define L_left   -71.5
  #define L_right   71.5
  #define L_ceiling 3633.75

  #define L_MOST_x0   10.0
  #define L_MOST_xmax 200.0

  // MOST ON for CS (Parameters.cuh의 default 1 그대로 사용)
  //  - geostrophic body force가 도메인 전체 momentum 유지
  //  - heat flux는 conduction OFF라 무영향

  // CS는 building block이 핵심 → IBM 활성화 보장
  #undef  DEACTIVATE_IBM
  #define DEACTIVATE_IBM 0

// ------------------------------------------------------------
// Leipzig (dx=2.5m)
// ------------------------------------------------------------
#elif CASE_ID == 1001

  #define L1  -900.0
  #define L2  -817.5
  #define L3  2200.0
  #define L4  2210.0

  #define L_left   -691.25
  #define L_right   691.25
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// Leipzig (Point surface input, dx=2.5m)
// ------------------------------------------------------------  
#elif CASE_ID == 1002

  // Leipzig_surface / Leipzig_surface_full particles span
  //   x[0, 1315]  y[-657.5, 657.5]  z[-8.75, 371.25]   (527 x 527 x 146, dx=2.5)
  // with h=3.75 and kappa=2, i.e. a kernel support radius of 7.5 m = 3 planes.
  // Every buffer band below must contain real particles: the matching BC kernel
  // early-returns on buffer_type, so an empty band leaves that face running as
  // plain interior fluid with truncated support (2026-09-10: that is what pinned
  // the x=1315 and y=+657.5 planes at the 20 m/s velocity limit).
  #define L1  0.0
  #define L2  8.75      // inlet buffer   x < L2       : 4 planes, 0 ... 7.5
  #define L3  1307.5    // outlet buffer  L3 <= x < L4 : 4 planes, 1307.5 ... 1315
  #define L4  1317.5

  // y-periodic: period = L_right - L_left = 1297.5 = 519 planes x 2.5, so the
  // buffer planes +-650 ... +-657.5 map exactly onto interior planes -+640 ... -+647.5.
  // KERNEL_periodic_boundary3D copies by exact position (tol 1e-3) and zeroes the
  // particle when no image is found, so the period must stay an exact multiple of dx.
  #define L_left   -648.75
  #define L_right   648.75
  #define L_ceiling 371.25

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// CEDVAL Complexity 2 (BL3-2-1 엇갈림 큐브 배열, dx=2.5m) — 확산 검증
//   input/CEDVAL_C2 (validation/cedval/make_input.py)
//   x[0,815] y[-367.5,367.5] z[-8.75,216.25]; 버퍼·주기 규칙은 1002 와 같음
// ------------------------------------------------------------
#elif CASE_ID == 1020
  #define L1  0.0
  #define L2  8.75
  #define L3  807.5
  #define L4  817.5
  #define L_left   -358.75     // 주기 717.5 = 287 dx
  #define L_right   358.75
  #define L_ceiling 216.25
  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// CEDVAL Complexity 2, dx=1.5m — 배열을 y 주기 432 m 로 무한 반복, 도메인 축소
//   input/CEDVAL_C2_dx15 (validation/cedval/make_input.py CEDVAL_C2_dx15)
//   x[0,669] y[-221.25,221.25] z[-5.25,195.75]; SPH x = CEDVAL x + 341.25
// ------------------------------------------------------------
#elif CASE_ID == 1021
  #define L1  0.0
  #define L2  5.25      // 유입 버퍼 4 면: 0 … 4.5
  #define L3  663.75    // 유출 버퍼 4 면: 664.5 … 669
  #define L4  670.5
  #define L_left   -216.0      // 주기 432 = 288 dx, 버퍼 ±216.75 … ±221.25
  #define L_right   216.0
  #define L_ceiling 195.75
  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// Leipzig (Point input, dx=2.5m, small)
// ------------------------------------------------------------
#elif CASE_ID == 1003

  #define L1  0.0
  #define L2  8.75
  #define L3  2200.0
  #define L4  2210.0

  #define L_left   -248.75
  #define L_right   248.75
  #define L_ceiling 15463.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// Leipzig (Interior input, dx=2.5m)
// ------------------------------------------------------------
  #elif CASE_ID == 1004

  #define L1  0.0
  #define L2  8.75
  #define L3  2200.0
  #define L4  2210.0

  #define L_left   -648.75
  #define L_right   648.75
  #define L_ceiling 15463.75

  #define L_MOST_x0   15.0
  #define L_MOST_xmax 5990.0

// ------------------------------------------------------------
// Leipzig (Interior input, small, dx=2.5m)
// ------------------------------------------------------------
  #elif CASE_ID == 1005

  #define L1  0.0
  #define L2  8.75
  #define L3  2200.0
  #define L4  2210.0

  #define L_left   -243.75
  #define L_right   243.75
  #define L_ceiling 15463.75

  #define L_MOST_x0   15.0
  #define L_MOST_xmax 1990.0

// ------------------------------------------------------------
// Gori (Surface input, dx=5.0m)
// ------------------------------------------------------------
  #elif CASE_ID == 1010

  #define L1  0.0
  #define L2  17.5
  #define L3  2855.0
  #define L4  2865.0

  #define L_left   -1272.5
  #define L_right   1272.5
  #define L_ceiling 15463.75

  #define L_MOST_x0   15.0
  #define L_MOST_xmax 1990.0

  // Kori는 건물(IBM)이 핵심 → IBM 활성화 보장
  #undef  DEACTIVATE_IBM
  #define DEACTIVATE_IBM 0

// ------------------------------------------------------------
// Scalar Transport test (dx=0.01m)
// ------------------------------------------------------------
  #elif CASE_ID == 2001   // Diffusion only

  #define L1  0.0
  #define L2  0.35
  #define L3  10.0
  #define L4  100.0

  #define L_left   -4.65
  #define L_right   4.65
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0


  #elif CASE_ID == 2002   // Advection-Diffusion

  #define L1  0.0
  #define L2  0.35
  #define L3  100.0
  #define L4  100.0

  #define L_left   -2.65
  #define L_right   2.65
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

  #elif CASE_ID == 3001   // MOST Validation (flat terrain, dx=2.5m)

  #define L1  0.0
  #define L2  8.75
  #define L3  1000.0
  #define L4  1000.0

  #define L_left   -53.75
  #define L_right   53.75
  #define L_ceiling 363.75

  #define L_MOST_x0   10.0      // Past inlet buffer (L2=8.75)
  #define L_MOST_xmax 500.0     // Fluid domain end

#else
  #error "Unknown CASE_ID. Please define a valid CASE_ID in case_config.h or via -DCASE_ID=..."
#endif

// MOST measurement height fallback (override per-case above if needed)
#ifndef L_MOST_z_lo
  #define L_MOST_z_lo 0.4
#endif
#ifndef L_MOST_z_hi
  #define L_MOST_z_hi 0.6
#endif
