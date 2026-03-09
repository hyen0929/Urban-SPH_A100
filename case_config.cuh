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
// Case A_n16 (dx=0.01m)
// ------------------------------------------------------------
  #define L1  0.0
  #define L2  0.0125
  #define L3  10.0
  #define L4  100.0

  #define L_left   -0.4275
  #define L_right   0.4275
  #define L_ceiling 363.75

  #define L_MOST_x0   0.0
  #define L_MOST_xmax 1990.0

#elif CASE_ID == 201
// ------------------------------------------------------------
// Case AS, n=10, Scale up
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

  #define L1  0.0
  #define L2  8.75
  #define L3  2200.0
  #define L4  2210.0

  #define L_left   -678.75
  #define L_right   678.75
  #define L_ceiling 15463.75

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
  #define L3  3200.0
  #define L4  3210.0

  #define L_left   -1272.5
  #define L_right   1272.5
  #define L_ceiling 15463.75

  #define L_MOST_x0   15.0
  #define L_MOST_xmax 1990.0

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
