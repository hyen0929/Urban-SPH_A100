#define		cu_memset				0xffffffff
#define 	vii_size				64
#define 	vif_size				32
#define		Max_GPU					10

#define 	Wcsph 					0
#define 	Isph 					  1

#define 	FLUID					  1
#define   AIR             2
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
#define   C_pres        0.2       // Pressure_smoothing_coefficient

// viscous model
#define     Cleary                  1
#define     Monaghan                2
#define     Violeau                  3

#define 	Lm						0.01
#define 	Laminar					0
#define 	K_LM					1
#define 	K_E						2
#define 	HB      			3
#define   SSM           4
#define   DDF           5

#define 	DIFF_DENSITY		    0
#define   Correction_Matrix_Size  3

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

#define Boundary_Velocity 1.0 // No-slip condition for moving boundary

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

#define ncell_init 2
#define Cell_division_factor 4
#define IBM_length 0.5
#define Pressure_length 1.0
#define XSPH_length 0.66667

/// MOST ///
#define MAX_ITER 200
#define TOLERANCE 1e-5
#define k_vonKarman 0.4
#define Businger_const 16.0   // For MOST psi function
#define MOST_z0 1.0           // Roughness length for city (>=2)
#define MOST_z0h 0.1          // Roughness length for heat
#define us_ini 0.2            // Friction velocity (m/s)
#define ABL_h 5.0
#define Pr_t 0.71             // Turbulent pr number

#define Solid_cp 900.0
#define Air_cp 1005         // J/kgK
#define Solid_k 1.0
#define Air_k 0.025
#define Solid_rho 2300.0

#define T_ini 300.0       // Initial T
#define T_con 300.0       // Converge T
#define T_start_t 0   // Start time
#define T_con_t 10.0    // Converge time
#define Pb 0.0
#define tempb 300.0     // Temp0 for boussinesq

#define Ck_sgs 0.1
#define Ce_sgs 0.93
#define Cs_sgs 0.15
#define Prt_sgs 0.7
#define Sigk_sgs 1.0

extern __device__ float us_update;
