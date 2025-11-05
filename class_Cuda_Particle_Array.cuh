// Cparticle Class Declaration
// Cparticle class contains particle information.
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
////////////////////////////////////////////////////////////////////////
typedef struct particles_array_1{
	uint_t i_type;													// Inner or Outer
	uint_t buffer_type;												// buffer_type (0: active , 1: inlet, 2: outlet)
	uint_t MOST_buffer;
	uint_t p_type;													// particle type: FLUID or BOUNDARY

	Real x,y,z;														// (Predicted) positions [m] ( Predictor_Corrector : Predicted position / Euler : Real Position )
	Real ux,uy,uz;													// (Predicted) velocity [m/s] ( Predictor_Corrector : Predicted velocity / Euler : Real Velocity )
	Real m;																	// mass [kg]
	Real h;																	// kernel distance
	Real temp;														// temperature [K]
	Real pres;														// pressure [Pa]
	Real pres0;
	Real rho;														// density [kg/m3]	( Predictor_Corrector : Predicted density / Euler : Real density )
	Real flt_s;														// Shepard filter
	Real enthalpy;
	Real concn;
	Real concentration;

	// turbulence (by esk)
	Real k_turb,e_turb;												// turbulence kinetic energy,dissipation rate --> check unit

	// ALE (by hsy)
 	Real eli;
	Real vol;

	// APR (by dhk)
	// int_t apr_cond;
	// int_t merge_flag;
	Real m_ref;															// reference mass [kg] (For APR)
	Real h_ref;															// reference smoothing length [m] (for APR)
	// Real M_num;
	// int_t blood_index;

	// LES (by TSC, KDH)
	Real vis_t;															// nu
	// Real a11,a12,a13,a21,a22,a23,a31,a32,a33;
	// Real st11,st12,st22,st23,st31,st33;
	// Real ux_f,uy_f,uz_f;										// u', v', w' (for LES by KDH)

	// for PPE
	// Real fpx, fpy, fpz;
	Real fbx, fby, fbz;
	Real fmx, fmy, hmz;
	// Real PPE1, PPE2;

	// Real fp_x,fp_y,fp_z;
	// Real fv_x,fv_y,fv_z;
	// Real fe_x,fe_y,fe_z;
	// Real fb_x,fb_y,fb_z;

	Real float1, float2, float3;

	int_t ncell;
	int_t source;

}part1;
////////////////////////////////////////////////////////////////////////
typedef struct particles_array_2{
	Real rho_ref;

	//// turbulence (by esk)
	Real SR;																// strain rate (2S:S)

	Real x0,y0,z0;													// Initial positions [m]
	Real ux0,uy0,uz0;												// Initial velocity [m/s]
	Real rho0;															// Initial density [kg/m3]
	Real drho0;															// Error compensation: divergence error

	Real temp0;

	// psh: concentration diffusion
	Real concn0;														// concentration
	Real enthalpy0;													// enthalpy [J/kg]
}part2;
////////////////////////////////////////////////////////////////////////
typedef struct particles_array_3{
	Real drho;															// Time Derivative of density [kg/m3 s]
	Real dconcn;														// concentration time derivative
	Real denthalpy;
	Real ftotalx,ftotaly,ftotalz;						// total force [m/s2]
	Real ftotal;
	Real fpx,fpy,fpz,fbody;
	Real dtemp;

	Real omx,omy,omz;             // Vorticity (kdh)

	Real dk_turb, de_trub;        // Turbulent kinetic energy, dissipation rate(for Deardorff model,KDH)
	// Real virial,radius;

	Real cc;
	Real nx, ny, nz;																// color code
	// Real lambda;

	Real Cm[Correction_Matrix_Size][Correction_Matrix_Size];
	// Real A[Correction_Matrix_Size][Correction_Matrix_Size];
	Real inv_cm_xx,inv_cm_yy,inv_cm_zz;
	Real inv_cm_xy,inv_cm_yz,inv_cm_zx;
}part3;
///////////////////////////////////////////////////////////////////////////
typedef struct p2p_particles_array_3{
	Real drho;															// Time Derivative of density [kg/m3 s]
	Real dconcn;														// concentration time derivative
	Real denthalpy;
	Real ftotalx,ftotaly,ftotalz;						// total force [m/s2]
	Real ftotal;
}p2p_part3;

typedef struct LDM_particles_array_1{
	int_t i_type;													// Inner or Outer
	uint_t buffer_type;											// buffer_type (0: active , 1: inlet, 2: outlet)
	uint_t p_type;													// particle type: FLUID or BOUNDARY
	uint_t ct_boundary;												// const temperatrue particle index (const temp particle :1 / else : 0)

	Real x,y,z;															// (Predicted) positions [m] ( Predictor_Corrector : Predicted position / Euler : Real Position )
	Real ux,uy,uz;													// (Predicted) velocity [m/s] ( Predictor_Corrector : Predicted velocity / Euler : Real Velocity )
	Real m;																	// mass [kg]
	Real m_ref;
	Real h;																	// kernel distance
	Real h_ref;
	Real temp;															// temperature [K]
	Real pres;															// pressure [Pa]
	//Real pres_ave;													// averaged pressure
	Real rho;																// density [kg/m3]	( Predictor_Corrector : Predicted density / Euler : Real density )

	Real flt_s;															// Shepard filter
	Real w_dx;															// w(dx) for particle shifting
	Real enthalpy;
	Real concn;
	Real sigma;

	Real grad_rhox,grad_rhoy,grad_rhoz;				// density gradient - only mass kernel using

	// turbulence (by esk)
	Real k_turb, e_turb;												// turbulence kinetic energy,dissipation rate --> check unit
	Real fx_cylinder, fy_cylinder;

	Real uxr, uyr, uzr;
	Real xwind, ywind, zwind;
	Real kver, khor;
	Real sigu, sigv, sigw, dsigw;
	Real disp_rho, disp_drho;

	Real concentration;

	Real yplus;
	Real Tscale;
	Real Sigvel;
	unsigned long long seed;
	// Real dx_pst, dy_pst, dz_pst; 						// particle shifting displacement
	uint_t ncell;
	// Real nx_s,ny_s,nz_s;									// surface normal vector --> not used
}L_part1;
////////////////////////////////////////////////////////////////////////
typedef struct LDM_particles_array_2{
	Real rho_ref;

	//// turbulence (by esk)
	Real SR;																	// strain rate (2S:S)

	Real x0,y0,z0;													// Initial positions [m]
	Real ux0,uy0,uz0;												// Initial velocity [m/s]
	Real rho0;															// Initial density [kg/m3]
	Real drho0;														// Error compensation: divergence error
	Real sigma0;
	// psh: concentration diffusion
	Real concn0;														// concentration
	Real enthalpy0;													// enthalpy [J/kg]

	Real uxr0,uyr0,uzr0;
}L_part2;
////////////////////////////////////////////////////////////////////////
typedef struct LDM_particles_array_3{
	Real drho;															// Time Derivative of density [kg/m3 s]
	Real dconcn;														// concentration time derivative
	Real denthalpy;
	Real ftotalx,ftotaly,ftotalz;						// total force [m/s2]
	Real ftotal;

	// turbulence (by esk)
	Real vis_t;																// turbulence viscosity
	Real Sxx,Sxy,Sxz,Syy,Syz,Szz;							// strain tensor... for SPH model
	Real dk_turb,de_turb;											// turbulence kinetic energy,dissipation rate --> check unit

	Real lbl_surf;													// surface label
	Real cc;																// color code

	Real Cm[Correction_Matrix_Size][Correction_Matrix_Size];
	// Real inv_cm_xx,inv_cm_yy,inv_cm_zz;
	// Real inv_cm_xy,inv_cm_yz,inv_cm_zx;

	Real nx,ny,nz;														// color code gradient for surface tension force (2016.09.02 jyb)
	Real nx_c,ny_c,nz_c;											// color code gradient for surface tension force (2016.09.02 jyb)
	Real nmag;																// 2017.04.20 jyb
	Real nmag_c;															// 2017.04.20 jyb
	/*
	// psh:: ISPH info
	Real fpx,fpy,fpz;													// pressure force [m/s2]
	Real x_adv,y_adv,z_adv;										// predicted position by advection forces
	Real rho_err;															// difference between predicted density and reference density
	Real stiffness;														// stiffness parameter
	//----------------------------------------

	uint_t hf_boundary;											// heat flux particle index for heat transfer (heat flux particle : 1 / else : 0)
	Real cm_d; 															//--- not use
	Real p001;															// extra data	//--- not use
	Real curv; 															//--- not use
	Real num_density;												// particle number density [1/m3] //--- not use

	Real nd_ref;
	Real cm_xx,cm_yy,cm_zz;
	Real cm_xy,cm_yz,cm_zx;
	//*/


}L_part3;
///////////////////////////////////////////////////////////////////////////
typedef struct LDM_p2p_particles_array_1{
	Real drho;															// Time Derivative of density [kg/m3 s]
	Real dconcn;														// concentration time derivative
	Real denthalpy;
	Real ftotalx,ftotaly,ftotalz;						// total force [m/s2]
	Real ftotal;
}L_p2p_part1;

typedef struct cfd_mesh{
	Real conc;
	Real e_turb, k_turb;
	Real u, v, w;
	Real yplus;
	int_t n;
}cfd_mesh;