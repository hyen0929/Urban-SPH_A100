// //////////////////////////////////////////////////////////////////////////////////
//
// // stability correction
// __device__ Real psi_m(Real zeta) {
//   if (zeta < 0.0f) {
//     Real x = powf(1.0f - Businger_const * zeta, 0.25f);
//     return 2.0f * logf((1.0f + x) * 0.5f);             // Simple Businger-Dyer form
//     // return 2.0f * logf((1.0f + x) * 0.5f) + 0.5f * logf((1.0f + x * x) * 0.5f) - 2.0f * atanf(chi) + 0.5f * M_PI;      // Original form
//   }
//   else {
//     return -5.0f * zeta;
//   }
// }
//
// // stability correction
// __device__ Real psi_h(Real zeta) {
//   if (zeta < 0.0f) {
//     Real y = powf(1.0f - Businger_const * zeta, 0.5f);
//     return 2.0f * logf((1.0f + y * y)*0.5f);
//   }
//   else {
//     return -5.0f * zeta;
//   }
// }
//
// __device__ Real calc_blend_alpha(Real z, Real z1, Real z2){
//   if (z <= z1)        return 1.0f;
//   else if (z >= z2)   return 0.0f;
//   else                return (z2 - z)/(z2 - z1);
// }
//
// // Neutral stability (L → ∞)
// __device__ Real calc_u_star_neutral(Real U, Real z, Real z0) {
//     const Real kappa = 0.4;
//     Real epsilon = 1e-6;
//     Real z_eff = max(z - z0, epsilon);  // avoid log(0)
//     return kappa * U / log(z_eff / z0 + epsilon);
// }
//
// __global__ void KERNEL_MOST_boundary3D(int_t* g_str, int_t* g_end, part1* P1) {
//     uint_t i = threadIdx.x + blockIdx.x * blockDim.x;
//     if (i >= k_num_part2) return;
//     if (P1[i].buffer_type!=6) return; // Only MOST boundary particles
//
//     Real z0 = 0.01; // 조도 길이 [m]
//     Real zmo = P1[i].z;
//
//     // 수평 유속 크기
//     Real uxy = sqrt(P1[i].ux * P1[i].ux + P1[i].uy * P1[i].uy);
//
//     // 마찰 속도 u_star 계산
//     Real u_star = calc_u_star_neutral(uxy, zmo, z0);
//
//     // 마찰력: tau = A*u*^2
//     Real tau = 25.0*u_star*u_star;
//
//     // 방향별 마찰 항 계산
//     Real fx = -tau*(P1[i].ux/(uxy+1e-6));
//     Real fy = -tau*(P1[i].uy/(uxy+1e-6));
//     Real fz = 0.0;
//
//     // 힘에 추가 (fb는 바디포스)
//     P1[i].fbx += fx;
//     P1[i].fby += fy;
//     P1[i].fbz += fz;
//
//     // 디버깅 (선택)
//     // P1[i].PPE1 = u_star;
// }
//
//
// #include <cmath>
// #include <cstdio>
//
// const Real kappa = 0.4f;
// const Real g = 9.81f;
// const Real EPS = 1e-6f;
// const Real PI = 3.14159265f;
//
// // 안정도 보정 함수: 불안정 조건 (z/L < 0)
// Real psi_m_unstable(Real zL) {
//     Real x = powf(1.0f - 16.0f * zL, 0.25f);
//     return 2.0f * logf((1.0f + x) / 2.0f)
//          + logf((1.0f + x * x) / 2.0f)
//          - 2.0f * atanf(x)
//          + PI / 2.0f;
// }
//
// Real psi_h_unstable(Real zL) {
//     Real x = powf(1.0f - 16.0f * zL, 0.25f);
//     return 2.0f * logf((1.0f + x * x) / 2.0f);
// }
//
// // 안정 조건 (z/L > 0)
// Real psi_m_stable(Real zL) {
//     return -5.0f * zL;
// }
//
// Real psi_h_stable(Real zL) {
//     return -5.0f * zL;
// }
//
// // 통합 similarity 함수 (log - psi)
// Real phi_m(Real z, Real z0, Real zL) {
//     Real psi = (zL < 0.0f) ? psi_m_unstable(zL) : psi_m_stable(zL);
//     return logf(z / z0) - psi;
// }
//
// Real phi_h(Real z, Real z0h, Real zL) {
//     Real psi = (zL < 0.0f) ? psi_h_unstable(zL) : psi_h_stable(zL);
//     return logf(z / z0h) - psi;
// }
//
// // 뉴턴 반복으로 Obukhov 길이 계산 (Dirichlet 조건 기준)
// Real compute_L_newton_dirichlet(Real u_h, Real theta_v_mo, Real theta_v_0, Real theta_ref, Real z_mo,
//                                  Real z0, Real z0h, int max_iter = 20) {
//     Real Ri_b = (g * z_mo * (theta_v_mo - theta_v_0)) / (theta_ref * (u_h * u_h) + EPS);
//     Real L = z_mo / (Ri_b + EPS);  // 초기 추정
//     Real L_prev;
//
//     for (int i = 0; i < max_iter; ++i) {
//         Real zL = z_mo / (L + EPS);
//         Real phiM = phi_m(z_mo, z0, zL);
//         Real phiH = phi_h(z_mo, z0h, zL);
//
//         Real f = Ri_b - (z_mo / L) * (phiH / (phiM * phiM));
//
//         // 수치 미분으로 f' 근사
//         Real delta = 0.01f * fabsf(L);
//         Real L_eps = L + ((delta > EPS) ? delta : EPS);
//         Real zL_eps = z_mo / L_eps;
//         Real phiM_eps = phi_m(z_mo, z0, zL_eps);
//         Real phiH_eps = phi_h(z_mo, z0h, zL_eps);
//         Real f_eps = Ri_b - (z_mo / L_eps) * (phiH_eps / (phiM_eps * phiM_eps));
//
//         Real df = (f_eps - f) / (L_eps - L + EPS);
//
//         L_prev = L;
//         L = L - f / (df + EPS);
//
//         if (fabsf(L - L_prev) < 1e-3f) break;
//     }
//
//     return L;
// }
////////////////////////////////////////////////////////////////////////////////
__device__ float clamp(float x, float minVal, float maxVal) {
    return fmaxf(minVal, fminf(x, maxVal));
}
// Stability functions (Psi functions) for momentum and heat
__device__ double psi_m(Real zeta) {
    if (zeta < 0.0f) {
        Real x = sqrtf(sqrtf(1.0f - 16.0f * zeta));
        return M_PI * 0.5f - 2.0f * atanf(x) + logf((1.0f + x) * (1.0f + x) * (1.0f + x * x) * 0.125f);
    }
    return -5.0f * zeta;
}

__device__ double psi_h(Real zeta) {
    if (zeta < 0.0f) {
        return 2.0f * logf(0.5f * (1.0f + sqrtf(1.0f - 16.0f * zeta)));
    }
    return -5.0f * zeta;
}

// Integrated stability functions (Phi functions) for momentum and heat
__device__ double Phi_m(Real L, Real z, Real z0) {
    return logf(z/z0)-psi_m(z/L)+psi_m(z0/L);
}

__device__ double Phi_h(Real L, Real z, Real z0h) {
    return logf(z/z0h)-psi_h(z/L)+psi_h(z0h/L);
}

__global__ void KERNEL_MOST_boundary3D(int_t*g_str, int_t*g_end, part1*P1, part3*P3, int_t tcount, Real tdt) {
    uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
    if(i>=k_num_part2) return;
    if(P1[i].MOST_buffer!=1) return;
    // if(P1[i].buffer_type!=0) return;      // Inlet/Outlet second layer 제외

    // Particle properties
    int_t icell,jcell,kcell;
    Real xi = P1[i].x;
    Real yi = P1[i].y;
    Real zi = P1[i].z;
    Real uxi = P1[i].ux;
    Real uyi = P1[i].uy;
    Real uzi = P1[i].uz;
    Real tempi = P1[i].temp;
    Real rhoi = P1[i].rho;
    Real mi = P1[i].m;
    Real space = P1[i].h/h_coeff;

    // Initialize MOST variables
    Real z0 = MOST_z0;
    Real z0h = MOST_z0h;
    Real u_star = 0.0;
    Real theta_star = 0.0;
    Real OB_length;
    int iter=0;
    int skip_flag=0;

    Real theta_surface;
    if(tdt*tcount>=T_start_t && tdt*tcount<T_con_t) theta_surface=T_ini+tdt*tcount*(T_con-T_ini)/(T_con_t-T_start_t);
    else if(tdt*tcount<T_start_t) theta_surface=T_ini;
    else theta_surface=T_con;
    // Calculate wind speed
    Real wind_speed = sqrtf(uxi*uxi+uyi*uyi);               // Horizontal wind speed
    wind_speed = clamp(wind_speed, 0.001, 100.0);

    Real ol_m, ol_l, ol_u;
    Real damping = 1.0;
    ol_m = ol_l = ol_u = 0.0;

    Real Ri_b = (Gravitational_CONST/tempi)*(tempi-theta_surface)*zi / (wind_speed*wind_speed+1e-10);
    // Initialize Obukhov length with a small value
    if(abs(Ri_b)>1e-3){
      OB_length = zi / Ri_b;  // Ri = z/L → L ≈ z / Ri
      OB_length = clamp(OB_length, -500.0, 500.0);  // 실측 기반 제한

      //if (P1[i].x>19&&P1[i].x<21&&P1[i].y>4&&P1[i].y<6) printf("Ri_b=%f, z0=%f, z0h=%f, OB_length=%f, wind=%f\n",Ri_b,z0,z0h,OB_length,wind_speed);
      // Newton iteration for OB_length using Rib
      if (tcount%10==0 && (P1[i].x>5&&P1[i].x<10&&P1[i].y>4&&P1[i].y<6)) printf("Initial L=%f Wind_speed=%f\n",OB_length,wind_speed);
      for (iter=0; iter<MAX_ITER; iter++) {
        if(skip_flag==1) break;
        ol_m = OB_length;
        ol_l = 0.999*ol_m;          // Perturbation: L-dL
        ol_u = 1.001*ol_m;          // Perturbation: L+dL

        Real fL = Ri_b-(zi/ol_m)*Phi_h(ol_m,zi,z0h) / powf(Phi_m(ol_m,zi,z0),2);              //  f(L) = Ri_b - Ri_MOST
        Real dfL = ((-(zi/ol_u)*Phi_h(ol_u,zi,z0h)/(Phi_m(ol_u,zi,z0)*Phi_m(ol_u,zi,z0))) - (-(zi/ol_l)*Phi_h(ol_l,zi,z0h)/(Phi_m(ol_l,zi,z0)*Phi_m(ol_l,zi,z0))) ) / (ol_u-ol_l);              // f'(L) = (f(L+dL) - f(L-dL))/2dL

        OB_length -= damping * fL/dfL;
        // if ((P1[i].x>19&&P1[i].x<21&&P1[i].y>4&&P1[i].y<6)&&iter==0 ) printf("ol_m=%f, ol_l=%f, ol_u=%f, fl=%f, dfl=%f\n",ol_m, ol_l, ol_u, fL, dfL);

        if (fabsf(OB_length) > 1e-6){
          if (fabsf((OB_length - ol_m) / OB_length) < TOLERANCE) break;
          else{
            if (fabsf(OB_length - ol_m) < TOLERANCE) break;
          }
        }
      }
      if (tcount%10==0 && (P1[i].x>0&&P1[i].x<10&&P1[i].y>4&&P1[i].y<6)) printf("Ri=%f, Iter=%d, L=%f\n",Ri_b,iter,OB_length);

      // Calculate friction velocity (u*) using updated OB_length, Calculate temperature scale (theta*)
      u_star = k_vonKarman*wind_speed / Phi_m(OB_length,zi,z0);
      theta_star = k_vonKarman*(tempi-theta_surface) / Phi_h(OB_length,zi,z0h);
    }
    else{           // 중립조건
      OB_length=1e6;
      u_star = k_vonKarman*wind_speed / logf(zi/z0);
      theta_star = k_vonKarman*(tempi-theta_surface) / logf(zi/z0h);
    }

    // Calculate momentum flux (tau) using wind speed
    float momentum_x = -rhoi*u_star*(uxi*u_star/wind_speed);        // N/m2 (kg/ms2)
    float momentum_y = -rhoi*u_star*(uyi*u_star/wind_speed);        // N/m2 (kg/ms2)

    // Calculate surface fluxes
    Real sensible_heat_flux = -rhoi*Air_cp*u_star*theta_star*100.0;      // W/m2(J/m2s)

    // calculate ghost node
    Real xgnode_i, ygnode_i, zgnode_i;
  	xgnode_i=xi;
  	ygnode_i=yi;
  	zgnode_i=zi+5.0;

    // calculate I,J,K in cell
  	if((k_x_max==k_x_min)){icell=0;}
  	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
  	if((k_y_max==k_y_min)){jcell=0;}
  	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
  	if((k_z_max==k_z_min)){kcell=0;}
  	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}

  	// out-of-range handling
  	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

  	for(int_t z=-1;z<=1;z++){
  		for(int_t y=-1;y<=1;y++){
  			for(int_t x=-1;x<=1;x++){
  				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
  				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
  				if(g_str[k]!=cu_memset){
  					int_t fend=g_end[k];
  					for(int_t j=g_str[k];j<fend;j++){
  						Real xj,yj,zj,tdist;

              xj=P1[j].x;
              yj=P1[j].y;
              zj=P1[j].z;

              tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
              if(tdist<0.1){
                // P1[j].fmx = 0.0*momentum_x*space*space/mi;       // *A/m (m/s2)
                // P1[j].fmy = 0.0*momentum_y*space*space/mi;
                // P1[j].vis_t = 0.4*k_vonKarman*u_star*(zi-0.0);
                P3[j].dtemp += 0.5*sensible_heat_flux/(space*rhoi*Air_cp);
                P1[j].hmz = 0.5*sensible_heat_flux/(space*rhoi*Air_cp);;
              }
            }
          }
        }
      }
    }

    P1[i].fmx = 1.0*momentum_x*space*space/mi;       // *A/m (m/s2)
    P1[i].fmy = 1.0*momentum_y*space*space/mi;

    //P1[i].vis_t += k_vonKarman*u_star*(zi-0.0)*rhoi;

    P3[i].dtemp += 0.5*sensible_heat_flux/(space*rhoi*Air_cp);

    P1[i].hmz = 0.5*sensible_heat_flux/(space*rhoi*Air_cp);

    if((xi-600.0)*(xi-600.0)+yi*yi<1e-3) us_update=u_star;

    // Store results in particle
    P1[i].float1 = u_star;
    P1[i].float2 = theta_star;
    P1[i].float3 = OB_length;
}

__global__ void KERNEL_turbulent_viscosity3D(int_t*g_str, int_t*g_end, part1*P1, part3*P3, int_t tcount, Real tdt) {
    uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
    if(i>=k_num_part2) return;
    if(P1[i].p_type!=1) return;      // Inlet/Outlet second layer 제외

    // Particle properties
    int_t icell,jcell,kcell;
    Real xi = P1[i].x;
    Real yi = P1[i].y;
    Real zi = P1[i].z;
    Real rhoi = P1[i].rho;

    Real u_star=0.0;
    // calculate ghost node
    Real xgnode_i, ygnode_i, zgnode_i;
  	xgnode_i=600.0;
  	ygnode_i=0.0;
  	zgnode_i=2.5;

    // calculate I,J,K in cell
  	if((k_x_max==k_x_min)){icell=0;}
  	else{icell=min(floor((xgnode_i-k_x_min)/k_dcell),k_NI-1);}
  	if((k_y_max==k_y_min)){jcell=0;}
  	else{jcell=min(floor((ygnode_i-k_y_min)/k_dcell),k_NJ-1);}
  	if((k_z_max==k_z_min)){kcell=0;}
  	else{kcell=min(floor((zgnode_i-k_z_min)/k_dcell),k_NK-1);}

  	// out-of-range handling
  	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

  	for(int_t z=-1;z<=1;z++){
  		for(int_t y=-1;y<=1;y++){
  			for(int_t x=-1;x<=1;x++){
  				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
  				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
  				if(g_str[k]!=cu_memset){
  					int_t fend=g_end[k];
  					for(int_t j=g_str[k];j<fend;j++){
  						Real xj,yj,zj,tdist;

              xj=P1[j].x;
              yj=P1[j].y;
              zj=P1[j].z;

              tdist=sqrt((xgnode_i-xj)*(xgnode_i-xj)+(ygnode_i-yj)*(ygnode_i-yj)+(zgnode_i-zj)*(zgnode_i-zj))-1e-20;
              if(tdist<0.1){
                u_star=P1[j].float1;
              }
            }
          }
        }
      }
    }
    //P1[i].vis_t += k_vonKarman*u_star*(zi-0.0)*rhoi*(1-zi/ABL_h)*(1-zi/ABL_h)*(zi<ABL_h);
    //P1[i].vis_t += k_vonKarman*u_star*(zi-0.0)*rhoi*(zi<ABL_h);
}
