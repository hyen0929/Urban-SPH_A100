// ====================================================================================
// Cell-based APR
// ====================================================================================

__global__ void KERNEL_APR_condition2D_cell(part1*P1,int_t*apr_cell,int_t*g_str,int_t*g_end,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //cell number

  if (i>k_num_cells-1) return;
  if (g_str[i]==-1) return;

  // Real x_cell=0.0;
  // Real y_cell=0.0;
  // Real vmag_cell=0.0;
  // Real apr_level=0.0;
  int_t flag=0;
  int_t nop=g_end[i]-g_str[i];
  Real M_val=0.0;
  int_t merge_flag=0;

  for(int_t j=g_str[i];j<g_end[i];j++){
    if(P1[j].p_type==1) break;
    if(j==g_end[i]-1) return;
  }

  Real mass_ratio=0.0;
  int_t merge_index=0;
  int_t blood=0;
  for(int_t j=g_str[i];j<g_end[i];j++){
    //if(P1[j].p_type!=1) continue;
    mass_ratio=log(P1[j].m_ref/P1[j].m)/log(4.0);
    if(P1[j].M_num-mass_ratio >= 0.75){
      P1[j].apr_cond=1;
      apr_cell[i]=1;
      return;
    }
    // else if(P1[j].M_num-mass_ratio >= 0.3){
    //   P1[j].apr_cond=2;
    //   apr_cell[i]=1;
    //   return;
    // }
    else if(P1[j].M_num-mass_ratio < -0.3){
      merge_index=j;
      blood=P1[j].blood_index;
      merge_flag=1;
    }
  }

  if(merge_flag==1){                   // Setting blood_index
    Real xi,yi,hi;
    int_t icell,jcell;
    Real search_range;

    xi=P1[merge_index].x;
    yi=P1[merge_index].y;
    hi=P1[merge_index].h;

    if((k_x_max==k_x_min)){icell=0;}
    else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
    if((k_y_max==k_y_min)){jcell=0;}
    else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}

    // out-of-range handling
    if(icell<0) icell=0;	if(jcell<0) jcell=0;

    search_range=k_search_kappa*hi;

    int_t merge_count=0;
    for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
  		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
    		int_t k=idx_cell(icell+x,jcell+y,0);
        if (k<0) continue;
    		if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;

        if(g_str[k]!=cu_memset){
    			int_t fend=g_end[k];
    			for(int_t j=g_str[k];j<fend;j++){
            Real xj,yj,tdist;
    				xj=P1[j].x;
    				yj=P1[j].y;

    				tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;
            if(tdist<search_range){
              if(P1[j].i_type==1 && P1[j].blood_index==blood){
                if(P1[j].M_num==1) return;
                merge_count ++;
                if(j>merge_index) {
                  merge_index=j;
                }
              }
            }
          }
        }
      }
    }
    if(merge_count==4) P1[merge_index].apr_cond=-1;         // merge_index is the largest particle index in the blood group
    apr_cell[i]=-1;
    return;
  }

  apr_cell[i]=0;
  return;

}

__global__ void KERNEL_cell_APR_2D(part1*P1,part2*P2,int_t*g_str,int_t*g_end,int_t*apr_cell,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //cell number

  if (i>k_num_cells-1) return;
  if (g_str[i]==-1) return;

  int_t idx_insert=0;
  int_t fend=0;
  int_t n=0;
  int_t bflag=0;
  Real msum=0.0;
  int_t pnumber[2]={0,0};
  int_t merge_index=0;

  if(apr_cell[i]<0){          // Cell merging for ESPH
    int_t blood=0;
    idx_insert=k_num_part2+apr_cell[i];

    for(int_t j=g_str[i];j<g_end[i];j++){
      if(P1[j].apr_cond==-1){
        merge_index=j;
        blood=P1[j].blood_index;
        bflag=1;
        break;
      }
    }
    if(blood==0) return;
    if(bflag==0) return;

    int_t icell,jcell;
    Real xi,yi;

    xi=P1[merge_index].x;
    yi=P1[merge_index].y;

    if((k_x_max==k_x_min)){icell=0;}
  	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  	if((k_y_max==k_y_min)){jcell=0;}
  	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}

    // out-of-range handling
  	if(icell<0) icell=0;	if(jcell<0) jcell=0;

    Real m_sum,vol_sum,x_avg,y_avg,ux_avg,uy_avg,pres_avg;
    m_sum=vol_sum=x_avg=y_avg=ux_avg=uy_avg=pres_avg=0.0;

    for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
  		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
  			int_t k=idx_cell(icell+x,jcell+y,0);

  			if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
  			if(g_str[k]!=cu_memset){
  				int_t fend=g_end[k];
  				for(int_t j=g_str[k];j<fend;j++){
            if(P1[j].i_type==1 && P1[j].blood_index==blood){
              m_sum+=P1[j].m;
              vol_sum+=P1[j].vol;
              x_avg+=0.25*P1[j].x;
              y_avg+=0.25*P1[j].y;
              ux_avg+=0.25*P1[j].ux;
              uy_avg+=0.25*P1[j].uy;
              pres_avg+=0.25*P1[j].pres;
              P1[j].i_type=3;
            }
          }
        }
      }
    }
    P1[idx_insert]=P1[merge_index];
    P1[idx_insert].i_type=1;
    P1[idx_insert].m=m_sum;
    P1[idx_insert].vol=vol_sum;
    P1[idx_insert].x=x_avg;
    P1[idx_insert].y=y_avg;
    P1[idx_insert].ux=ux_avg;
    P1[idx_insert].uy=uy_avg;
    P1[idx_insert].pres=pres_avg;
    if(k_h_change) P1[idx_insert].h=P1[idx_insert].h_ref*sqrt(m_sum/P1[idx_insert].m_ref);

    P1[idx_insert].apr_cond=0;

    return;
  }

  else if(apr_cell[i]>0){          // Cell splitting

    for(int_t j=g_str[i];j<g_end[i];j++){
      Real mass_ratio=log(P1[j].m_ref/P1[j].m)/log(4.0);

      if(P1[j].apr_cond==1){

        idx_insert=k_num_part+*num_buffer+4*apr_cell[i];

        Real xm=P1[j].x;
        Real ym=P1[j].y;
        Real Vm=P1[j].m/P1[j].rho;

        int_t k=0;
        while(k<4){
          P1[idx_insert+k]=P1[j];
          if (k==0){
            P1[idx_insert+k].x=xm-0.25*sqrt(Vm);
            P1[idx_insert+k].y=ym+0.25*sqrt(Vm);
            P1[idx_insert+k].pres=P1[j].pres+1.0e-6;
          }
          if (k==1){
            P1[idx_insert+k].x=xm+0.25*sqrt(Vm);
            P1[idx_insert+k].y=ym+0.25*sqrt(Vm);
            P1[idx_insert+k].pres=P1[j].pres-1.0e-6;
          }
          if(k==2){
            P1[idx_insert+k].x=xm-0.25*sqrt(Vm);
            P1[idx_insert+k].y=ym-0.25*sqrt(Vm);
            P1[idx_insert+k].pres=P1[j].pres-2.0e-6;
          }
          if(k==3){
            P1[idx_insert+k].x=xm+0.25*sqrt(Vm);
            P1[idx_insert+k].y=ym-0.25*sqrt(Vm);
            P1[idx_insert+k].pres=P1[j].pres+2.0e-6;
          }
          P1[idx_insert+k].m=0.25*P1[j].m;
          P1[idx_insert+k].vol=0.25*P1[j].vol;
          if(k_h_change) P1[idx_insert+k].h=0.5*P1[j].h;
          P1[idx_insert+k].apr_cond=0;
          k++;
        }

        P1[j].i_type=3;
        return;
      }

      // if(apr_cell[i]>0){          // Cell merging

      // else if(P1[j].apr_cond==2){
      //
      //   idx_insert=k_num_part+*num_buffer+4*apr_cell[i];
      //
      //   Real xm=P1[j].x;
      //   Real ym=P1[j].y;
      //   Real Vm=P1[j].m/P1[j].rho;
      //
      //   k=0;
      //   while(k<2){
      //     P1[idx_insert+k]=P1[j];
      //     if (k==0){
      //       P1[idx_insert+k].x=xm-0.25*sqrt(Vm);
      //       P1[idx_insert+k].y=ym+0.25*sqrt(Vm);
      //     }
      //     if (k==1){
      //       P1[idx_insert+k].x=xm+0.25*sqrt(Vm);
      //       P1[idx_insert+k].y=ym-0.25*sqrt(Vm);
      //     }
      //     P1[idx_insert+k].m=0.5*P1[j].m;
      //     if(k_h_change) P1[idx_insert+k].h=0.5*P1[j].h;
      //     P1[idx_insert+k].apr_cond=0;
      //     k++;
      //   }
      //
      //   P1[j].i_type=3;
      //   return;
      // }
    }
  }
}

__global__ void KERNEL_APR_condition3D_cell(part1*P1,int_t*apr_cell,int_t*g_str,int_t*g_end,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //cell number

  if (i>k_num_cells-1) return;
  if (g_str[i]==-1) return;

  Real x_cell=0.0;
  Real y_cell=0.0;
  Real z_cell=0.0;
  Real vmag_cell=0.0;
  Real apr_level=0.0;
  int_t flag=0;
  int_t merge_flag=0;
  int_t nop=g_end[i]-g_str[i];
  Real M_val=0.0;

  // for(int_t j=g_str[i];j<g_end[i];j++){      // 경계 입자 제
  //   if(P1[j].p_type==1) break;
  //   if(j==g_end[i]-1) return;
  // }


  // for(int_t j=g_str[i];j<g_end[i];j++){
  //   M_val+=float(P1[j].M_num)/nop;              // Cell의 평균 M_val 계산
  // }
  //
  // if(M_val>=0.55){      // splitting
  //   for(int_t j=g_str[i];j<g_end[i];j++){
  //     apr_level=P1[j].m/P1[j].m_ref;
  //     if(P1[j].p_type==1 && apr_level>=0.4) {
  //       apr_cell[i]=1;
  //       return;
  //     }
  //   }
  // }
  // else if(M_val<0.45){
  //   for(int_t j=g_str[i];j<g_end[i];j++){     // merging
  //     apr_level=P1[j].m/P1[j].m_ref;
  //     if(P1[j].p_type==1 && apr_level<=0.75) {
  //       apr_cell[i]=-1;
  //       return;
  //     }
  //   }
  // }

  Real mass_ratio=0.0;
  int_t merge_index=0;
  int_t blood=0;
  for(int_t j=g_str[i];j<g_end[i];j++){
    //if(P1[j].p_type!=1) continue;
    mass_ratio=log(P1[j].m_ref/P1[j].m)/log(8.0);
    if(P1[j].M_num-mass_ratio >= 0.75){
      P1[j].apr_cond=1;
      apr_cell[i]=1;
      return;
    }
    // else if(P1[j].M_num-mass_ratio >= 0.3){
    //   P1[j].apr_cond=2;
    //   apr_cell[i]=1;
    //   return;
    // }
    else if(P1[j].M_num-mass_ratio < -0.3){
      merge_index=j;
      blood=P1[j].blood_index;
      merge_flag=1;
    }
  }

  if (merge_flag==1){                   // Setting blood_index
    Real xi,yi,zi,hi;
    int_t icell,jcell,kcell;
    Real search_range;

    xi=P1[merge_index].x;
    yi=P1[merge_index].y;
    zi=P1[merge_index].z;
    hi=P1[merge_index].h;

    if((k_x_max==k_x_min)){icell=0;}
  	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  	if((k_y_max==k_y_min)){jcell=0;}
  	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  	if((k_z_max==k_z_min)){kcell=0;}
  	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
  	// out-of-range handling
  	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;

  	int_t cell_range=ceil(ncell_init*P1[merge_index].h/P1[merge_index].h_ref);
    search_range=k_search_kappa*hi;

  	for(int_t z=-cell_range;z<=cell_range;z++){
  		for(int_t y=-cell_range;y<=cell_range;y++){
  			for(int_t x=-cell_range;x<=cell_range;x++){
  				int_t k=idx_cell(icell+x,jcell+y,kcell+z);

  				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
  				if(g_str[k]!=cu_memset){
  					int_t fend=g_end[k];
  					for(int_t j=g_str[k];j<fend;j++){

              Real xj,yj,zj,tdist;
  						xj=P1[j].x;
  						yj=P1[j].y;
  						zj=P1[j].z;

  						tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;
              if(tdist<search_range){
                if(P1[j].i_type==1 && P1[j].blood_index==blood){
                  if(P1[j].M_num==1) return;
                  if(j>merge_index) merge_index=j;
                }
              }
            }
          }
        }
      }
    }
    P1[merge_index].apr_cond=-1;         // merge_index is the largest particle index in the blood group
    apr_cell[i]=-1;
    return;
  }

  apr_cell[i]=0;
  return;

  // Spatial condition
  // if(x_cell<0.015 || x_cell>0.085 || y_cell<0.015 || y_cell>0.085) apr_cell[i]=1;     // Splitting condition
  // if(x_cell>0.017 && x_cell<0.083 && y_cell>0.017 && ytim_cell<0.083) apr_cell[i]=-1;     // Merging condition

  // Velocity condition
  // if(vmag_cell > 0.01 && cell_level < 0.9) apr_cell[i]=1;     // Splitting condition
  // if(vmag_cell < 0.01 && cell_level > 0.1) apr_cell[i]=-1;     // Merging condition
  // if(M_val > 0.9 && cell_level < 0.5) {apr_cell[i]=1; printf("Cell num:%d M_val:%f Cell_lv:%f \n",i,M_val,cell_level);}     // Splitting condition
  // if(M_val < 0.1 && cell_level > 0.9) apr_cell[i]=-1;     // Merging condition
}

__global__ void KERNEL_cell_APR_3D(part1*P1,part2*P2,int_t*g_str,int_t*g_end,int_t*apr_cell,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //cell number

  if (i>k_num_cells-1) return;
  if (g_str[i]==-1) return;

  int_t j,k,nop,idx_insert;
  int_t n=0;
  int_t bflag=0;
  int_t pnumber[2]={0,0};

  // if(apr_cell[i]<0){          // Cell merging
  //   for(j=g_str[i];j<g_end[i];j++){
  //     if(P1[j].m<2.0*P1[j].m_ref && P1[j].p_type==1){
  //       pnumber[n]=j;
  //       n++;
  //       if(n==2) break;
  //     }
  //   }
  //   if (n!=2) return;
  //
  //   idx_insert=k_num_part2+apr_cell[i];
  //   j=pnumber[0];
  //   k=pnumber[1];
  //   Real msum, mj, mk;
  //   mj=P1[j].m;
  //   mk=P1[k].m;
  //   msum=mj+mk;
  //
  //   if (msum>P1[j].m_ref*2.0) return;
  //
  //   P1[idx_insert]=P1[j];
  //
  //   P1[idx_insert].m=msum;
  //   P1[idx_insert].x=(mj*P1[j].x+mk*P1[k].x)/msum;
  //   P1[idx_insert].y=(mj*P1[j].y+mk*P1[k].y)/msum;
  //   P1[idx_insert].ux=(mj*P1[j].ux+mk*P1[k].ux)/msum;
  //   P1[idx_insert].uy=(mj*P1[j].uy+mk*P1[k].uy)/msum;
  //   P1[idx_insert].pres=(mj*P1[j].pres+mk*P1[k].pres)/msum;
  //   P1[idx_insert].rho=(mj*P1[j].rho+mk*P1[k].rho)/msum;
  //   if(k_h_change) P1[idx_insert].h=P1[idx_insert].h_ref*cbrt(msum/P1[idx_insert].m_ref);
  //
  //   P1[idx_insert].apr_cond=0;
  //
  //   P1[j].i_type=3;
  //   P1[k].i_type=3;
  //
  //   return;
  // }

  // if(apr_cell[i]<0){          // Cell merging for ESPH
  //   int_t blood=0;
  //   idx_insert=k_num_part2+apr_cell[i];
  //
  //   for(j=g_str[i];j<g_end[i];j++){
  //     if(P1[j].apr_cond==-1){
  //       blood=P1[j].blood_index;
  //       if(blood==0) return;
  //       P1[idx_insert]=P1[j];
  //       break;
  //     }
  //   }
  //
  //   int_t icell,jcell,kcell;
  //   Real xi,yi,zi;
  //
  //   xi=P1[j].x;
  //   yi=P1[j].y;
  //   zi=P1[j].z;
  //
  //   if((k_x_max==k_x_min)){icell=0;}
  // 	else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  // 	if((k_y_max==k_y_min)){jcell=0;}
  // 	else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  // 	if((k_z_max==k_z_min)){kcell=0;}
  // 	else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
  //
  //   // out-of-range handling
  // 	if(icell<0) icell=0;	if(jcell<0) jcell=0;	if(kcell<0) kcell=0;
  //
  //   Real m_sum,vol_sum,x_avg,y_avg,z_avg,ux_avg,uy_avg,uz_avg,pres_avg;
  //   m_sum=vol_sum=x_avg=y_avg=z_avg=ux_avg=uy_avg=uz_avg=pres_avg=0.0;
  //
  // 	int_t cell_range=ceil(ncell_init*P1[i].h/P1[i].h_ref);  //caution
  //
  // 	for(int_t z=-cell_range;z<=cell_range;z++){
  // 		for(int_t y=-cell_range;y<=cell_range;y++){
  // 			for(int_t x=-cell_range;x<=cell_range;x++){
  // 				int_t k=idx_cell(icell+x,jcell+y,kcell+z);
  //
  // 				if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
  // 				if(g_str[k]!=cu_memset){
  // 					int_t fend=g_end[k];
  // 					for(int_t j=g_str[k];j<fend;j++){
  //             if(P1[j].i_type==1 && P1[j].blood_index==blood){
  //               m_sum+=P1[j].m;
  //               vol_sum+=P1[j].vol;
  //               x_avg+=0.125*P1[j].x;
  //               y_avg+=0.125*P1[j].y;
  //               z_avg+=0.125*P1[j].z;
  //               ux_avg+=0.125*P1[j].ux;
  //               uy_avg+=0.125*P1[j].uy;
  //               uz_avg+=0.125*P1[j].uz;
  //               pres_avg+=0.125*P1[j].pres;
  //               P1[j].i_type=3;
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  //   if(m_sum>1e-5) printf("m_sum=%f, blood=%d\n",m_sum,blood);
  //   P1[idx_insert].m=m_sum;
  //   P1[idx_insert].vol=vol_sum;
  //   P1[idx_insert].x=x_avg;
  //   P1[idx_insert].y=y_avg;
  //   P1[idx_insert].z=z_avg;
  //   P1[idx_insert].ux=ux_avg;
  //   P1[idx_insert].uy=uy_avg;
  //   P1[idx_insert].uz=uz_avg;
  //   P1[idx_insert].pres=pres_avg;
  //   P1[idx_insert].rho=m_sum/vol_sum;
  //   if(k_h_change) P1[idx_insert].h=P1[idx_insert].h_ref*cbrt(m_sum/P1[idx_insert].m_ref);
  //
  //   P1[idx_insert].apr_cond=0;
  //
  //   return;
  // }

  if(apr_cell[i]>0){          // Cell splitting

    for(j=g_str[i];j<g_end[i];j++){
      if (P1[j].p_type!=1 && P1[j].M_num>2.0) continue;
      Real mass_ratio=log(P1[j].m_ref/P1[j].m)/log(8.0);

      if(P1[j].apr_cond==1){

        idx_insert=k_num_part+*num_buffer+8*apr_cell[i];

        Real xm=P1[j].x;
        Real ym=P1[j].y;
        Real zm=P1[j].z;
        Real Vm=P1[j].m/P1[j].rho;

        k=0;
        while(k<8){
          P1[idx_insert+k]=P1[j];
          P2[idx_insert+k]=P2[j];

          if(k==0){
            P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym+0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
          }
          else if(k==1){
            P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym+0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm-0.25*cbrt(Vm);
          }
          else if(k==2){
            P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym-0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
          }
          else if(k==3){
            P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym-0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm-0.25*cbrt(Vm);
          }
          else if(k==4){
            P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym+0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
          }
          else if(k==5){
            P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym+0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm-0.25*cbrt(Vm);
          }
          else if(k==6){
            P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym-0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
          }
          else if(k==7){
            P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
            P1[idx_insert+k].y=ym-0.25*cbrt(Vm);
            P1[idx_insert+k].z=zm-0.25*cbrt(Vm);
          }

          P1[idx_insert+k].m=0.125*P1[j].m;
          P1[idx_insert+k].vol=0.125*P1[j].vol;
          P1[idx_insert+k].blood_index=P1[j].blood_index;
          if(k_h_change) P1[idx_insert+k].h=0.5*P1[j].h;
          P1[idx_insert+k].apr_cond=0;

          k++;
        }
        P1[j].i_type=3;
        return;
      }

      // else if(P1[j].p_type==1 && P1[j].m>0.4*P1[j].m_ref){
      //   idx_insert=k_num_part+*num_buffer+8*apr_cell[i];
      //
      //   Real xm=P1[j].x;
      //   Real ym=P1[j].y;
      //   Real zm=P1[j].z;
      //   Real Vm=P1[j].m/P1[j].rho;
      //
      //   k=0;
      //   while(k<4){
      //     P1[idx_insert+k]=P1[j];
      //
      //     if(k==0){
      //       P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
      //       P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
      //     }
      //     if(k==1){
      //       P1[idx_insert+k].x=xm+0.25*cbrt(Vm);
      //       P1[idx_insert+k].z=zm-0.25*cbrt(Vm);
      //     }
      //     if(k==1){
      //       P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
      //       P1[idx_insert+k].z=zm+0.25*cbrt(Vm);
      //     }
      //     else{
      //       P1[idx_insert+k].x=xm-0.25*cbrt(Vm);
      //       P1[idx_insert+k].z=ym-0.25*cbrt(Vm);
      //     }
      //
      //     P1[idx_insert+k].m=0.25*P1[j].m;
      //     if(k_h_change) P1[idx_insert+k].h=sqrt(0.5)*P1[j].h;
      //     else P1[idx_insert+k].h=P1[j].h;
      //
      //     P1[idx_insert+k].apr_cond=0;
      //
      //     k++;
      //   }
      //   P1[j].i_type=3;
      //   return;
      // }
    }
  }
}

// ====================================================================================
// BAPR
// ====================================================================================

__global__ void KERNEL_APR_condition2D_BAPR(part1*P1,int_t*b_idx,int_t*b_str,int_t*b_end,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //block number
  if (i>k_num_blocks-1) return;
  if (b_str[i]==-1) return;

  Real ux,uy,vel_mag_j;
  Real apr_cri=0.0;
  Real M_num_avg=0.0;
  int_t M_num=0;
  int_t flag=0;
  int_t nop=b_end[i]-b_str[i];

  for(int_t j=b_str[i];j<b_end[i];j++){
    Real vel_mag_j=sqrt(P1[j].ux*P1[j].ux+P1[j].uy*P1[j].uy);
    apr_cri+=vel_mag_j/nop;
    M_num_avg+=P1[j].M_num/nop;
    if(vel_mag_j >= 0.2) flag=1;
    // if(vel_mag_j >= 0.25) flag=2;
    // if(vel_mag_j >= 0.5) flag=3;
  }
  M_num=floor(M_num_avg+0.5);

  if (apr_cri > 0.15 || flag>=1){
    for(int_t j=b_str[i];j<b_end[i];j++){
      P1[j].M_num=1.0;
    }
  }
  else if (apr_cri < 0.12 && flag==0){
    for(int_t j=b_str[i];j<b_end[i];j++){
      P1[j].M_num=0.0;
    }
  }
  // else if (M_num==1 && (variable > 0.2 || flag>=2)){
  //   for(int_t z=-1;z<=1;z++){
  // 		for(int_t y=-1;y<=1;y++){
  // 			for(int_t x=-1;x<=1;x++){
  // 				int_t k=idx_block(i,x,y,z);
  //         if (k<=0 || k>=k_NBI*k_NBJ*k_NBK) continue;
  //         for (int_t j=b_str[k];j<b_end[k];j++){
  //           if(P1[j].M_num<1) return;
  //         }
  //       }
  //     }
  //   }
  //   for(int_t j=b_str[i];j<b_end[i];j++){
  //     P1[j].M_num=2.0;
  //   }
  // }
  // else if (M_num==2 && (variable > 0.4 || flag>=3)){
  //   for(int_t z=-1;z<=1;z++){
  // 		for(int_t y=-1;y<=1;y++){
  // 			for(int_t x=-1;x<=1;x++){
  // 				int_t k=idx_block(i,x,y,z);
  //         if (k<=0 || k>=k_NBI*k_NBJ*k_NBK) continue;
  //         for (int_t j=b_str[k];j<b_end[k];j++){
  //           if(P1[j].M_num<2) return;
  //         }
  //       }
  //     }
  //   }
  //   for(int_t j=b_str[i];j<b_end[i];j++){
  //     P1[j].M_num=3.0;
  //   }
  // }
}

__global__ void KERNEL_APR_condition3D_BAPR(part1*P1,int_t*b_idx,int_t*b_str,int_t*b_end,int_t tcount){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;    //block number
  if (i>k_num_blocks-1) return;
  if (b_str[i]==-1) return;

  Real vel_mag_j;
  Real variable=0.0;
  Real M_num_avg=0.0;
  int_t M_num=0;
  int_t flag=0;
  int_t nop=b_end[i]-b_str[i];

  for(int_t j=b_str[i];j<b_end[i];j++){
    Real vel_mag_j=sqrt(P1[j].ux*P1[j].ux+P1[j].uy*P1[j].uy+P1[j].uz*P1[j].uz);
    variable+=abs(vel_mag_j)/nop;
    M_num_avg+=P1[j].M_num/nop;
    if(vel_mag_j >= 0.15) flag=1;
    if(vel_mag_j >= 0.25) flag=2;
    if(vel_mag_j >= 0.5) flag=3;
  }
  M_num=floor(M_num_avg);

  if (M_num==0 && (variable > 0.1 || flag>=1)){
    for(int_t j=b_str[i];j<b_end[i];j++){
      P1[j].M_num=1.0;
    }
  }
  else if (M_num==1 && (variable > 0.2 || flag>=2)){
    for(int_t z=-1;z<=1;z++){
  		for(int_t y=-1;y<=1;y++){
  			for(int_t x=-1;x<=1;x++){
  				int_t k=idx_block(i,x,y,z);
          if (k<=0 || k>=k_NBI*k_NBJ*k_NBK) continue;
          for (int_t j=b_str[k];j<b_end[k];j++){
            if(P1[j].M_num<1) return;
          }
        }
      }
    }
    for(int_t j=b_str[i];j<b_end[i];j++){
      P1[j].M_num=2.0;
    }
  }
  else if (M_num==2 && (variable > 0.4 || flag>=3)){
    for(int_t z=-1;z<=1;z++){
  		for(int_t y=-1;y<=1;y++){
  			for(int_t x=-1;x<=1;x++){
  				int_t k=idx_block(i,x,y,z);
          if (k<=0 || k>=k_NBI*k_NBJ*k_NBK) continue;
          for (int_t j=b_str[k];j<b_end[k];j++){
            if(P1[j].M_num<2) return;
          }
        }
      }
    }
    for(int_t j=b_str[i];j<b_end[i];j++){
      P1[j].M_num=3.0;
    }
  }
  // else{
  //   for(int_t j=b_str[i];j<b_end[i];j++){
  //     P1[j].M_num=0.0;
  //   }
  // }
}


__global__ void KERNEL_variable_smoothing_length2D(int_t*g_str,int_t*g_end,part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if(i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  if(P1[i].p_type<1) return;

  int_t icell,jcell;
  Real xi,yi,hi,h_ref_i;
  Real search_range,tmp_A,tmp_R,tmp_flt,tmp_RR;

  tmp_R=tmp_flt=0.0;

  hi=P1[i].h;
  h_ref_i=P1[i].h_ref;
  tmp_A=calc_tmpA(hi);
  //search_range=k_search_kappa*hi;	// search range
  search_range=k_kappa*hi;	// search range

  xi=P1[i].x;
  yi=P1[i].y;

  // calculate I,J,K in cell
  if((k_x_max==k_x_min)){icell=0;}
  else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  if((k_y_max==k_y_min)){jcell=0;}
  else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  // out-of-range handling
  if(icell<0) icell=0;	if(jcell<0) jcell=0;

  tmp_RR=0.0;
  for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
      // int_t k=(icell+x)+k_NI*(jcell+y);
      int_t k=idx_cell(icell+x,jcell+y,0);
      if (k<0) continue;
      if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
      if(g_str[k]!=cu_memset){
        int_t fend=g_end[k];
        for(int_t j=g_str[k];j<fend;j++){
          Real xj,yj,mj,hj,rhoj,tdist,ptypej;

          xj=P1[j].x;
          yj=P1[j].y;
          mj=P1[j].m;
          hj=P1[j].h;
          rhoj=P1[j].rho;
          ptypej=P1[j].p_type;
          tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;

          //search_range=k_search_kappa*fmax(hi,hj);

          if(tdist<search_range){
            Real twij;

            twij=calc_kernel_wij(tmp_A,hi,tdist);
            tmp_R+=twij;
            tmp_RR+=hj*mj/rhoj*twij;
            tmp_flt+=mj/rhoj*twij;
          }
        }
      }
    }
  }
  //if (tcount%k_freq_h_reset==0)P1[i].h=1.5*(sqrt(1/tmp_R));
  //else P1[i].h=tmp_RR/tmp_flt;
  P1[i].h=tmp_RR/tmp_flt;
}
