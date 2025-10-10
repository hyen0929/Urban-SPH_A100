// ====================================================================================
// etc
// ====================================================================================

void Cell_index(int_t*APR_cell){
  int_t j=1;
  int_t k=-1;

  for (int_t i=0; i<=num_cells; i++){
    if(APR_cell[i]>0){
      APR_cell[i]=j;
      j++;
    }
    if(APR_cell[i]<0){
      APR_cell[i]=k;
      k--;
    }
  }
  return;
}

__global__ void KERNEL_set_p_type(part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>=k_num_part2) return;
  Real xi,yi;

  xi=P1[i].x;
  yi=P1[i].y;

  if((yi-0.5)*(yi-0.5)+xi*xi<0.26*0.26){

    if((yi-0.5)*(yi-0.5)+xi*xi<0.25*0.25) {        // bubble
      P1[i].p_type=2;
      P1[i].rho=100.0;
      P1[i].m=0.01;
      P1[i].m_ref=0.04;
    }
    else if ((P1[i].p_type!=0)&&((yi-0.5)*(yi-0.5)+xi*xi>=0.25*0.25)) {       // water
      P1[i].p_type=1;
      P1[i].rho=1000.0;
      P1[i].m=0.1;
      P1[i].m_ref=0.4;
    }
  }
}

__global__ void Kernel_change_ptype(part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if (i>k_num_part2-1) return;

  Real xi,yi,xcom,ycom,r;

  xi=P1[i].x;
  yi=P1[i].y;

  xcom=0.05;
  ycom=0.05;
  r=0.0049;

  if (((ycom-yi)*(ycom-yi)+(xcom-xi)*(xcom-xi))<r*r) P1[i].p_type=-2;
  else if (((ycom-yi)*(ycom-yi)+(xcom-xi)*(xcom-xi))<=1.1*r*r) P1[i].p_type=1;

  return;
}

__global__ void KERNEL_remove_boundary(part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>=k_num_part2) return;
  if(P1[i].p_type==1) return;
  Real xi,yi,zi;

  xi=P1[i].x;
  yi=P1[i].y;
  zi=P1[i].z;

  //if (xi < -0.025 || yi < -0.025 || xi > 3.225 ) P1[i].i_type=3;   // Dambreak
  if (xi < -0.051 || yi < -0.051 || xi > 1.051 || yi > 1.051) P1[i].i_type=3;     // Lid-driven cavity
  if (zi < -0.051 || zi > 1.051) P1[i].i_type=3;     // For 3D
  return;
}

__global__ void KERNEL_remember_blood(part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>=k_num_part2) return;
  P1[i].blood_index=i+1;
}

__global__ void KERNEL_reset_APR_variables(part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if(i>=k_num_part2) return;
  if(P1[i].i_type==3) return;

  P1[i].apr_cond=0;
  P1[i].merge_flag=0;
  return;
}

__global__ void KERNEL_find_miny(part1*P1,Real*min_y){

  uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>=k_num_part2) return;
  if((P1[i].p_type==2)||(P1[i].p_type==9)) min_y[i]=P1[i].y;
  else min_y[i]=10.0;

  return;
}

__global__ void KERNEL_store_min(part1*P1,Real*x_min_storage,Real*y_min_storage){

  uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>k_num_part2) return;

  if((P1[i].p_type==2)&&(P1[i].i_type==1)){
    x_min_storage[i]=P1[i].x;
    y_min_storage[i]=P1[i].y;
  }
  else {
    x_min_storage[i]=100.0;
    y_min_storage[i]=100.0;
  }
  return;
}

__global__ void KERNEL_store_max(part1*P1,Real*x_max_storage,Real*y_max_storage){

  uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
  if(i>k_num_part2) return;
  if(P1[i].i_type==3) return;

  if(P1[i].p_type==2){
    x_max_storage[i]=P1[i].x;
    y_max_storage[i]=P1[i].y;
  }
  else {
    x_max_storage[i]=-100.0;
    y_max_storage[i]=-100.0;
  }
  return;
}

__global__ void KERNEL_Smooth_M(int_t*g_str,int_t*g_end,part1*P1,part1*SP1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if(i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  int_t icell,jcell;
  Real xi,yi,hi;
  Real search_range,tmp_A,tmp_R,tmp_flt,tmp_RR;

  tmp_R=tmp_flt=0.0;

  hi=P1[i].h;
  tmp_A=calc_tmpA(hi);
  search_range=k_search_kappa*hi;	// search range

  xi=P1[i].x;
  yi=P1[i].y;

  // calculate I,J,K in cell
  if((k_x_max==k_x_min)){icell=0;}
  else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  if((k_y_max==k_y_min)){jcell=0;}
  else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  // out-of-range handling
  if(icell<0) icell=0;	if(jcell<0) jcell=0;

  for(int_t y=-P1[i].ncell;y<=P1[i].ncell;y++){
		for(int_t x=-P1[i].ncell;x<=P1[i].ncell;x++){
      // int_t k=(icell+x)+k_NI*(jcell+y);
      int_t k=idx_cell(icell+x,jcell+y,0);
      if (k<0) continue;
      if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
      if(g_str[k]!=cu_memset){
        int_t fend=g_end[k];
        for(int_t j=g_str[k];j<fend;j++){
          Real xj,yj,mj,hj,Mj,rhoj,tdist;

          xj=P1[j].x;
          yj=P1[j].y;
          mj=P1[j].m;
          hj=P1[j].h;
          Mj=SP1[j].M_num;
          rhoj=P1[j].rho;
          tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));

          if(k_apr_solv) search_range=k_search_kappa*fmax(hi,hj);

          if(tdist<search_range){
            Real twij;

            twij=calc_kernel_wij(tmp_A,P1[j].h_ref,tdist);
            tmp_R+=Mj*mj/rhoj*twij;
            tmp_flt+=mj/rhoj*twij;
          }
        }
      }
    }
  }
  P1[i].M_num=tmp_R/tmp_flt;
}

// ====================================================================================
// Smoothing length treatment
// ====================================================================================

__global__ void Kernel_variable_smoothing_length2D(int_t*g_str,int_t*g_end,part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if(i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  if(P1[i].p_type<1) return;

  int_t icell,jcell;
  Real xi,yi,hi;
  Real search_range,tmp_A,tmp_R,tmp_flt,tmp_RR;

  tmp_R=tmp_flt=0.0;

  hi=5.0*P1[i].h;
  tmp_A=calc_tmpA(hi);
  //search_range=k_search_kappa*hi;	// search range
  search_range=k_search_kappa*hi;	// search range

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
  for(int_t y=-5*P1[i].ncell;y<=5*P1[i].ncell;y++){
		for(int_t x=-5*P1[i].ncell;x<=5*P1[i].ncell;x++){
      // int_t k=(icell+x)+k_NI*(jcell+y);
      int_t k=idx_cell(icell+x,jcell+y,0);
      if (k<0) continue;
      if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))) continue;
      if(g_str[k]!=cu_memset){
        int_t fend=g_end[k];
        for(int_t j=g_str[k];j<fend;j++){
          Real xj,yj,mj,hj,rhoj,tdist;

          xj=P1[j].x;
          yj=P1[j].y;
          mj=P1[j].m;
          hj=P1[j].h;
          rhoj=P1[j].rho;
          tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj))+1e-20;

          search_range=k_search_kappa*fmax(hi,hj);

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

__global__ void KERNEL_variable_smoothing_length3D(int_t*g_str,int_t*g_end,part1*P1){
  int_t i=threadIdx.x+blockIdx.x*blockDim.x;

  if (i>k_num_part2-1) return;
  if(P1[i].i_type>i_type_crt) return;
  if(P1[i].p_type<1) return;

  int_t icell,jcell,kcell;
  Real xi,yi,zi,hi;
  Real search_range,tmp_A,tmp_R,tmp_RR,tmp_flt;

  tmp_R=tmp_RR=tmp_flt=0.0;
  Real smooth_range=3.0;
  hi=smooth_range*P1[i].h;
  tmp_A=calc_tmpA(hi);
  search_range=k_search_kappa*hi;	// search range

  xi=P1[i].x;
  yi=P1[i].y;
  zi=P1[i].z;

  // calculate I,J,K in cell
  if((k_x_max==k_x_min)){icell=0;}
  else{icell=min(floor((xi-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
  if((k_y_max==k_y_min)){jcell=0;}
  else{jcell=min(floor((yi-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
  if((k_z_max==k_z_min)){kcell=0;}
  else{kcell=min(floor((zi-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
  // out-of-range handling
  if(icell<0) icell=0;	if(jcell<0) jcell=0;  if(kcell<0) kcell=0;

  int_t ncell=smooth_range*P1[i].ncell;

  for(int_t z=-ncell;z<=ncell;z++){
    for(int_t y=-ncell;y<=ncell;y++){
      for(int_t x=-ncell;x<=ncell;x++){
        // int_t k=(icell+x)+k_NI*(jcell+y);
        int_t k=idx_cell(icell+x,jcell+y,kcell+z);
        if (k<0) continue;
        if(((icell+x)<0)||((icell+x)>(k_NI-1))||((jcell+y)<0)||((jcell+y)>(k_NJ-1))||((kcell+z)<0)||((kcell+z)>(k_NK-1))) continue;
        if(g_str[k]!=cu_memset){
          int_t fend=g_end[k];
          for(int_t j=g_str[k];j<fend;j++){
            Real xj,yj,zj,hj,mj,rhoj,tdist;

            xj=P1[j].x;
            yj=P1[j].y;
            zj=P1[j].z;
            mj=P1[j].m;
            hj=P1[j].h;
            rhoj=P1[j].rho;

            tdist=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))+1e-20;

            if(k_apr_solv) search_range=k_search_kappa*fmax(hi,hj);

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
  }
  //if(tcount%k_freq_h_reset==0)P1[i].h=1.5*(cbrt(1/tmp_R));
  // if(tcount%k_freq_h_reset==0)P1[i].h=P1[i].h_ref*cbrt(P1[i].m/P1[i].m_ref);
  // else P1[i].h=tmp_RR/tmp_flt;
  P1[i].h=tmp_RR/tmp_flt;
}

// ====================================================================================
// CPU part
// ====================================================================================

void Particle_refinement_cell(part1*dev_P1,part2*dev_P2,int_t*APR_cell,int_t*dev_APR_cell,int_t*g_str,int_t*g_end,int_t tcount){
  dim3 b,t;

  b.x=(num_cells-1)/t.x+1;
  if(dim==2) KERNEL_APR_condition2D_cell<<<b,t>>>(dev_P1,dev_APR_cell,g_str,g_end,count);
  if(dim==3) KERNEL_APR_condition3D_cell<<<b,t>>>(dev_P1,dev_APR_cell,g_str,g_end,count);
  cudaDeviceSynchronize();

  cudaMemcpy(APR_cell,dev_APR_cell,num_cells*sizeof(int_t),cudaMemcpyDeviceToHost);

  Cell_index(APR_cell);

  cudaMemcpy(dev_APR_cell,APR_cell,num_cells*sizeof(int_t),cudaMemcpyHostToDevice);

  //clock_t time_st=clock();

  t.x=128;
  b.x=(num_cells-1)/t.x+1;
  if(dim==2) KERNEL_cell_APR_2D<<<b,t>>>(dev_P1,dev_P2,g_str,g_end,dev_APR_cell,tcount);
  if(dim==3) KERNEL_cell_APR_3D<<<b,t>>>(dev_P1,dev_P2,g_str,g_end,dev_APR_cell,tcount);
  cudaDeviceSynchronize();

  // clock_t time_ed=clock();
	// double calc_time_step=(float)(time_ed-time_st)/CLOCKS_PER_SEC;
  // printf("Time: %f\n",calc_time_step);

}

void Set_block_condition(part1*dev_P1,part2*dev_P2,part1*dev_SP1,part2*dev_SP2,int_t*b_idx,int_t*b_idx_in,
  int_t*p_idx,int_t*p_idx_in,void*dev_sort_storage,size_t*sort_storage_bytes,int_t*b_str,int_t*b_end){
  dim3 b,t;
  t.x=128;
	int s=sizeof(int)*(t.x+1);

  // b_str을 리셋
  cudaMemset(b_str,cu_memset,sizeof(int_t)*num_blocks);

  // 입자의 블럭 번호 계산
  b.x=(num_part2-1)/t.x+1;
  KERNEL_index_particle_to_block<<<b,t>>>(b_idx_in,p_idx_in,dev_P1);
  cudaDeviceSynchronize();

  // 블럭 번호를 바탕으로 정렬
  cub::DeviceRadixSort::SortPairs(dev_sort_storage,*sort_storage_bytes,b_idx_in,b_idx,p_idx_in,p_idx,num_part2);
  cudaDeviceSynchronize();

  // 정렬한 입자를 재배치
  b.x=(num_part2-1)/t.x+1;
  KERNEL_reorder_single<<<b,t,s>>>(b_idx,p_idx,b_str,b_end,dev_P1,dev_P2,dev_SP1,dev_SP2);
  cudaDeviceSynchronize();

  cudaMemcpy(dev_P1,dev_SP1,sizeof(part1)*num_part2,cudaMemcpyDeviceToDevice);
  cudaMemcpy(dev_P2,dev_SP2,sizeof(part2)*num_part2,cudaMemcpyDeviceToDevice);
  pthread_barrier_wait(&barrier);

  // BAPR condition
  b.x=(num_part2-1)/t.x+1;
  if(dim==2)KERNEL_APR_condition2D_BAPR<<<b,t>>>(dev_SP1,b_idx,b_str,b_end,count);
  if(dim==3)KERNEL_APR_condition3D_BAPR<<<b,t>>>(dev_SP1,b_idx,b_str,b_end,count);
  cudaDeviceSynchronize();
}

void Particle_refinement_block(part1*dev_P1,part1*dev_SP1,part2*dev_P2,int_t*APR_cell,int_t*dev_APR_cell,int_t*g_str,int_t*g_end,int_t tcount){
  dim3 b,t;

  // t.x=128;
  // b.x=(num_part2-1)/t.x+1;
  // if(dim==2) KERNEL_Smooth_M<<<b,t>>>(g_str,g_end,dev_P1,dev_SP1);
  // cudaDeviceSynchronize();

  b.x=(num_cells-1)/t.x+1;
  if(dim==2) KERNEL_APR_condition2D_cell<<<b,t>>>(dev_P1,dev_APR_cell,g_str,g_end,count);
  if(dim==3) KERNEL_APR_condition3D_cell<<<b,t>>>(dev_P1,dev_APR_cell,g_str,g_end,count);
  cudaDeviceSynchronize();

  cudaMemcpy(APR_cell,dev_APR_cell,num_cells*sizeof(int_t),cudaMemcpyDeviceToHost);

  Cell_index(APR_cell);

  cudaMemcpy(dev_APR_cell,APR_cell,num_cells*sizeof(int_t),cudaMemcpyHostToDevice);

  t.x=128;
  b.x=(num_cells-1)/t.x+1;
  if(dim==2) KERNEL_cell_APR_2D<<<b,t>>>(dev_P1,dev_P2,g_str,g_end,dev_APR_cell,tcount);
  if(dim==3) KERNEL_cell_APR_3D<<<b,t>>>(dev_P1,dev_P2,g_str,g_end,dev_APR_cell,tcount);
  cudaDeviceSynchronize();

}

// ====================================================================================
// Particle counting
// ====================================================================================

void count_APR_buffer(part1*P1,int_t*count_buffer){
  int_t i,nop;
	nop=num_part2;
	int_t Nparticle=0;
	for(i=0;i<nop;i++) if(P1[i].i_type<3) Nparticle++;
  count_buffer[3]=Nparticle-num_part;

  return;
}

void count_buffer_particles(int_t*count_buffer, int_t*dev_count_buffer, int_t*num_buffer_temp, int_t tcount){
  int_t h_num_buffer[1], h_num_buffer_temp[1];

  if(count_buffer[3]!=0) *h_num_buffer=count_buffer[3];
  else{
    cudaMemcpy(count_buffer,dev_count_buffer,4*sizeof(int_t),cudaMemcpyDeviceToHost);
    cudaMemcpy(h_num_buffer_temp,num_buffer_temp,sizeof(int_t),cudaMemcpyDeviceToHost);

    *h_num_buffer=*h_num_buffer_temp;
    *h_num_buffer+=count_buffer[0];
    *h_num_buffer-=count_buffer[1];
    *h_num_buffer-=count_buffer[2];

  }
  cudaMemcpyToSymbol(num_buffer,h_num_buffer,sizeof(int_t),0,cudaMemcpyHostToDevice);
  if((tcount%freq_output)==0) printf("Num_part+buffer : %d\n\n",num_part+*h_num_buffer);
}
