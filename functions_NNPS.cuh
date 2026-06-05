void c_initial_inner_outer_particle_single(part1*HP1,part1*DHP1,int_t tid){

	int_t i,c_count;
	Real xi0;
	Real maxb,minb;
	c_count=0;

	for(i=0;i<num_part;i++){

		Real xi, yi, zi, scale_factor;

		xi=HP1[i].x;
		yi=HP1[i].y;
		zi=HP1[i].z;

		HP1[i].i_type=1;
		HP1[i].buffer_type=0;

		HP1[i].m_ref=HP1[i].m;
		HP1[i].ux=0.0;
		HP1[i].uy=0.0;
		HP1[i].uz=0.0;
		HP1[i].temp=300.0;
		HP1[i].concn=0.0;
		HP1[i].pres=0.0;
		HP1[i].pres0=0.0;   // incremental PPE anchor

		// Time-averaging accumulators zero-init
		HP1[i].ux_avg=0.0;
		HP1[i].uy_avg=0.0;
		HP1[i].uz_avg=0.0;
		HP1[i].k_turb_avg=0.0;
		HP1[i].e_turb_avg=0.0;

		// Initial condition selector
		if(CASE_ID>100 && CASE_ID<200) {
			scale_factor=1.0;
			HP1[i].ux=Interpolate_ux_caseA(HP1[i].z,scale_factor);
			HP1[i].k_turb=Interpolate_k_caseA(HP1[i].z,scale_factor);
		}
		if(CASE_ID>200 && CASE_ID<210) {
			scale_factor=0.08;
			HP1[i].ux=Interpolate_ux_caseA(HP1[i].z,scale_factor);
			//HP1[i].k_turb=Interpolate_k_caseC(HP1[i].z,scale_factor);
		}
		if(CASE_ID>210 && CASE_ID<220) {	// Stability comparison (unstable/neutral/stable)
			HP1[i].ux  = Interpolate_ux_MOST_Neutral(HP1[i].z);
			HP1[i].k_turb = Interpolate_k_MOST_Neutral(HP1[i].z);
			HP1[i].e_turb = Interpolate_e_MOST_Neutral(HP1[i].z);  // MOST eps profile for k-eps
			HP1[i].temp   = Interpolate_T_Linear_Stable(HP1[i].z);  // linear-init diagnostic
		}
		if(CASE_ID>500 && CASE_ID<600) {	// Case CS (single building block) — AS-style log-law init
			HP1[i].ux     = Interpolate_ux_MOST_Neutral(HP1[i].z);
			HP1[i].k_turb = Interpolate_k_MOST_Neutral(HP1[i].z);
			HP1[i].e_turb = Interpolate_e_MOST_Neutral(HP1[i].z);
		}
		if(CASE_ID>300 && CASE_ID<400) {
			scale_factor=1.0;
			HP1[i].ux=Interpolate_ux_caseC(HP1[i].z,scale_factor);
			//HP1[i].k_turb=Interpolate_k_caseC(HP1[i].z,scale_factor);
		}
		if(CASE_ID>1000 && CASE_ID<1100) {
			HP1[i].ux=Interpolate_ux_ABL(HP1[i].z);
			//HP1[i].ux=1.5*Interpolate_ux_caseA(HP1[i].z,0.08);
			HP1[i].k_turb=0.0;
		}
		if(CASE_ID==2001) HP1[i].ux=0.0;	// Scalar transport test
		if(CASE_ID==2002) HP1[i].ux=0.3;	// Scalar transport test with advection
		if(CASE_ID==3001) {	// MOST validation test
			HP1[i].ux=Interpolate_ux_MOST_Neutral(HP1[i].z);
			HP1[i].k_turb=Interpolate_k_MOST_Neutral(HP1[i].z);
			HP1[i].e_turb=Interpolate_e_MOST_Neutral(HP1[i].z);
		}
		if(CASE_ID==1500) {	// RANS/LES channel validation: band-limited random seed → bypass transition
			// uy/uz = sum of 6 Fourier modes (random phases, wavelengths >= ~20dx).
			// Smooth (no grid-scale 2dx content) yet broadband → triggers transition
			// without exciting the numerical checkerboard mode.
			Real xx=HP1[i].x, yy=HP1[i].y, zz=HP1[i].z;
			Real Lx=200.0, Ly=80.0, Lz=80.0;
			Real twoPI=2.0*PI;
			Real py = sin(twoPI*(2.0*xx/Lx +            1.0*zz/Lz)+0.7)
			        + sin(twoPI*(1.0*xx/Lx +1.0*yy/Ly + 2.0*zz/Lz)+2.1)
			        + sin(twoPI*(3.0*xx/Lx +2.0*yy/Ly           )+4.0)
			        + sin(twoPI*(1.0*xx/Lx +            3.0*zz/Lz)+5.3)
			        + sin(twoPI*(2.0*xx/Lx +3.0*yy/Ly + 1.0*zz/Lz)+1.1)
			        + sin(twoPI*(3.0*xx/Lx +1.0*yy/Ly + 2.0*zz/Lz)+3.8);
			Real pz = sin(twoPI*(2.0*xx/Lx +            2.0*zz/Lz)+1.4)
			        + sin(twoPI*(1.0*xx/Lx +1.0*yy/Ly + 1.0*zz/Lz)+3.2)
			        + sin(twoPI*(3.0*xx/Lx +2.0*yy/Ly + 1.0*zz/Lz)+0.3)
			        + sin(twoPI*(2.0*xx/Lx +4.0*yy/Ly           )+5.9)
			        + sin(twoPI*(1.0*xx/Lx +2.0*yy/Ly + 3.0*zz/Lz)+2.6)
			        + sin(twoPI*(3.0*xx/Lx +            2.0*zz/Lz)+4.7);
			HP1[i].ux=CHANNEL_U0;
			HP1[i].uy=(1.0/6.0)*py;   // peak ~±1 m/s, RMS ~0.3
			HP1[i].uz=(1.0/6.0)*pz;
			HP1[i].k_turb=0.0;
			HP1[i].e_turb=0.0;
		}

		// IBM building surface: initialize to wall temperature (T_con)
		if(HP1[i].p_type >= 1000) {
			HP1[i].temp = T_con;
			HP1[i].ux = 0.0;
			HP1[i].uy = 0.0;
			HP1[i].uz = 0.0;
		}


		// // Inner building particle ux=0.0
		// if(HP1[i].p_type==1 && zi<0.2){
		// 	if(xi>1.0 && xi<1.2 && yi>-0.5 && yi<-0.3) HP1[i].ux=0.0;
		// 	else if(xi>1.4 && xi<1.6 && yi>-0.5 && yi<-0.3) HP1[i].ux=0.0;
		// 	else if(xi>1.8 && xi<2.0 && yi>-0.5 && yi<-0.3) HP1[i].ux=0.0;
		// 	else if(xi>1.0 && xi<1.2 && yi>-0.1 && yi<0.1) HP1[i].ux=0.0;
		// 	else if(xi>1.4 && xi<1.6 && yi>-0.1 && yi<0.1) HP1[i].ux=0.0;
		// 	else if(xi>1.8 && xi<2.0 && yi>-0.1 && yi<0.1) HP1[i].ux=0.0;
		// 	else if(xi>1.0 && xi<1.2 && yi>0.3 && yi<0.5) HP1[i].ux=0.0;
		// 	else if(xi>1.4 && xi<1.6 && yi>0.3 && yi<0.5) HP1[i].ux=0.0;
		// 	else if(xi>1.8 && xi<2.0 && yi>0.3 && yi<0.5) HP1[i].ux=0.0;
		// }
		// else if(HP1[i].p_type==1 && zi<0.4){
		// 	if(xi>1.4 && xi<1.6 && yi>-0.1 && yi<0.1) HP1[i].ux=0.0;
		// }

		if(HP1[i].p_type==3) {
			//HP1[i].i_type=2;
			HP1[i].p_type=1;
		}

		// Interior fluid flag: p_type=99 from Python preprocessing → buffer_type=99
		if(HP1[i].p_type==99) {
			HP1[i].buffer_type=99;
			HP1[i].p_type=1;
			// HP1[i].ux=0.0;
			// HP1[i].uy=0.0;
			// HP1[i].uz=0.0;
		}

		if(HP1[i].p_type==0||HP1[i].p_type==1000) {
			//HP1[i].i_type=2;
			HP1[i].rho=1.225;
			HP1[i].rho=1.225;
		}


		if (open_boundary>0)			// (CAUTION)
		{
			// if (HP1[i].p_type>=0) HP1[i].ux=Inlet_Velocity;
			if(HP1[i].p_type==1){

				if(HP1[i].y<L_left) {
					HP1[i].buffer_type=3;
					HP1[i].i_type=2;
					// HP1[i].p_type=-1;
				}
				else if(HP1[i].y>L_right) {
					HP1[i].buffer_type=4;
					HP1[i].i_type=2;
					// HP1[i].p_type=-1;
				}
				else if((HP1[i].x>=L3-1e-6)&&(HP1[i].x<L4-1e-6)) {
					HP1[i].buffer_type=2;
					HP1[i].i_type=2;
				}
				else if(HP1[i].x<L2-1e-6) {
					HP1[i].buffer_type=1;
					HP1[i].i_type=2;
					//if(HP1[i].x<10.0) HP1[i].buffer_type=-1;

				}

				if((HP1[i].z>=L_ceiling)&&(HP1[i].y<L_left)&&(HP1[i].y>L_right)){
					HP1[i].buffer_type=5;
				}

				if(HP1[i].p_type==1 && HP1[i].x>L_MOST_x0 && HP1[i].x<L_MOST_xmax && HP1[i].z>L_MOST_z_lo && HP1[i].z<L_MOST_z_hi){
					HP1[i].MOST_buffer=1;
				}

				// else if((HP1[i].x>=L2-1e-6)&&(HP1[i].x<L3)) {   // Only buffer 1 for ESPH
				// 	HP1[i].buffer_type=0;
				// 	HP1[i].i_type=1;
				// }
				// else if ((HP1[i].x>(L3+1e-6))&&(HP1[i].x<=(1.001*L4))) {
				// 	HP1[i].buffer_type=2;
				// 	HP1[i].i_type=2;
				//
				// }
			}

			if(HP1[i].p_type<=0){
				// if(HP1[i].y<1e-6+L_left) {
				// 	HP1[i].buffer_type=3;
				// 	HP1[i].i_type=2;
				// }
				// else if(HP1[i].y>L_right-1e-6) {
				// 	HP1[i].buffer_type=4;
				// 	HP1[i].i_type=2;
				// }
			}
		}
		DHP1[c_count]=HP1[i];
		c_count++;
	}
}

void c_initial_inner_outer_particle_single_LDM(L_part1*HLP1,L_part1*DHLP1,int_t tid){

	int_t i,c_count;
	Real xi0;
	Real maxb,minb;

	c_count=0;

	for(i=0;i<num_part_LDM;i++){
		HLP1[i].i_type=3;
		// HP1[i].h=1.5;
		// HP1[i].h_ref=HP1[i].h;
		HLP1[i].m_ref=HLP1[i].m;
		HLP1[i].ux=0.0;
		HLP1[i].uy=0.0;
		HLP1[i].uz=0.0;

		HLP1[i].ncell=ncell_init;    	
		DHLP1[c_count]=HLP1[i];
		c_count++;
	}
}

//////////////////////////////////////////////////////////////////////////////
void c_initial_inner_outer_particle(part1*HP1,part1*DHP1,int_t tid){

	int_t i,c_count;
	Real xi0;
	Real maxb,minb;

	// 각 GPU의 x 축 왼쪽 좌표(minb)와 오른쪽 좌표(maxb) 계산
	minb=x_min+(Real)(calc_area*tid)*dcell;
	maxb=x_min+(Real)(calc_area*(tid+1))*dcell;

	// 입자의 개수(혹은 위치) 변수
	int_t cpucount;
	cpucount=0;

	// Host 입자정보(HP1)를 GPU별로 분할(DHP1)
	if(tid==0){		// 첫번째 GPU 라면
		c_count=0;
		for(i=0;i<num_part;i++){
			xi0=HP1[i].x;

			if(xi0<maxb&&xi0>=x_min) cpucount++;	// 입자의 위치가 영역안에 들어오는지 판별

			if(xi0>=maxb+dcell*2.0||xi0<x_min) continue;	// 입자의 위치가 계산영역(버퍼 포함)을 벗어나는지 판별
			if(xi0>=maxb&&xi0<maxb+dcell*2.0){	// 입자의 위치가 outer buffer 영역에 포함되는지 판별
				HP1[i].i_type=2;		// 버퍼영역이면 i_type=2
				DHP1[c_count]=HP1[i];	// 정보 복사
				c_count++;
				continue;
			}
			if((xi0<(maxb-dcell*2.0))&&(xi0>=(x_min+dcell*2.0))){		// 입자의 위치가 inner 영역에 있는지 판별
				HP1[i].i_type=0;
				DHP1[c_count]=HP1[i];
				c_count++;
			}else{		// 그 외에는 inner buffer 영역으로 간주
				HP1[i].i_type=1;
				DHP1[c_count]=HP1[i];
				c_count++;
			}
		}
	}else if(tid==ngpu-1){		// 마지막 GPU 라면
		c_count=0;
		for(i=0;i<num_part;i++){
			xi0=HP1[i].x;

			if(xi0<=x_max&&xi0>=minb) cpucount++;

			if(xi0>x_max||xi0<minb-dcell*2.0) continue;
			if(xi0<minb&&xi0>=minb-dcell*2.0){
				HP1[i].i_type=2;
				DHP1[c_count]=HP1[i];
				c_count++;
				continue;
			}
			if((xi0<(x_max-dcell*2.0))&&(xi0>=(minb+dcell*2.0))){
				HP1[i].i_type=0;
				DHP1[c_count]=HP1[i];
				c_count++;
			}else{
				HP1[i].i_type=1;
				DHP1[c_count]=HP1[i];
				c_count++;
			}
		}
	}else{	    // 그 밖의 GPU에 대해서
		c_count=0;
		for(i=0;i<num_part;i++){
			xi0=HP1[i].x;

			if(xi0<maxb&&xi0>=minb) cpucount++;

			if(xi0>=maxb+dcell*2.0||xi0<minb-dcell*2.0) continue;
			if(xi0>=maxb&&xi0<maxb+dcell*2.0){
				HP1[i].i_type=2;
				DHP1[c_count]=HP1[i];
				c_count++;
				continue;
			}
			if(xi0<minb&&xi0>=minb-dcell*2.0){
				HP1[i].i_type=2;
				DHP1[c_count]=HP1[i];
				c_count++;
				continue;
			}
			if((xi0<(maxb-dcell*2.0))&&(xi0>=(minb+dcell*2.0))){
				HP1[i].i_type=0;
				DHP1[c_count]=HP1[i];
				c_count++;
			}else{
				HP1[i].i_type=1;
				DHP1[c_count]=HP1[i];
				c_count++;
			}
		}
	}

	// if inner else outer
	//if((xi0<=(maxb-dcell*1.0))&&(xi0>=(minb+dcell*1.0))) HP1[i].i_type=1;
	//else HP1[i].i_type=2;
	printf("%d Particle Number : %d\n",tid,cpucount);
	printf("%d Ref+Real Number : %d\n",tid,c_count);

}
//////////////////////////////////////////////////////////////////////////////
// 필요한 cell의 개수
int_t clc_num_cells() {

	int_t result;
	int_t NI_max;

	NI_max=max(max(NI,NJ),NK);

	//printf("NI_max = %d\n", NI_max);

	result=NI*NJ*NK+1;

	return result;
}
//////////////////////////////////////////////////////////////////////////////
// cell의 index
__host__ __device__ int_t idx_block(int_t idx, int_t Ix, int_t Iy, int_t Iz) {

	int_t result;

	result=idx+(Ix)+k_NBI*(Iy)+k_NBI*k_NBJ*(Iz);
	return result;
}

__host__ __device__ int_t idx_cell(int_t Ix, int_t Iy, int_t Iz) {

	int_t result;

	result=(Ix)+k_NI*(Iy)+k_NI*k_NJ*(Iz);

	return result;
}
//////////////////////////////////////////////////////////////////////////////
__global__ void right_send_particle(int_t*p2p_af,int_t*p2p_idx,part1*P1,part3*P3,part1*SP1,p2p_part3*SP3,int_t tid){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type!=1){		// 입자가 처음부터 영역 1 에 있지 않았으면 (영역 0 또는 영역 2, 3)
		p2p_af[i]=2;					// 표식의 우선순위를 낮추고
		p2p_idx[i]=i;
		SP1[i].i_type=3;			// 더미 입자로 설정
		return;
	}

	Real maxb,xi0;
	xi0=P1[i].x;

	maxb=k_x_min+k_calc_area*k_dcell*(tid+1);

	if(xi0>=maxb-k_dcell*2.0){
		p2p_af[i]=1;
		p2p_idx[i]=i;

		SP1[i]=P1[i];
		SP3[i].drho=P3[i].drho;
		SP3[i].dconcn=P3[i].dconcn;
		SP3[i].denthalpy=P3[i].denthalpy;
		SP3[i].ftotalx=P3[i].ftotalx;
		SP3[i].ftotaly=P3[i].ftotaly;
		SP3[i].ftotalz=P3[i].ftotalz;
		SP3[i].ftotal=P3[i].ftotal;

		if(xi0>=maxb) SP1[i].i_type=1;
		else SP1[i].i_type=2;
	}else{
		p2p_af[i]=2;
		p2p_idx[i]=i;
		SP1[i].i_type=3;
	}
}
//////////////////////////////////////////////////////////////////////////////
__global__ void left_send_particle(int_t*p2p_af,int_t*p2p_idx,part1*P1,part3*P3,part1*SP1,p2p_part3*SP3,int_t tid){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type!=1){
		p2p_af[i]=2;
		p2p_idx[i]=i;
		SP1[i].i_type=3;
		return;
	}

	Real minb,xi0;
	xi0=P1[i].x;

	minb=k_x_min+k_calc_area*k_dcell*tid;

	if(xi0<minb+k_dcell*2.0){
		p2p_af[i]=1;
		p2p_idx[i]=i;

		SP1[i]=P1[i];
		SP3[i].drho=P3[i].drho;
		SP3[i].dconcn=P3[i].dconcn;
		SP3[i].denthalpy=P3[i].denthalpy;
		SP3[i].ftotalx=P3[i].ftotalx;
		SP3[i].ftotaly=P3[i].ftotaly;
		SP3[i].ftotalz=P3[i].ftotalz;
		SP3[i].ftotal=P3[i].ftotal;

		if(xi0<minb) SP1[i].i_type=1;
		else SP1[i].i_type=2;

	}else{
		p2p_af[i]=2;
		p2p_idx[i]=i;
		SP1[i].i_type=3;
	}
}
//////////////////////////////////////////////////////////////////////////////
__global__ void reorder_data_p2p(int_t*p2p_af,int_t*p2p_idx,part1*SP1,part1*s_SP1,p2p_part3*SP3,p2p_part3*s_SP3){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_p2p) return;
	//if(p2p_af[i]==2) return;

	int sortedIndex=p2p_idx[i];

	s_SP1[i]=SP1[sortedIndex];
	s_SP3[i]=SP3[sortedIndex];
}
//////////////////////////////////////////////////////////////////////////////
__global__ void p2p_copyData(part1*SP1,part3*SP3,part1*rP1,p2p_part3*rP3,int_t tid){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_p2p) return;
	//if(rP1[i].i_type!=1) return;

	int ps=k_num_part2-i-1;
	int ri=rP1[i].i_type;

	// if(SP1[ps].i_type==1) return;
	if(ri==3){
		SP1[ps].i_type=3;		// 들어온 입자가 더미이면 더미로 기록
	}else{								// 아니면 정보 복사
		SP1[ps]=rP1[i];
		SP3[ps].drho=rP3[i].drho;
		SP3[ps].dconcn=rP3[i].dconcn;
		SP3[ps].denthalpy=rP3[i].denthalpy;
		SP3[ps].ftotalx=rP3[i].ftotalx;
		SP3[ps].ftotaly=rP3[i].ftotaly;
		SP3[ps].ftotalz=rP3[i].ftotalz;
		SP3[ps].ftotal=rP3[i].ftotal;
	}
	rP1[i].i_type=3;			// 작업이 끝나면 더미로 설정
}
//////////////////////////////////////////////////////////////////////////////
__global__ void init_Recv(part1*rP1){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_p2p) return;

	rP1[i].i_type=3;
}
//////////////////////////////////////////////////////////////////////////////
__global__ void initial_particle(part1*P1,int_t tid){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;

	if(P1[i].i_type==3) return;
	if(P1[i].i_type==2){
		P1[i].i_type=3;
		return;
	}

	Real xi0=P1[i].x;

	Real maxb,minb;

	minb=k_x_min+(Real)(k_calc_area*tid)*k_dcell;
	maxb=k_x_min+(Real)(k_calc_area*(tid+1))*k_dcell;

	if(tid==0){
		if(xi0>=maxb+k_dcell*2.0||xi0<k_x_min){
			P1[i].i_type=3;
			return;
		}
		if(xi0>=maxb&&xi0<maxb+k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(maxb-k_dcell*2.0))&&(xi0>=(k_x_min+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}else if(tid==k_ngpu-1){
		if(xi0>k_x_max||xi0<minb-k_dcell*2.0){
			P1[i].i_type=3;
			return;
		}
		if(xi0<minb&&xi0>=minb-k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(k_x_max-k_dcell*2.0))&&(xi0>=(minb+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}else{
		if(xi0>=maxb+k_dcell*2.0||xi0<minb-k_dcell*2.0){
			P1[i].i_type=3;
			return;
		}
		if(xi0>=maxb&&xi0<maxb+k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if(xi0<minb&&xi0>=minb-k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(maxb-k_dcell*2.0))&&(xi0>=(minb+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}

}
//////////////////////////////////////////////////////////////////////////////
__global__ void inner_outer_particle(part1*P1,int_t tid){
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	// if(P1[i].i_type==2||P1[i].i_type==3) return;
	if(P1[i].i_type==3) return;

	Real xi0=P1[i].x;

	Real maxb,minb;

	minb=k_x_min+(Real)(k_calc_area*tid)*k_dcell;
	maxb=k_x_min+(Real)(k_calc_area*(tid+1))*k_dcell;

	if(tid==0){
		if(xi0>=maxb+k_dcell*2.0||xi0<k_x_min){
			P1[i].i_type=3;
			return;
		}
		if(xi0>=maxb&&xi0<maxb+k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(maxb-k_dcell*2.0))&&(xi0>=(k_x_min+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}else if(tid==k_ngpu-1){
		if(xi0>k_x_max||xi0<minb-k_dcell*2.0){
			P1[i].i_type=3;
			return;
		}
		if(xi0<minb&&xi0>=minb-k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(k_x_max-k_dcell*2.0))&&(xi0>=(minb+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}else{
		if(xi0>=maxb+k_dcell*2.0||xi0<minb-k_dcell*2.0){
			P1[i].i_type=3;
			return;
		}
		if(xi0>=maxb&&xi0<maxb+k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if(xi0<minb&&xi0>=minb-k_dcell*2.0){
			P1[i].i_type=2;
			return;
		}
		if((xi0<(maxb-k_dcell*2.0))&&(xi0>=(minb+k_dcell*2.0))){
			P1[i].i_type=0;
		}else{
			P1[i].i_type=1;
		}
	}

}
////////////////////////////////////////////////////////////////////////
// Deactivate IBM (building) particles for a building-free run.
// Sets i_type=3 on every p_type==1000 particle, which makes them:
//   - skipped as target i by every kernel  (i_type>i_type_crt guard)
//   - dropped into the dummy cell k_num_cells by KERNEL_index_particle_to_cell
//     -> never visited as a neighbour j either.
// One-shot: launch once at init. i_type=3 is preserved afterwards
// (inner_outer_particle / initial_particle early-return on i_type==3).
// Toggle: Parameters.cuh  #define DEACTIVATE_IBM
__global__ void KERNEL_deactivate_IBM(part1*P1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type==1000) P1[i].i_type=3;
}
////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_reorder(int_t*g_idx,int_t*p_idx,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part1*SP1,part2*SP2)
{
	extern __shared__ int sharedHash[];

	int idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part2) return;
	int hash;

	hash=g_idx[idx];

	//if (hash>=k_num_cells && P1[idx].i_type!=3) printf("itype=%d hash %d k_num_cells %d\n", P1[idx].i_type, hash, k_num_cells);

	sharedHash[threadIdx.x+1]=hash;
	if(idx>0&&threadIdx.x==0){
		/*save the end of the previous block g_idx*/
		sharedHash[0]=g_idx[idx-1];
	}
	__syncthreads();		// for sorting and reorder particle property
	if(idx==0||hash!=sharedHash[threadIdx.x]){
		//if(hash<=k_num_cells) {
			g_str[hash]=idx;
			if(idx>0) g_end[sharedHash[threadIdx.x]]=idx;
		//}
	}
	// if((idx==k_num_part2-1)&&(hash<k_num_cells)) g_end[hash]=idx+1;
	if((idx==k_num_part2-1)) g_end[hash]=idx+1;
	/*reorder data*/
	int sortedIndex=p_idx[idx];

	SP1[idx]=P1[sortedIndex];
	SP2[idx]=P2[sortedIndex];
}
__global__ void KERNEL_reorder_LDM(int_t*g_idx,int_t*p_idx,int_t*g_str,int_t*g_end,L_part1*P1,L_part2*P2,L_part1*SP1,L_part2*SP2)
{
	extern __shared__ int sharedHash[];

	int idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part_LDM) return;
	int hash;

	hash=g_idx[idx];

	//if (hash>=k_num_cells && P1[idx].i_type!=3) printf("itype=%d hash %d k_num_cells %d\n", P1[idx].i_type, hash, k_num_cells);

	sharedHash[threadIdx.x+1]=hash;
	if(idx>0&&threadIdx.x==0){
		/*save the end of the previous block g_idx*/
		sharedHash[0]=g_idx[idx-1];
	}
	__syncthreads();		// for sorting and reorder particle property
	if(idx==0||hash!=sharedHash[threadIdx.x]){
		//if(hash<=k_num_cells) {
			g_str[hash]=idx;
			if(idx>0) g_end[sharedHash[threadIdx.x]]=idx;
		//}
	}
	// if((idx==k_num_part_LDM-1)&&(hash<k_num_cells)) g_end[hash]=idx+1;
	if((idx==k_num_part_LDM-1)) g_end[hash]=idx+1;
	/*reorder data*/
	int sortedIndex=p_idx[idx];

	SP1[idx]=P1[sortedIndex];
	SP2[idx]=P2[sortedIndex];
}
////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_reorder_single(int_t*b_idx,int_t*p_idx,int_t*b_str,int_t*b_end,part1*P1,part2*P2,part1*SP1,part2*SP2)
{
	extern __shared__ int sharedHash[];

	int idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part2) return;
	int hash;

	hash=b_idx[idx];

	//if (hash>=k_num_cells && P1[idx].i_type!=3) printf("itype=%d hash %d k_num_cells %d\n", P1[idx].i_type, hash, k_num_cells);

	sharedHash[threadIdx.x+1]=hash;
	if(idx>0&&threadIdx.x==0){
		/*save the end of the previous block g_idx*/
		sharedHash[0]=b_idx[idx-1];
	}
	__syncthreads();		// for sorting and reorder particle property
	if(idx==0||hash!=sharedHash[threadIdx.x]){
		//if(hash<=k_num_cells) {
			b_str[hash]=idx;
			if(idx>0) b_end[sharedHash[threadIdx.x]]=idx;
		//}
	}
	// if((idx==k_num_part2-1)&&(hash<k_num_cells)) g_end[hash]=idx+1;
	if((idx==k_num_part2-1)) b_end[hash]=idx+1;
	/*reorder data*/
	int sortedIndex=p_idx[idx];

	SP1[idx]=P1[sortedIndex];
	SP2[idx]=P2[sortedIndex];
}
////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_index_particle_to_cell(int_t*g_idx,int_t*p_idx,part1*P1)
{
	int_t idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part2) return;
	// Not in Each GPUs Domain
	if(P1[idx].i_type==3){
		g_idx[idx]=k_num_cells-1;
		p_idx[idx]=idx;
	}
	else {
		int_t icell,jcell,kcell;
		// calculate I,J,K in cell
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((P1[idx].x-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((P1[idx].y-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((P1[idx].z-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
		// out-of-range handling
		if(icell<0) icell=0;
		if(jcell<0) jcell=0;
		if(kcell<0) kcell=0;
		// calculate cell index from I,J,K
		p_idx[idx]=idx;
		g_idx[idx]=idx_cell(icell,jcell,kcell);
	}

	if(g_idx[idx]==k_num_cells-1 && P1[idx].i_type!=3) printf("itype %d \n", P1[idx].i_type);
}

__global__ void KERNEL_index_particle_to_cell_LDM(int_t*g_idx,int_t*p_idx,L_part1*P1)
{
	int_t idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part_LDM) return;
	// Not in Each GPUs Domain
	if(P1[idx].i_type==3){
		g_idx[idx]=k_num_cells-1;
		p_idx[idx]=idx;
	}
	else {
		int_t icell,jcell,kcell;
		// calculate I,J,K in cell
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((P1[idx].x-k_x_min)/(k_x_max-k_x_min)*k_NI),k_NI-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((P1[idx].y-k_y_min)/(k_y_max-k_y_min)*k_NJ),k_NJ-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((P1[idx].z-k_z_min)/(k_z_max-k_z_min)*k_NK),k_NK-1);}
		// out-of-range handling
		if(icell<0) icell=0;
		if(jcell<0) jcell=0;
		if(kcell<0) kcell=0;
		// calculate cell index from I,J,K
		p_idx[idx]=idx;
		g_idx[idx]=idx_cell(icell,jcell,kcell);
	}

	if(g_idx[idx]==k_num_cells-1 && P1[idx].i_type!=3) printf("itype %d \n", P1[idx].i_type);
}

__global__ void KERNEL_index_particle_to_block(int_t*b_idx,int_t*p_idx,part1*P1)
{
	int_t idx=threadIdx.x+blockIdx.x*blockDim.x;
	if(idx>=k_num_part2) return;
	// Not in Each GPUs Domain
	if(P1[idx].i_type==3){
		b_idx[idx]=k_num_blocks;
		p_idx[idx]=idx;
	}
	else {
		int_t icell,jcell,kcell;
		int ki=k_NBI;
		int kj=k_NBJ;
		int kk=k_NBK;

		// calculate I,J,K in cell
		if((k_x_max==k_x_min)){icell=0;}
		else{icell=min(floor((P1[idx].x-k_x_min)/(k_x_max-k_x_min)*ki),ki-1);}
		if((k_y_max==k_y_min)){jcell=0;}
		else{jcell=min(floor((P1[idx].y-k_y_min)/(k_y_max-k_y_min)*kj),kj-1);}
		if((k_z_max==k_z_min)){kcell=0;}
		else{kcell=min(floor((P1[idx].z-k_z_min)/(k_z_max-k_z_min)*kk),kk-1);}

		// out-of-range handling
		if(icell<0) icell=0;
		if(jcell<0) jcell=0;
		if(kcell<0) kcell=0;
		// calculate cell index from I,J,K
		p_idx[idx]=idx;
		b_idx[idx]=(icell)+ki*(jcell)+ki*kj*(kcell);
	}
}

////////////////////////////////////////////////////////////////////////
__global__ void KERNEL_set_ncell(part1*P1)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].i_type==3) return;

	Real h0=P1[i].h_ref;
	Real h=P1[i].h;

	// P1[i].ncell = ceil(float(ncell_init)*(h/h0));
	P1[i].ncell = 1.0*ceil(float(ncell_init)/float(h0/h));
	//if(P1[i].y>2.4) P1[i].ncell = 8;
}
////////////////////////////////////////////////////////////////////////