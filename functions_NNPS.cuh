
//////////////////////////////////////////////////////////////////////////////
// Single building (Case H, n=10)
#define L1 	0   // min(x)+1.5*initial_particle_spacing
#define L2	0.028 // min(x)+3.5*initial_particle_spacing
#define L3	13.304 // max(x)-3.5*initial_particle_spacing
#define L4	2210 // max(x)+1.0*initial_particle_spacing

#define L_left -0.412// min(y)+3.5*initial_particle_spacing
#define L_right 0.412 // max(y)-3.5*initial_particle_spacing
#define L_ceiling 1.910

#define L_MOST_x0 0.0
#define L_MOST_xmax 1990.0
//
// // Single building (Case H, n=12)
// #define L1 	0   // min(x)+1.5*initial_particle_spacing
// #define L2	0.02333 // min(x)+3.5*initial_particle_spacing
// #define L3	13.304 // max(x)-3.5*initial_particle_spacing
// #define L4	2210 // max(x)+1.0*initial_particle_spacing
//
// #define L_left -0.416667// min(y)+3.5*initial_particle_spacing
// #define L_right 0.416667 // max(y)-3.5*initial_particle_spacing
// #define L_ceiling 1.910
//
// #define L_MOST_x0 0.0
// #define L_MOST_xmax 1990.0

// // Single building (Case H, n=10, Scale up)
// #define L1 	0   // min(x)+1.5*initial_particle_spacing
// #define L2	3.5 // min(x)+3.5*initial_particle_spacing
// #define L3	2013.304 // max(x)-3.5*initial_particle_spacing
// #define L4	2210 // max(x)+1.0*initial_particle_spacing

// #define L_left -62.5// min(y)+3.5*initial_particle_spacing
// #define L_right 62.5 // max(y)-3.5*initial_particle_spacing
// #define L_ceiling 121.0

// #define L_MOST_x0 0.0
// #define L_MOST_xmax 1990.0

// // Single building (Case K, n=10, Scale up)
// #define L1 	0   // min(x)+1.5*initial_particle_spacing
// #define L2	3.5 // min(x)+3.5*initial_particle_spacing
// #define L3	2013.304 // max(x)-3.5*initial_particle_spacing
// #define L4	2210 // max(x)+1.0*initial_particle_spacing
//
// #define L_left -92.91666666  // min(y)+3.5*initial_particle_spacing
// #define L_right 92.91666666  // max(y)-3.5*initial_particle_spacing
// #define L_ceiling 121.0
//
// #define L_MOST_x0 0.0
// #define L_MOST_xmax 1990.0


__host__ float hash_random(int i, int seed) {
    uint32_t x = i * 73856093 ^ seed;
    x = (x >> 13) ^ x;
    return ((x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 2147483647.0f;
}

__host__ float box_muller(int i, int seed_a, int seed_b) {
    float u1 = hash_random(i, seed_a);
    float u2 = hash_random(i, seed_b);

    return sqrtf(-2.0f * logf(u1 + 1e-6f)) * cosf(2 * M_PI * u2);
}
// Single building (dx=2.5m, boundary width = 4)
// #define L1 	0.0   // min(x)
// #define L2	8.75 // min(x)+3.5*initial_particle_spacing
// #define L3	1148.25 // max(x)-3.5*initial_particle_spacing
// #define L4	1150.5 // max(x)+1.0*initial_particle_spacing
//
// #define L_left -241.25// min(y)+3.5*initial_particle_spacing
// #define L_right 241.25 // max(y)-3.5*initial_particle_spacing
// #define L_ceiling 301.25
//////////////////////////////////////////////////////////////////////////////
// Single building (dx=2.5m, boundary width = 4)
// #define L1 	0.0   // min(x)
// #define L2	8.75 // min(x)+3.5*initial_particle_spacing
// #define L3	1148.25 // max(x)-3.5*initial_particle_spacing
// #define L4	1150.5 // max(x)+1.0*initial_particle_spacing
//
// #define L_left 8.75  // min(y)+3.5*initial_particle_spacing
// #define L_right 491.25  // max(y)-3.5*initial_particle_spacing

// // Flat terrain (dx=0.9m, boundary width = 4)
// #define L1 	0.0   // min(x)
// #define L2	3.15 // min(x)+3.5*initial_particle_spacing
// #define L3	196.5 // max(x)-3.5*initial_particle_spacing
// #define L4	1101.0 // max(x)+1.0*initial_particle_spacing
//
// #define L_left 3.15 // min(y)+3.5*initial_particle_spacing
// #define L_right 46.35 // max(y)-3.5*initial_particle_spacing
//////////////////////////////////////////////////////////////////////////////

// // Flat terrain (dx=1.1m, boundary width = 4)
// #define L1 	0.0   // min(x)
// #define L2	3.85 // min(x)+3.5*initial_particle_spacing
// #define L3	196.5 // max(x)-3.5*initial_particle_spacing
// #define L4	1101.0 // max(x)+1.0*initial_particle_spacing
//
// #define L_left 3.85 // min(y)+3.5*initial_particle_spacing
// #define L_right 45.65 // max(y)-3.5*initial_particle_spacing
//////////////////////////////////////////////////////////////////////////////

double Interpolate_ux_caseH(Real z_inp, Real scale){
  Real scale_factor=scale;             // (0.16 / Buidling height)
  Real z=scale_factor*z_inp;
  static const double z_data[] = {
        0.0050000, 0.0100000, 0.0200000, 0.0400000, 0.0600000,
        0.0800000, 0.1000000, 0.1200000, 0.1400000, 0.1600000,
        0.1650000, 0.1700000, 0.1800000, 0.1900000, 0.2100000,
        0.2300000, 0.2500000, 0.3000000, 0.3500000, 0.4000000,
        0.4500000, 0.5000000, 0.5500000, 0.6000000, 0.6528736,
        0.7724136, 0.8390808, 0.8689656, 0.8804600, 0.8919544
    };
    static const double u_data[] = {
        2.745, 2.935, 3.175, 3.435, 3.627,
        3.824, 4.021, 4.210, 4.362, 4.491,
        4.502, 4.586, 4.606, 4.712, 4.854,
        4.993, 5.132, 5.449, 5.782, 6.077,
        6.338, 6.588, 6.693, 6.751, 6.751,
        // 6.751, 6.316, 5.566, 5.132, 4.263
		    6.751, 6.751, 6.751, 6.751, 6.751
    };
    static const std::size_t N = sizeof(z_data)/sizeof(z_data[0]);

    // 아래/위 한 점으로 외삽
    if (z <= z_data[0]) {
        double dz = z_data[1] - z_data[0];
        double du = u_data[1] - u_data[0];
        return u_data[0] + ( (dz == 0.0) ? 0.0 : (z - z_data[0]) * (du/dz) );
    }
    if (z >= z_data[N-1]) {
        double dz = z_data[N-1] - z_data[N-2];
        double du = u_data[N-1] - u_data[N-2];
        return u_data[N-2] + ( (dz == 0.0) ? 0.0 : (z - z_data[N-2]) * (du/dz) );
    }

    // 구간 탐색 (z_data[i-1] <= z < z_data[i])
    const double* it = std::lower_bound(z_data, z_data + N, z);
    std::size_t i = static_cast<std::size_t>(it - z_data);
    // 안전장치 (이론상 도달 X)
    if (i == 0) i = 1;

    double z0 = z_data[i-1], z1 = z_data[i];
    double u0 = u_data[i-1], u1 = u_data[i];
    double t  = (z - z0) / (z1 - z0);       // 0~1
    return (1.0 - t) * u0 + t * u1;         // 선형 보간
}

double Interpolate_k_caseH(Real z_inp, Real scale){
  Real scale_factor=scale;             // (0.16 / Buidling height)
  Real z=scale_factor*z_inp;
  static const double z_data[] = {
        0.00393,  0.01178,  0.02356,  0.04123,  0.06086,
        0.08245,  0.10209,  0.11975,  0.14135,  0.16294,
        0.21202,  0.23166,  0.25325,  0.30233,  0.35141,
        0.40442,  0.4535,   0.50061,  0.55362,  0.60466,
        0.70282,  0.74209
    };
    static const double u_data[] = {
        0.3687,   0.39304,  0.4313,   0.49043,  0.53913,
        0.56696,  0.58087,  0.62261,  0.62957,  0.64696,
        0.66087,  0.64696,  0.62261,  0.5913,   0.51826,
        0.4313,   0.3513,   0.25043,  0.17043,  0.10087,
        0.0,      0.0
    };
    static const std::size_t N = sizeof(z_data)/sizeof(z_data[0]);

    // 아래/위 한 점으로 외삽
    if (z <= z_data[0]) {
        double dz = z_data[1] - z_data[0];
        double du = u_data[1] - u_data[0];
        return u_data[0] + ( (dz == 0.0) ? 0.0 : (z - z_data[0]) * (du/dz) );
    }
    if (z >= z_data[N-1]) {
        double dz = z_data[N-1] - z_data[N-2];
        double du = u_data[N-1] - u_data[N-2];
        return u_data[N-2] + ( (dz == 0.0) ? 0.0 : (z - z_data[N-2]) * (du/dz) );
    }

    // 구간 탐색 (z_data[i-1] <= z < z_data[i])
    const double* it = std::lower_bound(z_data, z_data + N, z);
    std::size_t i = static_cast<std::size_t>(it - z_data);
    // 안전장치 (이론상 도달 X)
    if (i == 0) i = 1;

    double z0 = z_data[i-1], z1 = z_data[i];
    double u0 = u_data[i-1], u1 = u_data[i];
    double t  = (z - z0) / (z1 - z0);       // 0~1
    return (1.0 - t) * u0 + t * u1;         // 선형 보간
}

void c_initial_inner_outer_particle_single(part1*HP1,part1*DHP1,int_t tid){

	int_t i,c_count;
	Real xi0;
	Real maxb,minb;

	c_count=0;

	for(i=0;i<num_part;i++){
		HP1[i].i_type=1;
		HP1[i].buffer_type=0;
		// HP1[i].h=1.5;
		// HP1[i].h_ref=HP1[i].h;
		HP1[i].m_ref=HP1[i].m;
		HP1[i].uy=0.0;
		HP1[i].uz=0.0;
		HP1[i].temp=300.0;
    Real scale_factor=0.16/0.16;       // 0.16/building_height
    HP1[i].ux=Interpolate_ux_caseH(HP1[i].z,scale_factor);
    HP1[i].k_turb=Interpolate_k_caseH(HP1[i].z,scale_factor);
    // if((HP1[i].x-765.0)*(HP1[i].x-765.0)+(HP1[i].y-0.0)*(HP1[i].y-0.0)+(HP1[i].z-2.5)*(HP1[i].z-2.5)<1e-3) {HP1[i].source=1;}
    // else if((HP1[i].x-425.0)*(HP1[i].x-425.0)+(HP1[i].y-5.0)*(HP1[i].y-5.0)+(HP1[i].z-2.5)*(HP1[i].z-2.5)<1e-3) {HP1[i].source=1;}
    // else if((HP1[i].x-435.0)*(HP1[i].x-435.0)+(HP1[i].y-155.0)*(HP1[i].y-155.0)+(HP1[i].z-12.5)*(HP1[i].z-12.5)<1e-3) {HP1[i].source=1;}
    // else if((HP1[i].x-500.0)*(HP1[i].x-500.0)+(HP1[i].y-5.0)*(HP1[i].y-5.0)+(HP1[i].z-57.5)*(HP1[i].z-57.5)<1e-3) {HP1[i].source=1;}
    //HP1[i].temp=310.0-10.0*(log(HP1[i].z+1.0)/log(301.0));
		//HP1[i].pres=HP1[i].rho*Gravitational_CONST*(L_ceiling+20.0-HP1[i].z);
    //HP1[i].pres0=HP1[i].pres;
		// //HP1[i].blood_index=0;
		// //HP1[i].M_num=0;
		// if(HP1[i].p_type!=1){
		// 	if(HP1[i].z<0) HP1[i].p_type=0;
		// 	else if(HP1[i].y<0||HP1[i].y>50||HP1[i].z>30) {
		// 		HP1[i].p_type=-1;
		// 		HP1[i].ux=1.0;
		// 	}
		// }
		//HP1[i].ux=0.0;
		//if(HP1[i].p_type==1000) HP1[i].i_type=3;
		// if(HP1[i].p_type==3) {    // for SNU_IBM (KDH)
		// 	if(HP1[i].y>261 || HP1[i].y<16) HP1[i].p_type=-1;
		// 	else if(HP1[i].z<4) HP1[i].i_type=3;
		// 	else HP1[i].p_type=1;
		// }

		if(HP1[i].p_type==3) {
			//HP1[i].i_type=2;
			HP1[i].p_type=1;
		}

		if(HP1[i].p_type==0||HP1[i].p_type==1000) {
			//HP1[i].i_type=2;
			HP1[i].rho=1.225;
			HP1[i].rho=1.225;
		}

		if (open_boundary>0)			// (CAUTION)
		{
			// if (HP1[i].p_type>=0) HP1[i].ux=Inlet_Velocity;
			if(HP1[i].p_type==1||HP1[i].p_type==-1){

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

				// if(HP1[i].buffer_type<3 && HP1[i].x>L_MOST_x0 && HP1[i].x<L_MOST_xmax){
				// 	if(HP1[i].z<=5.0) HP1[i].MOST_buffer=1;
				// 	else if(HP1[i].z<=10.0) HP1[i].MOST_buffer=2;
				// }

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
		if(dim==3){
			if (HP1[i].p_type==1 || HP1[i].p_type==-1){
				//if (HP1[i].buffer_type==1){

					Real a=0.27;
					Real z_ref=0.0;
					Real u_ref=5.0;
					Real z=0.0;

					switch(terrain_type){
						case 0:
						a=0.14;
						z_ref=2.0;
						printf("terrain type: flat\n");
						break;
						case 1:
						a=0.19;
						z_ref=350.0;
						printf("terrain type: Grass\n");
						break;
						case 2:
						a=0.19;
						z_ref=350.0;
						printf("terrain type: Forest\n");
						break;
						case 3:
						a=0.19;
						z_ref=350.0;
						//printf("terrain type: City\n");
						break;
						case 4:
						a=0.19;
						z_ref=350.0;
						printf("terrain type: Ocean\n");
						break;
						default:
						a=0.19;
						z_ref=350.0;
						printf("terrain type: flat\n");
						break;
					}
					if (HP1[i].z>=0) z=HP1[i].z;
					else z=0.0;
					//HP1[i].ux=u_ref*pow(z/z_ref,a);
					//HP1[i].ux=u_ref;
					//HP1[i].ux=(us_ini/k_vonKarman)*log(HP1[i].z/MOST_z0);

					if (HP1[i].z<0) HP1[i].ux=0.0;

					// // C++/CUDA 식 예시
					// float rand1 = box_muller(i, 1234, 5678);
					// float rand2 = box_muller(i, 2345, 6789);
					// float rand3 = box_muller(i, 3456, 7890);
					//
					// float delta = 0.02;  // perturbation scale
					//
					// HP1[i].ux += delta * rand1;
					// HP1[i].uy += delta * rand2;
					// HP1[i].uz += delta * rand3;
			}
			else HP1[i].ux=0.0;
		}
		DHP1[c_count]=HP1[i];
		c_count++;
	}
}

void c_initial_velocity(part1*HP1,part1*DHP1,int_t tid){

	int_t c_count=0;
	Real a=0.0;
	Real z_ref=0.0;

	switch(terrain_type){
		case 0:
			a=0.14;
			z_ref=300.0;
			printf("terrain type: flat\n");
			break;
		case 1:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: Other\n");
			break;
		case 2:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		case 3:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		case 4:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: other\n");
			break;
		default:
			a=0.19;
			z_ref=350.0;
			printf("terrain type: flat\n");
			break;
	}
	Real u_ref=10.0;

	for(int_t i=0;i<num_part;i++){
		// Real z=max(HP1[i].y,0.0);
		// if(HP1[i].buffer_type==1) HP1[i].ux=u_ref*pow(z/z_ref,a);
		// else if(HP1[i].p_type==0) HP1[i].ux=0.0;//u_ref*pow(z/z_ref,a);
		// else HP1[i].ux=u_ref*pow(z/z_ref,a);//u_ref*pow(z/z_ref,a);

		if(HP1[i].p_type==1) HP1[i].ux=1.0;
		else HP1[i].ux=0.0;

		DHP1[c_count]=HP1[i];
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
	}else{	// 그 밖의 GPU에 대해서
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
		g_idx[idx]=k_num_cells;
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

	if(g_idx[idx]==k_num_cells && P1[idx].i_type!=3) printf("itype %d \n", P1[idx].i_type);
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

// __global__ void KERNEL_check(int_t*g_idx,int_t*p_idx,int_t*g_str,int_t*g_end,part1*P1,part2*P2,part1*SP1,part2*SP2)
// {
// 	extern __shared__ int sharedHash[];
//
// 	int idx=threadIdx.x+blockIdx.x*blockDim.x;
// 	if(idx>=k_num_part2) return;
//   Real xi, yi, zi;
// 	xi=P1[idx].x;
// 	yi=P1[idx].y;
// 	zi=P1[idx].z;
// 	if (xi>220 && xi<230 && yi<30 && yi>20 && zi> 25 && zi<30) printf("ptype:%d p_idx:%d g_idx:%d\n",P1[idx].p_type,p_idx[idx],g_idx[idx]);
// 	if (xi>215 && xi<225 && yi<35 && yi>25 && zi> 25 && zi<30) printf("ptype:%d p_idx:%d g_idx:%d\n",P1[idx].p_type,p_idx[idx],g_idx[idx]);
// }
