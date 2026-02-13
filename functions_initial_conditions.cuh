double Interpolate_ux_caseC(Real z_inp, Real scale){
  Real z_scale=scale;             // (Building height / 0.08)
  Real z=z_scale*z_inp;
  Real u_scale=1.0;
  static const double z_data[] = {
        0.01, 0.02, 0.03, 0.05, 0.07,
		0.1, 0.15, 0.2, 0.25, 0.3,
		0.35, 0.4, 0.6, 0.8, 1.0, 1.2
    };
    static const double u_data[] = {
        2.372, 2.434, 2.548, 2.912, 3.042,
		3.392, 3.483, 3.654, 3.82, 4.019,
		4.205, 4.343, 4.985, 5.713, 6.12, 6.201
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
    return u_scale*((1.0 - t) * u0 + t * u1);         // 선형 보간
}

double Interpolate_ux_caseA(Real z_inp, Real scale){
  Real z_scale=scale;             // (Building height / 0.16)
  Real z=z_scale*z_inp;
  Real u_scale=1.0;
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
        //6.751, 6.0, 4.0, 2.0, 0.0
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
    return u_scale*((1.0 - t) * u0 + t * u1);         // 선형 보간
}

double Interpolate_k_caseA(Real z_inp, Real scale){
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

double Interpolate_k_caseAS(Real z_inp){
  Real z_scale = 0.16/20.0;             // (0.16 / Buidling height)
  Real u_scale = 1.5;
  Real z=z_scale*z_inp;
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
    return ((1.0 - t) * u0 + t * u1)*u_scale;         // 선형 보간
}

double Interpolate_ux_ABL(Real z_inp){

	Real a=0.19;
	Real u_ref=10.0;
	Real z_ref=100.0;
	
	Real ux=0.0;
	// switch(terrain_type){
	// 	case 0:
	// 		a=0.14;
	// 		z_ref=300.0;
	// 		printf("terrain type: flat\n");
	// 		break;
	// 	case 1:
	// 		a=0.19;
	// 		z_ref=350.0;
	// 		printf("terrain type: Other\n");
	// 		break;
	// 	case 2:
	// 		a=0.19;
	// 		z_ref=350.0;
	// 		printf("terrain type: other\n");
	// 		break;
	// 	case 3:
	// 		a=0.19;
	// 		z_ref=350.0;
	// 		printf("terrain type: other\n");
	// 		break;
	// 	case 4:
	// 		a=0.19;
	// 		z_ref=350.0;
	// 		printf("terrain type: other\n");
	// 		break;
	// 	default:
	// 		a=0.19;
	// 		z_ref=350.0;
	// 		printf("terrain type: flat\n");
	// 		break;
	// }

	ux=u_ref*pow(z_inp/z_ref,a);
	if(z_inp<0.0) ux=0.0;

	return ux;
}