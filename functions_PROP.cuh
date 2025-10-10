////////////////////////////////////////////////////////////////////////
#define NU0_HB		1.0			//Herschel-Bulkey model parameter
#define TAU0_HB		18.24		//Herschel-Bulkey model parameter (for lava flow)
#define K0_HB		1.90		//Herschel-Bulkey model parameter (for lava flow)
#define	N0_HB		0.53		//Herschel-Bulkey model parameter (for lava flow)
////////////////////////////////////////////////////////////////////////

__host__ __device__ Real interp2(Real *x_data,Real *y_data,int size,Real x)
{

	Real y;
	int i;
	int end_idx=size-1;

	if(x_data[end_idx]<x){
		//y=y_data[end_idx]+(y_data[end_idx]-y_data[end_idx-1])/(x_data[end_idx]-x_data[end_idx-1])*(x-x_data[end_idx]);
		y=y_data[end_idx]-y_data[end_idx-1];
		y/=(x_data[end_idx]-x_data[end_idx-1]);
		y*=(x-x_data[end_idx]);
		y+=y_data[end_idx];
	}else if(x<=x_data[0]){
		//y=y_data[0]+(y_data[1]-y_data[0])/(x_data[1]-x_data[0])*(x-x_data[0]);
		y=y_data[1]-y_data[0];
		y/=(x_data[1]-x_data[0]);
		y*=(x-x_data[0]);
		y+=y_data[0];
	}else{
		for(i=0;i<size;i++){
			if((x_data[i]<x)&(x<=x_data[i+1])){
				//y=y_data[i]+(y_data[i+1]-y_data[i])/(x_data[i+1]-x_data[i])*(x-x_data[i]);
				y=y_data[i+1]-y_data[i];
				y/=(x_data[i+1]-x_data[i]);
				y*=(x-x_data[i]);
				y+=y_data[i];
				break;
			}
		}
	}

	return y;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__ Real htoT(Real enthalpy,uint_t p_type)
{
	// corium data
	Real x_data1[5]={-554932,215789,757894,905263,1268421};
	Real y_data1[5]={280,1537,2300,2450,2650};

	Real y;
	int index0,table_size0;

	p_type=abs(p_type);

	if (k_prop_table)
	{
		index0=k_table_index[p_type];
		table_size0=k_table_size[p_type];

		y=interp2(&k_Tab_h[index0],&k_Tab_T[index0],table_size0,enthalpy);
	}
	else
	{
		y=interp2(x_data1,y_data1,5,enthalpy);
	}

	return y;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__ Real Ttoh(Real temp,uint_t p_type)
{
	// corium data
	Real y_data1[5]={-554932,215789,757894,905263,1268421};
	Real x_data1[5]={280,1537,2300,2450,2650};

	Real y;
	int index0,table_size0;

	p_type=abs(p_type);

	if (k_prop_table)
	{
		index0=k_table_index[p_type];
		table_size0=k_table_size[p_type];

		y=interp2(&k_Tab_T[index0],&k_Tab_h[index0],table_size0,temp);
	}
	else
	{
		y=interp2(x_data1,y_data1,5,temp);
	}

	return y;
}

////////////////////////////////////////////////////////////////////////
__host__ __device__ Real viscosity(Real temp,uint_t p_type)
{
	Real vis;
	int index0,table_size0;

	p_type=abs(p_type);

	if (k_prop_table)
	{
		index0=k_table_index[p_type];
		table_size0=k_table_size[p_type];

		vis=interp2(&k_Tab_T[index0],&k_Tab_vis[index0],table_size0,temp);

		if(p_type==1000) printf("index0:%d table_size0:%d vis=%f p_type=%d\n",index0, table_size0, vis, p_type);
	}
	else
	{
			vis=1.09e-3;
	}

	return vis;
}


__host__ __device__ Real viscosity2(Real temp,uint_t p_type)
{
	Real vis;

	if(p_type==1) vis=1.8e-5;
	else if(p_type==0) vis=1.0e-2;
	else if(p_type==1000) vis=1.0e-2;
	else vis=1.8e-5;

	return vis;
}
////////////////////////////////////////////////////////////////////////
__device__ Real specific_heat(Real temp, uint_t p_type)
{
	Real cpi;
	int index0,table_size0;

	p_type=abs(p_type);

	cpi=0.0;

	if(k_prop_table){
		index0=k_table_index[p_type];
		table_size0=k_table_size[p_type];
		printf("index0=%d\n",index0);

		cpi=interp2(&k_Tab_T[index0],&k_Tab_cp[index0],table_size0,temp);
	}
	// else{
	// 	//printf("NO");
	// 		cpi=4200;
	// }

	return cpi;
}

////////////////////////////////////////////////////////////////////////
__device__ Real conductivity(Real temp,uint_t p_type)
{
	Real cond;
	int index0,table_size0;

	p_type=abs(p_type);

	cond=0.0;

	if (k_prop_table)
	{
		index0=k_table_index[p_type];
		table_size0=k_table_size[p_type];

		cond=interp2(&k_Tab_T[index0],&k_Tab_k[index0],table_size0,temp);
	}
	// else {
	// 	//printf("NO");
	// 		cond=1.65*200;
	// }

	return cond;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__ Real sigma(Real temp,uint_t p_type)
{
	Real y;

	p_type=abs(p_type);

	y = 24.5;
	return y;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__ Real diffusion_coefficient(Real temp,uint_t p_type)
{
	Real y;

	p_type=abs(p_type);

	y=2.0;

	return y;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__ Real thermal_expansion(Real temp,uint_t p_type)
{
	Real y;

	p_type=abs(p_type);

	y=1.0/T_ini;

	return y;
}
////////////////////////////////////////////////////////////////////////
__host__ __device__  Real K_to_eta(Real tK_stiff)
{
	// corium data
	Real x_data1[5]={0,1000,5000,25000,100000};
	Real y_data1[5]={1.0,1.8,2.3,2.5,2.5};

	Real y=1;

	y=interp2(x_data1,y_data1,5,tK_stiff);

	return y;
}
