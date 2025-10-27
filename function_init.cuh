// int_t type variable define
#define		solver_type							vii[0]		// solver type: WCSPH/ISPH
#define		dim									vii[1]		// dimension
#define		prop_table							vii[2]		// property_table(Yes/NO)
#define		kernel_type							vii[3]		// kernel type
#define		time_type							vii[4]		// time stepping type
#define		open_boundary						vii[5]		// open boundary (YES / NO)
#define		flag_timestep_update				vii[6]		// flag for varying timestep update
#define		nb_cell_type						vii[7]		// neighbor cell type
#define		freq_filt							vii[8]		// filtering frequency
#define		freq_mass_sum						vii[9]		// mass summation frequency
#define		freq_temp							vii[10]		// temperature filtering frequency
#define		freq_output							vii[11]		// output frequency
#define		fp_solve							vii[12]		// solve pressure force ?
#define		fv_solve							vii[13]		// sovle viscous force ?
#define		fg_solve							vii[14]		// solve gravity force ?
#define		con_solve							vii[15]		// solve conduction?
#define		boussinesq_solve					vii[16]		// solve boussinesq approximation based natural convection?
#define		kgc_solve							vii[17]		// solve kernel gradient correction ?
#define		turbulence_model					vii[18]		// turbulence model	(by esk)
#define		concn_solve							vii[19]		// concentration diffusion model(PSH)
//
//psh:: ISPH input
#define		minIteration						vii[20]		// minimum number of PCISPH iteration
#define		maxIteration						vii[21]		// maximum number of PCISPH iteration
//
//solution variables
#define		nb_cell_number						vii[22]
#define		num_part							vii[23]
#define		NI									vii[24]
#define		NJ									vii[25]
#define		NK									vii[26]
#define		count								vii[27]
#define		num_cells							vii[28]
#define		ngpu								vii[29]
#define		calc_area							vii[30]
#define		num_p2p								vii[31]
#define		num_part2							vii[32]
#define		switch_ptype						vii[33]
#define   noslip_bc								vii[34]
#define   apr_solv                				vii[35]
#define 	h_change							vii[36]
#define   freq_apr								vii[37]
#define		freq_h_reset						vii[38]
#define		num_blocks							vii[39]
#define		freq_bapr							vii[40]
#define		NBI									vii[41]
#define		NBJ									vii[42]
#define		NBK									vii[43]
#define   structure_solv          				vii[44]
#define   atmos_stab		          			vii[45]  	// AeroSPHere (KDH)
#define   terrain_type	          				vii[46]  	// AeroSPHere (KDH)
#define	  num_part_LDM   						vii[47]  	// LDM frequency (KDH)
#define	  freq_LDM   							vii[48]  	// LDM frequency (KDH)

// Real type variable define
#define		kappa								vif[0]		// k in k*h
#define		p_ref								vif[1]		// reference pressure(for EOS)
#define		dt									vif[2]		// time-step(s)
#define		time								vif[3]		// time(s)
#define		time_end							vif[4]
// margin for simulation range
#define		u_limit								vif[5]

//psh:: ISPH input
#define		drho_th								vif[6]		// density convergence criterion
#define		dp_th								vif[7]		// pressure convergence criterion
#define		p_relaxation						vif[8]		// relaxation factor for PCISPH pressure
//
//solution variables
#define		x_min								vif[9]
#define		x_max								vif[10]
#define		y_min								vif[11]
#define		y_max								vif[12]
#define		z_min								vif[13]
#define		z_max								vif[14]
#define		nd_ref								vif[15]
//
#define		search_incr_factor					vif[16]
#define		search_kappa						vif[17]
#define		dcell								vif[18]
#define		buffersize							vif[19]
#define		dt_structure						vif[20]
////////////////////////////////////////////////////////////////////////
// int_t type variable define
#define		k_solver_type						k_vii[0]		// solver type: WCSPH/ISPH
#define		k_dim								k_vii[1]		// dimension
#define		k_prop_table						k_vii[2]		// property table (YES/NO)
#define		k_kernel_type						k_vii[3]		// kernel type
#define		k_time_type							k_vii[4]		// time stepping type
#define		k_open_boundary						k_vii[5]		// open boundary (YES / NO)
#define		k_flag_timestep_update				k_vii[6]		// flag for varying timestep update
#define		k_nb_cell_type						k_vii[7]		// neighbor cell type
#define		k_freq_filt							k_vii[8]		// filtering frequency
#define		k_freq_mass_sum						k_vii[9]		// mass summation frequency
#define		k_freq_temp							k_vii[10]		// temperature filtering frequency
#define		k_freq_output						k_vii[11]		// output frequency
#define		k_fp_solve							k_vii[12]		// solve pressure force ?
#define		k_fv_solve							k_vii[13]		// sovle viscous force ?
#define		k_fg_solve							k_vii[14]		// solve gravity force ?
#define		k_con_solve							k_vii[15]		// solve conduction?
#define		k_boussinesq_solve					k_vii[16]		// solve boussinesq approximation based natural convection?
#define		k_kgc_solve							k_vii[17]		// solve kernel gradient correction ?
#define		k_turbulence_model					k_vii[18]		// turbulence model	(by esk)
#define		k_concn_solve						k_vii[19]		// concentration diffusion model(PSH)
//
//psh:: ISPH input
#define		k_minIteration						k_vii[20]		// minimum number of PCISPH iteration
#define		k_maxIteration						k_vii[21]		// maximum number of PCISPH iteration
//
//solution variables
#define		k_nb_cell_number					k_vii[22]
#define		k_num_part							k_vii[23]
#define		k_NI								k_vii[24]
#define		k_NJ								k_vii[25]
#define		k_NK								k_vii[26]
#define		k_count								k_vii[27]
#define		k_num_cells							k_vii[28]
#define		k_ngpu								k_vii[29]
#define		k_calc_area							k_vii[30]
#define		k_num_p2p							k_vii[31]
#define		k_num_part2							k_vii[32]
#define		k_switch_ptype						k_vii[33]
#define   k_noslip_bc							k_vii[34]
#define   k_apr_solv                			k_vii[35]
#define  	k_h_change							k_vii[36]
#define   k_freq_apr							k_vii[37]
#define		k_freq_h_reset						k_vii[38]
#define		k_num_blocks						k_vii[39]
#define		freq_bapr							k_vii[40]
#define		k_NBI								k_vii[41]
#define		k_NBJ								k_vii[42]
#define		k_NBK								k_vii[43]
#define   	k_structure_solv          			k_vii[44]
#define   	k_atmos_stab		        		k_vii[45]	 	// AeroSPHere (KDH)
#define   	k_terrain_type	          			k_vii[46]  		// AeroSPHere (KDH)
#define 	k_num_part_LDM						k_vii[47]
#define		k_freq_LDM   						k_vii[48]  		// LDM frequency (KDH)
//
// Real type variable define
#define		k_kappa								k_vif[0]		// k in k*h
#define		k_p_ref								k_vif[1]		// reference pressure(for EOS)
#define		k_dt								k_vif[2]		// time-step(s)
#define		k_time								k_vif[3]		// time(s)
#define		k_time_end							k_vif[4]
// margin for simulation range
#define		k_u_limit							k_vif[5]

//psh:: ISPH input
#define		k_drho_th							k_vif[6]		// density convergence criterion
#define		k_dp_th								k_vif[7]		// pressure convergence criterion
#define		k_p_relaxation						k_vif[8]		// relaxation factor for PCISPH pressure
//
//solution variables
#define		k_x_min								k_vif[9]
#define		k_x_max								k_vif[10]
#define		k_y_min								k_vif[11]
#define		k_y_max								k_vif[12]
#define		k_z_min								k_vif[13]
#define		k_z_max								k_vif[14]
#define		k_nd_ref							k_vif[15]
//
#define		k_search_incr_factor				k_vif[16]
#define		k_search_kappa						k_vif[17]
#define		k_dcell								k_vif[18]
#define		k_buffersize						k_vif[19]
#define		k_dt_structure						k_vif[20]

////////////////////////////////////////////////////////////////////////
void read_solv_input(int_t*vii,Real*vif,const char*FileName)
{
	solver_type=Wcsph;

	fp_solve=0;
	fv_solve=0;
	fg_solve=0;
	con_solve=0;
	boussinesq_solve=0;
	kgc_solve=0;
	turbulence_model=0;
	concn_solve=0;
	apr_solv=0;
	structure_solv=0;
	atmos_stab=0;
	terrain_type=0;

	char inputString[1000];

	//inFile.open("../Result/output.txt");
	FILE*fd;
	fd=fopen(FileName,"r");

	int end;

	while(1){
		end=fscanf(fd,"%s",&inputString);		// reading one data from cc
		if(strcmp(inputString,"solver_type")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"WCSPH")==0) solver_type=Wcsph;
			if(strcmp(inputString,"ISPH")==0) solver_type=Isph;
		}
		if(strcmp(inputString,"dimension(1/2/3):")==0){
			fscanf(fd,"%s",&inputString);
			dim=atoi(inputString);
		}
		if(strcmp(inputString,"property_table(YES/NO):")==0){
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) prop_table=1;
			if(strcmp(inputString,"NO")==0) prop_table=0;
		}
		if(strcmp(inputString,"atmospheric-stability")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Unstable")==0) atmos_stab=0;
			if(strcmp(inputString,"Neutral")==0) atmos_stab=1;
			if(strcmp(inputString,"Stable")==0) atmos_stab=2;
		}
		if(strcmp(inputString,"terrain-type")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Flat")==0) terrain_type=0;
			if(strcmp(inputString,"Grass")==0) terrain_type=1;
			if(strcmp(inputString,"Forest")==0) terrain_type=2;
			if(strcmp(inputString,"City")==0) terrain_type=3;
			if(strcmp(inputString,"Ocean")==0) terrain_type=4;
		}
		if(strcmp(inputString,"kernel")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Quartic")==0) kernel_type=Quartic;
			if(strcmp(inputString,"Gaussian")==0)kernel_type=Gaussian;
			if(strcmp(inputString,"Quintic")==0) kernel_type=Quintic;
			if(strcmp(inputString,"Wendland2")==0) kernel_type=Wendland2;
			if(strcmp(inputString,"Wendland4")==0) kernel_type=Wendland4;
			if(strcmp(inputString,"Wendland6")==0) kernel_type=Wendland6;
		}
		if(strcmp(inputString,"time")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Euler")==0) time_type=Euler;
			if(strcmp(inputString,"Predictor_Corrector")==0) time_type=Pre_Cor;
		}
		if(strcmp(inputString,"structure")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) structure_solv=1;
			if(strcmp(inputString,"NO")==0) structure_solv=0;
		}
		if(strcmp(inputString,"open")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"OPEN")==0) open_boundary=1;
			if(strcmp(inputString,"PERIODIC")==0) open_boundary=2;
			if(strcmp(inputString,"NO")==0) open_boundary=0;
		}
		if(strcmp(inputString,"pressure-force")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) fp_solve=1;
			if(strcmp(inputString,"NO")==0) fp_solve=0;
		}
		if(strcmp(inputString,"viscous-force")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"NO")==0) fv_solve=0;
			if(strcmp(inputString,"Cleary")==0) fv_solve=1;
			if(strcmp(inputString,"Monaghan")==0) fv_solve=2;
			if(strcmp(inputString,"Violeau")==0) fv_solve=3;
		}
		if(strcmp(inputString,"turbulence-model")==0){		//(by esk)
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Laminar")==0) turbulence_model=0;
			if(strcmp(inputString,"k-lm")==0) turbulence_model=1;
			if(strcmp(inputString,"k-e")==0) turbulence_model=2;
			if(strcmp(inputString,"HB")==0) turbulence_model=3;
			if(strcmp(inputString,"SSM")==0) turbulence_model=4;
			if(strcmp(inputString,"DDF")==0) turbulence_model=5;
		}
		if(strcmp(inputString,"gravitational-force")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) fg_solve=1;
			if(strcmp(inputString,"NO")==0) fg_solve=0;
		}
		if(strcmp(inputString,"Conduction(YES/NO):")==0){
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) con_solve=1;
			if(strcmp(inputString,"NO")==0) con_solve=0;
		}
		if(strcmp(inputString,"Boussinesq-natural-convection")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) boussinesq_solve=1;
			if(strcmp(inputString,"NO")==0) boussinesq_solve=0;
		}
		if(strcmp(inputString,"Concentration-diffusion(YES/NO):")==0){
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) concn_solve=1;
			if(strcmp(inputString,"NO")==0) concn_solve=0;
		}
		if(strcmp(inputString,"kernel-gradient-correction")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"NO")==0) kgc_solve=0;
			if(strcmp(inputString,"KGC")==0) kgc_solve=1;
			if(strcmp(inputString,"FPM")==0) kgc_solve=2;
			if(strcmp(inputString,"DFPM")==0) kgc_solve=3;
			if(strcmp(inputString,"KGF")==0) kgc_solve=4;
		}
		if(strcmp(inputString,"switch_ptype")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) switch_ptype=1;
			if(strcmp(inputString,"NO")==0) switch_ptype=0;
		}
		if(strcmp(inputString,"noslip_boundary")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) noslip_bc=1;
			if(strcmp(inputString,"NO")==0) noslip_bc=0;
		}
		if(strcmp(inputString,"APR_solve")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Cell")==0) apr_solv=1;
			if(strcmp(inputString,"Switch")==0) apr_solv=2;
			if(strcmp(inputString,"Block")==0) apr_solv=3;
			if(strcmp(inputString,"NO")==0) apr_solv=0;
		}
		if(strcmp(inputString,"Block_configuration")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			NBI=atoi(inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			NBJ=atoi(inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			NBK=atoi(inputString);
		}
		if(strcmp(inputString,"Variable")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"Abrupt")==0) h_change=1;
			if(strcmp(inputString,"Smooth")==0) h_change=2;
			if(strcmp(inputString,"NO")==0) h_change=0;
		}
		if(strcmp(inputString,"Buffersize:")==0){
			fscanf(fd,"%s",&inputString);
			buffersize=atof(inputString);
		}
		if(strcmp(inputString,"reference-pressure")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			p_ref=atof(inputString);
		}
		if(strcmp(inputString,"kappa:")==0){
			fscanf(fd,"%s",&inputString);
			kappa=atof(inputString);
		}
		if(strcmp(inputString,"minimum-iteration:")==0){
			fscanf(fd,"%s",&inputString);
			minIteration=atoi(inputString);
		}
		if(strcmp(inputString,"maximum-iteration:")==0){
			fscanf(fd,"%s",&inputString);
			maxIteration=atoi(inputString);
		}
		if(strcmp(inputString,"pressure-convergence")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			dp_th=atof(inputString);
		}
		if(strcmp(inputString,"density-convergence")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			drho_th=atof(inputString);
		}
		if(strcmp(inputString,"pressure-relaxation")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			p_relaxation=atof(inputString);
		}
		if(strcmp(inputString,"time-step")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			dt=atof(inputString);
		}
		if(strcmp(inputString,"start-time")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			time=atof(inputString);
		}
		if(strcmp(inputString,"end-time")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			time_end=atof(inputString);
		}
		if(strcmp(inputString,"neighbor")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"3X3")==0) nb_cell_type=0;
			if(strcmp(inputString,"5X5")==0) nb_cell_type=1;
		}
		if(strcmp(inputString,"timestep")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			if(strcmp(inputString,"YES")==0) flag_timestep_update=1;
			if(strcmp(inputString,"NO")==0) flag_timestep_update=0;
		}
		if(strcmp(inputString,"structure-timestep")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			dt_structure=atof(inputString);
		}
		if(strcmp(inputString,"filtering")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			freq_filt=atoi(inputString);
		}
		if(strcmp(inputString,"temperature-filtering")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			freq_temp=atoi(inputString);
		}
		if(strcmp(inputString,"plot-output")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			freq_output=atoi(inputString);
		}
		if(strcmp(inputString,"LDM-frequency")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			freq_LDM=atoi(inputString);
		}
		if(strcmp(inputString,"APR_frequency:")==0){
			fscanf(fd,"%s",&inputString);
			freq_apr=atoi(inputString);
		}
		if(strcmp(inputString,"h-reset")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			freq_h_reset=atoi(inputString);
		}
		if(strcmp(inputString,"BAPR_frequency:")==0){
			fscanf(fd,"%s",&inputString);
			freq_bapr=atoi(inputString);
		}
		if(strcmp(inputString,"plot-variables:")==0){
			fscanf(fd,"%s",&inputString);
			num_plot_data=atoi(inputString);
			for(int ccount=0;ccount<num_plot_data;ccount++) {
				fscanf(fd,"%s",plot_data[ccount]);
			}
		}
		if(strcmp(inputString,"velocity-limit")==0){
			fscanf(fd,"%s",&inputString);
			fscanf(fd,"%s",&inputString);
			u_limit=atof(inputString);
		}
		if(end==-1) break;
	}
	fclose(fd);
}
////////////////////////////////////////////////////////////////////////
// function calculating number of particles from input file
int_t gpu_count_particle_numbers(const char*FileName)
{
	int_t idx=0;
	int_t tmp,end;

	FILE*inFile;
	inFile=fopen(FileName,"r");

	while(1){
		end=fscanf(inFile,"%d\n",&tmp);
		if(end==-1) break;
		idx+=1;
	}
	fclose(inFile);
	return idx;
}
////////////////////////////////////////////////////////////////////////
// function calculating number of particles from input file
int_t gpu_count_particle_numbers2(const char*FileName)
{
	int_t idx=0;
	int_t nop;
	char buffer[1024];

	FILE*inFile;

	inFile=fopen(FileName,"r");
	while(fgets(buffer,1024,inFile)!=NULL) idx+=1;
	fclose(inFile);

	nop=idx-1;

	return nop;
}
////////////////////////////////////////////////////////////////////////
// function calculating number of boundary particles from input file
int_t gpu_count_boundary_numbers(const char*FileName)
{
	FILE*fd;
	fd=fopen("./input/p_type.txt","r");
	int_t end;
	int_t nb=0;		// number of boundary particles

	int tmp;

	// calculation of number of boundary particles
	while(1){
		end=fscanf(fd,"%d\n",&tmp);
		if(end==-1) break;
		else if(tmp==0) nb+=1;
	}
	nb=nb;				// actual number of bondary particles

	return nb;
}
/////////////////////////////////////////////////////////////
#define inp_x   1
#define inp_y   2
#define inp_z   3
#define inp_ux  4
#define inp_uy  5
#define inp_uz  6
#define inp_m   7
#define inp_ptype 8
#define inp_h 9
#define inp_temp  10
#define inp_pres  11
#define inp_rho 12
#define inp_rhoref  13
#define inp_ftotal  14
#define inp_concn 15
#define inp_cc 16
#define inp_vist  17
#define inp_ct_boundary 18
#define inp_hf_boundary 19
#define inp_lbl_surf  20
#define inp_drho  21
#define inp_denthalpy 22
#define inp_dconcn  23
#define inp_dk  24
#define inp_de  25
#define inp_h0  26
#define inp_m0  27
///////////////////////////////////////////////////////////////////////
void read_input(part1*Pa1,char*INPUT)
{
  char FileName[256];
	//strcpy(FileName,"./input/CaseH_n5.txt");
	strcpy(FileName,INPUT);
  char buffer[1024];
  char *tok;    //token of string

  int j,end,tmp,nov,nop;   //number of data, number of variables, number of partices
  int lbl_var[100];
  nov=nop=0;

  FILE*inFile;
  inFile=fopen(FileName,"r");

  // count number of variables
  fgets(buffer,1024-1,inFile);				// read first line
  tok=strtok(buffer,"\t");						// line segmentation
  while(tok!=NULL){
	  tmp=atoi(tok);
	  lbl_var[nov]=tmp;
	  nov++;														// count number of segments(variables)
	  tok=strtok(NULL,"\t");
  }

  // read data
  while(1){
		for(j=0;j<nov;j++){
			end=fscanf(inFile,"%s\n",buffer);
			if(end==-1) break;
			switch(lbl_var[j]){
			case inp_x:
				Pa1[nop].x=atof(buffer);
				//Pa2[nop].x0=atof(buffer);
				break;
			case inp_y:
				Pa1[nop].y=atof(buffer);
				//Pa2[nop].y0=atof(buffer);
				break;
			case inp_z:
				Pa1[nop].z=atof(buffer);
				//Pa2[nop].z0=atof(buffer);
				break;
			case inp_ux:
				Pa1[nop].ux=atof(buffer);
				//Pa2[nop].ux0=atof(buffer);
				break;
			case inp_uy:
				Pa1[nop].uy=atof(buffer);
				//Pa2[nop].uy0=atof(buffer);
				break;
			case inp_uz:
				Pa1[nop].uz=atof(buffer);
				//Pa2[nop].uz0=atof(buffer);
				break;
			case inp_m:
				Pa1[nop].m=atof(buffer);
				break;
			case inp_ptype:
				Pa1[nop].p_type=atoi(buffer);
				break;
			case inp_h:
				Pa1[nop].h=atof(buffer);
				break;
			case inp_h0:
				Pa1[nop].h_ref=atof(buffer);
				break;
			case inp_m0:
				Pa1[nop].m_ref=atof(buffer);
				break;
			case inp_temp:
				Pa1[nop].temp=atof(buffer);
				break;
			case inp_pres:
				Pa1[nop].pres=atof(buffer);
				break;
			case inp_rho:
				Pa1[nop].rho=atof(buffer);
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_rhoref:
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_ftotal:
				//Pa3[nop].ftotal=atof(buffer);
				break;
			case inp_concn:
				Pa1[nop].concn=atof(buffer);
				//Pa2[nop].concn0=atof(buffer);
				break;
			case inp_cc:
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_ct_boundary:
				break;
			case inp_vist:
				//Pa3[nop].vis_t=atof(buffer);
				break;
			case inp_lbl_surf:
				break;
			case inp_drho:
				//Pa3[nop].drho=atof(buffer);
				break;
			case inp_denthalpy:
				//Pa3[nop].denthalpy=atof(buffer);
				break;
			case inp_dconcn:
				//Pa3[nop].dconcn=atof(buffer);
				break;
			case inp_dk:
				//Pa3[nop].dk_turb=atof(buffer);
				break;
			case inp_de:
				//Pa3[nop].de_turb=atof(buffer);
				break;
			default:
				printf("undefined variable name");
				break;
			}
		}
		if(end==-1) break;
		Pa1[nop].i_type=1;
		nop++;
  }
  fclose(inFile);
  printf("Input Files have been sucessfully read!!\n");
}
void read_input_LDM(L_part1*Pa1,char*INPUT)
{
  char FileName[256];
	//strcpy(FileName,"./input/CaseH_n5.txt");
	strcpy(FileName,INPUT);
  char buffer[1024];
  char *tok;    //token of string

  int j,end,tmp,nov,nop;   //number of data, number of variables, number of partices
  int lbl_var[100];
  nov=nop=0;

  FILE*inFile;
  inFile=fopen(FileName,"r");

  // count number of variables
  fgets(buffer,1024-1,inFile);				// read first line
  tok=strtok(buffer,"\t");						// line segmentation
  while(tok!=NULL){
	  tmp=atoi(tok);
	  lbl_var[nov]=tmp;
	  nov++;														// count number of segments(variables)
	  tok=strtok(NULL,"\t");
  }

  // read data
  while(1){
		for(j=0;j<nov;j++){
			end=fscanf(inFile,"%s\n",buffer);
			if(end==-1) break;
			switch(lbl_var[j]){
			case inp_x:
				Pa1[nop].x=atof(buffer);
				//Pa2[nop].x0=atof(buffer);
				break;
			case inp_y:
				Pa1[nop].y=atof(buffer);
				//Pa2[nop].y0=atof(buffer);
				break;
			case inp_z:
				Pa1[nop].z=atof(buffer);
				//Pa2[nop].z0=atof(buffer);
				break;
			case inp_ux:
				Pa1[nop].ux=atof(buffer);
				//Pa2[nop].ux0=atof(buffer);
				break;
			case inp_uy:
				Pa1[nop].uy=atof(buffer);
				//Pa2[nop].uy0=atof(buffer);
				break;
			case inp_uz:
				Pa1[nop].uz=atof(buffer);
				//Pa2[nop].uz0=atof(buffer);
				break;
			case inp_m:
				Pa1[nop].m=atof(buffer);
				break;
			case inp_ptype:
				Pa1[nop].p_type=atoi(buffer);
				break;
			case inp_h:
				Pa1[nop].h=atof(buffer);
				break;
			case inp_h0:
				Pa1[nop].h_ref=atof(buffer);
				break;
			case inp_m0:
				Pa1[nop].m_ref=atof(buffer);
				break;
			case inp_temp:
				Pa1[nop].temp=atof(buffer);
				break;
			case inp_pres:
				Pa1[nop].pres=atof(buffer);
				break;
			case inp_rho:
				Pa1[nop].rho=atof(buffer);
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_rhoref:
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_ftotal:
				//Pa3[nop].ftotal=atof(buffer);
				break;
			case inp_concn:
				Pa1[nop].concn=atof(buffer);
				//Pa2[nop].concn0=atof(buffer);
				break;
			case inp_cc:
				//Pa2[nop].rho0=atof(buffer);
				break;
			case inp_ct_boundary:
				break;
			case inp_vist:
				//Pa3[nop].vis_t=atof(buffer);
				break;
			case inp_lbl_surf:
				break;
			case inp_drho:
				//Pa3[nop].drho=atof(buffer);
				break;
			case inp_denthalpy:
				//Pa3[nop].denthalpy=atof(buffer);
				break;
			case inp_dconcn:
				//Pa3[nop].dconcn=atof(buffer);
				break;
			case inp_dk:
				//Pa3[nop].dk_turb=atof(buffer);
				break;
			case inp_de:
				//Pa3[nop].de_turb=atof(buffer);
				break;
			default:
				printf("undefined variable name");
				break;
			}
		}
		if(end==-1) break;
		Pa1[nop].i_type=1;
		nop++;
  }
  fclose(inFile);
  printf("Input Files have been sucessfully read!!\n");
}
////////////////////////////////////////////////////////////////////////
void find_minmax(int_t*vii,Real*vif,part1*Pa1)
{
	int_t i;
	int_t end=num_part-1;

	Real min_x=Pa1[0].x;	Real max_x=Pa1[0].x;
	Real min_y=Pa1[0].y;	Real max_y=Pa1[0].y;
	Real min_z=Pa1[0].z;	Real max_z=Pa1[0].z;
	//
	Real tmp_x,tmp_y,tmp_z;

	for(i=0;i<end;i++){
		tmp_x=Pa1[i].x;
		tmp_y=Pa1[i].y;
		tmp_z=Pa1[i].z;

		if(tmp_x<min_x) min_x=tmp_x;
		if(tmp_x>max_x) max_x=tmp_x;
		if(tmp_y<min_y) min_y=tmp_y;
		if(tmp_y>max_y) max_y=tmp_y;
		if(tmp_z<min_z) min_z=tmp_z;
		if(tmp_z>max_z) max_z=tmp_z;
	}

	x_min=min_x;
	x_max=max_x;
	y_min=min_y;
	y_max=max_y;
	z_min=min_z;
	z_max=max_z;
}


void test_array(Real* datax, int size)
{
	printf("Test Array ---------------- \n");

	for(int i=0;i<size;i++)
	{
		printf("datax[%d]=%f\n",i,datax[i]);
	}

	printf("/n");

}

void read_table(const char*FileName)
{
	char inputString[1000];
	float data;
	char ch[2000];
	char *ptr;

	int data_numbers=0;
	int count; count=0;
	int count2;
	int number_of_tables=0;
	int N_data[2000];

	//inFile.open("../Result/output.txt");
	FILE*fd;
	fd=fopen(FileName,"r");

	while(1)
	{

		// 한줄 읽기
		if(fgets(inputString,sizeof(inputString),fd)==NULL) break;
		//printf("%s\n",inputString);

		// 첫마디 읽기
		ptr=strtok(inputString, " ");
		//printf("%s\n",inputString);

		if (!strncmp(ptr,"#p_type",7)) {
		 		ptr=strtok(NULL," ");
				int number=atoi(ptr);
				//printf("%d\n",number);

				number_of_tables++;

				while(1)
				{
					fgets(inputString,sizeof(inputString),fd);
					//printf("%s\n",inputString);

					if (!strncmp(inputString,"#end",4)) break;

					ptr=strtok(inputString, " ");
					//printf("%s\n",ptr);

					if (!strncmp(ptr,"//T",3)) {
						//printf("OK!\n");
						fgets(inputString,sizeof(inputString),fd);
						ptr=strtok(inputString, ",");
						Real data=atof(ptr);
						//printf("%f\n",data);
						host_Tab_T[data_numbers+count]=data;
						// printf("data[%d]=%f\n",data_numbers+count,data);
						count++;
						while(ptr!=NULL)
						{
							ptr=strtok(NULL, ",");
							if(ptr==NULL) break;
							data=atof(ptr);
							host_Tab_T[data_numbers+count]=data;
							// printf("data[%d]=%f\n",data_numbers+count,data);
							count++;
						}
						count2=count;
						count=0;
						//fgets(inputString,sizeof(inputString),fd);
						//ptr=inputString;
					}
					else if (!strncmp(ptr,"//h",3)) {
						//printf("OK!\n");
						fgets(inputString,sizeof(inputString),fd);
						ptr=strtok(inputString, ",");
						Real data=atof(ptr);
						//printf("%f\n",data);
						host_Tab_h[data_numbers+count]=data;
						// printf("data[%d]=%f\n",data_numbers+count,data);
						count++;
						while(ptr!=NULL)
						{
							ptr=strtok(NULL, ",");
							if(ptr==NULL) break;
							data=atof(ptr);
							host_Tab_h[data_numbers+count]=data;
							// printf("data[%d]=%f\n",data_numbers+count,data);
							count++;
						}
						count2=count;
						count=0;
						//fgets(inputString,sizeof(inputString),fd);
						//ptr=inputString;
					}
					else	if (!strncmp(ptr,"//k",3)) {
						//printf("OK!\n");
						fgets(inputString,sizeof(inputString),fd);
						ptr=strtok(inputString, ",");
						Real data=atof(ptr);
						//printf("%f\n",data);
						host_Tab_k[data_numbers+count]=data;
						// printf("data[%d]=%f\n",data_numbers+count,data);
						count++;
						while(ptr!=NULL)
						{
							ptr=strtok(NULL, ",");
							if(ptr==NULL) break;
							data=atof(ptr);
							host_Tab_k[data_numbers+count]=data;
							// printf("data[%d]=%f\n",data_numbers+count,data);
							count++;
						}
						count2=count;
						count=0;
						//fgets(inputString,sizeof(inputString),fd);
						//ptr=inputString;
					}
					else	if (!strncmp(ptr,"//cp",3)) {
						//printf("OK!\n");
						fgets(inputString,sizeof(inputString),fd);
						ptr=strtok(inputString, ",");
						Real data=atof(ptr);
						//printf("%f\n",data);
						host_Tab_cp[data_numbers+count]=data;
						// printf("data[%d]=%f\n",data_numbers+count,data);
						count++;
						while(ptr!=NULL)
						{
							ptr=strtok(NULL, ",");
							if(ptr==NULL) break;
							data=atof(ptr);
							host_Tab_cp[data_numbers+count]=data;
							// printf("data[%d]=%f\n",data_numbers+count,data);
							count++;
						}
						count2=count;
						count=0;
						//fgets(inputString,sizeof(inputString),fd);
						//ptr=inputString;
					}
					else	if (!strncmp(ptr,"//mu",3)) {
						//printf("OK!\n");
						fgets(inputString,sizeof(inputString),fd);
						ptr=strtok(inputString, ",");
						Real data=atof(ptr);
						//printf("%f\n",data);
						host_Tab_vis[data_numbers+count]=data;
						// printf("data[%d]=%f\n",data_numbers+count,data);
						count++;
						while(ptr!=NULL)
						{
							ptr=strtok(NULL, ",");
							if(ptr==NULL) break;
							data=atof(ptr);
							host_Tab_vis[data_numbers+count]=data;
							// printf("data[%d]=%f\n",data_numbers+count,data);
							count++;
						}
						count2=count;
						count=0;
						//fgets(inputString,sizeof(inputString),fd);
						//ptr=inputString;
					}
				}

				data_numbers=data_numbers+count2;
				// printf("data_numbers=%d\n",data_numbers);
				N_data[number]=data_numbers;
		}
	}

		fclose(fd);

		// 추가계산
		// int Table_size[10];
		// int Table_index[10];

		host_table_index[0]=0;
		host_table_size[0]=N_data[0];

		for (int i=1;i<number_of_tables;i++)
		{
			host_table_index[i]=N_data[i-1];
			host_table_size[i]=N_data[i]-N_data[i-1];
		}


		// print for check
		printf("\nsaved results-----------\n\n");
		printf("number of tables=%d\n\n",number_of_tables);
		printf("data numbers=%d\n\n",data_numbers);

		for (int i=0;i<data_numbers;i++)
		{
			printf("T[%d]=%f\n",i,host_Tab_T[i]);
		}

		printf("\n");

		for (int i=0;i<data_numbers;i++)
		{
			printf("h[%d]=%f\n",i,host_Tab_h[i]);
		}

		printf("\n");

		for (int i=0;i<data_numbers;i++)
		{
			printf("k[%d]=%f\n",i,host_Tab_k[i]);
		}

		printf("\n");

		for (int i=0;i<data_numbers;i++)
		{
			printf("cp[%d]=%f\n",i,host_Tab_cp[i]);
		}

		printf("\n");

		for (int i=0;i<data_numbers;i++)
		{
			printf("mu[%d]=%f\n",i,host_Tab_vis[i]);
		}

		printf("\n");

		for (int i=0;i<number_of_tables;i++)
		{
			printf("N_data[%d]=%d\n",i,N_data[i]);
		}

		printf("\n");

		for (int i=0;i<number_of_tables;i++)
		{
			printf("Table_size[%d]=%d\n",i,host_table_size[i]);
		}

		printf("\n");

		for (int i=0;i<number_of_tables;i++)
		{
			printf("Table_index[%d]=%d\n",i,host_table_index[i]);
		}

		printf("\n");

		//test_array(&host_Tab_T[5],5);
}
