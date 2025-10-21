double cclssw[256],cclgsw[256];  // COtimers
////////////////////////////////////////////////////////////////////////
double setsw(int seq)
{
	struct timespec AA;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&AA);
	cclssw[seq]=(double)AA.tv_sec+(double)AA.tv_nsec/1000000000.0;
	return(cclssw[seq]);
}
////////////////////////////////////////////////////////////////////////
double getsw(int seq)
{
	struct timespec AA;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&AA);
	cclgsw[seq]=(double)AA.tv_sec+(double)AA.tv_nsec/1000000000.0-cclssw[seq];
	return(cclgsw[seq]);
}
////////////////////////////////////////////////////////////////////////
__global__ void kernel_copy_max(part1*P1,part2*P2,part3*P3,Real*mrho,Real*mft,Real*mu)
{
	uint_t i=threadIdx.x+blockIdx.x*blockDim.x;
	if(i>=k_num_part2) return;
	if(P1[i].p_type==0||P1[i].i_type>=i_type_crt){
		mu[i]=0;
		mrho[i]=0;
		mft[i]=0;
		return;
	}
	mu[i]=sqrt(P1[i].ux*P1[i].ux+P1[i].uy*P1[i].uy+P1[i].uz*P1[i].uz);
	// mrho[i]=P1[i].rho;
	mrho[i]=(P1[i].rho-P2[i].rho_ref)/P2[i].rho_ref*100.0;
	mft[i]=P3[i].ftotal;
}
////////////////////////////////////////////////////////////////////////
float FloatSwap( float f )
{
   union
   {
      float f;
      unsigned char b[4];
      //unsigned char b[8];
   } dat1,dat2;

   dat1.f=f;
   dat2.b[0]=dat1.b[3];
   dat2.b[1]=dat1.b[2];
   dat2.b[2]=dat1.b[1];
   dat2.b[3]=dat1.b[0];
	 /*
   dat2.b[0]=dat1.b[7];
   dat2.b[1]=dat1.b[6];
   dat2.b[2]=dat1.b[5];
   dat2.b[3]=dat1.b[4];
   dat2.b[4]=dat1.b[3];
   dat2.b[5]=dat1.b[2];
   dat2.b[6]=dat1.b[1];
   dat2.b[7]=dat1.b[0];
	 //*/

   return dat2.f;
}

int IntSwap( int d )
{
   union
   {
      int d;
      unsigned char b[4];
      //unsigned char b[8];
   } dat1,dat2;

   dat1.d=d;
   dat2.b[0]=dat1.b[3];
   dat2.b[1]=dat1.b[2];
   dat2.b[2]=dat1.b[1];
   dat2.b[3]=dat1.b[0];
	 /*
   dat2.b[0]=dat1.b[7];
   dat2.b[1]=dat1.b[6];
   dat2.b[2]=dat1.b[5];
   dat2.b[3]=dat1.b[4];
   dat2.b[4]=dat1.b[3];
   dat2.b[5]=dat1.b[2];
   dat2.b[6]=dat1.b[1];
   dat2.b[7]=dat1.b[0];
	 //*/
   return dat2.d;
}
////////////////////////////////////////////////////////////////////////
void save_computing_time(Real*computing_time)
{
	int_t i;
	char FileName[256];
	sprintf(FileName,"./plotdata/time/computing_time_%dstp.txt",count);
	FILE*outFile;
	outFile=fopen(FileName,"w");

	fprintf(outFile,"Total computing time at %d stp is %f s\n\n",count,computing_time[0]);
	fprintf(outFile,"Computing time before NNPS : %f s\n",computing_time[1]);
	fprintf(outFile,"Computing time for NNPS : %f s\n",computing_time[2]);
	fprintf(outFile,"Computing time for PREP : %f s\n",computing_time[3]);
	fprintf(outFile,"Computing time for MASS : %f s\n",computing_time[4]);
	fprintf(outFile,"Computing time for EOS : %f s\n",computing_time[5]);
	fprintf(outFile,"Computing time for INTERACTION : %f s\n",computing_time[6]);
	fprintf(outFile,"Computing time for Time update : %f s\n\n",computing_time[7]);
	fprintf(outFile,"Computing time for PST : %f s\n\n",computing_time[8]);
	fprintf(outFile,"Computing time for APR : %f s\n\n",computing_time[9]);

	// fprintf(outFile,"Computing time until NNPS : %f s\n",computing_time[11]);
	// fprintf(outFile,"Computing time for Prep : %f s\n",computing_time[12]);
	// fprintf(outFile,"Computing time for BC : %f s\n",computing_time[13]);
	// fprintf(outFile,"Computing time for Density : %f s\n",computing_time[14]);
	// fprintf(outFile,"Computing time for EOS & Interaction : %f s\n",computing_time[15]);
	// fprintf(outFile,"Computing time for Time update : %f s\n",computing_time[16]);
	// fprintf(outFile,"Computing time for PST : %f s\n",computing_time[17]);

	fclose(outFile);
}
////////////////////////////////////////////////////////////////////////

void save_vtk_bin_single_flag(part1*P1,part2*P2,part3*P3)
{
	int_t i,nop;//,nob;
	nop=num_part2;
	// nob=number_of_boundaries;
	int_t Nparticle=0;							// number of fluid particles (x>0.00) for 3D PGSFR calculation
	int_t *plot_flag;
	plot_flag = (int*)malloc(sizeof(int)*nop);
	//for(i=0;i<nop;i++) if(P1[i].x>0) Nparticle++;
	//for(i=0;i<nop;i++) if(((P1[i].i_type==1)||(P1[i].i_type==2))&&P1[i].y>=0) {
	for(i=0;i<nop;i++) if((P1[i].i_type==1)||(P1[i].i_type==2)) {
		Nparticle++;
		plot_flag[i]=1;
	}
	printf("Number of Particles = %d\n\n",nop);

	float val;
	int valt;

	// Filename: It should be series of frame numbers(nameXXX.vtk) for the sake of auto-reading in PARAVIEW.
	char FileName_vtk[256];
	char DirName[256]="CaseH_n10_DDF(WL6,dt=2.5e-5,KGC,IBM_WM,u=3.0,b=0.3)";
	// sprintf(FileName_vtk,"./900X500X300(dx=5.0,dt=0.1,mu=1e-4)/fluid_%dstp.vtk",count);
	sprintf(FileName_vtk,"./%s/fluid_%dstp.vtk",DirName, count);
	printf("Saving plot to '%s', count=%d\n",DirName,count);
	// If the file already exists,its contents are discarded and create the new one.
	FILE*outFile_vtk;
	outFile_vtk=fopen(FileName_vtk,"w");

	fprintf(outFile_vtk,"# vtk DataFile Version 3.0\n");					// version & identifier: it must be shown.(ver 1.0/2.0/3.0)
	fprintf(outFile_vtk,"Print out results in vtk format\n");			// header: description of file,it never exceeds 256 characters
	fprintf(outFile_vtk,"BINARY\n");														// format of data (ACSII / BINARY)
	fprintf(outFile_vtk,"DATASET POLYDATA\n");										// define DATASET format: 'POLYDATA' is proper to represent SPH particles

	//Define SPH particles---------------------------------------------------------------
	fprintf(outFile_vtk,"POINTS\t%d\tfloat\n",Nparticle);					// define particles position as POINTS
	for(i=0;i<nop;i++){							// print out (x,y,z) coordinates of particles
		//if(P1[i].x>0){
		if(plot_flag[i]==1){
			val=FloatSwap(P1[i].x);
			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
			val=FloatSwap(P1[i].y);
			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
			val=FloatSwap(P1[i].z);
			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		}
	}

	fprintf(outFile_vtk,"POINT_DATA\t%d\n",Nparticle);

	fprintf(outFile_vtk,"FIELD FieldData\t%d\n",num_plot_data);

	for (int ccount=0;ccount<num_plot_data;ccount++)
	{
		char data_label[20];
		strcpy(data_label,plot_data[ccount]);

		// buffer_type
		if (!strncmp(data_label,"buffer_type",3)) {
			fprintf(outFile_vtk,"buffer_type\t1\t%d\tint\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				if(plot_flag[i]==1){
					valt=IntSwap(P1[i].buffer_type);
					fwrite((void*)&valt,sizeof(int),1,outFile_vtk);
				}
			}
		}

		// i_type
		if (!strncmp(data_label,"i_type",3)) {
			fprintf(outFile_vtk,"i_type\t1\t%d\tint\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				if(plot_flag[i]==1){
					valt=IntSwap(P1[i].i_type);
					fwrite((void*)&valt,sizeof(int),1,outFile_vtk);
				}
			}
		}

		// p_type
		if (!strncmp(data_label,"p_type",3)) {
			fprintf(outFile_vtk,"p_type\t1\t%d\tint\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				if(plot_flag[i]==1){
					valt=IntSwap(P1[i].p_type);
					fwrite((void*)&valt,sizeof(int),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"source",5)) {
			fprintf(outFile_vtk,"source\t1\t%d\tint\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				if(plot_flag[i]==1){
					valt=IntSwap(P1[i].source);
					fwrite((void*)&valt,sizeof(int),1,outFile_vtk);
				}
			}
		}

		// density
		if (!strncmp(data_label,"rho",3)) {
			fprintf(outFile_vtk,"rho\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].rho);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"vol",3)) {
			fprintf(outFile_vtk,"vol\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].vol);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"velocity",4)) {
			fprintf(outFile_vtk,"velocity\t3\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].ux);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].uy);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].uz);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"normal",5)) {
			fprintf(outFile_vtk,"normal\t3\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P3[i].nx);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P3[i].ny);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P3[i].nz);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// if (!strncmp(data_label,"v_force",5)) {
		// 	fprintf(outFile_vtk,"fv\t3\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].fv_x);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fv_y);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fv_z);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }
		//
		// if (!strncmp(data_label,"p_force",5)) {
		// 	fprintf(outFile_vtk,"fp\t3\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].fp_x);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fp_y);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fp_z);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }
		//
		// if (!strncmp(data_label,"e_force",5)) {
		// 	fprintf(outFile_vtk,"fe\t3\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].fe_x);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fe_y);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].fe_z);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }
		//
		if (!strncmp(data_label,"b_force",5)) {
			fprintf(outFile_vtk,"fb\t3\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].fbx);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].fby);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].fbz);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"MOST_force",7)) {
			fprintf(outFile_vtk,"fm\t2\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].fmx);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].fmy);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"MOST_heat",7)) {
			fprintf(outFile_vtk,"hm\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].hmz);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// if (!strncmp(data_label,"vel_f",5)) {
		// 	fprintf(outFile_vtk,"vel_f\t3\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].ux_f);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].uy_f);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 			val=FloatSwap(P1[i].uz_f);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }

		// pressure
		if (!strncmp(data_label,"pressure",3)) {
			fprintf(outFile_vtk,"pressure\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].pres);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					//printf("write pres=%f\n",P1[i].pres);

				}
			}
		}

		// vis_t
		if (!strncmp(data_label,"vis_t",5)) {
			fprintf(outFile_vtk,"vis_t\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].vis_t);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// temp
		if (!strncmp(data_label,"temp",3)) {
			fprintf(outFile_vtk,"temp\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].temp);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// flt_s
		if (!strncmp(data_label,"flt_s",3)) {
			fprintf(outFile_vtk,"flt_s\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].flt_s);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// concn
		if (!strncmp(data_label,"concn",5)) {
			fprintf(outFile_vtk,"concn\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].concn);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// concn
		if (!strncmp(data_label,"concentration",6)) {
			fprintf(outFile_vtk,"concentration\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].concentration);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		// concn
		if (!strncmp(data_label,"enthalpy",4)) {
			fprintf(outFile_vtk,"enthalpy\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].enthalpy);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"mass",4)) {
			fprintf(outFile_vtk,"mass\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].m);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"eli",3)) {
			fprintf(outFile_vtk,"eli\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].eli);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"smoothing",5)) {
			fprintf(outFile_vtk,"h\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].h);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"monitor",5)) {
			fprintf(outFile_vtk,"monitor\t3\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].float1);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].float2);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
					val=FloatSwap(P1[i].float3);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}
		if (!strncmp(data_label,"ncell",5)) {
			fprintf(outFile_vtk,"ncell\t1\t%d\tint\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					valt=IntSwap(P1[i].ncell);
					fwrite((void*)&valt,sizeof(int),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"rho_ref",5)) {
			fprintf(outFile_vtk,"rho_ref\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P2[i].rho_ref);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"k_turb",4)) {
			fprintf(outFile_vtk,"k_turb\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].k_turb);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}

		if (!strncmp(data_label,"e_turb",4)) {
			fprintf(outFile_vtk,"e_turb\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P1[i].e_turb);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}
		// if (!strncmp(data_label,"ppe1",4)) {
		// 	fprintf(outFile_vtk,"ppe1\t1\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].PPE1);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }
		// if (!strncmp(data_label,"ppe2",4)) {
		// 	fprintf(outFile_vtk,"ppe2\t1\t%d\tfloat\n",Nparticle);
		// 	for(i=0;i<nop;i++){
		// 		//if(P1[i].x>0){
		// 		// if((P1[i].p_type==2)|(P1[i].p_type==9)){
		// 		if(plot_flag[i]==1){
		// 			val=FloatSwap(P1[i].PPE2);
		// 			fwrite((void*)&val,sizeof(float),1,outFile_vtk);
		// 		}
		// 	}
		// }
		if (!strncmp(data_label,"dtemp",5)) {
			fprintf(outFile_vtk,"dtemp\t1\t%d\tfloat\n",Nparticle);
			for(i=0;i<nop;i++){
				//if(P1[i].x>0){
				// if((P1[i].p_type==2)|(P1[i].p_type==9)){
				if(plot_flag[i]==1){
					val=FloatSwap(P3[i].dtemp);
					fwrite((void*)&val,sizeof(float),1,outFile_vtk);
				}
			}
		}
	}


	fclose(outFile_vtk);
}
