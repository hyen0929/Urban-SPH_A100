__device__ Real calc_kernel_wij_ipf(Real tH,Real rr){

	Real tR,wij_ipf,tA;
	tR=wij_ipf=0.0;
	tA=1.0;
	Real eps;
	eps=tH/3.5;

	tR=rr/tH*0.5;
	wij_ipf=(tR<1)*tA*(1-tR)*(1-tR)*(1-tR)*(1-tR)*(1+4*tR);

	return wij_ipf;
}

__device__ Real calc_kernel_wij_half(Real tH,Real rr){

	Real tR,wij_half,tA;
	tR=wij_half=0.0;
	tA=1.0;
	Real eps0;
	eps0=tH/3.5*0.4;

	tR=rr/tH*0.5;
	wij_half=((tR*2)<1)*tA*(1-(tR*2))*(1-(tR*2))*(1-(tR*2))*(1-(tR*2))*(1+4*(tR*2));

	return wij_half;
}
