/*******************************************************
* define the pore diameters of the membrane(s)
*******************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

FILE *fppore;
//double maxx, maxy;
//int factorial(int n){
//	return n < 2 ? 1 : factorial(n-1) * n;
//}
//double power(double x, int n){
//	return !n ? 1 : power(x, n-1) * x;
//}
//double Binomial(int i, int n, double p){
//	return factorial(n) / (factorial(n-i)*factorial(i)) * power(p, i) * power(1-p, n-i);
//}

void MakeMembrane(){
	Area1 = 0; Area2 = 0;
	double d;
	if(!MemType){
		fppore = fopen("pore_size.dat","w");
		for(int j=0; j<my; j++){
			for(int i=0; i<mx; i++){
				Area1 += atan(1)*SQR(D);//-sigmaA
				Area2 += atan(1)*SQR(D);//-sigmaB
				fprintf(fppore,"%lf\n", D);
	 		}
 		}
 		fprintf(fppore,"#%lf\t%lf\t%lf\n", atan(1), Area1, Area2);
 		fclose(fppore);
	}
	else{
		fppore = fopen("pore_size.dat","r");
		for(int i = 0; i < M; i++){
			fscanf(fppore,"%lf",&d);
			Area1 += atan(1)*SQR(d);//-sigmaA
			Area2 += atan(1)*SQR(d);//-sigmaB
 		}
 		fclose(fppore);
	}
	return;
}
