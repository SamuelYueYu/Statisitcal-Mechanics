/***************************************************
* check overlaps between every pair
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void MemOverlapCheck(void){
	bool overlap;
	double tor = 9.99*10.E-3;
	for(int i=0; i<NumberOfParticles; i++){
		if((position[i].z <= Z1 - sigma[i]/2 + tor && position[i].z >= sigma[i]/2 - tor) || (position[i].z <= Z3 - sigma[i]/2 + tor && position[i].z >= Z2 + sigma[i]/2 - tor)){overlap = false;}
		else if((position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)){
			overlap = sigma[i] < enter_pore[i].d && sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) <= tor + (enter_pore[i].d - sigma[i]) / 2 ? false : true;
			if(overlap){
				printf("coln = %d, overlap occurs for particle %d at (%lf,%lf,%lf,%lf) with pore (%lf,%lf,%lf,%lf), sqrt(SQR(position[i].x - enter_porex) + SQR(position[i].y - enter_porey)) - (enter_pored - sigma[i]) / 2 = %lf\n",coln,i,position[i].x,position[i].y,position[i].z,sigma[i], enter_pore[i].x, enter_pore[i].y, enter_pore[i].z, enter_pore[i].d, sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) - (enter_pore[i].d - sigma[i]) / 2); exit(1);
				return;
			}
		}
		else{
			if(enter_pore[i].d == 0.){overlap = (position[i].z > Z1 - sigma[i]/2 && position[i].z < Z1) || (position[i].z > Z2 && position[i].z < Z2 + sigma[i]/2) || (position[i].z > Z3 - sigma[i]/2 && position[i].z < Z3) || position[i].z < sigma[i]/2 ? true : false;}
			
			else{
				overlap = sigma[i]/2 - sqrt(SQR(position[i].z - enter_pore[i].z) + SQR(enter_pore[i].d/2 - sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)))) > tor ? true : false;
			}
		}
		if(overlap){printf("overlap occurs for particle %d at coln = %d, (%lf,%lf,%lf,%lf) with pore (%lf,%lf,%lf,%lf) and velocity (%lf,%lf,%lf)\n",i,coln,position[i].x,position[i].y,position[i].z,sigma[i],enter_pore[i].x,enter_pore[i].y,enter_pore[i].z, enter_pore[i].d, velocity[i].x, velocity[i].y, velocity[i].z); exit(1);}	
	}
	return;
}

void OverlapCheck(void) //
{
 int i,j;
 double rij2;
 double sigmaij2;
 double tor;

 tor = 10.E-4;

 for(i=0;i<NumberOfParticles-1;i++)
 for(j=i+1;j<NumberOfParticles;j++)
 {
  rij2 = Distance(i,j);
  sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));

  if(rij2 - sigmaij2 < -tor)
  {
   printf("overlap occurs between %d and %d, rij = %.30lf\tsigmaij = %.30lf\n",i,j,sqrt(rij2),sqrt(sigmaij2));
   exit(1);
  }
 }

 return;
}
