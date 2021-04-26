#include <stdio.h>
#include <math.h>
#include "system.h"

void Jlist(double t){
	int i,X,Y;
	for(i = 0; i < NA; i++){
		if(position[i].z >= Z1 && position[i].z < Z2){
			if(position[i].z >= (Z1 + Z2)/2 && position[i].z - velocity[i].z*t < (Z1 + Z2)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JA[X + mx*Y] += 4/(M_PI*SQR(enter_pore[i].d));
			}
			else if(position[i].z < (Z1 + Z2)/2 && position[i].z - velocity[i].z*t >= (Z1 + Z2)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JA[X + mx*Y] -= 4/(M_PI*SQR(enter_pore[i].d));
			}
		}
		else if(position[i].z >= Z3 && position[i].z < Lz){
			if(position[i].z < (Z3 + Lz)/2 && position[i].z - velocity[i].z*t >= (Z3 + Lz)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JA[M + X + mx*Y] += 4/(M_PI*SQR(enter_pore[i].d));
			}
			else if(position[i].z >= (Z3 + Lz)/2 && position[i].z - velocity[i].z*t < (Z3 + Lz)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JA[M + X + mx*Y] -= 4/(M_PI*SQR(enter_pore[i].d));
			}
		}
	}
	for(i = NA; i < NumberOfParticles; i++){
		if(position[i].z >= Z1 && position[i].z < Z2){
			if(position[i].z >= (Z1 + Z2)/2 && position[i].z - velocity[i].z*t < (Z1 + Z2)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JB[X + mx*Y] += 4/(M_PI*SQR(enter_pore[i].d));
			}
			else if(position[i].z < (Z1 + Z2)/2 && position[i].z - velocity[i].z*t >= (Z1 + Z2)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JB[X + mx*Y] -= 4/(M_PI*SQR(enter_pore[i].d));
			}
		}
		else if(position[i].z >= Z3 && position[i].z < Lz){
			if(position[i].z < (Z3 + Lz)/2 && position[i].z - velocity[i].z*t >= (Z3 + Lz)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JB[M + X + mx*Y] += 4/(M_PI*SQR(enter_pore[i].d));
			}
			else if(position[i].z >= (Z3 + Lz)/2 && position[i].z - velocity[i].z*t < (Z3 + Lz)/2){
				X = (int) (enter_pore[i].x*mx/Lx - 0.5);
				Y = (int) (enter_pore[i].y*my/Ly - 0.5);
				JB[M + X + mx*Y] -= 4/(M_PI*SQR(enter_pore[i].d));
			}
		}
	}
	return;
}

void PrintJlist(void){
	char filename[20];
	FILE *fpJlist;
	
	sprintf(filename, "JAlist.txt");
	fpJlist = fopen(filename, "w");
	for(int i = 0; i < 2*M; i++){
		fprintf(fpJlist, "%.8lf\n", JA[i]);
	}
	fclose(fpJlist);
	
	sprintf(filename, "JBlist.txt");
	fpJlist = fopen(filename, "w");
	for(int i = 0; i < 2*M; i++){
		fprintf(fpJlist, "%.8lf\n", JB[i]);
	}
	fclose(fpJlist);
	return;
}
