/***************************************************
* collision dyanmics involving particle i and j
* i and j are in contact at the moment
* virial is calculated
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void MemCollision(int i){
	VECTOR rij, v_par;

	if((position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)){
		rij.x = sigma[i] * (enter_pore[i].x - position[i].x) / (enter_pore[i].d - sigma[i]);
		rij.y = sigma[i] * (enter_pore[i].y - position[i].y) / (enter_pore[i].d - sigma[i]);

		v_par.x = (rij.x*velocity[i].x + rij.y*velocity[i].y)/(SQR(rij.x) + SQR(rij.y)) * rij.x;
		v_par.y = (rij.x*velocity[i].x + rij.y*velocity[i].y)/(SQR(rij.x) + SQR(rij.y)) * rij.y;

		velocity[i].x -= 2*v_par.x;
		velocity[i].y -= 2*v_par.y;
	}
	else{		
		if(enter_pore[i].d == 0./* || (sigma[i] > enter_pore[i].d && sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) < tolerance*1000)*/){velocity[i].z *= -1;}
		else{
			rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
			rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
			
			if(MIN(MIN(fabs(position[i].z), fabs(position[i].z - Z1)), MIN(fabs(position[i].z - Z2), fabs(position[i].z - Z3))) == fabs(position[i].z)){rij.z = position[i].z;}
			else if(MIN(MIN(fabs(position[i].z), fabs(position[i].z - Z1)), MIN(fabs(position[i].z - Z2), fabs(position[i].z - Z3))) == fabs(position[i].z - Z1)){rij.z = position[i].z - Z1;}
			else if(MIN(MIN(fabs(position[i].z), fabs(position[i].z - Z1)), MIN(fabs(position[i].z - Z2), fabs(position[i].z - Z3))) == fabs(position[i].z - Z2)){rij.z = position[i].z - Z2;}
			else{rij.z = position[i].z - Z3;}
		
			v_par.x = (rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z) / (SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) * rij.x;
			v_par.y = (rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z) / (SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) * rij.y;
			v_par.z = (rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z) / (SQR(rij.x) + SQR(rij.y) + SQR(rij.z)) * rij.z;
			
			velocity[i].x -= 2*v_par.x;
			velocity[i].y -= 2*v_par.y;
			velocity[i].z -= 2*v_par.z;
		}
	}
	return;
}

void Collision(int i, int j) //
{
 VECTOR rij,vij,pij;
 double sigmaij2,bij;
 double mu;

 sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
 mu = mass[i]*mass[j]/(mass[i]+mass[j]);

 rij.x = position[i].x - position[j].x;
 rij.y = position[i].y - position[j].y;
 rij.z = position[i].z - position[j].z;
 //pbc
 MinimumImage(&(rij.x),Lx);
 MinimumImage(&(rij.y),Ly);
 MinimumImage(&(rij.z),Lz);
 
 vij.x = velocity[i].x - velocity[j].x;
 vij.y = velocity[i].y - velocity[j].y;
 vij.z = velocity[i].z - velocity[j].z;

 bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z;

 pij.x = -rij.x*bij/sigmaij2*2.*mu;
 pij.y = -rij.y*bij/sigmaij2*2.*mu;
 pij.z = -rij.z*bij/sigmaij2*2.*mu;

 velocity[i].x = velocity[i].x + pij.x/mass[i];
 velocity[i].y = velocity[i].y + pij.y/mass[i];
 velocity[i].z = velocity[i].z + pij.z/mass[i];
 
 velocity[j].x = velocity[j].x - pij.x/mass[j];
 velocity[j].y = velocity[j].y - pij.y/mass[j];
 velocity[j].z = velocity[j].z - pij.z/mass[j];

 Virial = pij.x*rij.x + pij.y*rij.y + pij.z*rij.z;
 if(position[i].z >= 0. && position[i].z < Z1 && position[j].z >= 0. && position[j].z < Z1){
 	VirialLeft = Virial; VirialRight = 0; VirialLeftMem = 0; VirialRightMem = 0;
 }
 else if(position[i].z >= Z1 && position[i].z < Z2 && position[j].z >= Z1 && position[j].z < Z2){
 	VirialLeft = 0; VirialRight = 0; VirialLeftMem = Virial; VirialRightMem = 0;
 }
 else if(position[i].z >= Z2 && position[i].z < Z3 && position[j].z >= Z2 && position[j].z < Z3){
 	VirialLeft = 0; VirialRight = Virial; VirialLeftMem = 0; VirialRightMem = 0;
 }
 else if(position[i].z >= Z3 && position[i].z < Lz && position[j].z >= Z3 && position[j].z < Lz){
 	VirialLeft = 0; VirialRight = 0; VirialLeftMem = 0; VirialRightMem = Virial;
 }
 else{VirialLeft = 0; VirialRight = 0; VirialLeftMem = 0; VirialRightMem = 0;}

 return;
}
