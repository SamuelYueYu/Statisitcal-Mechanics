#include <stdio.h>
#include <math.h>
#include "system.h"

void VanishInfo(int i){
	VECTOR rij;
 	double bij; // rij*vij
 	double rij2,vij2; // square
 	double discr; // discriminate b^2 - 4ac
 	double tvanish;
	VanishTime[i] = TimeBig;
	
	rij.z = position[i].z - (Z2 + Z3)/2;
 	bij = rij.z*velocity[i].z;

  	if(bij < 0.){
   		rij2 = SQR(rij.z);
   		vij2 = SQR(velocity[i].z); 
   		discr = SQR(bij) - vij2*(rij2-SQR((Z3 - Z2) / 2));
		
   		if(discr > 0.){
			tvanish = (-bij - sqrt(discr))/vij2;
			
    			if(tvanish < VanishTime[i]){
    				VanishTime[i] = tvanish;
   			}
   		} // endif discriminant > 0
  	} //endif bij < 0
	return;
}

void VanishUpdate(int i){
	int ii,j,k;
	int cello1, cello2, cello3, celln;
	
	if(identity[i] == 1){
		j = NA - 1; k = NumberOfParticles-1;
		
		for(int l = 0; l < MAX(N_0_1, N_Z2_1); l++){
			if(n_0_1[l] == j){
				n_0_1[l] = i; break;
			}
		}
		for(int l = 0; l < MAX(N_0_2, N_Z2_2); l++){
			if(n_0_2[l] == k){
				n_0_2[l] = j; break;
			}
		}
			
		Mcom -= mass[i];
   		Vcom.x -= mass[i]*velocity[i].x;
   		Vcom.y -= mass[i]*velocity[i].y;
   		Vcom.z -= mass[i]*velocity[i].z;
   			
   		NumberOfParticles--; NA--;
   		fA = 1.*NA/NumberOfParticles;
 		fB = 1.*NB/NumberOfParticles;
		Nf -= 3;
			   
		Kinstant -= 0.5*mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
   		Tinstant = 2.0*Kinstant/Nf/kB;
			
		sigma[j] = sigmaB;
		identity[j] = 2;
		mass[j] = massB;
			
		if(i != j){
			position[i].x = position[j].x;
   			position[i].y = position[j].y;
   			position[i].z = position[j].z;
   			
   			position[j].x = position[k].x;
   			position[j].y = position[k].y;
	   		position[j].z = position[k].z;
	   		
	   		position[k].x = -nan("");
   			position[k].y = -nan("");
	   		position[k].z = -nan("");
   		
   			velocity[i].x = velocity[j].x;
   			velocity[i].y = velocity[j].y;
	   		velocity[i].z = velocity[j].z;
	   		
   			velocity[j].x = velocity[k].x;
   			velocity[j].y = velocity[k].y;
	   		velocity[j].z = velocity[k].z;
	   		
	   		velocity[k].x = -nan("");
   			velocity[k].y = -nan("");
	   		velocity[k].z = -nan("");
	   		
	   		enter_pore[i].x = enter_pore[j].x; enter_pore[i].y = enter_pore[j].y;
    		enter_pore[i].z = enter_pore[j].z; enter_pore[i].d = enter_pore[j].d;
    		
    		enter_pore[j].x = enter_pore[k].x; enter_pore[j].y = enter_pore[k].y;
    		enter_pore[j].z = enter_pore[k].z; enter_pore[j].d = enter_pore[k].d;
    		
    		enter_pore[k].x = 0.; enter_pore[k].y = 0.;
    		enter_pore[k].z = 0.; enter_pore[k].d = 0.;
    		
    		VanishTime[i] = VanishTime[j];
    		VanishTime[j] = VanishTime[k];
    		VanishTime[k] = TimeBig;
		}
		
		else{
			position[i].x = position[k].x;
   			position[i].y = position[k].y;
	   		position[i].z = position[k].z;
	   		
	   		position[k].x = -nan("");
   			position[k].y = -nan("");
	   		position[k].z = -nan("");
	   		
	   		velocity[i].x = velocity[k].x;
   			velocity[i].y = velocity[k].y;
	   		velocity[i].z = velocity[k].z;
	   		
	   		velocity[k].x = -nan("");
   			velocity[k].y = -nan("");
	   		velocity[k].z = -nan("");
	   		
	   		enter_pore[i].x = enter_pore[k].x; enter_pore[i].y = enter_pore[k].y;
    		enter_pore[i].z = enter_pore[k].z; enter_pore[i].d = enter_pore[k].d;
    		
    		enter_pore[k].x = 0.; enter_pore[k].y = 0.;
    		enter_pore[k].z = 0.; enter_pore[k].d = 0.;
    		
    		VanishTime[i] = VanishTime[k];
    		VanishTime[k] = TimeBig;
		}
   			
   		if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells){
			cello1 = CellTrack[i].WhichCell; RemoveFromCell(i, cello1);
			if(i != j){
				cello2 = CellTrack[j].WhichCell; RemoveFromCell(j, cello2);
				CellTrack[i].WhichCell = cello2; AddToCell(i, cello2); CellTrack[i].Prev = -1;
				
				cello3 = CellTrack[k].WhichCell; RemoveFromCell(k, cello3);
				CellTrack[j].WhichCell = cello3; AddToCell(j, cello3); CellTrack[j].Prev = -1;
				
				CellTrack[k].WhichCell = -1; CellTrack[k].Next = -1; CellTrack[k].Prev = -1;
				
				EscapeInfo(i); EscapeInfo(j);
			}
			else{
				cello2 = CellTrack[k].WhichCell; RemoveFromCell(k, cello2);
				CellTrack[i].WhichCell = cello2; AddToCell(i, cello2); CellTrack[i].Prev = -1;
				
				CellTrack[k].WhichCell = -1; CellTrack[k].Next = -1; CellTrack[k].Prev = -1;
				
				EscapeInfo(i);
			}
		}
   			
   		for(int l = 0; l<NumberOfParticles; l++){InsdelCollisionUpdate(l, CellTrack[l].WhichCell); MemCollisionInfo(l); VanishInfo(l);}
   			
   		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
	}
	else{
		j = NumberOfParticles-1;
		for(int l = 0; l < MAX(N_0_2, N_Z2_2); l++){
			if(n_0_2[l] == j){
				n_0_2[l] = i; break;
			}
		}
		
		Mcom -= mass[i];
   		Vcom.x -= mass[i]*velocity[i].x;
   		Vcom.y -= mass[i]*velocity[i].y;
   		Vcom.z -= mass[i]*velocity[i].z;
   		
   		NumberOfParticles--; NB--;
   		fA = 1.*NA/NumberOfParticles;
 		fB = 1.*NB/NumberOfParticles;
		Nf -= 3;
		
		Kinstant -= 0.5*mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
   		Tinstant = 2.0*Kinstant/Nf/kB;
		
   		if(i != j){
   			position[i].x = position[j].x;
   			position[i].y = position[j].y;
   			position[i].z = position[j].z;
   			
   			position[j].x = -nan("");
   			position[j].y = -nan("");
 	  		position[j].z = -nan("");
 	  		
 	  		velocity[i].x = velocity[j].x;
			velocity[i].y = velocity[j].y;
			velocity[i].z = velocity[j].z;
			
			velocity[j].x = -nan("");
   			velocity[j].y = -nan("");
   			velocity[j].z = -nan("");
   			
   			enter_pore[i].x = enter_pore[j].x; enter_pore[i].y = enter_pore[j].y;
    		enter_pore[i].z = enter_pore[j].z; enter_pore[i].d = enter_pore[j].d;
    		
    		enter_pore[j].x = 0.; enter_pore[j].y = 0.;
    		enter_pore[j].z = 0.; enter_pore[j].d = 0.;
    		
    		VanishTime[i] = VanishTime[j];
    		VanishTime[j] = TimeBig;
		}
		else{
			position[i].x = -nan("");
   			position[i].y = -nan("");
 	  		position[i].z = -nan("");
 	  		
 	  		velocity[i].x = -nan("");
   			velocity[i].y = -nan("");
   			velocity[i].z = -nan("");
   			
   			enter_pore[i].x = 0.; enter_pore[i].y = 0.;
    		enter_pore[i].z = 0.; enter_pore[i].d = 0.;
    		
    		VanishTime[i] = TimeBig;
   		}
   		
   		if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells){
   			cello1 = CellTrack[i].WhichCell; RemoveFromCell(i, cello1);
			if(i != j){
				cello2 = CellTrack[j].WhichCell; RemoveFromCell(j, cello2);
				CellTrack[i].WhichCell = cello2; AddToCell(i, cello2); CellTrack[i].Prev = -1;
				
				CellTrack[j].WhichCell = -1; CellTrack[j].Next = -1; CellTrack[j].Prev = -1;
				
				EscapeInfo(i);
			}
			else{
				CellTrack[i].WhichCell = -1; CellTrack[i].Next = -1; CellTrack[i].Prev = -1;
			}
		}
 		
		for(int l = 0; l<NumberOfParticles; l++){InsdelCollisionUpdate(l, CellTrack[l].WhichCell); MemCollisionInfo(l); VanishInfo(l);}
   			
		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
	}
	return;
}
