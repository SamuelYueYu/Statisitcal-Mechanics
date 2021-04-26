/********** insertion and/or deletion of particles at any moment of time **********/

#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Initialize(void){
	N_0_1 = 0; N_Z1_1 = 0; N_Z2_1 = 0; N_Z3_1 = 0;
	N_0_2 = 0; N_Z1_2 = 0; N_Z2_2 = 0; N_Z3_2 = 0;
   
    for(int i=0; i<MAX_NumberOfParticles; i++){
   	 n_0_1[i] = -1;
   	 n_Z2_1[i] = -1;
   	 
   	 n_0_2[i] = -1;
   	 n_Z2_2[i] = -1;
    }
    for(int i=0; i<NumberOfParticles; i++){
     if(identity[i] == 1){
     	if(position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2){N_0_1++; n_0_1[N_0_1-1] = i;}
   	 	else if(position[i].z > Z1 - sigma[i]/2 && position[i].z < Z2 + sigma[i]/2){N_Z1_1++;}
   	 	else if(position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2){N_Z2_1++; n_Z2_1[N_Z2_1-1] = i;}
   	 	else{N_Z3_1++;}
	 }
   	 else{
   	 	if(position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2){N_0_2++; n_0_2[N_0_2-1] = i;}
   	 	else if(position[i].z > Z1 - sigma[i]/2 && position[i].z < Z2 + sigma[i]/2){N_Z1_2++;}
   	 	else if(position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2){N_Z2_2++; n_Z2_2[N_Z2_2-1] = i;}
   	 	else{N_Z3_2++;}
	 }
    }
	return;
}

bool Overlap(double x, double y, double z, double d, double vx, double vy, double vz){
	VECTOR rij;
	double xtemp, ytemp, ztemp;
	int ncell,nx,ny,nz;
	
	if(NumberOfCells <= NumberOfNeighborCells || CellSwitch == 0){
		for(int j=0; j<NumberOfParticles; j++){
			rij.d = (d+sigma[j])/2;
			
			for(int l=-1;l<=1;l++)
    		for(int m=-1;m<=1;m++)
    		for(int n=-1;n<=1;n++){
				xtemp = position[j].x + l*Lx;
   				ytemp = position[j].y + m*Ly;
   				ztemp = position[j].z + n*Lz;
				
   				rij.x = x - xtemp;
   				rij.y = y - ytemp;
   				rij.z = z - ztemp;
   				
   				if(SQR(rij.x) + SQR(rij.y) + SQR(rij.z) < SQR(rij.d /*+ 10E-6*/)){return true;}
  			}// end loop periodic images
 		}//end loop j
	}
	else{
		nx = (int)(x/rcellx); ny = (int)(y/rcelly); nz = (int)(z/rcellz);
		
 		if(fabs(x/rcellx - (double)(nx)) < tolerance && vx < 0.){nx = (nx-1+LCellx)%LCellx;}
 		else if(fabs(x/rcellx - (double)(nx+1)) < tolerance && vx > 0.){nx = (nx+1)%LCellx;}
 		
 		if(fabs(y/rcelly - (double)(ny)) < tolerance && vy < 0.){ny = (ny-1+LCelly)%LCelly;}
 		else if(fabs(y/rcelly - (double)(ny+1)) < tolerance && vy > 0.){ny = (ny+1)%LCelly;}
 		
 		if(fabs(z/rcellz - (double)(nz)) < tolerance && vz < 0.){nz = (nz-1+LCellz)%LCellz;}
 		else if(fabs(z/rcellz - (double)(nz+1)) < tolerance && vz > 0.){nz = (nz+1)%LCellz;}
 		
 		ncell = (nz*LCellz+ny)*LCelly+nx;
 		
 		for(int cell=0;cell<NumberOfNeighborCells;cell++){
  			int jCell = NeighborCellList[ncell][cell];
  			
  			for(int j = HeadOfChain[jCell]; j != -1; j = CellTrack[j].Next){
  				rij.x = x - position[j].x;
  				rij.y = y - position[j].y;
  				rij.z = z - position[j].z;
  				rij.d = (d+sigma[j])/2;
  				//pbc
  				MinimumImage(&(rij.x),Lx);
  				MinimumImage(&(rij.y),Ly);
  				MinimumImage(&(rij.z),Lz);

    			if(SQR(rij.x) + SQR(rij.y) + SQR(rij.z) < SQR(rij.d)){return true;}
  			}//end while j != -1
 		}//end loop cell
	}
 	return false;
}

void Insert(double z_small, double z_great, int Identity, double t){
	int i, j, N_space_id;
	int cello1;
	double mu, prob, scale;
	double PosX, PosY, PosZ, Sigma;
	double vx, vy, vz;
	
	if(z_small == 0. && Identity == 1){N_space_id = N_0_1; N_insert_trial_0_1++;}
	else if(z_small == 0. && Identity == 2){N_space_id = N_0_2; N_insert_trial_0_2++;}
	else if(z_small != 0. && Identity == 1){N_space_id = N_Z2_1; N_insert_trial_Z2_1++;}
	else{N_space_id = N_Z2_2; N_insert_trial_Z2_2++;}
	
	Sigma = Identity == 1 ? sigmaA : sigmaB;
	PosX = RandomNumber() * Lx;
	PosY = RandomNumber() * Ly;
	PosZ = z_small + Sigma/2 + (z_great - z_small - Sigma) * RandomNumber();
	vx = BoxMuller(0.,1.); //*sqrt(T/Tinstant)
	vy = BoxMuller(0.,1.);
	vz = BoxMuller(0.,1.);
	
	mu = z_small == 0. ? muA : muB;
	prob = Overlap(PosX, PosY, PosZ, Sigma, vx, vy, vz) ? 0. : Lx*Ly*(z_great - z_small - Sigma) / (lambda * (N_space_id + 1)) * exp(beta*mu);
	
	bool yes = RandomNumber() < prob;
	if(yes && Identity == 1){
		if(z_small == 0.){
			N_insert_0_1++; if(t >= SCT){M1_1++;}
			N_0_1++; n_0_1[N_0_1-1] = NA;
		}
		else{
			N_insert_Z2_1++; if(t >= SCT){M2_1++;}
			N_Z2_1++; n_Z2_1[N_Z2_1-1] = NA;
		}
		
		NumberOfParticles++; NA++;
		i = NA-1; j = NumberOfParticles-1;
		sigma[j] = sigmaB; sigma[i] = sigmaA;
		identity[j] = 2; identity[i] = 1;
		mass[j] = massB; mass[i] = massA;
		
		for(int k = 0; k < N_0_2; k++){
			if(n_0_2[k] == i){
				n_0_2[k] = j; break;
			}
		}
		
		fA = 1.*NA/NumberOfParticles;
 		fB = 1.*NB/NumberOfParticles;
		Nf += 3;
		
		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
		
		// begin velocity assignment
		velocity[j].x = velocity[i].x;
   		velocity[j].y = velocity[i].y;
   		velocity[j].z = velocity[i].z;
		
		velocity[i].x = vx;
   		velocity[i].y = vy;
   		velocity[i].z = vz;
   	
   		Mcom += mass[i];

  		Vcom.x += mass[i]*velocity[i].x;
  		Vcom.y += mass[i]*velocity[i].y;
  		Vcom.z += mass[i]*velocity[i].z;
   		
   		//Vcom.x /= Mcom;
   		//Vcom.y /= Mcom;
   		//Vcom.z /= Mcom;
   		
   		Kinstant += 0.5*mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
   		Tinstant = 2.0*Kinstant/Nf/kB;
    	// end of velocity assignment
    	
    	// begin position assignment
    	position[j].x = position[i].x;
   		position[j].y = position[i].y;
   		position[j].z = position[i].z;
   		position[i].x = PosX;
		position[i].y = PosY;
		position[i].z = PosZ;
   		
   		if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells){
			cello1 = CellTrack[i].WhichCell; RemoveFromCell(i, cello1);
			CellTrack[i].WhichCell = CellDetermine(i); CellTrack[i].Next = HeadOfChain[CellTrack[i].WhichCell];
   			if(CellTrack[i].Next != -1){CellTrack[CellTrack[i].Next].Prev = i;} // avoid array[-1] which is segmentation  fault
   			HeadOfChain[CellTrack[i].WhichCell] = i; CellTrack[i].Prev = -1;
			
			CellTrack[j].WhichCell = CellDetermine(j); CellTrack[j].Next = HeadOfChain[CellTrack[j].WhichCell];
   			if(CellTrack[j].Next != -1){CellTrack[CellTrack[j].Next].Prev = j;} // avoid array[-1] which is segmentation  fault
   			HeadOfChain[CellTrack[j].WhichCell] = j; CellTrack[j].Prev = -1;
			
			EscapeInfo(i); EscapeInfo(j);
		}
    	// end of position assignment
    	
    	enter_pore[j].x = enter_pore[i].x; enter_pore[j].y = enter_pore[i].y;
    	enter_pore[j].z = enter_pore[i].z; enter_pore[j].d = enter_pore[i].d;
    	
    	for(int k = 0; k<NumberOfParticles; k++){
			InsdelCollisionUpdate(k, CellTrack[k].WhichCell);
			MemCollisionInfo(k); VanishInfo(k);
		}
	}
	else if(yes){
		if(z_small == 0.){
			N_insert_0_2++; if(t >= SCT){M1_2++;}
			N_0_2++; n_0_2[N_0_2-1] = NumberOfParticles;
		}
		else{
			N_insert_Z2_2++; if(t >= SCT){M2_2++;}
			N_Z2_2++; n_Z2_2[N_Z2_2-1] = NumberOfParticles;
		}
		
		NumberOfParticles++; NB++; i = NumberOfParticles-1;
		sigma[i] = sigmaB;
		identity[i] = 2;
		mass[i] = massB;
		
		fA = 1.*NA/NumberOfParticles;
 		fB = 1.*NB/NumberOfParticles;
		Nf += 3;
		
		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
			
		// begin velocity assignment
		velocity[i].x = vx;
   		velocity[i].y = vy;
   		velocity[i].z = vz;
   		
   		Mcom += mass[i];

   		Vcom.x += mass[i]*velocity[i].x;
   		Vcom.y += mass[i]*velocity[i].y;
   		Vcom.z += mass[i]*velocity[i].z;
   		
   		Vcom.x /= Mcom;
   		Vcom.y /= Mcom;
   		Vcom.z /= Mcom;
    		
    	Kinstant += 0.5*mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
   		Tinstant = 2.0*Kinstant/Nf/kB;
    	// end of velocity assignment
    		
    	// begin position assignment
    	position[i].x = PosX;
		position[i].y = PosY;
		position[i].z = PosZ;
			
		if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells){
			CellTrack[i].WhichCell = CellDetermine(i); CellTrack[i].Next = HeadOfChain[CellTrack[i].WhichCell];
   			if(CellTrack[i].Next != -1){CellTrack[CellTrack[i].Next].Prev = i;} // avoid array[-1] which is segmentation  fault
   			HeadOfChain[CellTrack[i].WhichCell] = i; CellTrack[i].Prev = -1;
   			
   			EscapeInfo(i);
		}
    	// end of position assignment
    	
    	for(int k = 0; k<NumberOfParticles; k++){
			InsdelCollisionUpdate(k, CellTrack[k].WhichCell);
			MemCollisionInfo(k); VanishInfo(k);
		}
	}
	return;
}

void Delete(double z_small, double z_great, int Identity, double t){
	int i, ii, j, k, N_space_id;
	int cello1, cello2, cello3, celln;
	double mu, prob, scale;
	
	if(z_small == 0. && Identity == 1){
		N_space_id = N_0_1;
		if(N_0_1){
			N_delete_trial_0_1++;
			ii = (int) (RandomNumber() * N_0_1); i = n_0_1[ii];
		}
	}
	else if(z_small == 0. && Identity == 2){
		N_space_id = N_0_2;
		if(N_0_2){
			N_delete_trial_0_2++; 
			ii = (int) (RandomNumber() * N_0_2); i = n_0_2[ii];
		}
	}
	else if(z_small == Z2 && Identity == 1){
		N_space_id = N_Z2_1;
		if(N_Z2_1){
			N_delete_trial_Z2_1++; 
			ii = (int) (RandomNumber() * N_Z2_1); i = n_Z2_1[ii];
		}
	}
	else{
		N_space_id = N_Z2_2;
		if(N_Z2_2){
			N_delete_trial_Z2_2++; 
			ii = (int) (RandomNumber() * N_Z2_2); i = n_Z2_2[ii];
		}
	}
	
	mu = z_small == 0. ? muA : muB;
	prob = N_space_id ? N_space_id*lambda / (Lx*Ly*(z_great - z_small - sigma[i])) * exp(-beta*mu) : 0.;
		
	bool yes = RandomNumber() < prob;
	if(yes && Identity == 1){
		if(z_small == 0.){
			N_delete_0_1++; if(t >= SCT){M1_1--;}
			n_0_1[ii] = n_0_1[N_0_1-1]; n_0_1[N_0_1-1] = -1;
			N_0_1--;
		}
		else{
			N_delete_Z2_1++; if(t >= SCT){M2_1--;}
			n_Z2_1[ii] = n_Z2_1[N_Z2_1-1]; n_Z2_1[N_Z2_1-1] = -1;
			N_Z2_1--;
		}
		
		j = NA - 1; k = NumberOfParticles-1;
		for(int l = 0; l < N_0_1; l++){
			if(n_0_1[l] == j){
				n_0_1[l] = i; break;
			}
		}
		for(int l = 0; l < N_0_2; l++){
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
   			
   		for(int l = 0; l<NumberOfParticles; l++){
		   InsdelCollisionUpdate(l, CellTrack[l].WhichCell);
		   MemCollisionInfo(l); VanishInfo(l);
		}
   		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
	}
	else if(yes){
		if(z_small == 0.){
			N_delete_0_2++; if(t >= SCT){M1_2--;}
			n_0_2[ii] = n_0_2[N_0_2-1]; n_0_2[N_0_2-1] = -1;
			N_0_2--;
		}
		else{
			N_delete_Z2_2++; if(t >= SCT){M2_2--;}
			n_Z2_2[ii] = n_Z2_2[N_Z2_2-1]; n_Z2_2[N_Z2_2-1] = -1;
			N_Z2_2--;
		}
		
		j = NumberOfParticles-1;
		for(int l = 0; l < N_0_2; l++){
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
 		
		for(int l = 0; l<NumberOfParticles; l++){
			InsdelCollisionUpdate(l, CellTrack[l].WhichCell);
			MemCollisionInfo(l); VanishInfo(l);
		}
   			
		rho = NumberOfParticles/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
		rhoA = NA/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
		rhoB = NB/RECTANGULAR(Lx,Ly,Lz)*RECTANGULAR(sigmaB,sigmaB,sigmaB);
		packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
	}
	return;
}

void InsDel(double t){
	Initialize();
	
	if(RandomNumber() < .5){Delete(0., Z1, 1, t);}
	else{Insert(0., Z1, 1, t);}
	
//	if(RandomNumber() < .5){Delete(Z2, Z3, 1, t);}
//	else{Insert(Z2, Z3, 1, t);}
	
	if(RandomNumber() < .5){Delete(0., Z1, 2, t);}
	else{Insert(0., Z1, 2, t);}
	
//	if(RandomNumber() < .5){Delete(Z2, Z3, 2, t);}
//	else{Insert(Z2, Z3, 2, t);}
	
	return;
}
