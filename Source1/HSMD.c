/***************************************************
* 3D Molecular Dynamics (MD) simulation of
* binary mixtures: particle A B
* of hard spheres
* dual control volume 
* flux
* Kai Zhang, Yue Yu, Duke Kunshan University, 2019
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

//int main(void)
int main(int argc, char *argv[])
{
 int event, Mevent; // 0: particle collision 1: transfer 3: membrane collision
 double time, RN;
 double sampletime, sampletime_L, sampletime_R; // sampling time
 int Nleft = 0, Nright = 0, NMemleft = 0, NMemright = 0;
 int NLeftList[MAX_NumberOfParticles],NRightList[MAX_NumberOfParticles];
 int NLeftMemList[MAX_NumberOfParticles],NRightMemList[MAX_NumberOfParticles];
 double rhoLeft,rhoRight,rhoMemLeft,rhoMemRight;
 double tij; // collision time
 double tescape; // time needs to escape a cell
 double t_mem; //membrane collision time
 int I,J; // collision particle index
 int cello,celln;
 int i,k,n=0;
 char filename[20];
 FILE *fp;
 FILE *fpmovie;
 FILE *fpvelocity;
 FILE *fpsample;
 FILE *fpr; // position
 FILE *fpv; // velocity
 FILE *fpinsdel_1; //particle numbers for type 1
 FILE *fpinsdel_2; //particle numbers for type 2
 FILE *fpfreqinsdel_1; //insdel frequency for type 1
 FILE *fpfreqinsdel_2; //insdel frequency for type 2
 int onboundx[MAX_NumberOfParticles], onboundy[MAX_NumberOfParticles], onboundz[MAX_NumberOfParticles];
 int rim1[MAX_NumberOfParticles], rim2[MAX_NumberOfParticles];
 int rim3[MAX_NumberOfParticles], rim4[MAX_NumberOfParticles], rim5[MAX_NumberOfParticles];
 VECTOR rij;

 double Sum_K,Sum_U,Sum_E,Sum_W,Sum_T;
 double Sum_WLeft,Sum_WRight,Sum_WMemLeft,Sum_WMemRight;
 double Sum_RhoL,Sum_RhoR,Sum_RhoMemL,Sum_RhoMemR;
 double Count, pCount;

 double maxL;

 sscanf(argv[1],"%d",&JobIndex);
 
 printf("**************** 3D Hard Sphere Molecular Dynamics simulation ********************");
 printf("\n");

 ReadInput();
 MakeMembrane();
 Tkcommand();

 sprintf(filename,"insdel_1.txt");
 fpinsdel_1=fopen(filename,"w");
 fprintf(fpinsdel_1,"#time\tN_0_1\tN_Z1_1\tN_Z2_1\tN_Z3_1\trho_0_1\trho_Z1_1\trho_Z2_1\trho_Z3_1\n");
 fclose(fpinsdel_1);
 	 
 sprintf(filename,"insdel_2.txt");
 fpinsdel_2=fopen(filename,"w");
 fprintf(fpinsdel_2,"#time\tN_0_2\tN_Z1_2\tN_Z2_2\tN_Z3_2\trho_0_2\trho_Z1_2\trho_Z2_2\trho_Z3_2\n");
 fclose(fpinsdel_2);
 
 sprintf(filename,"freqinsdel_1.txt");
 fpfreqinsdel_1=fopen(filename,"w");
 if(N_insert_trial_0_1 && N_insert_trial_Z2_1 && N_delete_trial_0_1 && N_delete_trial_Z2_1){
 	fprintf(fpfreqinsdel_1,"#time\tinsertN_0_1\tinsertN_Z2_1\tdeleteN_0_1\tdeleteN_Z2_1\tinsert_0_1\tinsert_Z2_1\tdelete_0_1\tdelete_Z2_1\n");
 }
 fclose(fpfreqinsdel_1);
  	
 sprintf(filename,"freqinsdel_2.txt");
 fpfreqinsdel_2=fopen(filename,"w");
 if(N_insert_trial_0_2 && N_insert_trial_Z2_2 && N_delete_trial_0_2 && N_delete_trial_Z2_2){
 	fprintf(fpfreqinsdel_2,"#time\tinsertN_0_1\tinsertN_Z2_1\tdeleteN_0_1\tdeleteN_Z2_1\tinsert_0_1\tinsert_Z2_1\tdelete_0_1\tdelete_Z2_1\n");
 }
 fclose(fpfreqinsdel_2);

 if(randomseed == 0.)  randomseed = (double)(JobIndex);
 printf("randomseed = %lf\n",randomseed);
 InitializeRandomNumberGenerator(randomseed); //time(0L)
 Initialization(JobIndex);

 OverlapCheck();
 MakeCell();
 fflush(stdout);

for(i=0;i<NumberOfParticles;i++)
{
 CollisionTime[i] = TimeBig;
 MemCollisionTime[i] = TimeBig;
 CollisionPartner[i] = NumberOfParticles;
 EscapeTime[i] = TimeBig;
 if(!(position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2) && !(position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2)){
 	for(int j = 0; j < M; j++){
		if(SQR(position[i].x - pore_position[j].x)+SQR(position[i].y - pore_position[j].y) < SQR(pore_position[j].d/2)){
			enter_pore[i].x = pore_position[j].x;
			enter_pore[i].y = pore_position[j].y;
			enter_pore[i].d = pore_position[j].d;
			break;
		}
	}
 }
}
tescape = TimeBig;

//initialize collision and cell escaping time
for(i=0;i<NumberOfParticles;i++)
{
 CollisionInfo(i); // return collision time and partner
 MemCollisionInfo(i);

 if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)//use Cell list
 EscapeInfo(i); // return excape time
}

//initialize density values at different positions
for(i=0;i<MAX_NumberOfCells;i++){rho_1[i] = 0; rho_2[i] = 0;}
for(i=0;i<MAX_M;i++){rhoMem_left_1[i] = 0; rhoMem_right_1[i] = 0; rhoMem_left_2[i] = 0; rhoMem_right_2[i] = 0;}

n0_1 = 0., nZ1_1 = 0., nZ2_1 = 0., nZ3_1 = 0.;
n0_2 = 0., nZ1_2 = 0., nZ2_2 = 0., nZ3_2 = 0.;
//initialize left and right chamber particle lists
//for(i=0;i<MAX_NumberOfParticles;i++){NLeftList[i] = -1; NRightList[i] = -1;}

fflush(stdout);

 /***********************************************************
 sprintf(filename,"position0_%d.dat",JobIndex);
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity0_%d.dat",JobIndex);
 fpv=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%lf\t%lf\t%lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%lf\t%lf\t%lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 ***********************************************************/

 printf("start MD loop ......\n");
 printf("\n");

collisionfrequency = 0.0;  

 if(MovieMultiplier < 10000) 
 {
  sprintf(filename,"movie_%d.xyz",JobIndex);
  fpmovie=fopen(filename,"w");
  Writemovie(fpmovie);
  fclose(fpmovie);
 }

time = 0.;
tij = 0.;
t_mem = 0.;
sampletime = 0., sampletime_L = 0., sampletime_R = 0.;

 sprintf(filename,"sample_%d.dat",JobIndex);
 fpsample=fopen(filename,"w");
   fprintf(fpsample,"collision = %d\tdt = %lf\tt = %lf\tK = %lf\tVirial = %lf\tTinstant = %lf\tV = %lf\tLx = %lf\tLy = %lf\tLz = %lf\tX = %lf\tY = %lf\tZ = %lf\tphi = %lf\tncell = %d\n", \
   coln,tij,time,Kinstant/NumberOfParticles,Virial/V/3.,Tinstant,V,Lx,Ly,Lz,X,Y,Z,packingfraction,NumberOfCells);
 fclose(fpsample);

 Sum_K = 0.;
 Sum_U = 0.;
 Sum_E = 0.;
 Sum_W = 0.;
 Sum_WLeft = 0.; Sum_WRight = 0.;
 Sum_RhoL = 0.; Sum_RhoR = 0.;
 Sum_T = 0.;
 Count = 0.;
 pCount = 0.;
 
time = 0.;
tij = TimeBig;
t_mem = TimeBig;
sampletime = 0.;
event = 100;
k = 1; stream(time);

while(1) // relaxation loop
{
  // next collision
  /*********************************/
   MemOverlapCheck();
   if(RandomNumber() < 1.0/2){InsDel(time);}
   if(time >= SCT){Jlist(MIN(MIN(tij,tescape), t_mem));}

   if(n == 999){
	   Initialize(); InstantBulkDensities(time);
	   if(time >= SCT){Rho();}
	   stream(time); n = 0;
	   sprintf(filename,"position.dat");
	   fpr=fopen(filename,"w");
	   sprintf(filename,"velocity.dat");
	   fpv=fopen(filename,"w");
	   for(i=0;i<NumberOfParticles;i++){
		   fprintf(fpr,"%.12lf\t%.12lf\t%.12lf\n",position[i].x,position[i].y,position[i].z);
		   fprintf(fpv,"%.12lf\t%.12lf\t%.12lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
	   }
	   fclose(fpr);
	   fclose(fpv);
	   PrintInput();
   }
   else{n++;}

   tij = TimeBig;
   t_mem = TimeBig;

   for(i=0;i<NumberOfParticles;i++)
   {
   	onboundx[i] = 0; onboundy[i] = 0; onboundz[i] = 0;
   	rim1[i] = 0; rim2[i] = 0;
	rim3[i] = 0; rim4[i] = 0; rim5[i] = 0;
   	
    if(CollisionTime[i] < tij)
    {
     tij = CollisionTime[i];
     I = i; // collision particle 1
    }
    if(MemCollisionTime[i] < t_mem){
    	t_mem = MemCollisionTime[i];
	}
   }
   J = CollisionPartner[I]; // collision particle 2
  /*********************************/
  //printf("Here after collision and memcolllision time check\n");
  // next transfer
  /*********************************/
  if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells) // if use cell list 
  {
   tescape = TimeBig;
   for(i=0;i<NumberOfParticles;i++)
   {
    if(EscapeTime[i] < tescape)
    {
     tescape = EscapeTime[i];
     EscapeID = i;
    }
   }
  }
  if(time >= SCT){for(int i=0;i<NumberOfParticles;i++){position_old[i] = position[i];}}
/*********************************/
if(time+MIN(MIN(tij,tescape),t_mem) > timewindow) // if time exceeds timewindow, evolve dt/2, break
{
  event = 2;
  if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells) tij = MIN(MIN(tij,tescape),t_mem);
  time += tij/2.; // evolve dt/2 so that not particle on contact
  collisionfrequency = coln/time; // collision frequency in timewindow

  for(i=0;i<NumberOfParticles;i++)
  {
    position[i].x += velocity[i].x * tij/2.;
    position[i].y += velocity[i].y * tij/2.;
    position[i].z += velocity[i].z * tij/2.;
    PBC(&(position[i].x),Lx); // periodic boundary condition 
    PBC(&(position[i].y),Ly); // periodic boundary condition
    PBC(&(position[i].z),Lz); // periodic boundary condition
    CollisionTime[i] -= tij/2.;
    MemCollisionTime[i] -= tij/2.;
 	
    if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
     EscapeTime[i] -= tij/2.;
  }
  break; // getting out of relaxation loop and no particles are at contact
}
/************  transfer or collision ********************/
 if(tij == MIN(MIN(tij,tescape),t_mem)) // collision
  {
   event = 0;
   time += tij; // increase time
   coln++;
 
   for(i=0;i<NumberOfParticles;i++)
   {
	CollisionTime[i] -= tij; MemCollisionTime[i] -= tij;
    position[i].x += velocity[i].x * tij;
    position[i].y += velocity[i].y * tij;
    position[i].z += velocity[i].z * tij;
    PBC(&(position[i].x),Lx); // periodic boundary condition
    PBC(&(position[i].y),Ly); // periodic boundary condition
    PBC(&(position[i].z),Lz); // periodic boundary condition
    if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells){EscapeTime[i] -= tij;}
    
    if(position[i].z > 0. && position[i].z < sigma[i]/2 \
	&& f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
	    rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		rij.z = position[i].z;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim1[i] = 1;}
	}
    else if(position[i].z > Z1-sigma[i]/2 && position[i].z < Z1 \
	&& f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		rij.z = position[i].z - Z1;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim2[i] = 1;}
	}
    else if(position[i].z > Z2 && position[i].z < Z2+sigma[i]/2 \
	&& f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
	  	rij.z = position[i].z - Z2;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim3[i] = 1;}
	}
    else if(position[i].z > Z3-sigma[i]/2 && position[i].z < Z3 \
	&& f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
	  	rij.z = position[i].z - Z3;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim4[i] = 1;}
	}
	else if((!position[i].z || (position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)) \
	&& sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) - (enter_pore[i].d - sigma[i]) / 2 > -1.0E-3){
		rij.x = sigma[i] * (enter_pore[i].x - position[i].x) / (enter_pore[i].d - sigma[i]);
		rij.y = sigma[i] * (enter_pore[i].y - position[i].y) / (enter_pore[i].d - sigma[i]);
		if(rij.x*velocity[i].x + rij.y*velocity[i].y < 0.){rim5[i] = 1;}
	}
	if(rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){MemCollision(i);}
   }
   Collision(I,J); // collision dynamics

   if(coln>=NumberOfInitialSteps)
   {
    if(time > SCT){
    	Sum_W += Virial; sampletime += tij;
    	if(VirialLeft != 0.){Sum_WLeft += VirialLeft; sampletime_L += tij;}
    	if(VirialLeftMem != 0.){Sum_WMemLeft += VirialLeftMem;}
    	if(VirialRight != 0.){Sum_WRight += VirialRight; sampletime_R += tij;}
    	if(VirialRightMem != 0.){Sum_WMemRight += VirialRightMem;}
	}
   }

  if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells) //use Cell list
  {
  	for(i=0;i<NumberOfParticles;i++){
	  if(i == I || CollisionPartner[i] == I || i == J || CollisionPartner[i] == J){CollisionInfo(i);}
	  
	  if(rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){
       	MemCollisionInfo(i); EscapeInfo(i);
		for(int j = 0; j<NumberOfParticles; j++){if(j == i || CollisionPartner[j] == i){CollisionInfo(j);}}
      }
	}
	EscapeInfo(I); EscapeInfo(J);
    MemCollisionInfo(I); MemCollisionInfo(J);
  }
  
  	if(time > SCT && !n){
  		for(int i=0;i<NumberOfParticles;i++){
   			if(position[i].z >= 0. && position[i].z < Z1){Nleft++; NLeftList[Nleft-1] = i;}
   			else if(position[i].z >= Z1 && position[i].z < Z2){NMemleft++; NLeftMemList[NMemleft-1] = i;}
   			else if(position[i].z >= Z2 && position[i].z < Z3){Nright++; NRightList[Nright-1] = i;}
   			else{NMemright++; NRightMemList[NMemright-1] = i;}
      	}
      	rhoLeft = Nleft/(Lx*Ly*Z1); rhoRight = Nright/(Lx*Ly*(Z3-Z2));
      	rhoMemLeft = NMemleft/(Area1*(Z2-Z1)); rhoMemRight = NMemright/(Area1*(Lz-Z3));
   	  	
   	  	Sum_RhoL += rhoLeft; Sum_RhoMemL += rhoMemLeft;
   	  	Sum_RhoR += rhoRight; Sum_RhoMemR += rhoMemRight;
   	  	
   	  	pLeftGr += rhoLeft*kB*T * (1 + gTerm(Nleft, NLeftList, rhoLeft));
   	  	pRightGr += rhoRight*kB*T * (1 + gTerm(Nright, NRightList, rhoRight));
   	  	
   	  	pLeftCS += rhoLeft*T*Compressibility(M_PI*rhoLeft/6);
	  	pRightCS += rhoRight*T*Compressibility(M_PI*rhoRight/6);
	  	pCount++; k++;
	  	Nleft = 0.; Nright = 0.;
	  	NMemleft = 0.; NMemright = 0.;
  		for(int i=0;i<NumberOfParticles;i++){
		  NLeftList[i] = -1; NRightList[i] = -1;
		  NLeftMemList[NMemleft-1] = -1; NRightList[Nright-1] = -1;
		}
  	}
}// end if collision
 else if(tescape == MIN(MIN(tij,tescape),t_mem)) // if transfer 
 {
 	event = 1;
    tij = tescape;
    if(coln>=NumberOfInitialSteps) sampletime += tij;
    time += tij; // increase time
   
    if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells){
	 for(i=0;i<NumberOfParticles;i++){
   	  CollisionTime[i] -= tij; MemCollisionTime[i] -= tij;
   	  position[i].x += velocity[i].x * tij;
      position[i].y += velocity[i].y * tij;
      position[i].z += velocity[i].z * tij;
      PBC(&(position[i].x),Lx); // periodic boundary condition
      PBC(&(position[i].y),Ly); // periodic boundary condition
      PBC(&(position[i].z),Lz); // periodic boundary condition
    
      if((fabs(position[i].x/rcellx - (double)((int)(position[i].x/rcellx))) < tolerance) && (velocity[i].x < 0.) \
       ||(fabs(position[i].x/rcellx - (double)((int)(position[i].x/rcellx))+1) < tolerance) && (velocity[i].x > 0.))
      onboundx[i] = 1;
    
      if((fabs(position[i].y/rcelly - (double)((int)(position[i].y/rcelly))) < tolerance) && (velocity[i].y < 0.) \
       ||(fabs(position[i].y/rcelly - (double)((int)(position[i].y/rcelly))+1) < tolerance) && (velocity[i].y > 0.))
      onboundy[i] = 1;
    
      if((fabs(position[i].z/rcellz - (double)((int)(position[i].z/rcellz))) < tolerance) && (velocity[i].z < 0.) \
       ||(fabs(position[i].z/rcellz - (double)((int)(position[i].z/rcellz))+1) < tolerance) &&(velocity[i].z > 0.))
      onboundz[i] = 1;
    
      if(onboundx[i] == 1) position[i].x = round(position[i].x/rcellx)*rcellx;
      if(onboundy[i] == 1) position[i].y = round(position[i].y/rcelly)*rcelly;
      if(onboundz[i] == 1) position[i].z = round(position[i].z/rcellz)*rcellz;

      if(position[i].z > 0. && position[i].z < sigma[i]/2 \
	  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
	  	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		rij.z = position[i].z;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim1[i] = 1;}
	  }
      else if(position[i].z > Z1-sigma[i]/2 && position[i].z < Z1 \
	  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		rij.z = position[i].z - Z1;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim2[i] = 1;}
	  }
      else if(position[i].z > Z2 && position[i].z < Z2+sigma[i]/2 \
	  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
	  	rij.z = position[i].z - Z2;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim3[i] = 1;}
	  }
      else if(position[i].z > Z3-sigma[i]/2 && position[i].z < Z3 \
	  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
       	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
	    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
	  	rij.z = position[i].z - Z3;
		if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim4[i] = 1;}
	  }
	  else if((!position[i].z || (position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)) \
	  && sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) - (enter_pore[i].d - sigma[i]) / 2 > -1.0E-3){
		rij.x = sigma[i] * (enter_pore[i].x - position[i].x) / (enter_pore[i].d - sigma[i]);
		rij.y = sigma[i] * (enter_pore[i].y - position[i].y) / (enter_pore[i].d - sigma[i]);
		if(rij.x*velocity[i].x + rij.y*velocity[i].y < 0.){rim5[i] = 1;}
	  }
	  if(rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){MemCollision(i);}
     }
     for(i=0;i<NumberOfParticles;i++)
     {
      if(i == EscapeID || onboundx[i]==1 || onboundy[i]==1 || onboundz[i]==1) 
      {
	   cello = CellTrack[i].WhichCell;
	   celln = CellDetermine(i);

       RemoveFromCell(i,cello);
       AddToCell(i,celln);
       CellTrack[i].WhichCell = celln;

       EscapeInfo(i); CollisionUpdate(i,celln); MemCollisionInfo(i);
      }// end if particle on boundary
      else{
       EscapeTime[i] -= tij; //only for particles not on boundary
  	  }
  	  if(rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){
       	MemCollisionInfo(i);
		for(int j = 0; j<NumberOfParticles; j++){if(j == i || CollisionPartner[j] == i){CollisionInfo(j);}}
		
		if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells){EscapeInfo(i);}
      }
     }//end loop i
     
    if(time > SCT && !n){
  		for(int i=0;i<NumberOfParticles;i++){
   			if(position[i].z >= 0. && position[i].z < Z1){Nleft++; NLeftList[Nleft-1] = i;}
   			else if(position[i].z >= Z1 && position[i].z < Z2){NMemleft++; NLeftMemList[NMemleft-1] = i;}
   			else if(position[i].z >= Z2 && position[i].z < Z3){Nright++; NRightList[Nright-1] = i;}
   			else{NMemright++; NRightMemList[NMemright-1] = i;}
      	}
      	rhoLeft = Nleft/(Lx*Ly*Z1); rhoRight = Nright/(Lx*Ly*(Z3-Z2));
      	rhoMemLeft = NMemleft/(Area1*(Z2-Z1)); rhoMemRight = NMemright/(Area1*(Lz-Z3));
   	  	
   	  	Sum_RhoL += rhoLeft; Sum_RhoMemL += rhoMemLeft;
   	  	Sum_RhoR += rhoRight; Sum_RhoMemR += rhoMemRight;
   	  	
   	  	pLeftGr += rhoLeft*kB*T * (1 + gTerm(Nleft, NLeftList, rhoLeft));
   	  	pRightGr += rhoRight*kB*T * (1 + gTerm(Nright, NRightList, rhoRight));
   	  	
   	  	pLeftCS += rhoLeft*T*Compressibility(M_PI*rhoLeft/6);
	  	pRightCS += rhoRight*T*Compressibility(M_PI*rhoRight/6);
  		pCount++; k++;
	  	Nleft = 0.; Nright = 0.;
	  	NMemleft = 0.; NMemright = 0.;
  		for(int i=0;i<NumberOfParticles;i++){
		  NLeftList[i] = -1; NRightList[i] = -1;
		  NLeftMemList[NMemleft-1] = -1; NRightList[Nright-1] = -1;
		}
  	}
   }
}// end if transfer
 
 else if(t_mem == MIN(MIN(tij,tescape),t_mem)){ // membrane collision
	  event = 3;
	  tij = t_mem;
	  time += tij;

	  for(i=0;i<NumberOfParticles;i++){
	  	  CollisionTime[i] -= tij; MemCollisionTime[i] -= tij;
		  if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells){EscapeTime[i] -= tij;}
		  position[i].x += velocity[i].x * tij;
		  position[i].y += velocity[i].y * tij;
		  position[i].z += velocity[i].z * tij;
		  PBC(&(position[i].x),Lx);
		  PBC(&(position[i].y),Ly);
		  PBC(&(position[i].z),Lz);
		  if(MemCollisionTime[i] == 0.){MemCollision(i);}
		  
		  if(position[i].z > 0. && position[i].z < sigma[i]/2 \
	  	  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
		  	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
		    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
			rij.z = position[i].z;
			if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim1[i] = 1;}
		  }
      	  else if(position[i].z > Z1-sigma[i]/2 && position[i].z < Z1 \
		  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
      	  	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
		    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
			rij.z = position[i].z - Z1;
			if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim2[i] = 1;}
		  }
          else if(position[i].z > Z2 && position[i].z < Z2+sigma[i]/2 \
		  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
          	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
		    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		  	rij.z = position[i].z - Z2;
			if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim3[i] = 1;}
		  }
          else if(position[i].z > Z3-sigma[i]/2 && position[i].z < Z3 \
		  && f(enter_pore[i].d, position[i].x, enter_pore[i].x, position[i].y, enter_pore[i].y, position[i].z, enter_pore[i].z, sigma[i]) < 1.0E-6){
          	rij.x = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].x - position[i].x);
		    rij.y = (.5 * enter_pore[i].d / sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) - 1) * (enter_pore[i].y - position[i].y);
		  	rij.z = position[i].z - Z3;
			if(rij.x*velocity[i].x + rij.y*velocity[i].y + rij.z*velocity[i].z < 0.){rim4[i] = 1;}
		  }
    	  else if((!position[i].z || (position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)) \
		  && sqrt(SQR(position[i].x - enter_pore[i].x) + SQR(position[i].y - enter_pore[i].y)) - (enter_pore[i].d - sigma[i]) / 2 > -1.0E-3){
		  	rij.x = sigma[i] * (enter_pore[i].x - position[i].x) / (enter_pore[i].d - sigma[i]);
		    rij.y = sigma[i] * (enter_pore[i].y - position[i].y) / (enter_pore[i].d - sigma[i]);
		    if(rij.x*velocity[i].x + rij.y*velocity[i].y < 0.){rim5[i] = 1;}
		  }
		  if(rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){MemCollision(i);}
	  }
	  for(i=0;i<NumberOfParticles;i++){
	  	if(MemCollisionTime[i] == 0. || rim1[i] || rim2[i] || rim3[i] || rim4[i] || rim5[i]){
	  		MemCollisionInfo(i);
	  		for(int j = 0; j<NumberOfParticles; j++){if(j == i || CollisionPartner[j] == i){CollisionInfo(j);}}
			 
			if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells){EscapeInfo(i);}
		}
	  }
	  
	  if(time > SCT && !n){
  		for(int i=0;i<NumberOfParticles;i++){
   			if(position[i].z >= 0. && position[i].z < Z1){Nleft++; NLeftList[Nleft-1] = i;}
   			else if(position[i].z >= Z1 && position[i].z < Z2){NMemleft++; NLeftMemList[NMemleft-1] = i;}
   			else if(position[i].z >= Z2 && position[i].z < Z3){Nright++; NRightList[Nright-1] = i;}
   			else{NMemright++; NRightMemList[NMemright-1] = i;}
      	}
      	rhoLeft = Nleft/(Lx*Ly*Z1); rhoRight = Nright/(Lx*Ly*(Z3-Z2));
      	rhoMemLeft = NMemleft/(Area1*(Z2-Z1)); rhoMemRight = NMemright/(Area1*(Lz-Z3));
   	  	
   	  	Sum_RhoL += rhoLeft; Sum_RhoMemL += rhoMemLeft;
   	  	Sum_RhoR += rhoRight; Sum_RhoMemR += rhoMemRight;
   	  	
   	  	pLeftGr += rhoLeft*kB*T * (1 + gTerm(Nleft, NLeftList, rhoLeft));
   	  	pRightGr += rhoRight*kB*T * (1 + gTerm(Nright, NRightList, rhoRight));
   		
   	  	pLeftCS += rhoLeft*T*Compressibility(M_PI*rhoLeft/6);
	  	pRightCS += rhoRight*T*Compressibility(M_PI*rhoRight/6);
	  	pCount++; k++;
	  	Nleft = 0.; Nright = 0.;
	  	NMemleft = 0.; NMemright = 0.;
  		for(int i=0;i<NumberOfParticles;i++){
		  NLeftList[i] = -1; NRightList[i] = -1;
		  NLeftMemList[NMemleft-1] = -1; NRightList[Nright-1] = -1;
		}
  	  }
  } // end if membrane collision
  
  fflush(stdout);
  
 /************** sample *******************/
  if(event == 0 && (int) coln%SampleMultiplier==0)
  {
   Kinetic();
   
   sprintf(filename,"sample_%d.dat",JobIndex);
   fpsample=fopen(filename,"a+");
   fprintf(fpsample,"collision = %d\tdt = %.12lf\tt = %lf\tK = %lf\tVirial = %lf\tTinstant = %lf\tV = %lf\tLx = %lf\tLy = %lf\tLz = %lf\tX = %f\tY = %f\tZ = %f\tphi = %lf\tncell = %d\n",\
   coln,tij,time,Kinstant/NumberOfParticles,Virial/V/3.,Tinstant,V,Lx,Ly,Lz,X,Y,Z,packingfraction,NumberOfCells);
   fclose(fpsample);
   printf("collision %d\n",coln);
   printf("collision time %lf\n",time);
   printf("Number of cells %d\n",NumberOfCells);
   fflush(stdout);
  if(coln>=NumberOfInitialSteps)
  {
   Sum_K += Kinstant;
   Sum_T += Tinstant;
   Count += 1.;
  }//endif after equilibration
 }//endif every ? steps
  /************** end of sample *******************/

 if(event == 0 && coln%MovieMultiplier==0) 
 {
  sprintf(filename,"movie_%d.xyz",JobIndex);
  fpmovie=fopen(filename,"a+");
  Writemovie(fpmovie);
  fclose(fpmovie);
 }
} // end relaxation while loop

 Sum_K /= Count;
 Sum_U /= Count;
 Sum_E /= Count;
 Sum_T /= Count;
 Sum_W /= sampletime;
 Sum_W /= 3.;
 Sum_WLeft /= time-SCT; Sum_WRight /= time-SCT;
 Sum_WMemLeft /= time-SCT; Sum_WMemRight /= time-SCT;
 Sum_WLeft /= 3.; Sum_WRight /= 3.;
 Sum_WMemLeft /= 3.; Sum_WMemRight /= 3.;
 Sum_RhoL /= pCount; Sum_RhoR /= pCount;
 Sum_RhoMemL /= pCount; Sum_RhoMemR /= pCount;
 pLeftVirial = Sum_RhoL*kB*T + Sum_WLeft/(Lx*Ly*Z1);
 pRightVirial = Sum_RhoR*kB*T + Sum_WRight/(Lx*Ly*(Z3-Z2));
 pLeftMemVirial = Sum_RhoMemL*kB*T + Sum_WMemLeft/(Lx*Ly*(Z2-Z1));
 pRightMemVirial = Sum_RhoMemR*kB*T + Sum_WMemRight/(Lx*Ly*(Lz-Z3));
 pLeftGr /= pCount; pRightGr /= pCount;
 pLeftCS /= pCount; pRightCS /= pCount;
 
 J_1 = (fabs(M1_1) + fabs(M2_1)) / (4*Area1*(time-SCT));
 J_2 = (fabs(M1_2) + fabs(M2_2)) / (4*Area2*(time-SCT));
 JAtot = 0; JBtot = 0;
 for(i = 0; i < 2*M; i++){
 	JA[i] /= time-SCT; JB[i] /= time-SCT;
 	JAtot += JA[i]*M_PI*SQR(pore_position[i % M].d/2);
	JBtot += JB[i]*M_PI*SQR(pore_position[i % M].d/2);
 }
 JAtot /= 2*Area1; JBtot /= 2*Area2;
 
 if(averaged_over){
 	for(double j = 0.; j < Lz; j+=.5){rho_1[(int) (2*j)] /= averaged_over; rho_2[(int) (2*j)] /= averaged_over;}
 	for(i = 0; i < M; i++){
 		rhoMem_left_1[i] /= averaged_over; rhoMem_right_1[i] /= averaged_over;
		rhoMem_left_2[i] /= averaged_over; rhoMem_right_2[i] /= averaged_over;
	}
 	n0_1 /= averaged_over; nZ1_1 /= averaged_over;
 	nZ2_1 /= averaged_over; nZ3_1 /= averaged_over;
 	
 	n0_2 /= averaged_over; nZ1_2 /= averaged_over;
 	nZ2_2 /= averaged_over; nZ3_2 /= averaged_over;
 	
 	sprintf(filename,"insdel_1.txt");
 	fpinsdel_1=fopen(filename,"a+");
 	fprintf(fpinsdel_1,"//Chamber 1 total packing fraction : %lf\n", M_PI/6*(n0_1*CUBIC(sigmaA) + n0_2*CUBIC(sigmaB)));
 	fclose(fpinsdel_1);
 	 
 	sprintf(filename,"insdel_2.txt");
 	fpinsdel_2=fopen(filename,"a+");
 	fprintf(fpinsdel_1,"//Chamber 1 total packing fraction : %lf\n", M_PI/6*(n0_1*CUBIC(sigmaA) + n0_2*CUBIC(sigmaB)));
 	fclose(fpinsdel_2);
 	
 	xA = n0_1/(n0_1+n0_2); xB = n0_2/(n0_1+n0_2);
 	xAM = nZ1_1/(nZ1_1+nZ1_2); xBM = nZ1_2/(nZ1_1+nZ1_2);
	yA = nZ2_1/(nZ2_1+nZ2_2); yB = nZ2_2/(nZ2_1+nZ2_2);
	yAM = nZ3_1/(nZ3_1+nZ3_2); yBM = nZ3_2/(nZ3_1+nZ3_2);
	
	pAL = pLeftVirial*xA; pAR = pRightVirial*yA;
	pALM = pLeftMemVirial*xAM; pARM = pRightMemVirial*yAM;
	MFPLA_p = kB*T / (sqrt(2)*M_PI*pAL*SQR(sigmaA));
 	MFPRA_p = kB*T / (sqrt(2)*M_PI*pAR*SQR(sigmaA));
 	MFPLA_rho = 1 / (M_PI*n0_1*SQR(sigmaA));
 	MFPRA_rho = 1 / (M_PI*nZ2_1*SQR(sigmaA));
 	MFPLMA_p = kB*T / (sqrt(2)*M_PI*pALM*SQR(sigmaA));
 	MFPRMA_p = kB*T / (sqrt(2)*M_PI*pARM*SQR(sigmaA));
 	MFPLMA_rho = 1 / (M_PI*nZ1_1*SQR(sigmaA));
 	MFPRMA_rho = 1 / (M_PI*nZ3_1*SQR(sigmaA));
 	KnA_p[0] = MFPLA_p/D; KnA_p[2] = MFPRA_p/D;
 	KnA_p[1] = MFPLMA_p/D; KnA_p[3] = MFPRMA_p/D;
 	KnA_rho[0] = MFPLA_rho/D; KnA_rho[2] = MFPRA_rho/D;
 	KnA_rho[1] = MFPLMA_rho/D; KnA_rho[3] = MFPRMA_rho/D;
	
	pBL = pLeftVirial*xB; pBR = pRightVirial*yB;
	pBLM = pLeftMemVirial*xBM; pBRM = pRightMemVirial*yBM;
	MFPLB_p = kB*T / (sqrt(2)*M_PI*pBL*SQR(sigmaB));
 	MFPRB_p = kB*T / (sqrt(2)*M_PI*pBR*SQR(sigmaB));
 	MFPLB_rho = 1 / (M_PI*n0_2*SQR(sigmaB));
 	MFPRB_rho = 1 / (M_PI*nZ2_2*SQR(sigmaB));
 	MFPLMB_p = kB*T / (sqrt(2)*M_PI*pBLM*SQR(sigmaB));
 	MFPRMB_p = kB*T / (sqrt(2)*M_PI*pBRM*SQR(sigmaB));
 	MFPLMB_rho = 1 / (M_PI*nZ1_2*SQR(sigmaB));
 	MFPRMB_rho = 1 / (M_PI*nZ3_2*SQR(sigmaB));
 	KnB_p[0] = MFPLB_p/D; KnB_p[2] = MFPRB_p/D;
 	KnB_p[1] = MFPLMB_p/D; KnB_p[3] = MFPRMB_p/D;
 	KnB_rho[0] = MFPLB_rho/D; KnB_rho[2] = MFPRB_rho/D;
 	KnB_rho[1] = MFPLMB_rho/D; KnB_rho[3] = MFPRMB_rho/D;
 	
	P_1 = JAtot*(Z2-Z1)/(pAL-pAR); P_2 = JBtot*(Z2-Z1)/(pBL-pBR);
 	alphaABP = P_2/P_1; alphaABrho = (yB/yA)/(xB/xA);
 }
 BulkDensities(); PrintPlot(); JPalpha(); PrintJlist();
 
 sprintf(filename,"movie_%d.xyz",JobIndex);
 fpmovie=fopen(filename,"a+");
 Writemovie(fpmovie);
 fclose(fpmovie);

 sprintf(filename,"position.dat");
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity.dat");
 fpv=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%.12lf\t%.12lf\t%.12lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%.12lf\t%.12lf\t%.12lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 PrintInput();
 
 printf("\n");
 printf("MD simulation is finished, averaged over = %lf\n", averaged_over);
 printf("final packing fraction is %lf\n",packingfraction);
 
 printf("\n");

 printf("****************************** the end *******************************");

 return 0;
}
