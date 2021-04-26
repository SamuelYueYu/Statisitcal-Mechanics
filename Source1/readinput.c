/*********************************************************************
 * input simulation parameters from file "input"
 *********************************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void ReadInput(void)
{
 int i;
 FILE *fp;
 
 fp=fopen("input","r"); 
 fscanf(fp,"%d",&EnsembleType);
 fscanf(fp,"%d %d %d %d",&rInitialType,&PackType,&NumberOfLatticeSites,&vInitialType);
 fscanf(fp,"%d %d",&PotentialType,&MemType);
 fscanf(fp,"%d %d %d",&NumberOfParticles,&NA,&NB);
 fscanf(fp,"%lf %lf",&massA,&massB);
 fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaAB);
 fscanf(fp,"%lf",&T);
 fscanf(fp,"%lf",&V);
 fscanf(fp,"%lf %lf %lf",&Lx,&Ly,&Lz);
 fscanf(fp,"%lf %lf %lf",&X,&Y,&Z);
 fscanf(fp,"%lf %lf %lf",&Z1,&Z2,&Z3);
 fscanf(fp,"%d %d %d",&M,&mx,&my);
 fscanf(fp,"%lf %lf %lf",&D,&Sigma,&p);
 fscanf(fp,"%lf %lf",&mu1,&mu2);
 fscanf(fp,"%d %lf",&coln,&SCT);
 fscanf(fp,"%lf %lf",&N_insert_0_1,&N_insert_trial_0_1);
 fscanf(fp,"%lf %lf",&N_insert_Z2_1,&N_insert_trial_Z2_1);
 fscanf(fp,"%lf %lf",&N_delete_0_1,&N_delete_trial_0_1);
 fscanf(fp,"%lf %lf",&N_delete_Z2_1,&N_delete_trial_Z2_1);
 fscanf(fp,"%lf %lf",&N_insert_0_2,&N_insert_trial_0_2);
 fscanf(fp,"%lf %lf",&N_insert_Z2_2,&N_insert_trial_Z2_2);
 fscanf(fp,"%lf %lf",&N_delete_0_2,&N_delete_trial_0_2);
 fscanf(fp,"%lf %lf",&N_delete_Z2_2,&N_delete_trial_Z2_2);
 fscanf(fp,"%d %d %d",&LCellx,&LCelly,&LCellz);
 fscanf(fp,"%d",&NumberOfInitialSteps);
 fscanf(fp,"%d %d",&SampleMultiplier,&MovieMultiplier);
 fscanf(fp,"%lf",&randomseed);
 fscanf(fp,"%d %d",&CellSwitch,&RhoSwitch);
 fscanf(fp,"%lf",&timewindow);
 //fscanf(fp,"%d",&VType);
 fclose(fp);

 rc = MAX(sigmaA,sigmaB);
 rc = MAX(sigmaAB,rc);

 TimeBig = 1.0E10;
 tolerance = 1.0E-10;
 timetol = 1.0E-9;

 fA = 1.*NA/NumberOfParticles;
 fB = 1.*NB/NumberOfParticles;

 Nf = 3*NumberOfParticles-3;

 if(EnsembleType == 0) printf("NVE ensemble\n");

 T0 = T; // initial temperature
 kB = 1.; lambda = 1.;
 beta = 1.0 / T / kB;
 
 pLeftVirial = 0; pRightVirial = 0;
 pLeftGr = 0; pRightGr = 0;
 pLeftCS = 0; pRightCS = 0;

 rho = NumberOfParticles/RECTANGULAR(X,Y,Z)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 rhoA = NA/RECTANGULAR(X,Y,Z)*RECTANGULAR(sigmaA,sigmaA,sigmaA);
 rhoB = NB/RECTANGULAR(X,Y,Z)*RECTANGULAR(sigmaB,sigmaB,sigmaB);

 packingfraction = M_PI/6.*rho*(fA*RECTANGULAR(sigmaA,sigmaA,sigmaA)+fB*RECTANGULAR(sigmaB,sigmaB,sigmaB));
 
 M1_1 = 0; M2_1 = 0; M1_2 = 0; M2_2 = 0;
 averaged_over = 0;

 //V = NumberOfParticles / rho;

 printf("\n");
 printf("Number of equilibrium MD steps: %d\n",NumberOfInitialSteps);
 printf("Sample multiplier (sampling frequency): %d\n",SampleMultiplier);
 printf("Movie multiplier (draw snapshot frequency): %d\n",MovieMultiplier);
 printf("\n");

 printf("\n");
 printf("T = %lf\tbeta = %lf\n",T,beta);
 if(VType == 0)
 {
 printf("rho = %lf\n",rho);
 printf("packing fraction = %lf\n",packingfraction);
 printf("V = %lf\n",V);
 printf("Lx = %lf\tLy = %lf\tLz = %lf\n",Lx, Ly, Lz);
 printf("X = %lf\tY = %lf\tZ = %lf\n", X, Y, Z);
 printf("Z1 = %lf\tZ2 = %lf\tZ3 = %lf\n", Z1, Z2, Z3);
 }

 printf("\n");
 printf("N = %d\tNA:NB = %d:%d\tfA:fB = %lf:%lf\n",NumberOfParticles,NA,NB,fA,fB);
 printf("massA = %lf\tmassB = %lf\n",massA,massB);
 printf("sigmaA = %lf\tsigmaB = %lf\tsigmaAB = %lf\n",sigmaA,sigmaB,sigmaAB);
 printf("M = %d\tmx = %d\tmy = %d\tD = %lf\tMemType = %d\n", M, mx, my, D, MemType);
 printf("\n");

 printf("simulation time window = %lf\n",timewindow);
 printf("\n");
 
 return;
}
