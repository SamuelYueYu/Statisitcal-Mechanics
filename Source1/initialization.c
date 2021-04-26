/***************************************************
* initialize particle parameters
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Initialization(int job)
{ 
 char filename[20];
 int i;
 double I;
 int nx,ny,nz;
// double Mcom; // total mass
// VECTOR Vcom;
 double scale;
 double particlepercell;
 //double Kinstant,scale,Tinstant;

 double goldenratio;
 double icobond;
 double sideradius;
 double radiusside;
 int checkboard[10000];

 int ishift;

 FILE *fp;
 FILE *fppore;

 sideradius = 1.05146223;
 radiusside = 0.95105652;
 goldenratio = (sqrt(5.)+1.)/2.;
// icobond = radiusside*sigmaB/1.8;
 icobond = sigmaB/1.9022;
 
 if(PackType == 4) // AB13
 printf("icosahedron bond = %lf\n",icobond);
  

 for(i=0;i<NumberOfParticles;i++)
 {
  if(i<NA)
  {
  sigma[i] = sigmaA;
  mass[i] = massA;
  identity[i] = 1;
  }
  else
  {
  sigma[i] = sigmaB;
  mass[i] = massB;
  identity[i] = 2;
  }
 }

 fppore = fopen("pore_size.dat", "r");
 for(int j = 0; j < my; j++){
	 for(int i = 0; i < mx; i++){
		fscanf(fppore,"%lf",&pore_position[i + mx*j].d);
		pore_position[i + mx*j].x = (i+0.5)*Lx/mx;
		pore_position[i + mx*j].y = (j+0.5)*Ly/my;
		JA[i + mx*j] = 0; JB[i + mx*j] = 0;
		JA[M + i + mx*j] = 0; JB[M + i + mx*j] = 0;
	 }
 }
 fclose(fppore);

/****************  position **********/
 if(rInitialType == 0) // read from file
 {
  fp = fopen("position.dat","r");
//  sprintf(filename,"position_%d",job);
//  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++){
  	fscanf(fp,"%lf\t%lf\t%lf\n",&position[i].x,&position[i].y,&position[i].z);
  }
  fclose(fp);
 }
 else if(rInitialType == 1) // on lattice 
 {
  /***********************place particle on lattice**********************************/
  if(PackType == 1) // simple cubic
  {
   //latticeconst = pow(4.0*V/NumberOfLatticeSites,1.0/3.0);
   particlepercell = 1.;
   latticeconst = pow(X*Y*(Z - MAX(sigmaA, sigmaB)/2 - .01)/NumberOfParticles,1.0/3.0);
   printf("SC with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0; 
   do
   {
    for(nx=0;nx*latticeconst<X;nx++){
    	for(ny=0;ny*latticeconst<Y;ny++){
    		for(nz=0;nz*latticeconst<Z;nz++){
    			position[i].x = nx*latticeconst;
     			position[i].y = ny*latticeconst;
     			position[i].z = .01+MAX(sigmaA, sigmaB)/2+nz*latticeconst;
     			i++;	
     			if(i >= NumberOfParticles) break;
			}
			if(i >= NumberOfParticles) break;
		}
		if(i >= NumberOfParticles) break;
	}
   }//end do
   while(i<NumberOfParticles);
 }// end if sc
  if(PackType == 2) // bcc
  {
   particlepercell = 2.;
   latticeconst = pow(2.0/rho,1.0/3.0);
   printf("BCC with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0;	
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites/2,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/2,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/2,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst;
     position[i].y = ny*latticeconst;
     position[i].z = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
    for(nx=0;nx<pow(NumberOfLatticeSites/2,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/2,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/2,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst+0.5*latticeconst;
     position[i].y = ny*latticeconst+0.5*latticeconst;
     position[i].z = nz*latticeconst+0.5*latticeconst;
     i++;
     if(i >= NumberOfParticles) break;
    }
   } //end do
   while(i<NumberOfParticles);
  }//end if bcc
  if(PackType == 3) // fcc
  {
   //latticeconst = pow(4.0*V/NumberOfLatticeSites,1.0/3.0);
   
   particlepercell = 4.;
   latticeconst = pow(4.0/rho,1.0/3.0);

   printf("FCC with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0;
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites/4,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/4,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/4,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst;
     position[i].y = ny*latticeconst;
     position[i].z = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
   
    for(nx=0;nx<pow(NumberOfLatticeSites/4,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/4,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/4,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst+0.5*latticeconst;
     position[i].y = ny*latticeconst+0.5*latticeconst;
     position[i].z = nz*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
 
    for(nx=0;nx<pow(NumberOfLatticeSites/4,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/4,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/4,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst+0.5*latticeconst;
     position[i].y = ny*latticeconst;
     position[i].z = nz*latticeconst+0.5*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }

    for(nx=0;nx<pow(NumberOfLatticeSites/4,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites/4,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites/4,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst;
     position[i].y = ny*latticeconst+0.5*latticeconst;
     position[i].z = nz*latticeconst+0.5*latticeconst;
     i++;	
     if(i >= NumberOfParticles) break;
    }
   }//end do
   while(i<NumberOfParticles);
  }// end if fcc
  if(PackType == 4) // AB13
  {
   //latticeconst = pow(4.0*V/NumberOfLatticeSites,1.0/3.0);
   particlepercell = 1.;
   latticeconst = pow(1.0/rhoA,1.0/3.0);

   printf("AB13 with Nc = %d and a = %lf\n",NumberOfLatticeSites,latticeconst);
   i=0;
   do
   {
    for(nx=0;nx<pow(NumberOfLatticeSites,1.0/3.0)*X/pow(X*Y*Z,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites,1.0/3.0)*Y/pow(X*Y*Z,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites,1.0/3.0)*Z/pow(X*Y*Z,1.0/3.0);nz++) 	
    {
     position[i].x = nx*latticeconst;
     position[i].y = ny*latticeconst;
     position[i].z = nz*latticeconst;
     i++;	
     if(i >= NA) break;
    }
   }//end do
   while(i<NA);
   
   i=NA;
   do // put central B
   {
    for(nx=0;nx<pow(NumberOfLatticeSites,1.0/3.0)*X/pow(X*Y*Z,1.0/3.0);nx++) 	
    for(ny=0;ny<pow(NumberOfLatticeSites,1.0/3.0)*Y/pow(X*Y*Z,1.0/3.0);ny++) 	
    for(nz=0;nz<pow(NumberOfLatticeSites,1.0/3.0)*Z/pow(X*Y*Z,1.0/3.0);nz++) 	
    {
     if((nx+ny+nz)%2 == 0)
      checkboard[i] = 0;
     else
      checkboard[i] = 1;
      
     position[i].x = nx*latticeconst+0.5*latticeconst;
     position[i].y = ny*latticeconst+0.5*latticeconst;
     position[i].z = nz*latticeconst+0.5*latticeconst;
     i++;	
     if(i >= 2*NA) break;
    }
   }//end do
   while(i<2*NA);

   for(i=NA;i<2*NA;i++) // NA, NA+1, ... 2NA-1
   {
   ishift = 2*NA + (i-NA)*12;
   if(checkboard[i] == 0)
   {
   position[ishift+0].x = position[i].x + 0.;
   position[ishift+0].y = position[i].y + 1.*icobond;
   position[ishift+0].z = position[i].z + goldenratio*icobond;
   
   position[ishift+1].x = position[i].x + 0.;
   position[ishift+1].y = position[i].y - 1.*icobond;
   position[ishift+1].z = position[i].z + goldenratio*icobond;
   
   position[ishift+2].x = position[i].x + 0.;
   position[ishift+2].y = position[i].y + 1.*icobond;
   position[ishift+2].z = position[i].z - goldenratio*icobond;
   
   position[ishift+3].x = position[i].x + 0.;
   position[ishift+3].y = position[i].y - 1.*icobond;
   position[ishift+3].z = position[i].z - goldenratio*icobond;

   position[ishift+4].x = position[i].x + 1.*icobond;
   position[ishift+4].y = position[i].y + goldenratio*icobond;
   position[ishift+4].z = position[i].z + 0.;
   
   position[ishift+5].x = position[i].x - 1.*icobond;
   position[ishift+5].y = position[i].y + goldenratio*icobond;
   position[ishift+5].z = position[i].z + 0.;
   
   position[ishift+6].x = position[i].x + 1.*icobond;
   position[ishift+6].y = position[i].y - goldenratio*icobond;
   position[ishift+6].z = position[i].z + 0.;
   
   position[ishift+7].x = position[i].x - 1.*icobond;
   position[ishift+7].y = position[i].y - goldenratio*icobond;
   position[ishift+7].z = position[i].z + 0.;
   
   position[ishift+8].x = position[i].x + goldenratio*icobond;
   position[ishift+8].y = position[i].y + 0.;
   position[ishift+8].z = position[i].z + 1.*icobond;
   
   position[ishift+9].x = position[i].x + goldenratio*icobond;
   position[ishift+9].y = position[i].y + 0.;
   position[ishift+9].z = position[i].z - 1.*icobond;
   
   position[ishift+10].x = position[i].x - goldenratio*icobond;
   position[ishift+10].y = position[i].y + 0.;
   position[ishift+10].z = position[i].z + 1.*icobond;
   
   position[ishift+11].x = position[i].x - goldenratio*icobond;
   position[ishift+11].y = position[i].y + 0.;
   position[ishift+11].z = position[i].z - 1.*icobond;
   }
   else //rotate 90
   {
   position[ishift+0].y = position[i].y + 0.;
   position[ishift+0].x = position[i].x - 1.*icobond;
   position[ishift+0].z = position[i].z + goldenratio*icobond;
   
   position[ishift+1].y = position[i].y + 0.;
   position[ishift+1].x = position[i].x + 1.*icobond;
   position[ishift+1].z = position[i].z + goldenratio*icobond;
   
   position[ishift+2].y = position[i].y + 0.;
   position[ishift+2].x = position[i].x - 1.*icobond;
   position[ishift+2].z = position[i].z - goldenratio*icobond;
   
   position[ishift+3].y = position[i].y + 0.;
   position[ishift+3].x = position[i].x + 1.*icobond;
   position[ishift+3].z = position[i].z - goldenratio*icobond;

   position[ishift+4].y = position[i].y + 1.*icobond;
   position[ishift+4].x = position[i].x - goldenratio*icobond;
   position[ishift+4].z = position[i].z + 0.;
   
   position[ishift+5].y = position[i].y - 1.*icobond;
   position[ishift+5].x = position[i].x - goldenratio*icobond;
   position[ishift+5].z = position[i].z + 0.;
   
   position[ishift+6].y = position[i].y + 1.*icobond;
   position[ishift+6].x = position[i].x + goldenratio*icobond;
   position[ishift+6].z = position[i].z + 0.;
   
   position[ishift+7].y = position[i].y - 1.*icobond;
   position[ishift+7].x = position[i].x + goldenratio*icobond;
   position[ishift+7].z = position[i].z + 0.;
   
   position[ishift+8].y = position[i].y + goldenratio*icobond;
   position[ishift+8].x = position[i].x - 0.;
   position[ishift+8].z = position[i].z + 1.*icobond;
   
   position[ishift+9].y = position[i].y + goldenratio*icobond;
   position[ishift+9].x = position[i].x - 0.;
   position[ishift+9].z = position[i].z - 1.*icobond;
   
   position[ishift+10].y = position[i].y - goldenratio*icobond;
   position[ishift+10].x = position[i].x - 0.;
   position[ishift+10].z = position[i].z + 1.*icobond;
   
   position[ishift+11].y = position[i].y - goldenratio*icobond;
   position[ishift+11].x = position[i].x - 0.;
   position[ishift+11].z = position[i].z - 1.*icobond;
   }// end rotate 90
   } // end loop icosahedron

  
  }// end if AB13
  /**********************************************************************************/

  for(i=0;i<NumberOfParticles;i++) 
  {
 //  position_old[i].x = position[i].x - velocity[i].x*dt;
 //  position_old[i].y = position[i].y - velocity[i].y*dt;
 //  position_old[i].z = position[i].z - velocity[i].z*dt;
    
    //PBC
    if(position_old[i].x < 0.0)
     position_old[i].x += Lx;
    else if(position_old[i].x >= Lx)
     position_old[i].x -= Lx; 
    
    if(position_old[i].y < 0.0)
     position_old[i].y += Ly;
    else if(position_old[i].y >= Ly)
     position_old[i].y -= Ly; 
    
    if(position_old[i].z < 0.0)
     position_old[i].z += Lz;
    else if(position_old[i].z >= Lz)
     position_old[i].z -= Lz; 
  }
 }// end if on lattice 
  
  
/********************************** velocity ********************************************/

if(vInitialType == 0) // read from file
{
  fp = fopen("velocity.dat","r");
//  sprintf(filename,"velocity_%d",job);
//  fp = fopen(filename,"r");
  for(i=0;i<NumberOfParticles;i++){
  	fscanf(fp,"%lf\t%lf\t%lf\n",&velocity[i].x,&velocity[i].y,&velocity[i].z);
  }
  //fscanf(fp,"i = %*d\tvx = %lf\tvy = %lf\tvz = %lf\n",&velocity[i].x,&velocity[i].y,&velocity[i].z);
  fclose(fp);
  
  // velocity check
  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;
  Kinstant = 0.;
  fp = fopen("vmaxwell.txt","w");
  for(i=0;i<NumberOfParticles;i++) 
  {
   Mcom += mass[i];
   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;
   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
   fprintf(fp,"%lf\n",sqrt(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z)));
  }
  fclose(fp);

   Kinstant *= 0.5;
   Tinstant = 2.0*Kinstant/Nf/kB;
   printf("center of mass velosity (%lf, %lf, %lf) and instantenous T = %lf\n",Vcom.x/Mcom,Vcom.y/Mcom,Vcom.z/Mcom,Tinstant);
}
else if(vInitialType == 1)
{
  /********************random velocity**************************************/
  Mcom = 0.;
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;

  for(i=0;i<NumberOfParticles;i++)
  {
   velocity[i].x = BoxMuller(0.,1.);
   velocity[i].y = BoxMuller(0.,1.);
   velocity[i].z = BoxMuller(0.,1.);

   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;

   Mcom += mass[i];
  }

  Vcom.x /= Mcom;
  Vcom.y /= Mcom;
  Vcom.z /= Mcom;

  Kinstant = 0.;
  for(i=0;i<NumberOfParticles;i++) 
  {
   velocity[i].x -= Vcom.x;
   velocity[i].y -= Vcom.y;
   velocity[i].z -= Vcom.z;

   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
  }
 
   Kinstant *= 0.5; // instantenous kinetic energy

   // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
   Tinstant = 2.0*Kinstant/Nf/kB;

   scale = sqrt(T/Tinstant);

   for(i=0;i<NumberOfParticles;i++) 
   {
    velocity[i].x *= scale;
    velocity[i].y *= scale;
    velocity[i].z *= scale;
   }
  /**********************************************************/
 // check velocity
  Vcom.x = 0.;
  Vcom.y = 0.;
  Vcom.z = 0.;
  Kinstant = 0.;
  fp = fopen("vmaxwell.txt","w");
  for(i=0;i<NumberOfParticles;i++)
  {
   Vcom.x += mass[i]*velocity[i].x;
   Vcom.y += mass[i]*velocity[i].y;
   Vcom.z += mass[i]*velocity[i].z;
   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));

   fprintf(fp,"%lf\n",sqrt(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z)));
  }
  fclose(fp);

  Vcom.x /= Mcom;
  Vcom.y /= Mcom;
  Vcom.z /= Mcom;
  Kinstant *= 0.5;
  Tinstant = 2.0*Kinstant/Nf/kB;
  printf("center of mass velocity (%lf, %lf, %lf) and instantenous T = %lf\n",Vcom.x,Vcom.y,Vcom.z,Tinstant);
}//end if maxwell velocity

return;
}
