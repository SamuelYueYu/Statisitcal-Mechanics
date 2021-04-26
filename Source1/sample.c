/***************************************************
* sampling
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"


/*****************
void Sample(void) // now useless
{
 int i,j;
 int jbegin,jend,jlist;
 double rij;
 double wij;
 //double Kinstant;
 double w; //virial
 double u;

 FILE *fp;

  Kinstant = 0.;
  u = 0.;
  for(i=0;i<NumberOfParticles-1;i++) 
  {
   Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));

   jbegin = point[i];
   jend = point[i+1]-1;
   if(jbegin<=jend)
   {
    for(jlist=jbegin;jlist<=jend;jlist++)
    {
     j = VerletList[jlist];
     
     rij = sqrt(Distance(i,j)); // avoid sqrt if possible
     
     u += Potential(identity[i],identity[j],rij);
     
    }//end loop jlist

   }//endif jbegin <= jend
  }//end loop i
 
  Kinstant *= 0.5; // instantenous kinetic energy
  Upotential = u;

   // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
  Tinstant = 2.0*Kinstant/Nf/kB;

 return;
}
*****************/


void Kinetic(void)
{
 int i;

 Kinstant = 0.;
 for(i=0;i<NumberOfParticles;i++){
    Kinstant += mass[i]*(SQR(velocity[i].x) + SQR(velocity[i].y) + SQR(velocity[i].z));
    if(fabs(Kinstant) == nan("")){printf("i = %d, velocity[i].x = %lf, velocity[i].y = %lf, velocity[i].z = %lf\n", i, velocity[i].x, velocity[i].y, velocity[i].z);}
 }
 
 Kinstant *= 0.5; // instantenous kinetic energy
   // 0.5*kT*Nf = K = 0.5* sum_mv^2, Nf = 3N-3
 //Tinstant = 2.0*Kinstant/(3.*NumberOfParticles-3.)/kB;
 Tinstant = 2.0*Kinstant/Nf/kB;

 return;
}
