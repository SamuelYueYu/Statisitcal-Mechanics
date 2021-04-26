/***************************************************
* return distance r, u(r), du/dr, u_tail, w_tail
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

double Distance(int i,int j) // return rij^2
{
  double r2,dx,dy,dz;
  
  // miminum image convention under periodic boundary condition
  dx = position[i].x-position[j].x;
  dx = dx - Lx*round(dx/Lx);

  dy = position[i].y-position[j].y;
  dy = dy - Ly*round(dy/Ly);

  dz = position[i].z-position[j].z;
  dz = dz - Lz*round(dz/Lz);
 
  r2 = SQR(dx)+SQR(dy)+SQR(dz); // reduced unit
  //r = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

  return r2;
}

double SigmaIJ(int iID,int jID) // input identity[i]!
{
 double sigmaij;

 if(iID == jID)
 {
  if(iID == 1) // A-A
  sigmaij = sigmaA;
  else if (iID == 2) // B-B
  sigmaij = sigmaB;
 }
 else
  sigmaij = sigmaAB;

 return sigmaij;
}

