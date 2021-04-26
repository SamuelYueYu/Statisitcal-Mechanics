/***********************************************************
* determine which cell does particle i with
* coordinate x,y,z belong to
* be careful when particles are on the boundary of cells
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include <math.h>

int CellDetermine(int i)
{
 int ncell;
 int nx,ny,nz;

 nx = (int)(position[i].x/rcellx);
 ny = (int)(position[i].y/rcelly);
 nz = (int)(position[i].z/rcellz);
 
 //avoid n = LCell
 if(nx == -1) nx = 0;
 if(ny == -1) ny = 0;
 if(nz == -1) nz = 0;

 if(nx == LCellx) nx = LCellx-1;
 if(ny == LCelly) ny = LCelly-1;
 if(nz == LCellz) nz = LCellz-1;


 if(fabs(position[i].x/rcellx - (double)(nx)) < tolerance)
 {
  if(velocity[i].x < 0.)
    nx = (nx-1+LCellx)%LCellx;
 }
 else if(fabs(position[i].x/rcellx - (double)(nx+1)) < tolerance)
 {
  if(velocity[i].x > 0.)
    nx = (nx+1)%LCellx;
 }
 
 if(fabs(position[i].y/rcelly - (double)(ny)) < tolerance)
 {
  if(velocity[i].y < 0.)
    ny = (ny-1+LCelly)%LCelly;
 }
 else if(fabs(position[i].y/rcelly - (double)(ny+1)) < tolerance)
 {
  if(velocity[i].y > 0.)
    ny = (ny+1)%LCelly;
 }
 
 if(fabs(position[i].z/rcellz - (double)(nz)) < tolerance)
 {
  if(velocity[i].z < 0.)
    nz = (nz-1+LCellz)%LCellz;
 }
 else if(fabs(position[i].z/rcellz - (double)(nz+1)) < tolerance)
 {
  if(velocity[i].z > 0.)
    nz = (nz+1)%LCellz;
 }
 
 ncell = (nz*LCellz+ny)*LCelly+nx;

 return ncell;
}
