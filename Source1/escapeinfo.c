/***************************************************
* find escape time of each  particle i
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"


void EscapeInfo(int i) //
{
 double xold,yold,zold;

 int cello,celln;

 int nx,ny,nz;
 double xa,xb;
 double ya,yb;
 double za,zb;
 double tx,ty,tz;
 int onboundx,onboundy,onboundz;

if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
{
 //indicator of particles on boundary of cells
 onboundx = 0;
 onboundy = 0;
 onboundz = 0;

 nx = (int)(position[i].x/rcellx);
 ny = (int)(position[i].y/rcelly);
 nz = (int)(position[i].z/rcellz);

//avoid nx = LCellx, ny = LCelly, or nz = LCellz
 if(nx == -1) nx = 0;
 if(ny == -1) ny = 0;
 if(nz == -1) nz = 0;

 if(nx == LCellx) nx = LCellx-1;
 if(ny == LCelly) ny = LCelly-1;
 if(nz == LCellz) nz = LCellz-1;
 
 if(fabs(position[i].x/rcellx - (double)(nx)) < tolerance)
 {
  if(velocity[i].x < 0.)
  {
    nx = (nx-1+LCellx)%LCellx;
    onboundx = 1;
  }
 }
 else if(fabs(position[i].x/rcellx - (double)(nx+1)) < tolerance)
 {
  if(velocity[i].x > 0.)
  {
    nx = (nx+1)%LCellx;
    onboundx = 1;
  }
 }
 
 if(fabs(position[i].y/rcelly - (double)(ny)) < tolerance)
 {
  if(velocity[i].y < 0.)
  {
    ny = (ny-1+LCelly)%LCelly;
    onboundy = 1;
  }
 }
 else if(fabs(position[i].y/rcelly - (double)(ny+1)) < tolerance)
 {
  if(velocity[i].y > 0.)
  {
    ny = (ny+1)%LCelly;
    onboundy = 1;
  }
 }
 
 if(fabs(position[i].z/rcellz - (double)(nz)) < tolerance)
 {
  if(velocity[i].z < 0.)
  {
    nz = (nz-1+LCellz)%LCellz;
    onboundz = 1;
  }
 }
 else if(fabs(position[i].z/rcellz - (double)(nz+1)) < tolerance)
 {
  if(velocity[i].z > 0.)
  {
    nz = (nz+1)%LCellz;
    onboundz = 1;
  }
 }
 
 xa = nx*rcellx;
 xb = (nx+1)*rcellx;
 ya = ny*rcelly;
 yb = (ny+1)*rcelly;
 za = nz*rcellz;
 zb = (nz+1)*rcellz;
 
 if(velocity[i].x > 0.) //->
 {
   if(onboundx == 1)
    tx = rcellx/velocity[i].x;
   else
    tx = (xb - position[i].x)/velocity[i].x;
 }
 else if(velocity[i].x < 0.) // <-
 {
   if(onboundx == 1)
    tx = -rcellx/velocity[i].x;
   else
    tx = (xa - position[i].x)/velocity[i].x;
 }
 else{tx = TimeBig;}
 
 if(velocity[i].y > 0.) //->
 {
   if(onboundy == 1)
    ty = rcelly/velocity[i].y;
   else
    ty = (yb - position[i].y)/velocity[i].y;
 }
 else if(velocity[i].y < 0.)// <-
 {
   if(onboundy == 1)
    ty = -rcelly/velocity[i].y;
   else
    ty = (ya - position[i].y)/velocity[i].y;
 }
 else{ty = TimeBig;}
 
 if(velocity[i].z > 0.) //->
 {
   if(onboundz == 1)
    tz = rcellz/velocity[i].z;
   else
    tz = (zb - position[i].z)/velocity[i].z;
 }
 else if(velocity[i].z < 0.)// <-
 {
   if(onboundz == 1)
    tz = -rcellz/velocity[i].z;
   else
    tz = (za - position[i].z)/velocity[i].z;
 }
 else{tz = TimeBig;}

 EscapeTime[i] = MIN(tx,ty);
 EscapeTime[i] = MIN(EscapeTime[i],tz);
} // end if cell list is used

 return;
}
