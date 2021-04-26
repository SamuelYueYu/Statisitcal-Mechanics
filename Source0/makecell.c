/*******************************************************
* divide simulation box into small cells
* of size rcell (rcell > sigma)
*******************************************************/

#include <stdio.h>
#include <math.h>
#include "system.h"

FILE *fp; char filename[20];
double min(double a, double b, double c)
{
 if(a <= b && a <= c) return a;
 else if(b <= a && b <= c) return b;
 else return c;
}

void MakeCell()
{
 int i,j;
 rcellx = Lx/LCellx; rcelly = Ly/LCelly; rcellz = Lz/LCellz;

 if((rcellx<rc)||(rcelly<rc)||(rcellz<rc)) // rc is the max of sigmaAA, sigmaBB,sigmaCC, sigmaAB, sigmaBC,sigmaCA
 {
  if(rcellx == min(rcellx, rcelly, rcellz)){LCellx = (int)(Lx/rc); rcellx = Lx/LCellx; rcelly = Ly/LCellx; rcellz = Lz/LCellx;}
  else if(rcelly == min(rcellx, rcelly, rcellz)){LCelly = (int)(Ly/rc); rcellx = Lx/LCelly; rcelly = Ly/LCelly; rcellz = Lz/LCelly;}
  else{LCellz = (int)(Lz/rc); rcellx = Lx/LCellz; rcelly = Ly/LCellz; rcellz = Lz/LCellz;}
 }
 NumberOfNeighborCells = 27;
 NumberOfCells = RECTANGULAR(LCellx,LCelly,LCellz);
 
 printf("\n");
 printf("rcellx = %lf\trcelly = %lf\trcellz = %lf\tNCell = %d\n",rcellx,rcelly,rcellz,NumberOfCells);
 printf("\n");

 if(CellSwitch == 0)
 printf("cell list is NOT used\n");
 printf("\n");

 NeighborCell(); // assign each cell i of all its neighboring cells

 //group particles into cells

 for(i=0;i<(LCellz*(LCellz-1)+LCelly-1)*LCelly+LCellx;i++)
   HeadOfChain[i] = -1;  

 for(i=0;i<NumberOfParticles;i++){
   CellTrack[i].WhichCell = CellDetermine(i);

   CellTrack[i].Next = HeadOfChain[CellTrack[i].WhichCell];

   if(CellTrack[i].Next != -1) // avoid array[-1] which is segmentation fault
    CellTrack[CellTrack[i].Next].Prev = i;

   HeadOfChain[CellTrack[i].WhichCell] = i;
 }// end loop particle i

return;
} 
