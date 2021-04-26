/****************************************************************
* the array NeighborCellList[i][j], i = 0, ..., NumberOfCells-1, 
* j = 0, 1, ..., 26, stores the neighbors of 
* the cell i 
********************************************************************/
#include <stdio.h>
#include "system.h"

/******************************************************
* Cells are stored in a 1d array. this one 
* gets 3d coordinates from 1d chain on a SC lattice
*******************************************************/
void CoordinateTrans(int i)
{

 int quotient;
 int d;
 int LCell[2];
 LCell[0] = LCelly; LCell[1] = LCellz;

 quotient=i;

 for(d=0;d<3;d++){

	if(d == 2) CellCoordinate[i][d]=quotient;
	else{
			CellCoordinate[i][d]=quotient%LCell[d];
        		quotient=(quotient-CellCoordinate[i][d])/LCell[d];
	}
 }
//i = (x[2]*L+x[1])*L+x[0]
 return;
}

void NeighborCell(void) // Periodic Boundary Condition is used
{
        int i;

        for(i=0;i<(LCellz*(LCellz-1)+LCelly-1)*LCelly+LCellx;i++){
           CoordinateTrans(i); // get 3d coordinate of current cell in CellCoordinate[i][d]

           NeighborCellList[i][0]=i;                                                                            // (0, 0, 0)
           //NeighborCellList[i][0]=(CellCoordinate[i][2]*LCell+CellCoordinate[i][1])*LCell+CellCoordinate[i][0];

           // nearest neighbor nn
           NeighborCellList[i][1]=(CellCoordinate[i][2]*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]+1)%LCellx;   // (0, 0, +1)
           NeighborCellList[i][2]=(CellCoordinate[i][2]*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (0, 0, -1)
           NeighborCellList[i][3]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+CellCoordinate[i][0];   // (0, +1, 0)
           NeighborCellList[i][4]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+CellCoordinate[i][0]; // (0, -1, 0)
           NeighborCellList[i][5]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+CellCoordinate[i][0];   // (+1, 0, 0)
           NeighborCellList[i][6]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+CellCoordinate[i][0]; // (-1, 0, 0)

           // next nearest neighbor  nnn
           NeighborCellList[i][7]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx;  // (0, +1, +1)
           NeighborCellList[i][8]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx;  // (0, +1, -1)
           NeighborCellList[i][9]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx;  // (0, -1, +1)
           NeighborCellList[i][10]=(CellCoordinate[i][2]*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (0, -1, -1)
           NeighborCellList[i][11]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (+1, 0, +1)
           NeighborCellList[i][12]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (+1, 0, -1)
           NeighborCellList[i][13]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (-1, 0, +1)
           NeighborCellList[i][14]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+CellCoordinate[i][1])*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (-1, 0, -1)
           NeighborCellList[i][15]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+CellCoordinate[i][0]; // (+1, +1, 0)
           NeighborCellList[i][16]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+CellCoordinate[i][0]; // (+1, -1, 0)
           NeighborCellList[i][17]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+CellCoordinate[i][0]; // (-1, +1, 0)
           NeighborCellList[i][18]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+CellCoordinate[i][0]; // (-1, -1, 0)

           // next next nearest neighbor  nnnn    
           NeighborCellList[i][19]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (+1, +1, +1)
           NeighborCellList[i][20]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (+1, +1, -1)
           NeighborCellList[i][21]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (+1, -1, +1)
           NeighborCellList[i][22]=(((CellCoordinate[i][2]+1)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (+1, -1, -1)
           NeighborCellList[i][23]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (-1, +1, +1)
           NeighborCellList[i][24]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]+1)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (-1, +1, -1)
           NeighborCellList[i][25]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]+1)%LCellx; // (-1, -1, +1)
           NeighborCellList[i][26]=(((CellCoordinate[i][2]-1+LCellz)%LCellz)*LCellz+(CellCoordinate[i][1]-1+LCelly)%LCelly)*LCelly+(CellCoordinate[i][0]-1+LCellx)%LCellx; // (-1, -1, -1)
        }
 return;
}

