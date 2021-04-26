#include <stdio.h>
#include <math.h>
#include "system.h"

double gTerm(int N, int Nlist[], double rho){
	int i, j;
	int iCell, jCell;
	VECTOR dr;
	double rij, gAA = 0., gBB = 0., gAB = 0.;
	bool AA = false, BB = false, AB = false;
	
	if(NumberOfCells <= NumberOfNeighborCells || CellSwitch == 0){ // not use cell list
		for(i=0; i<N-1; i++)
		for(j=i+1; j<N; j++)
		{
			dr.x = position[Nlist[i]].x - position[Nlist[j]].x;
			MinimumImage(&(dr.x),Lx);
							 
			dr.y = position[Nlist[i]].y - position[Nlist[j]].y;
			MinimumImage(&(dr.y),Ly);
			
			dr.z = position[Nlist[i]].z - position[Nlist[j]].z;
							                                                                                    
			rij = sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
									 
			if(sigma[Nlist[i]] + sigma[Nlist[j]] == 2*sigmaA && rij >= sigmaA && rij < 1.01*sigmaA){gAA += 1.0/N; AA = true;}
			else if(sigma[Nlist[i]] + sigma[Nlist[j]] == 2*sigmaB && rij >= sigmaB && rij < 1.01*sigmaB){gBB += 1.0/N; BB = true;}
			else if(sigma[Nlist[i]] + sigma[Nlist[j]] == 2*sigmaAB && rij >= sigmaAB && rij < 1.01*sigmaAB){gAB += 1.0/N; AB = true;}
		}
	}
	else // use cell list when CellSwitch == 1 && NumberOfCells>27
	{
		for(i=0; i<N; i++){
			iCell = CellTrack[Nlist[i]].WhichCell;
			 
 			for(int cell=0;cell<NumberOfNeighborCells;cell++){
  				jCell = NeighborCellList[iCell][cell];
  				j=HeadOfChain[jCell];
  				
  				while(j != -1){
  					if(j != Nlist[i] && ((position[j].z >= 0. && position[j].z < Z1) || (position[j].z >= Z2 && position[j].z < Z3))){
  						dr.x = position[Nlist[i]].x - position[j].x;
						MinimumImage(&(dr.x),Lx);
								 
						dr.y = position[Nlist[i]].y - position[j].y;
						MinimumImage(&(dr.y),Ly);
						
						dr.z = position[Nlist[i]].z - position[j].z;
										                                                                                    
						rij = sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
						
						if(sigma[Nlist[i]] + sigma[j] == 2*sigmaA && rij >= sigmaA && rij < 1.01*sigmaA){gAA += 1.0/N; AA = true;}
						else if(sigma[Nlist[i]] + sigma[j] == 2*sigmaB && rij >= sigmaB && rij < 1.01*sigmaB){gBB += 1.0/N; BB = true;}
						else if(sigma[Nlist[i]] + sigma[j] == 2*sigmaAB && rij >= sigmaAB && rij < 1.01*sigmaAB){gAB += 1.0/N; AB = true;}
 					} //end if j is not i
					
    				j = CellTrack[j].Next;
  				}//end while j != -1
 			}//end loop cell
		}
	}//end cell list is used
	if(AA){gAA /= 2*(CUBIC(1.01)-1); return gAA;}
	if(BB){gBB /= 2*(CUBIC(1.01)-1); return gBB;}
	if(AB){gAB /= 2*(CUBIC(1.01)-1);}
	return gAB;
}
