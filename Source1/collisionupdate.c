/**********************************************************
* update collision info involving the transferred particle
* do not need to loop all the 27 cells
* the new 9 cells are enough, ideally
* need to update both i and j collision time
***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
void InsdelCollisionUpdate(int i, int iCell){
 int j, jCell;
 VECTOR rij,vij;
 double bij; // rij*vij
 double rij2,vij2,sigmaij2; // square
 double discr; // discriminate b^2 - 4ac
 double tij;
 int l,m,n; // pbc box index
 double xtemp,ytemp,ztemp;
 
 CollisionTime[i] = TimeBig;
 CollisionPartner[i] = MAX_NumberOfParticles;
 
 if(NumberOfCells <= NumberOfNeighborCells || CellSwitch == 0){
 	for(j=0;j<NumberOfParticles;j++) {
		if(j != i){
	 		vij.x = velocity[i].x - velocity[j].x;
  			vij.y = velocity[i].y - velocity[j].y;
  			vij.z = velocity[i].z - velocity[j].z;
  			
  			for(l=-1;l<=1;l++)
  			for(m=-1;m<=1;m++)
  			for(n=-1;n<=1;n++){
   				xtemp = position[j].x + l*Lx;
   				ytemp = position[j].y + m*Ly;
   				ztemp = position[j].z + n*Lz;
				
   				rij.x = position[i].x - xtemp;
   				rij.y = position[i].y - ytemp;
   				rij.z = position[i].z - ztemp;
  				
   				bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z;
   				
   				if(bij < 0.){
					rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   					vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   					sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   					discr = SQR(bij)- vij2*(rij2-sigmaij2);
  					
   					if(discr > 0.){
    					tij = (-bij - sqrt(discr))/vij2;
						
    					if(tij < CollisionTime[i]){
     						CollisionTime[i] = tij;
     						CollisionPartner[i] = j;
    					}
    					if(tij < CollisionTime[j]){ // update j simutaneously
     						CollisionTime[j] = tij;
     						CollisionPartner[j] = i;
    					}
   					} // endif discriminant > 0
 	 			} //endif bij < 0
  			}// end loop periodic images
		
 		} //end if j is not i
 	}//end loop j
 }
 else{
 	for(int cell=0;cell<NumberOfNeighborCells;cell++){
  		jCell = NeighborCellList[iCell][cell];
  		j=HeadOfChain[jCell]; //printf("Here about to loop over neighbor particles\n");
  		while(j != -1){
  			if(j != i){
  				rij.x = position[i].x - position[j].x;
  				rij.y = position[i].y - position[j].y;
  				rij.z = position[i].z - position[j].z;
  				//pbc
  				MinimumImage(&(rij.x),Lx);
  				MinimumImage(&(rij.y),Ly);
  				MinimumImage(&(rij.z),Lz);

  				vij.x = velocity[i].x - velocity[j].x;
  				vij.y = velocity[i].y - velocity[j].y;
  				vij.z = velocity[i].z - velocity[j].z;
				
  				bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z;

  				if(bij < 0.){
   					rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   					vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   					sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   					discr = SQR(bij) - vij2*(rij2-sigmaij2);
					
   					if(discr > 0.){
    					tij = (-bij - sqrt(discr))/vij2;
						
    					if(tij < CollisionTime[i]){
     						CollisionTime[i] = tij;
     						CollisionPartner[i] = j;
    					}
    					if(tij < CollisionTime[j]){
     						CollisionTime[j] = tij;
     						CollisionPartner[j] = i;
    					}
   					} // endif discriminant > 0
  				} //endif bij < 0
  			} //end if j != i
			
    		j = CellTrack[j].Next;
  		}//end while j != -1
 	}//end loop cell
 }
 if(CollisionTime[i] < -tolerance){printf("Particle %d has negative CollisionTime[i]*10E10 = %lf with partner %ld\n", i, CollisionTime[i]*10E10, CollisionPartner[i]); exit(1);}
}

void CollisionUpdate(int i, int iCell) //
{
 int j;
 VECTOR rij,vij;
 double bij; // rij*vij
 double rij2,vij2,sigmaij2; // square
 double discr; // discriminate b^2 - 4ac
 double tij;
 int cell,jCell;

 for(cell=0;cell<NumberOfNeighborCells;cell++)
 {
  jCell = NeighborCellList[iCell][cell];
  j=HeadOfChain[jCell]; //printf("Here about to loop over neighbor particles\n");
  while(j != -1)
  {
  if(j != i)
  {
  rij.x = position[i].x - position[j].x;
  rij.y = position[i].y - position[j].y;
  rij.z = position[i].z - position[j].z;
  //pbc
  MinimumImage(&(rij.x),Lx);
  MinimumImage(&(rij.y),Ly);
  MinimumImage(&(rij.z),Lz);

  vij.x = velocity[i].x - velocity[j].x;
  vij.y = velocity[i].y - velocity[j].y;
  vij.z = velocity[i].z - velocity[j].z;

  bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z;

  if(bij < 0.)
  {
   rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   discr = SQR(bij) - vij2*(rij2-sigmaij2);

   if(discr > 0.)
   {
    tij = (-bij - sqrt(discr))/vij2;

    if(tij < CollisionTime[i])
    {
     CollisionTime[i] = tij;
     CollisionPartner[i] = j;
    }
    if(tij < CollisionTime[j])
    {
     CollisionTime[j] = tij;
     CollisionPartner[j] = i;
    }
   } // endif discriminant > 0
  } //endif bij < 0
  } //end if j != i

    j = CellTrack[j].Next;
  }//end while j != -1
 }//end loop cell
 if(CollisionTime[i] < -tolerance){printf("Particle %d has negative CollisionTime[i]*10E6 = %lf with partner %ld\n", i, CollisionTime[i]*10E6, CollisionPartner[i]); exit(1);}
 return;
}
