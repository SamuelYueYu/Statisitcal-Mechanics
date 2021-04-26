/***************************************************
* find collision time and partner of particle i
* and particle j
* can use cell list
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
double tol = 0.00000001;

double other_side(double z){
	if(!z){return Z3 - Lz;}
	else if(z == Z1){return Z2;}
	else if(z == Z2){return Z1;}
	else if(z == Z3){return Lz;}
	return Z3;
}

void cube_function_roots(double a, double b, double c){
	double p = b - SQR(a)/3.0, q = 2*a*a*a/27.0 - a*b/3.0 + c;
	
	if(p*p*p/27 + SQR(q/2) > 0.){
		u[0] = cbrt(sqrt(SQR(q/2) + p*p*p/27) - q/2) - cbrt(sqrt(SQR(q/2) + p*p*p/27) + q/2) - a/3.0;
		u[1] = u[0]; u[2] = u[0];
	}
	else if(p*p*p/27 + SQR(q/2) == 0.){
		u[0] = 2.0*cbrt(-q/2) - a/3.0;
		u[1] = -cbrt(-q/2) - a/3.0; u[2] = u[1];
	}
	else{
		double theta = acos(-q/2 / sqrt(-p*p*p/27)) / 3.0;
		u[0] = 2.0*sqrt(-p/3)*cos(theta) - a/3.0;
		u[1] = 2.0*sqrt(-p/3)*cos(theta + 2.0*M_PI/3.0) - a/3.0;
		u[2] = 2.0*sqrt(-p/3)*cos(theta - 2.0*M_PI/3.0) - a/3.0;
	}
	return;
}

double MemCollisionT(double d, double x, double vx, double xp, double y, double vy, double yp, double z, double vz, double zp, double sig){
	double t = TimeBig, tp;
	double r_v = (x-xp)*vx + (y-yp)*vy + (z-zp)*vz, SQR_Delr = SQR(x-xp) + SQR(y-yp) + SQR(z-zp), SQR_v = SQR(vx) + SQR(vy) + SQR(vz);
	double A = 4*r_v / SQR_v, B = (4*SQR(r_v) - SQR(d)*(SQR(vx) + SQR(vy))) / SQR(SQR_v) + (4*SQR_Delr + SQR(d) - SQR(sig)) / (2*SQR_v);
	double C = (2*SQR(d)*(z-zp)*vz + r_v*(4*SQR_Delr - SQR(d) - SQR(sig))) / SQR(SQR_v);
	double D = (SQR(SQR_Delr + (SQR(d) - SQR(sig)) / 4) - SQR(d)*(SQR(x-xp) + SQR(y-yp))) / SQR(SQR_v);
	
	double p = B - .375*SQR(A), q = .125*A*A*A - .5*A*B + C, w = -3.0*SQR(SQR(A))/256.0 + SQR(A)*B/16.0 - A*C/4.0 + D;
	cube_function_roots(-p/2.0, -w, (4.0*w*p - q*q) / 8.0);
	for(int i = 0; i < 3; i++){
		if(q >= 0.){
			if(2.0*u[i] - p >= 0. && SQR(u[i]) - w >= 0. && 4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p >= 0.){
				tp = (-sqrt(2.0*u[i] - p) + sqrt(4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
				
				tp = (-sqrt(2.0*u[i] - p) - sqrt(4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
			}
			if(2.0*u[i] - p >= 0. && SQR(u[i]) - w >= 0. && -4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p >= 0.){
				tp = (sqrt(2.0*u[i] - p) + sqrt(-4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
				
				tp = (sqrt(2.0*u[i] - p) - sqrt(-4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
			}
		}
		else{
			if(2.0*u[i] - p >= 0. && SQR(u[i]) - w >= 0. && 4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p >= 0.){
				tp = (sqrt(2.0*u[i] - p) + sqrt(4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
				
				tp = (sqrt(2.0*u[i] - p) - sqrt(4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
			}
			if(2.0*u[i] - p >= 0. && SQR(u[i]) - w >= 0. && -4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p >= 0.){
				tp = (-sqrt(2.0*u[i] - p) + sqrt(-4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
				
				tp = (-sqrt(2.0*u[i] - p) - sqrt(-4.0*sqrt(SQR(u[i]) - w) - 2.0*u[i] - p))/2 - r_v/SQR_v;
				t = tp > tolerance ? MIN(t, tp) : t;
			}
		}
	}
	return t;
}

void MemCollisionInfo(int i)
{
 double tij, t;
 double COS1, SIN2;
 double PosX, PosY;
 MemCollisionTime[i] = TimeBig;
 
 if((position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2) || (position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2)){
 	enter_pore[i].x = 0.; enter_pore[i].y = 0.; enter_pore[i].d = 0.;
 	enter_pore[i].z = velocity[i].z > 0. ? (Z1 - sigma[i]/2 >= position[i].z ? Z1 : Z3) : (position[i].z >= Z2 + sigma[i]/2 ? Z2 : 0.);
 	tij = velocity[i].z > 0. ? (enter_pore[i].z - sigma[i]/2 - position[i].z) / velocity[i].z : (enter_pore[i].z + sigma[i]/2 - position[i].z) / velocity[i].z;
 	
 	PosX = position[i].x + tij * velocity[i].x; PosY = position[i].y + tij * velocity[i].y;
	for(int j = 0; j < M; j++){
		if(SQR(PosX - pore_position[j].x)+SQR(PosY - pore_position[j].y) < SQR(pore_position[j].d/2)){
			enter_pore[i].x = pore_position[j].x;
			enter_pore[i].y = pore_position[j].y;
			enter_pore[i].d = pore_position[j].d;
		}
	}
 }
 else{
 	for(int j = 0; j < M; j++){
		if(SQR(position[i].x - pore_position[j].x)+SQR(position[i].y - pore_position[j].y) < SQR(pore_position[j].d/2)){
			enter_pore[i].x = pore_position[j].x;
			enter_pore[i].y = pore_position[j].y;
			enter_pore[i].d = pore_position[j].d;
		}
	}
 }
 if(velocity[i].z > 0. && ((position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2) || (position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2))){
	if(enter_pore[i].d == 0.){MemCollisionTime[i] = tij;}

	else if(sigma[i] > enter_pore[i].d){
		MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
	}

	else{
		MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
		
		if(MemCollisionTime[i] >= TimeBig || position[i].z + velocity[i].z*MemCollisionTime[i] > enter_pore[i].z){
			COS1 = ((enter_pore[i].x-position[i].x)*velocity[i].x + (enter_pore[i].y-position[i].y)*velocity[i].y) / (sqrt(SQR(enter_pore[i].x-position[i].x) + SQR(enter_pore[i].y-position[i].y)) * sqrt(SQR(velocity[i].x) + SQR(velocity[i].y)));
			
			SIN2 = sqrt(1-SQR(COS1)) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) / ((enter_pore[i].d - sigma[i])/2);
				
			t = sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y) + SQR((enter_pore[i].d - sigma[i])/2) + (enter_pore[i].d - sigma[i]) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * (COS1*sqrt(1-SQR(SIN2)) - SIN2*sqrt(1-SQR(COS1)))) / sqrt(SQR(velocity[i].x) + SQR(velocity[i].y));
			
			if(position[i].z + velocity[i].z*t <= other_side(enter_pore[i].z)){
				MemCollisionTime[i] = t;
			}
			else{
				MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, other_side(enter_pore[i].z), sigma[i]);
			}
		}
	}
 }	 

 else if(velocity[i].z < 0. && ((position[i].z >= sigma[i]/2 && position[i].z <= Z1 - sigma[i]/2) || (position[i].z >= Z2 + sigma[i]/2 && position[i].z <= Z3 - sigma[i]/2))){
	 if(enter_pore[i].d == 0.){MemCollisionTime[i] = tij;}

	 else if(sigma[i] > enter_pore[i].d){
	 	MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
	 }

	 else{
	 	MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
		 
		if(MemCollisionTime[i] >= TimeBig || position[i].z + velocity[i].z*MemCollisionTime[i] < enter_pore[i].z){
			COS1 = ((enter_pore[i].x-position[i].x)*velocity[i].x + (enter_pore[i].y-position[i].y)*velocity[i].y) / (sqrt(SQR(enter_pore[i].x-position[i].x) + SQR(enter_pore[i].y-position[i].y)) * sqrt(SQR(velocity[i].x) + SQR(velocity[i].y)));
			
			SIN2 = sqrt(1-SQR(COS1)) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) / ((enter_pore[i].d - sigma[i])/2);
				
			t = sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y) + SQR((enter_pore[i].d - sigma[i])/2) + (enter_pore[i].d - sigma[i]) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * (COS1*sqrt(1-SQR(SIN2)) - SIN2*sqrt(1-SQR(COS1)))) / sqrt(SQR(velocity[i].x) + SQR(velocity[i].y));
			
			if(position[i].z + velocity[i].z*t >= other_side(enter_pore[i].z)){
				MemCollisionTime[i] = t;
			}
			else{
				MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, other_side(enter_pore[i].z), sigma[i]);
			}
		}
	}
 }

 else if((!position[i].z && velocity[i].z < 0.) || (position[i].z >= Z1 && position[i].z <= Z2) || (position[i].z >= Z3 && position[i].z <= Lz)){
	if(position[i].z >= Z1 && position[i].z <= Z2){
 		if(velocity[i].z > 0.){enter_pore[i].z = Z2;}
 		else{enter_pore[i].z = Z1;}
	}
	else if(position[i].z >= Z3 && position[i].z <= Lz){
		if(velocity[i].z > 0.){enter_pore[i].z = Lz;}
 		else{enter_pore[i].z = Z3;}
	}
	else{enter_pore[i].z = Z3 - Lz;}
	
	COS1 = ((enter_pore[i].x - position[i].x)*velocity[i].x + (enter_pore[i].y - position[i].y)*velocity[i].y) / (sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * sqrt(SQR(velocity[i].x) + SQR(velocity[i].y)));
	
	SIN2 = sqrt(1-SQR(COS1)) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) / ((enter_pore[i].d - sigma[i])/2);

	tij = sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y) + SQR((enter_pore[i].d - sigma[i])/2) + (enter_pore[i].d - sigma[i]) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * (COS1*sqrt(1-SQR(SIN2)) - SIN2*sqrt(1-SQR(COS1)))) / sqrt(SQR(velocity[i].x) + SQR(velocity[i].y));
	
	if(tij * fabs(velocity[i].z) <= fabs(enter_pore[i].z - position[i].z)){MemCollisionTime[i] = tij;}

	else{
		MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
	}
 }

 else if((position[i].z <= Z1 && position[i].z > Z1-sigma[i]/2) || (position[i].z <= Z3 && position[i].z > Z3-sigma[i]/2) || (position[i].z >= 0. && position[i].z < sigma[i]/2) || (position[i].z >= Z2 && position[i].z < Z2+sigma[i]/2)){
	if(position[i].z <= Z1 && position[i].z > Z1-sigma[i]/2){enter_pore[i].z = Z1;}
	else if(position[i].z <= Z3 && position[i].z > Z3-sigma[i]/2){enter_pore[i].z = Z3;}
	else if(position[i].z >= 0. && position[i].z < sigma[i]/2){enter_pore[i].z = 0.;}
	else{enter_pore[i].z = Z2;}
	
	if(sigma[i] > enter_pore[i].d){
		MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
	}
	else{
		if((enter_pore[i].z - position[i].z) * velocity[i].z > 0.){
			MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
			
			if(MemCollisionTime[i] >= TimeBig || (velocity[i].z > 0. && position[i].z + velocity[i].z*MemCollisionTime[i] > enter_pore[i].z) || (velocity[i].z < 0. && position[i].z + velocity[i].z*MemCollisionTime[i] < enter_pore[i].z)){
				COS1 = ((enter_pore[i].x - position[i].x)*velocity[i].x + (enter_pore[i].y - position[i].y)*velocity[i].y) / (sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * sqrt(SQR(velocity[i].x) + SQR(velocity[i].y)));
				
				SIN2 = sqrt(1-SQR(COS1)) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) / ((enter_pore[i].d - sigma[i])/2);
				
				tij = sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y) + SQR((enter_pore[i].d - sigma[i])/2) + (enter_pore[i].d - sigma[i]) * sqrt(SQR(enter_pore[i].x - position[i].x) + SQR(enter_pore[i].y - position[i].y)) * (COS1*sqrt(1-SQR(SIN2)) - SIN2*sqrt(1-SQR(COS1)))) / sqrt(SQR(velocity[i].x) + SQR(velocity[i].y));
				
				if(position[i].z + velocity[i].z*tij >= MIN(enter_pore[i].z, other_side(enter_pore[i].z)) && position[i].z + velocity[i].z*tij <= MAX(enter_pore[i].z, other_side(enter_pore[i].z))){
					MemCollisionTime[i] = tij;
				}
				else{
					MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, other_side(enter_pore[i].z), sigma[i]);
				}
			}
		}
		else{
			MemCollisionTime[i] = MemCollisionT(enter_pore[i].d, position[i].x, velocity[i].x, enter_pore[i].x, position[i].y, velocity[i].y, enter_pore[i].y, position[i].z, velocity[i].z, enter_pore[i].z, sigma[i]);
		}
	}
 }
 if(MemCollisionTime[i] < 0.){printf("Particle %d has negative MemCollisionTime[i] = %lf at (%lf %lf %lf)\n", i, MemCollisionTime[i], position[i].x, position[i].y, position[i].z);}
 return;
}

void CollisionInfo(int i) //
{
 int j;
 VECTOR rij,vij;
 double bij; // rij*vij
 double rij2,vij2,sigmaij2; // square
 double discr; // discriminate b^2 - 4ac
 double tij;
 int iCell,cell,jCell;
 int l,m,n; // pbc box index
 double xtemp,ytemp,ztemp;

 CollisionTime[i] = TimeBig;
 CollisionPartner[i] = MAX_NumberOfParticles;

if(NumberOfCells <= NumberOfNeighborCells || CellSwitch == 0) // not use cell list
{
 for(j=0;j<NumberOfParticles;j++) 
 {
 if(j != i)
 {
  vij.x = velocity[i].x - velocity[j].x;
  vij.y = velocity[i].y - velocity[j].y;
  vij.z = velocity[i].z - velocity[j].z;
  for(l=-1;l<=1;l++)
  for(m=-1;m<=1;m++)
  for(n=-1;n<=1;n++)
  {
   xtemp = position[j].x + l*Lx;
   ytemp = position[j].y + m*Ly;
   ztemp = position[j].z + n*Lz;

   rij.x = position[i].x - xtemp;
   rij.y = position[i].y - ytemp;
   rij.z = position[i].z - ztemp;
  
   bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z;
   
   if(bij < 0.)
   {
   rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   discr = SQR(bij)- vij2*(rij2-sigmaij2);
  
   if(discr > 0.)
   {
    tij = (-bij - sqrt(discr))/vij2;

    if(tij < CollisionTime[i])
    {
     CollisionTime[i] = tij;
     CollisionPartner[i] = j;
    }
    if(tij < CollisionTime[j]) // update j simutaneously
    {
     CollisionTime[j] = tij;
     CollisionPartner[j] = i;
    }
   } // endif discriminant > 0
  } //endif bij < 0
  }// end loop periodic images

 } //end if j is not i
 }//end loop j
}//end if cell list not used
else // use cell list when CellSwitch == 1 && NumberOfCells>27
{
 iCell = CellTrack[i].WhichCell;
 
 for(cell=0;cell<NumberOfNeighborCells;cell++)
 {
  jCell = NeighborCellList[iCell][cell];
  j=HeadOfChain[jCell];
  
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
 } //end if j is not i

    j = CellTrack[j].Next;
  }//end while j != -1
 }//end loop cell
}//end cell list is used
 return;
}
