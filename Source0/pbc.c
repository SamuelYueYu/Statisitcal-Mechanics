/***************************************
 * periodic boundary condition
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include <math.h>

/************************/
void PBC(double *xx, double l)
{
  if(*xx < 0.0)
   *xx += l;
  else if(*xx >= l)
   *xx -= l; 

 return;
}
/************************/

void MinimumImage(double *xx, double l)
{
 *xx = *xx - l*round(*xx/l);

 return;
}

/***************
void PBC(void)
{
    if(position[i].x < 0.0)
     position[i].x += L;
    else if(position[i].x >= L)
     position[i].x -= L; 
    
    if(position[i].y < 0.0)
     position[i].y += L;
    else if(position[i].y >= L)
     position[i].y -= L; 
    
    if(position[i].z < 0.0)
     position[i].z += L;
    else if(position[i].z >= L)
     position[i].z -= L; 

return;
}
***************/
