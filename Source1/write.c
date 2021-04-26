/***************************************
 * write movie snapshot in .xyz format
 *
 * number of particles
 * blank
 * name x y z
 * C 0.0000 0.0000 0.0000
 * S 0.0000 1.0000 0.2345
 *
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"

void Writemovie(FILE *FilePtr)
{
  int i;
  
  fprintf(FilePtr,"%d\n",NumberOfParticles);
  fprintf(FilePtr,"%lf %lf %lf %lf\n",Lx,Ly,Lz,T);

  for(i=0;i<NumberOfParticles;i++)
  {
     if(identity[i] == 1)	
      fprintf(FilePtr,"%s\t","N");
     else if(identity[i] == 2)	
      fprintf(FilePtr,"%s\t","S");
     else
      fprintf(FilePtr,"%s\t","O");
     

     fprintf(FilePtr,"%lf\t",position[i].x);
     fprintf(FilePtr,"%lf\t",position[i].y);
     fprintf(FilePtr,"%lf\n",position[i].z);
  }

}
