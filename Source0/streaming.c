/************ check for streaming ************/
#include <stdio.h>
#include "system.h"

void stream(double t){
	char filename[20];
	FILE *fpstream;
	
	sprintf(filename,"stream.txt");
 	if(t == 0.){
 		fpstream=fopen(filename,"w");
 		fprintf(fpstream,"time\tVcom.x\tVcom.y\tVcom.z\tTinstant\n");
	}
	else{
		fpstream=fopen(filename,"a+");
		fprintf(fpstream,"%lf\t%lf\t%lf\t%lf\t%lf\n", t, Vcom.x, Vcom.y, Vcom.z, Tinstant);
	}
 	fclose(fpstream);
	
	return;
}
