#include <stdio.h>
#include <math.h>
#include "system.h"

void Tkcommand(){
	FILE *fptk, *fppore;
	double d;
	
	fppore = fopen("pore_size.dat", "r");
	fptk = fopen("commands.dat","w");
	
	fprintf(fptk, "topo readvarxyz movie_%d.xyz\n", JobIndex);
	fprintf(fptk, "draw delete all\n");
	fprintf(fptk, "\n");
	
	fprintf(fptk, "[atomselect top \"all\"] set radius %lf\n[atomselect top \"name S\"] set radius %lf\n", 2*sigmaA, 2*sigmaB);
	fprintf(fptk, "\n");
	
	fprintf(fptk, "pbc set {%lf %lf %lf 90 90 90} -all\npbc box -on -color green\n", Lx, Ly, Lz);
	fprintf(fptk, "\n");
	
	fprintf(fptk, "draw color red\ndraw material Transparent\n");
	fprintf(fptk, "\n");
	
	for(int j=0; j<my; j++){
		for(int i=0; i<mx; i++){
			fscanf(fppore, "%lf", &d);
			fprintf(fptk,"draw cylinder {%lf %lf %lf} {%lf %lf %lf} radius %lf resolution 80\n", Lx*(i+0.5)/mx, Ly*(j+0.5)/my, Z1, Lx*(i+0.5)/mx, Ly*(j+0.5)/my, Z2, d/2);
			fprintf(fptk,"draw cylinder {%lf %lf %lf} {%lf %lf %lf} radius %lf resolution 80\n", Lx*(i+0.5)/mx, Ly*(j+0.5)/my, Z3, Lx*(i+0.5)/mx, Ly*(j+0.5)/my, Lz, d/2);
			fprintf(fptk,"draw cylinder {%lf %lf %lf} {%lf %lf %lf} radius %lf resolution 80\n", Lx*(i+0.5)/mx, Ly*(j+0.5)/my, 0., Lx*(i+0.5)/mx, Ly*(j+0.5)/my, .01, d/2);
	 	}
 	}
	
	fclose(fppore); fclose(fptk);
	return;
}
