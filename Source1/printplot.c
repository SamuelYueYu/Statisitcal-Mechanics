#include <stdio.h>
#include <math.h>
#include "system.h"

void InstantBulkDensities(double t){
	char filename[20];
	FILE *fpinsdel_1; //particle numbers for type 1
 	FILE *fpinsdel_2; //particle numbers for type 2
 	FILE *fpfreqinsdel_1; //insdel frequency for type 1
 	FILE *fpfreqinsdel_2; //insdel frequency for type 2
    
    sprintf(filename,"insdel_1.txt");
    fpinsdel_1=fopen(filename,"a+");
	fprintf(fpinsdel_1,"%lf\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",\
	t, N_0_1, N_Z1_1, N_Z2_1, N_Z3_1,\
	N_0_1/(Lx*Ly*(Z1-sigmaA)), N_Z1_1/(Area1*(Z2-Z1+sigmaA)), N_Z2_1/(Lx*Ly*(Z3-Z2-sigmaA)), N_Z3_1/(Area1*(Lz-Z3+sigmaA)));
   	fclose(fpinsdel_1);
    
	sprintf(filename,"insdel_2.txt");
    fpinsdel_2=fopen(filename,"a+");
	fprintf(fpinsdel_2,"%lf\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n",\
	t, N_0_2, N_Z1_2, N_Z2_2, N_Z3_2,\
	N_0_2/(Lx*Ly*(Z1-sigmaB)), N_Z1_2/(Area2*(Z2-Z1+sigmaB)), N_Z2_2/(Lx*Ly*(Z3-Z2-sigmaB)), N_Z3_2/(Area2*(Lz-Z3+sigmaB)));
    fclose(fpinsdel_2);
    
	sprintf(filename,"freqinsdel_1.txt");
    fpfreqinsdel_1=fopen(filename,"a+");
    if(N_insert_trial_0_1 && N_insert_trial_Z2_1 && N_delete_trial_0_1 && N_delete_trial_Z2_1){
    	fprintf(fpfreqinsdel_1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t,\
    	N_insert_0_1, N_insert_Z2_1, N_delete_0_1, N_delete_Z2_1,\
		N_insert_0_1*1./N_insert_trial_0_1, N_insert_Z2_1*1./N_insert_trial_Z2_1,\
		N_delete_0_1*1./N_delete_trial_0_1, N_delete_Z2_1*1./N_delete_trial_Z2_1);
	}
    fclose(fpfreqinsdel_1);
    
    sprintf(filename,"freqinsdel_2.txt");
    fpfreqinsdel_2=fopen(filename,"a+");
    if(N_insert_trial_0_2 && N_insert_trial_Z2_2 && N_delete_trial_0_2 && N_delete_trial_Z2_2){
    	fprintf(fpfreqinsdel_2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t,\
    	N_insert_0_2, N_insert_Z2_2, N_delete_0_2, N_delete_Z2_2,\
		N_insert_0_2/N_insert_trial_0_2, N_insert_Z2_2/N_insert_trial_Z2_2,\
		N_delete_0_2/N_delete_trial_0_2, N_delete_Z2_2/N_delete_trial_Z2_2);
    }
    fclose(fpfreqinsdel_2);
    
    return;
}

void Rho(void){
	char filename[20];
	FILE *fprhoz_1;
	FILE *fprhoz_2;
	
	averaged_over++;
   	for(int j = 0; j < NA; j++){
    	if((position[j].z >= sigmaA/2 && position[j].z < Z1 - sigmaA/2) || (position[j].z >= Z2 + sigmaA/2 && position[j].z < Z3 - sigmaA/2)){
    		rho_1[(int) (2*position[j].z)] += 1/(Lx*Ly*MAX(sigmaA, sigmaB)/2);
		}
		else{
			rho_1[(int) (2*position[j].z)] += 1/(Area1*MAX(sigmaA, sigmaB)/2);
			for(int i=0; i<M; i++){
				if(SQR(position[j].x - pore_position[i].x)+SQR(position[j].y - pore_position[i].y) <= SQR((pore_position[i].d - sigma[j])/2)){
					if(position[j].z >= Z1 && position[j].z < Z2){rhoMem_left_1[i] += 1/(Area1/M);}
					else{rhoMem_right_1[i] += 1/(Area1/M);}
				}
	 		}
		}
	}
	for(int j = NA; j < NumberOfParticles; j++){
		if((position[j].z >= sigmaB/2 && position[j].z < Z1 - sigmaB/2) || (position[j].z >= Z2 + sigmaB/2 && position[j].z < Z3 - sigmaB/2)){
    		rho_2[(int) (2*position[j].z)] += 1/(Lx*Ly*MAX(sigmaA, sigmaB)/2);
		}
		else{
			rho_2[(int) (2*position[j].z)] += 1/(Area2*MAX(sigmaA, sigmaB)/2);
			
			for(int i=0; i<M; i++){
				if(SQR(position[j].x - pore_position[i].x)+SQR(position[j].y - pore_position[i].y) <= SQR((pore_position[i].d - sigma[j])/2)){
					if(position[j].z >= Z1 && position[j].z < Z2){rhoMem_left_2[i] += 1/(Area2/M);}
					else{rhoMem_right_2[i] += 1/(Area2/M);}
				}
	 		}
		}
	}
	sprintf(filename,"rhoz_1");
 	fprhoz_1 = fopen(filename,"a+");
 	fprintf(fprhoz_1, "%d\t%lf\t%lf\n", coln, rho_2[60]/averaged_over, rho_2[220]/averaged_over);
 	fclose(fprhoz_1);
	
	sprintf(filename,"rhoz_2");
 	fprhoz_2 = fopen(filename,"a+");
 	fprintf(fprhoz_2, "%d\t%lf\t%lf\n", coln, rho_2[60]/averaged_over, rho_2[220]/averaged_over);
 	fclose(fprhoz_2);
	
	n0_1 += N_0_1/(Lx*Ly*(Z1-sigmaA)); nZ1_1 += N_Z1_1/(Area1*(Z2-Z1+sigmaA));
	nZ2_1 += N_Z2_1/(Lx*Ly*(Z3-Z2-sigmaA)); nZ3_1 += N_Z3_1/(Area1*(Lz-Z3+sigmaA));
	
	n0_2 += N_0_2/(Lx*Ly*(Z1-sigmaB)); nZ1_2 += N_Z1_2/(Area2*(Z2-Z1+sigmaB));
	nZ2_2 += N_Z2_2/(Lx*Ly*(Z3-Z2-sigmaB)); nZ3_2 += N_Z3_2/(Area2*(Lz-Z3+sigmaB));
	
	return;
}

void PrintPlot(void){
	char filename[20];
 	FILE *fpplotscript_1; //file for plot, type 1
 	FILE *fpplotscript_2; //file for plot, type 2
	FILE *fpdistribution_1; //distribution of particle A along z-axis
 	FILE *fpdistribution_2; //distribution of particle B along z-axis
 
 	sprintf(filename,"numplotscript_1");
 	fpplotscript_1 = fopen(filename,"w");
	fprintf(fpplotscript_1, "set xlabel \"lg(time)\"\n");
 	fprintf(fpplotscript_1, "set ylabel \"Particle numbers\"\n");
	fprintf(fpplotscript_1, "set title \"MC particle number flux, type 1\"\n");
 	fprintf(fpplotscript_1, "set grid\n");
 	fprintf(fpplotscript_1, "plot \"insdel_1.txt\" using 1:2 title \"0 to Z_1 number\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:3 title \"Z_1 to Z_2 number\" with lines lc \"yellow\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:4 title \"Z_2 to Z_3 number\" with lines lc \"green\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:5 title \"Z_3 to L_z number\" with lines lc \"blue\" lw 2\n");
 	fprintf(fpplotscript_1, "pause -1 \"Hit return to continue\"\n");
 	fclose(fpplotscript_1);
 	
	sprintf(filename,"plotscript_1");
 	fpplotscript_1 = fopen(filename,"w");
	fprintf(fpplotscript_1,"set term png\n");
	fprintf(fpplotscript_1,"set output \"figure1.png\"\n");

	fprintf(fpplotscript_1, "set xlabel \"time\"\n");
 	fprintf(fpplotscript_1, "set ylabel \"Particle densities\"\n");
	fprintf(fpplotscript_1, "set title \"MC particle density flux, type 1\"\n");
 	fprintf(fpplotscript_1, "set grid\n");
 	fprintf(fpplotscript_1, "plot \"insdel_1.txt\" using 1:6 title \"0 to Z_1 density\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:7 title \"Z_1 to Z_2 density\" with lines lc \"yellow\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:8 title \"Z_2 to Z_3 density\" with lines lc \"green\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "\"insdel_1.txt\" using 1:9 title \"Z_3 to L_z density\" with lines lc \"blue\" lw 2\n");
 	fprintf(fpplotscript_1, "pause -1 \"Hit return to continue\"\n");
 	fclose(fpplotscript_1);
 	
 	sprintf(filename,"rhozplot_1");
 	fpplotscript_1 = fopen(filename,"w");
 	fprintf(fpplotscript_1, "set xlabel \"coln\"\n");
 	fprintf(fpplotscript_1, "set ylabel \"rho_1\"\n");
 	fprintf(fpplotscript_1, "set title \"rho_1 flux, type 2\"\n");
 	fprintf(fpplotscript_1, "set grid\n");
 	fprintf(fpplotscript_1, "plot \"rhoz_1\" using 1:2 title \"rho_1[30]\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_1, "plot \"rhoz_1\" using 1:3 title \"rho_1[110]\" with lines lc \"blue\" lw 2\n");
 	fclose(fpplotscript_1);
 	
// 	sprintf(filename,"insdelnplotscript_1");
// 	fpplotscript_1 = fopen(filename,"w");
// 	fprintf(fpplotscript_1, "set xlabel \"lg(time)\"\n");
// 	fprintf(fpplotscript_1, "set ylabel \"Insertion/deletion numbers\"\n");
// 	fprintf(fpplotscript_1, "set title \"MC insertion/deletion numbers, type 1\"\n");
// 	fprintf(fpplotscript_1, "set grid\n");
// 	fprintf(fpplotscript_1, "plot \"freqinsdel_1.txt\" using 1:2 title \"0 to Z_1 insert\" with lines lc \"red\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:3 title \"Z_2 to Z_3 insert\" with lines lc \"blue\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:4 title \"0 to Z_1 delete\" with lines lc \"yellow\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:5 title \"Z_2 to Z_3 delete\" with lines lc \"green\" lw 2\n");
// 	fprintf(fpplotscript_1, "pause -1 \"Hit return to continue\"\n");
// 	fclose(fpplotscript_1);
// 	
// 	sprintf(filename,"freqplotscript_1");
// 	fpplotscript_1 = fopen(filename,"w");
// 	fprintf(fpplotscript_1, "set xlabel \"lg(time)\"\n");
// 	fprintf(fpplotscript_1, "set ylabel \"Frequencies\"\n");
// 	fprintf(fpplotscript_1, "set title \"MC frequencies, type 1\"\n");
// 	fprintf(fpplotscript_1, "set grid\n");
// 	fprintf(fpplotscript_1, "plot \"freqinsdel_1\" using 1:6 title \"0 to Z_1 insert\" with lines lc \"red\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:7 title \"Z_2 to Z_3 insert\" with lines lc \"blue\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:8 title \"0 to Z_1 delete\" with lines lc \"yellow\" lw 2,\\\n");
// 	fprintf(fpplotscript_1, "\"freqinsdel_1\" using 1:9 title \"Z_2 to Z_3 delete\" with lines lc \"green\" lw 2\n");
// 	fprintf(fpplotscript_1, "pause -1 \"Hit return to continue\"\n");
// 	fclose(fpplotscript_1);
//	
	sprintf(filename,"numplotscript_2");
 	fpplotscript_2 = fopen(filename,"w");
 	fprintf(fpplotscript_2, "set xlabel \"time\"\n");
 	fprintf(fpplotscript_2, "set ylabel \"Particle numbers\"\n");
 	fprintf(fpplotscript_2, "set title \"MC particle number flux, type 2\"\n");
 	fprintf(fpplotscript_2, "set grid\n");
 	fprintf(fpplotscript_2, "plot \"insdel_2.txt\" using 1:2 title \"0 to Z_1\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:3 title \"Z_1 to Z_2\" with lines lc \"yellow\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:4 title \"Z_2 to Z_3\" with lines lc \"green\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:5 title \"Z_3 to L_z\" with lines lc \"blue\" lw 2\n");
 	fprintf(fpplotscript_2, "pause -1 \"Hit return to continue\"\n");
 	fclose(fpplotscript_2);
	 
 	sprintf(filename,"plotscript_2");
 	fpplotscript_2 = fopen(filename,"w");
	fprintf(fpplotscript_2,"set term png\n");
	fprintf(fpplotscript_2,"set output \"figure2.png\"\n");

 	fprintf(fpplotscript_2, "set xlabel \"time\"\n");
 	fprintf(fpplotscript_2, "set ylabel \"Particle densities\"\n");
 	fprintf(fpplotscript_2, "set title \"MC particle density flux, type 2\"\n");
 	fprintf(fpplotscript_2, "set grid\n");
 	fprintf(fpplotscript_2, "plot \"insdel_2.txt\" using 1:6 title \"0 to Z_1\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:7 title \"Z_1 to Z_2\" with lines lc \"yellow\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:8 title \"Z_2 to Z_3\" with lines lc \"green\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "\"insdel_2.txt\" using 1:9 title \"Z_3 to L_z\" with lines lc \"blue\" lw 2\n");
 	fprintf(fpplotscript_2, "pause -1 \"Hit return to continue\"\n");
 	fclose(fpplotscript_2);
 	
 	sprintf(filename,"rhozplot_2");
 	fpplotscript_2 = fopen(filename,"w");
 	fprintf(fpplotscript_2, "set xlabel \"coln\"\n");
 	fprintf(fpplotscript_2, "set ylabel \"rho\"\n");
 	fprintf(fpplotscript_2, "set title \"rho flux, type 2\"\n");
 	fprintf(fpplotscript_2, "set grid\n");
 	fprintf(fpplotscript_2, "plot \"rhoz_2\" using 1:2 title \"rho_2[30]\" with lines lc \"red\" lw 2,\\\n");
 	fprintf(fpplotscript_2, "plot \"rhoz_2\" using 1:3 title \"rho_2[110]\" with lines lc \"blue\" lw 2\n");
 	fclose(fpplotscript_2);
 
 	sprintf(filename,"dis_1.txt");
 	fpdistribution_1 = fopen(filename,"w");
 	fprintf(fpdistribution_1,"z\trho\n");
 	for(double i=0.; i < Lz; i+=.5){
 		fprintf(fpdistribution_1,"%lf\t%lf\n", i, rho_1[(int) (2*i)]);
 	}
 	fclose(fpdistribution_1);
 	
 	sprintf(filename,"disMem_1.txt");
 	fpdistribution_1 = fopen(filename,"w");
 	fprintf(fpdistribution_1,"pore number\trho\n");
 	for(int i=0; i < M; i++){
 		fprintf(fpdistribution_1,"%d\t%lf\t%lf\n", i, rhoMem_left_1[i], rhoMem_right_1[i]);
 	}
 	fclose(fpdistribution_1);
 
 	sprintf(filename,"dis_2.txt");
 	fpdistribution_2 = fopen(filename,"w");
 	fprintf(fpdistribution_2,"z\trho\n");
 	for(double i=0.; i < Lz; i+=.5){
 		fprintf(fpdistribution_2,"%lf\t%lf\n", i, rho_2[(int) (2*i)]);
 	}
	fclose(fpdistribution_2);
	
	sprintf(filename,"disMem_2.txt");
 	fpdistribution_2 = fopen(filename,"w");
 	fprintf(fpdistribution_2,"pore number\trho\n");
 	for(int i=0; i < M; i++){
 		fprintf(fpdistribution_2,"%d\t%lf\t%lf\n", i, rhoMem_left_2[i], rhoMem_right_2[i]);
 	}
	fclose(fpdistribution_2);
	
 	return;
}

void PrintInput(void){
	char filename[20];
	FILE *fpinput; //input for next run
	
	sprintf(filename,"input");
 	fpinput = fopen(filename, "w");
 	fprintf(fpinput,"%d\n",EnsembleType);
 	fprintf(fpinput,"0 %d %d 0\n",PackType,NumberOfLatticeSites);
 	fprintf(fpinput,"%d %d\n",PotentialType, MemType);
 	fprintf(fpinput,"%d %d %d\n",NumberOfParticles,NA,NB);
 	fprintf(fpinput,"%lf %lf\n",massA,massB);
 	fprintf(fpinput,"%lf %lf %lf\n",sigmaA,sigmaB,sigmaAB);
 	fprintf(fpinput,"%lf\n",T);
 	fprintf(fpinput,"%lf\n",V);
 	fprintf(fpinput,"%lf %lf %lf\n",Lx,Ly,Lz);
 	fprintf(fpinput,"%lf %lf %lf\n",X,Y,Z);
 	fprintf(fpinput,"%lf %lf %lf\n",Z1,Z2,Z3);
 	fprintf(fpinput,"%d %d %d\n",M,mx,my);
 	fprintf(fpinput,"%lf %lf %lf\n",D,Sigma,p);
 	fprintf(fpinput,"%lf %lf\n",mu1,mu2);
 	fprintf(fpinput,"%d 0\n",coln);
 	fprintf(fpinput,"%lf %lf\n",N_insert_0_1,N_insert_trial_0_1);
 	fprintf(fpinput,"%lf %lf\n",N_insert_Z2_1,N_insert_trial_Z2_1);
 	fprintf(fpinput,"%lf %lf\n",N_delete_0_1,N_delete_trial_0_1);
 	fprintf(fpinput,"%lf %lf\n",N_delete_Z2_1,N_delete_trial_Z2_1);
 	fprintf(fpinput,"%lf %lf\n",N_insert_0_2,N_insert_trial_0_2);
 	fprintf(fpinput,"%lf %lf\n",N_insert_Z2_2,N_insert_trial_Z2_2);
 	fprintf(fpinput,"%lf %lf\n",N_delete_0_2,N_delete_trial_0_2);
 	fprintf(fpinput,"%lf %lf\n",N_delete_Z2_2,N_delete_trial_Z2_2);
 	fprintf(fpinput,"%d %d %d\n",LCellx,LCelly,LCellz);
 	fprintf(fpinput,"%d\n",NumberOfInitialSteps);
 	fprintf(fpinput,"%d %d\n",SampleMultiplier,MovieMultiplier);
 	fprintf(fpinput,"%lf\n",randomseed);
 	fprintf(fpinput,"%d 1\n",CellSwitch);
 	fprintf(fpinput,"%lf\n",timewindow);
 	fclose(fpinput);
 	
 	return;
}

void BulkDensities(void){
	char filename[30];
	FILE *fpinsdel_1; //particle numbers for type 1
 	FILE *fpinsdel_2; //particle numbers for type 2
    
    	sprintf(filename,"insdel_cumulative_1.txt");
 	fpinsdel_1=fopen(filename,"a+");
 	fprintf(fpinsdel_1,"%lf\t%lf\t%lf\t%lf\n", n0_1, nZ1_1, nZ2_1, nZ3_1);
 	fclose(fpinsdel_1);
 	
 	sprintf(filename,"insdel_1.txt");
 	fpinsdel_1=fopen(filename,"a+");
 	fprintf(fpinsdel_1,"//M1_1 = %lf\tM2_1 = %lf\tJ_1 = %lf\n", M1_1, M2_1, J_1);
 	fclose(fpinsdel_1);
 
 	sprintf(filename,"insdel_cumulative_2.txt");
 	fpinsdel_2=fopen(filename,"a+");
 	fprintf(fpinsdel_2,"%lf\t%lf\t%lf\t%lf\n", n0_2, nZ1_2, nZ2_2, nZ3_2);
 	fclose(fpinsdel_2);
 	
 	sprintf(filename,"insdel_2.txt");
 	fpinsdel_2=fopen(filename,"a+");
 	fprintf(fpinsdel_2,"//M1_2 = %lf\tM2_2 = %lf\tJ_2 = %lf\n", M1_2, M2_2, J_2);
 	fclose(fpinsdel_2);
}

void JPalpha(void){
	char filename[20];
	FILE *fpJ; // J value for types 1 and 2
	FILE *fpP; // pressure values
	FILE *fpKn; // mean free paths
 	FILE *fpalpha; //separation factor
 	
 	sprintf(filename,"J.txt");
 	fpJ=fopen(filename,"w");
 	fprintf(fpJ,"J_1\tJ_2\tJAtot\tJBtot\n");
 	fprintf(fpJ,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", J_1, J_2, JAtot, JBtot);
 	fclose(fpJ);
 	
 	sprintf(filename,"pressure.txt");
 	fpP=fopen(filename,"w");
 	fprintf(fpP,"PLeftVirial\tPLeftGr\tPLeftCS\tPRightVirial\tPRightGr\tPRightCS\tMFPL\tMFPR\n");
 	fprintf(fpP,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",pLeftVirial,pLeftGr,pLeftCS,pRightVirial,pRightGr,pRightCS);
 	fclose(fpP);
 	
 	sprintf(filename,"MFPA.txt");
 	fpKn=fopen(filename,"w");
 	fprintf(fpKn,"MFPLA_p\tMFPLMA_p\tMFPRA_p\tMFPRMA_p\tMFPLA_rho\tMFPLMA_rho\tMFPRA_rho\tMFPRMA_rho\n");
 	fprintf(fpKn,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", MFPLA_p, MFPLMA_p, MFPRA_p, MFPRMA_p, MFPLA_rho, MFPLMA_rho, MFPRA_rho, MFPRMA_rho);
 	fclose(fpKn);
 	
 	sprintf(filename,"KnA.txt");
 	fpKn=fopen(filename,"w");
 	fprintf(fpKn,"KnA_p[0]\tKnA_p[1]\tKnA_p[2]\tKnA_p[3]\tKnA_rho[0]\tKnA_rho[1]\tKnA_rho[2]\tKnA_rho[3]\n");
 	fprintf(fpKn,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", KnA_p[0], KnA_p[1], KnA_p[2], KnA_p[3], KnA_rho[0], KnA_rho[1], KnA_rho[2], KnA_rho[3]);
 	fclose(fpKn);
 	
 	sprintf(filename,"MFPB.txt");
 	fpKn=fopen(filename,"w");
 	fprintf(fpKn,"MFPLB_p\tMFPLMB_p\tMFPRB_p\tMFPRMB_p\tMFPLB_rho\tMFPLMB_rho\tMFPRB_rho\tMFPRMB_rho\n");
 	fprintf(fpKn,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", MFPLB_p, MFPLMB_p, MFPRB_p, MFPRMB_p, MFPLB_rho, MFPLMB_rho, MFPRB_rho, MFPRMB_rho);
 	fclose(fpKn);
 	
 	sprintf(filename,"KnB.txt");
 	fpKn=fopen(filename,"w");
 	fprintf(fpKn,"KnB_p[0]\tKnB_p[1]\tKnB_p[2]\tKnB_p[3]\tKnB_rho[0]\tKnB_rho[1]\tKnB_rho[2]\tKnB_rho[3]\n");
 	fprintf(fpKn,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", KnB_p[0], KnB_p[1], KnB_p[2], KnB_p[3], KnB_rho[0], KnB_rho[1], KnB_rho[2], KnB_rho[3]);
 	fclose(fpKn);
 	
 	sprintf(filename,"alpha.txt");
 	fpalpha=fopen(filename,"w");
 	fprintf(fpalpha,"P_A\tP_B\talphaABP\talphaABrho\n");
 	fprintf(fpalpha,"%lf\t%lf\t%lf\t%lf\n", P_1, P_2, alphaABP, alphaABrho);
 	fclose(fpalpha);
}
