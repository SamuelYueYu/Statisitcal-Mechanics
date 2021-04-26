#include "system.h"
int CellSwitch; // use or not use cell list
struct LinkedList CellTrack[MAX_NumberOfParticles];
int HeadOfChain[MAX_NumberOfCells];
double rcellx, rcelly, rcellz; // sizes of cell dimentions > sigma
int LCellx,LCelly,LCellz,NumberOfCells;
int NumberOfNeighborCells;
double rc; // rcutoff
int NeighborCellList[MAX_NumberOfCells][MAX_NumberOfNeighborCells];
int CellCoordinate[MAX_NumberOfCells][8];

double randomseed;

double collisionfrequency; 
double timewindow;

double tolerance,timetol;
double TimeBig; // collision time upper bound
double u[3];
double CollisionTime[MAX_NumberOfParticles];
double MemCollisionTime[MAX_NumberOfParticles];
long int CollisionPartner[MAX_NumberOfParticles];
double EscapeTime[MAX_NumberOfParticles]; // time needed to escape a cell
int EscapeID;

int JobIndex;
double rate; // cooling rate

int rInitialType,vInitialType;
int PackType; // SC, BCC or FCC
int PotentialType; //1: shifted-force L-J
int EnsembleType; // 0: NVE 1: isokinetic 2: NVT 3:NpT
int VType; // volume input type
int MemType; // pore type
int RhoSwitch;

int NumberOfParticles; //N
int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
int NumberOfLatticeSites;
int NA,NB;
int N_0_1,N_Z1_1,N_Z2_1,N_Z3_1;
int N_0_2,N_Z1_2,N_Z2_2,N_Z3_2;
double fA,fB; // fraction
int coln;
double SCT,averaged_over;
double N_insert_trial_0_1, N_insert_trial_0_2;
double N_insert_0_1, N_insert_0_2;
double N_insert_trial_Z2_1, N_insert_trial_Z2_2;
double N_insert_Z2_1, N_insert_Z2_2;
double N_delete_trial_0_1, N_delete_trial_0_2;
double N_delete_0_1, N_delete_0_2;
double N_delete_trial_Z2_1, N_delete_trial_Z2_2;
double N_delete_Z2_1, N_delete_Z2_2;
double M1_1, M2_1;
double M1_2, M2_2;
double J_1, J_2;
double JAtot, JBtot;
double P_1, P_2;
double rho_1[MAX_NumberOfCells], rho_2[MAX_NumberOfCells];
double rhoMem_left_1[MAX_M], rhoMem_right_1[MAX_M];
double rhoMem_left_2[MAX_M], rhoMem_right_2[MAX_M];
double JA[MAX_M], JB[MAX_M];
double n0_1, nZ1_1, nZ2_1, nZ3_1;
double n0_2, nZ1_2, nZ2_2, nZ3_2;
double xA,xB; // fraction of both types upstream
double yA,yB; // fraction of both types downstream
double xAM,xBM; // fraction of both types upstream membrane
double yAM,yBM; // fraction of both types downstream membrane
double alphaABP, alphaABrho; // separation factor through two formulae

int NumberOfInitialSteps;
int SampleMultiplier,MovieMultiplier;

double kB; // Boltzmann constant, set to 1
double lambda; // de Broglie wavelength, set to 1
double T,T0; // temperature
double Tinstant; // temperature
double Kinstant,Upotential;
double Virial,VirialLeft,VirialRight,VirialLeftMem,VirialRightMem; //Virial values
double Mcom; // total mass
double dT; // temperature increment
double beta; // 1/(kB*T)
double rho,rhoA,rhoB; // number density rho = N/V
double packingfraction; // phi = pi/6 rho in 3D
double V,Vo;  // Volume
double Lx,Ly,Lz,Lo; //simulation box length V = Lx*Ly*Lz
double X,Y,Z; //cut box short at x=X,y=Y,z=Z
double Z1,Z2,Z3; //membrane endpoints
int M,mx,my; //number of pores, M=mx*my
double D,Sigma,p; //diameter of pores
double Area1, Area2; //Porous area
double mu1,mu2; // chemical potentials
double latticeconst; // lattice constant
double pLeftVirial, pRightVirial; // pressure as measured with virial values
double pLeftMemVirial, pRightMemVirial;
double pLeftGr, pRightGr; // pressure as measured with radial distribution
double pLeftCS, pRightCS; // pressure as predicted by Carnahan and Starling formula
double pAL, pBL, pAR, pBR; // partial pressures
double pALM, pBLM, pARM, pBRM; // partial pressures
double MFPLA_p, MFPRA_p, MFPLB_p, MFPRB_p; // particle mean free path in left & right chambers
double MFPLA_rho, MFPRA_rho, MFPLB_rho, MFPRB_rho; // particle mean free path in left & right chambers
double MFPLMA_p, MFPRMA_p, MFPLMB_p, MFPRMB_p;
double MFPLMA_rho, MFPRMA_rho, MFPLMB_rho, MFPRMB_rho;
double KnA_p[4], KnB_p[4]; // Knudsen number
double KnA_rho[4], KnB_rho[4];

double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
double mass[MAX_NumberOfParticles]; // particle mass m_i
double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
double n_0_1[MAX_NumberOfParticles], n_Z2_1[MAX_NumberOfParticles];
double n_0_2[MAX_NumberOfParticles], n_Z2_2[MAX_NumberOfParticles];
double sigmaA,sigmaB,sigmaAB;
double massA,massB;

VECTOR Vcom;
VECTOR position[MAX_NumberOfParticles]; // t
VECTOR position_old[MAX_NumberOfParticles]; // t-dt
VECTOR position_new[MAX_NumberOfParticles]; // t+dt
VECTOR velocity[MAX_NumberOfParticles];
VECTOR velocity_old[MAX_NumberOfParticles];
VECTOR velocity_new[MAX_NumberOfParticles];
VECTOR pore_position[MAX_M];
VECTOR enter_pore[MAX_NumberOfParticles];
