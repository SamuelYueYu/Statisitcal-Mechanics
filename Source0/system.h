/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_NumberOfParticles 100000
#define MAX_drBins 10000 // for g(r)

#define MAX_NumberOfCells 1000000
#define MAX_NumberOfNeighborCells 30

#define MAX_M 10000

#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define RECTANGULAR(x,y,z) ((x)*(y)*(z))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define Normal(x, mu, s) (exp(-SQR((x-mu)/s) / 2))

#define f(d, x, xp, y, yp, z, zp, sig) (sqrt(SQR(d/2 - sqrt(SQR(x-xp) + SQR(y-yp))) + SQR(z - zp)) - sig / 2)
#define Compressibility(nu) ((1+nu+SQR(nu)-CUBIC(nu)) / CUBIC(1-nu))
/************************/
extern int CellSwitch; // use or not use cell list
struct LinkedList{
 int WhichCell;
 int Prev;
 int Next;
};
extern struct LinkedList CellTrack[MAX_NumberOfParticles];
extern int HeadOfChain[MAX_NumberOfCells];
extern double rcellx,rcelly,rcellz; // sizes of cell dimentions > sigma
extern int LCellx,LCelly,LCellz,NumberOfCells;
extern int NumberOfNeighborCells;
extern double rc;// rcutoff
extern int NeighborCellList[MAX_NumberOfCells][MAX_NumberOfNeighborCells];
extern int CellCoordinate[MAX_NumberOfCells][8];
/************************/
 
extern double randomseed;

extern double collisionfrequency; 
extern double timewindow;

extern double tolerance,timetol;
extern double TimeBig; // collision time upper bound
extern double u[3]; // cubic equation solutions
extern double CollisionTime[MAX_NumberOfParticles];
extern double MemCollisionTime[MAX_NumberOfParticles];
extern long int CollisionPartner[MAX_NumberOfParticles];
extern double EscapeTime[MAX_NumberOfParticles]; // time needed to escape a cell
extern double VanishTime[MAX_NumberOfParticles];
extern int EscapeID;
extern int VanishID;

extern int JobIndex;
extern double rate; // cooling rate

extern int rInitialType,vInitialType;
extern int PackType; // SC, BCC or FCC
extern int PotentialType; //1: shifted-force L-J
extern int EnsembleType; // 0: NVE 1: isokinetic 2: NVT 3:NpT
extern int VType; // volume input type
extern int MemType; // pore type
extern int RhoSwitch;

extern int NumberOfParticles; //N
extern int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
extern int NumberOfLatticeSites;
extern int NA,NB;
extern int N_0_1,N_Z1_1,N_Z2_1,N_Z3_1;
extern int N_0_2,N_Z1_2,N_Z2_2,N_Z3_2;
extern double fA,fB; // fraction
extern int coln;
extern double SCT,averaged_over;
extern double N_insert_trial_0_1, N_insert_trial_0_2;
extern double N_insert_0_1, N_insert_0_2;
extern double N_insert_trial_Z2_1, N_insert_trial_Z2_2;
extern double N_insert_Z2_1, N_insert_Z2_2;
extern double N_delete_trial_0_1, N_delete_trial_0_2;
extern double N_delete_0_1, N_delete_0_2;
extern double N_delete_trial_Z2_1, N_delete_trial_Z2_2;
extern double N_delete_Z2_1, N_delete_Z2_2;
extern double M1_1, M2_1;
extern double M1_2, M2_2;
extern double J_1, J_2;
extern double JAtot, JBtot;
extern double P_1, P_2;
extern double rho_1[MAX_NumberOfCells], rho_2[MAX_NumberOfCells];
extern double rhoMem_left_1[MAX_M], rhoMem_right_1[MAX_M];
extern double rhoMem_left_2[MAX_M], rhoMem_right_2[MAX_M];
extern double JA[MAX_M], JB[MAX_M];
extern double n0_1, nZ1_1, nZ2_1, nZ3_1;
extern double n0_2, nZ1_2, nZ2_2, nZ3_2;
extern double xA,xB; // fraction of both types upstream
extern double yA,yB; // fraction of both types downstream
extern double alphaABP, alphaABrho; // separation factor through two formulae

extern int NumberOfInitialSteps;
extern int SampleMultiplier,MovieMultiplier;

/**************************************************/
extern double kB; // Boltzmann constant, set to 1
extern double lambda; // de Broglie wavelength, set to 1
extern double T,T0; // temperature
extern double Tinstant; // temperature
extern double Kinstant,Upotential;
extern double Virial,VirialLeft,VirialRight; // Virial values
extern double Mcom; // total mass
extern double dT; // temperature increment
extern double beta; // 1/(kB*T)
extern double rho,rhoA,rhoB; // number density rho = N/V
extern double packingfraction; // phi = pi/6 rho in 3D
extern double V,Vo;  // Volume
extern double Lx,Ly,Lz,Lo; // simulation box length V = Lx*Ly*Lz
extern double X,Y,Z; // cut box short at x=X,y=Y,z=Z
extern double Z1,Z2,Z3; // membrane endpoints
extern int M,mx,my; //number of pores
extern double D,Sigma,p; //diameter of pores & related parameters
extern double Area1,Area2; //Porous area accounting for excluded volume
extern double muA,muB,sigmaV; // chemical potentials
extern double latticeconst; // lattice constant
extern double pLeftVirial, pRightVirial; // pressure as measured with virial values
extern double pLeftGr,pRightGr; // pressure as measured with radial distribution
extern double pLeftCS,pRightCS; // pressure as predicted by Carnahan and Starling formula
extern double pAL,pBL,pAR,pBR; // partial pressures
extern double MFPL,MFPR; // particle mean free path in left & right chambers
extern double Kn[4]; // Knudsen number
/**************************************************/

extern double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
extern double mass[MAX_NumberOfParticles]; // particle mass  m_i
extern double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
extern double n_0_1[MAX_NumberOfParticles], n_Z2_1[MAX_NumberOfParticles];
extern double n_0_2[MAX_NumberOfParticles], n_Z2_2[MAX_NumberOfParticles];
extern double sigmaA,sigmaB,sigmaAB;
extern double massA,massB;

typedef struct
{
	double x;
	double y;
	double z;
	double d;
} VECTOR;

extern VECTOR Vcom;
extern VECTOR position[MAX_NumberOfParticles]; // t
extern VECTOR position_old[MAX_NumberOfParticles]; // t-dt
extern VECTOR position_new[MAX_NumberOfParticles]; // t+dt
extern VECTOR velocity[MAX_NumberOfParticles];
extern VECTOR velocity_old[MAX_NumberOfParticles];
extern VECTOR velocity_new[MAX_NumberOfParticles];
extern VECTOR pore_position[MAX_M];
extern VECTOR enter_pore[MAX_NumberOfParticles];

void ReadInput(void);
void Initialization(int job);
void Initialize(void);
double BoxMuller(double mm, double ss);

double Distance(int i,int j); // return rij^2
double SigmaIJ(int iID,int jID);
void PBC(double *xx, double l);
void MinimumImage(double *xx, double l);

void Kinetic(void);
void CollisionInfo(int i);
void MemCollisionInfo(int i);
void Collision(int i,int j);
void MemCollision(int i);
void VanishInfo(int i);
void VanishUpdate(int i);
void Writemovie(FILE *FilePtr);
void MakeMembrane(void);
void Tkcommand(void);
void Insert(double z_small, double z_great, int Identity, double t);
void Delete(double z_small, double z_great, int Identity, double t);
void InsDel(double t);
double gTerm(int N, int Nlist[], double rho);
void InstantBulkDensities(double t);
void stream(double t);
void Rho(void);
void BulkDensities(void);
void PrintPlot(void);
void PrintInput(void);
void JPalpha(void);
void Jlist(double t);
void PrintJlist(void);

// cell list
int CellDetermine(int i);
void NeighborCell(void); // generate neighbor cell list
void CoordinateTrans(int i);
void MakeCell(void);
void AddToCell(int i,int iCell);
void RemoveFromCell(int i,int iCell);
void UpdateCell(int i);
void EscapeInfo(int i); 
void OverlapCheck(void);
void MemOverlapCheck(void);
void CollisionUpdate(int i, int iCell);
void InsdelCollisionUpdate(int i, int iCell);
