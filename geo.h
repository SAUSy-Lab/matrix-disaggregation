#ifndef GEO_H
#define GEO_H

#define MAXLENGTH  250
#define GEOVERSION "0.1"

// Max. Required Heap Memory:
// -------------------------
// Integer: SIK, DJK             :: 2KK(N+1)
//          SIK_S, DJK_S         :: 2N (Temporary) -> Not counted
// Double : XMOD                 :: KK(N+1)^2
//          XIJ, DIS, Q, CIJ     :: 4(N+1)^2
//          AIK, SUMI, SUMJ, BJK :: 4(N+1)
// --------------------------------------------------------------
// Tot: 8KK(N+1) + 8KK(N+1)^2 + 32(N+1)^2 + 32(N+1)
//      (N+1) [ 8KK(N+1) + 8KK + 32(N+1) + 32]
//      (N+1) [ (8KK+32)(N+1) +(8KK+32)]
//      8 (N+2)(N+1)(KK+4) Bytes
//

// Boolean Type
typedef enum{false,true} bool;

// Parse Input
typedef struct INPUT
{
   int N;                      // Number of Rows
   int KK;                     // Number of Columns
   double Beta;                // Damping Factor

   char fileOik[MAXLENGTH];    // Name Input File for SIK Matrix
   char fileDjk[MAXLENGTH];    // Name Input File for DJK Matrix
   char fileOd[MAXLENGTH];     // Name Input File for XIJ Matrix
   char fileDis[MAXLENGTH];    // Name Input File for DIS Matrix
   char fileOut[MAXLENGTH];    // Name Output File

   int maxIter;                // Max. #iterations         [99]
   double threshXij;           // Low. limit XIJ           [0.001]
   double bcDis;               // Boundary Value for DIS   [50.0]
   double initBjk;             // Start value for BJK      [1.00]
   double convI;               // Conv. devI               [50.0]
   double convJ;               // Conv. devJ               [50.0]
   double convK;               // Conv. devK               [100.0]
   
   int verbosity;              // Verbosity of the Output
} *Input;

// Header
void printHeader(void);

// Input
void printAllOptionsCmdLineInput();
void printInputParameters(Input inp);
Input parseCmdLineInput(int argc,char **argv);

// Init functions
int * createSIK(char * filename, int nrow, int ncol);
int * createSIKred(int *SIK, int nrow, int ncol);
int createSIKscal(int *SIKred, int nrow);

double *createXIJ(char *filename, int nrow, double THRESHOLD, int *SIK_S, int * DJK_S);
double * createDIS(char *filename, int nrow, double BOUNDARY_VALUE);
double * createBJK(int nrow, int ncol, double START_VALUE);
double *createQ(int nrow, double beta, double *DIS);
double *createCIJ(double *Q, int *SIK, int *DJK, int nrow, int ncol);

// Loop functions
void calcAIK(double *Q, double *BJK, int *DJK, double *CIJ, double *XIJ,
             int nrow, int ncol, double *AIK);
void calcBJK(double *Q, double *AIK, int *SIK, double *CIJ, double *XIJ,
             int nrow, int ncol, double *BJK);
void calcCIJ(double *Q, double *AIK, int *SIK, double *BJK, int *DJK,
             int nrow, int ncol, double *CIJ);
void calcSUMIJ_XMOD(double *Q, double *CIJ, double *XIJ,
                    double *AIK, int *SIK,
                    double *BJK, int *DJK,
                    int nrow, int ncol,
                    double *SUMI, double *SUMJ, double *XMOD);
double calcDevIJ(int *A, double *B, int nrow, int ncol);
double calcDevK(double *XMOD, double *XIJ, int nrow, int ncol);

// Dump functions
void dumpOutput(char *filename, double *XIJ, double *XMOD,
                int nrow, int ncol, double threshold);

#endif
