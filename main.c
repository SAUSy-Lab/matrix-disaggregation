#include "geo.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv)
{ 
    Input inp;
    bool converged;
    int itrap, cloop,NPLUS;
    double devi, devj, devk;

    // Vectors 
    int *SIK, *SIK_S;
    int *DJK, *DJK_S;
    double *XIJ, *DIS;
    double *Q;
    double *AIK,*BJK,*CIJ;
    double *SUMI, *SUMJ, *XMOD;

    // INPUT
    printHeader();
    inp=parseCmdLineInput(argc,argv);
    printInputParameters(inp);
    NPLUS=(inp->N+1); // adds the imaginary zone

    // A.INITIALIZATION::
    SIK=createSIK(inp->fileOik,NPLUS,inp->KK);                    // (int*) SIK    : (N+1) x KK
    SIK_S=createSIKred(SIK,NPLUS,inp->KK);                        // (int*) SIK_S  : (N) x 1
    DJK=createSIK(inp->fileDjk,NPLUS,inp->KK);                    // (int*) DJK    : (N+1) x KK
    DJK_S=createSIKred(DJK,NPLUS,inp->KK);                        // (int*) DJK_S  : (N) x 1
    XIJ=createXIJ(inp->fileOd,NPLUS,inp->threshXij,SIK_S,DJK_S);  // (double*) XIJ : (N+1) x (N+1)  
    free(SIK_S); free(DJK_S);
    DIS=createDIS(inp->fileDis,NPLUS,inp->bcDis);                 // (double*) DIS : (N+1) x (N+1) 
    BJK=createBJK(NPLUS,inp->KK,inp->initBjk);                    // (double*) BJK : (N+1) x (KK)
    Q=createQ(NPLUS,inp->Beta,DIS);                               // (double*) Q   : (N+1) x (N+1)
    CIJ=createCIJ(Q,SIK,DJK,NPLUS,inp->KK);                       // (double*) CIJ : (N+1) x (N+1)

    // B.LOOP::
    AIK=(double*)malloc(sizeof(double)*NPLUS*inp->KK);            // (double*) AIK : (N+1) x (KK) 
    SUMI=(double*)malloc(sizeof(double)*NPLUS*inp->KK);           // (double*) SUMI: (N+1) x (KK)
    SUMJ=(double*)malloc(sizeof(double)*NPLUS*inp->KK);           // (double*) SUMJ: (N+1) x (KK)
    XMOD=(double*)malloc(sizeof(double*)*NPLUS*NPLUS*inp->KK);    // (double*) XMOD: (N+1) x (N+1) x (KK)

    converged=false;
    cloop=0;
    itrap=0;
    while(!converged)
    {
       itrap++; 
       // Calculate AIK, BJK & CIJ
       calcAIK(Q,BJK,DJK,CIJ,XIJ,NPLUS,inp->KK,AIK);
       calcBJK(Q,AIK,SIK,CIJ,XIJ,NPLUS,inp->KK,BJK);
       calcCIJ(Q,AIK,SIK,BJK,DJK,NPLUS,inp->KK,CIJ);
       cloop++;

       // Calculate SUMI, SUMJ, XMOD
       calcSUMIJ_XMOD(Q,CIJ,XIJ,
                      AIK,SIK,BJK,DJK,
                      NPLUS,inp->KK,
                      SUMI,SUMJ,XMOD);
       devi=calcDevIJ(SIK,SUMI,NPLUS,inp->KK);
       devj=calcDevIJ(DJK,SUMJ,NPLUS,inp->KK);
       devk=calcDevK(XMOD,XIJ,NPLUS,inp->KK);

       printf("     iter=%4d  devi=%10.2lf  devj=%10.2lf  devk=%10.2lf\n", itrap, devi, devj, devk);
 
       if(itrap>inp->maxIter)
       {
          printf("Runaway loop for a and b -> exiting \n");
          exit(-1);
       }
       // check for convergence
       if((devi <= inp->convI) && (devj <= inp->convJ) && (devk <= inp->convK))
       {
          printf("     iter=%4d --> CONVERGED\n\n",itrap);
          converged=true;
       }
    }

    // C. WRITE TO FILE 
    dumpOutput(inp->fileOut, XIJ,XMOD,NPLUS,inp->KK,inp->threshXij);

    // D. RELEASE MEMORY FROM HEAP
    free(SIK); 
    free(DJK); 
    free(XIJ); free(DIS);
    free(Q);
    free(AIK); free(BJK); free(CIJ);
    free(SUMI); free(SUMJ); free(XMOD);  
    free(inp);
    return 0;
}
