#include "geo.h"
#include <math.h>
#include <stdio.h>

/**
 * Method to calculate the AIK matrix
 *   @parameter  (double*) Q  : (N+1) x (N+1)
 *   @parameter  (double*) BJK: (N+1) x (KK)
 *   @parameter  (int*)    DJK: (N+1) x (KK)
 *   @parameter  (double*) CIJ: (N+1) x (N+1)
 *   @parameter  (double*) XIJ: (N+1) x (N+1)
 *   @parameter  (double*) AIK: (N+1) x (KK)  <== RESULT     
 */ 
void calcAIK(double *Q, double *BJK, int *DJK, double *CIJ, double *XIJ, 
             int nrow, int ncol, double *AIK)
{
    int i,j,k,p;
    int N=nrow-1;
    double prefac;

    // Initialization
    // AIK(1:N,1:KK)=0.
    for(k = 0 ; k<N*ncol ; k++)
        AIK[k]=0.0; 

    // AIK (N+1,1:KK)=1.0
    for(k = 0 ; k<ncol ; k++)
        AIK[N*ncol+k]=1.0;


    // Calculate AIK
    for(i=0; i<N; i++)
    {
        for(j=0; j<nrow; j++)
        {
            p = i*nrow + j;
            prefac = Q[p] * CIJ[p] * XIJ[p] ; 
            for(k=0; k<ncol; k++)
                AIK[i*ncol + k] += prefac * BJK[j*ncol + k] * DJK[j*ncol + k];
        }
    }  

    // Invert AIK(1:N,1:KK)
    for(k=0; k<N*ncol; k++)
        AIK[k]=1.0/AIK[k];

    return;
}


/**
 * Method to calculate the BJK matrix
 *   @parameter (double*) Q : (N+1)x(N+1)
 *   @parameter (double*) AIK: (N+1) x (KK)
 *   @parameter (int*)    SIK: (N+1) x (KK)
 *   @parameter (double*) CIJ: (N+1) x (N+1)
 *   @parameter (double*) XIJ: (N+1) x (N+1)
 *   @parameter (double*) BJK: (N+1) x (KK)  <== RESULT
 */
void calcBJK(double *Q, double *AIK, int *SIK, double *CIJ, double *XIJ,
             int nrow, int ncol, double *BJK)
{
    int i,j,k,p;
    int N=nrow-1;
    double prefac;
 
    // BJK(1:N,1:KK)=0.00
    for(k=0; k<N*ncol; k++)
        BJK[k]=0.00;

    // BJK(N+1,1:KK)=1.0
    for(k=0; k<ncol; k++)
        BJK[N*ncol+k]=1.0;

    // Calculate BJK
    for(j=0; j<N; j++)
    {
        for(i=0; i<nrow; i++)
        {
            p = i*nrow + j;
            prefac = Q[p] * CIJ[p] * XIJ[p] ;
            for(k=0; k<ncol; k++)
                BJK[j*ncol + k] +=  (prefac * AIK[i*ncol + k] * SIK[i*ncol + k]);
        }
    }

    // Invert BJK
    for(k=0; k<N*ncol; k++)
        BJK[k]=1.0/BJK[k];
    
    return;
} 


/**
 * Method to calculate the CIJ matrix
 *   @parameter  (double*) Q : (N+1)x(N+1) 
 *   @parameter  (double*) AIK: (N+1) x (KK)
 *   @parameter  (int*)    SIK: (N+1) x (KK) 
 *   @parameter  (double*) BJK: (N+1) x (KK)
 *   @parameter  (int*)    DJK: (N+1) x (KK)
 *   @parameter  int       nrow: (N+1) 
 *   @parameter  int       ncol: KK
 *   ---------------------------
 *   @parameter  (double*) CIJ: (N+1) x (N+1)  => RESULT
 */
void calcCIJ(double *Q, double *AIK, int *SIK, double *BJK, int *DJK,
             int nrow, int ncol, double *CIJ)
{
    int i,j,k,p;

    // Initialization 
    for(k=0; k<nrow*nrow; k++)
        CIJ[k]=0.;

    // Calculate CIJ
    for(i = 0; i<nrow; i++)
    {
        for(j = 0; j<nrow; j++)
        {  
            p = i*nrow + j;
            for(k = 0; k<ncol; k++)
                CIJ[p] += (AIK[i*ncol+k] * SIK[i*ncol+k] * BJK[j*ncol+k] * DJK[j*ncol+k]);
            CIJ[p] *= Q[p];
            CIJ[p] = (1.0 / CIJ[p]);
        }
    }
    return;
}

/**
 * Method to calculate the matrices SUMI, SUMJ, XMOD 
 * at the same time
 * The matrices SUMI, SUMJ and XMOD have been allocated previously
 *
 * @parameter: double *Q(1:N+1,1:N+1)
 * @parameter: double *CIJ(1:N+1,1:N+1)
 * @parameter: double *XIJ(1:N+1,1:N+1)
 * @parameter: double *AIK(1:N+1,1:KK)
 * @parameter: int    *SIK(1:N+1,1:KK)
 * @parameter: double *BJK(1:N+1,1:KK)
 * @parameter: int    *DJK(1:N+1,1:KK)
 * @parameter: int nrow
 * @parameter: int ncol
 * ------------------------------------------
 * @parameter: double *SUMI(1:N+1,1:KK)
 * @parameter: double *SUMJ(1:N+1,1:KK)
 * @parameter: double *XMOD(1:N+1,1:N+1,1:KK)
 */
void calcSUMIJ_XMOD(double *Q, double *CIJ, double *XIJ,
                    double *AIK, int *SIK,
                    double *BJK, int *DJK, 
                    int nrow, int ncol,
                    double *SUMI, double *SUMJ, double *XMOD) 
{
   int i,j,k,p, ik,jk;
   double prefac, temp;

   // Initialize SUMI(1:N+1,1:KK), SUMJ(1:N+1,1:KK) to 0.00
   for(k=0;k<nrow*ncol; k++)
   {
       SUMI[k]=0.00;
       SUMJ[k]=0.00;
   }  

   // Calculate SUMI, SUMJ and XMOD
   for(i=0; i<nrow; i++)
   {
       for(j=0; j<nrow; j++)
       {
           p = i*nrow + j;
           prefac = Q[p] * CIJ[p] * XIJ[p]; 
           for(k=0; k<ncol; k++)
           {
               ik = i*ncol + k;
               jk = j*ncol + k;
               temp = prefac * AIK[ik] * SIK[ik] * BJK[jk] * DJK[jk];
               XMOD[i*nrow*ncol + jk] =  temp;
               SUMI[ik] += temp;
               SUMJ[jk] += temp;
          }              
       }
   }
   return;
}


/**
 * METHOD to calculate the L1-Norm of the 
 * difference of two matrices (EXCLUDING the last row)
 *   @parameter A(1:N+1,1:KK)
 *   @parameter B(1:N+1,1:KK)
 *   @parameter  nrow
 *   @parameter  ncol
 * 
 *   @return L1(A-B) without the LAST ROW!!!
 */
double calcDevIJ(int *A, double *B, int nrow, int ncol)
{
   int k;
   int N = nrow-1;
   double dev;
   dev = 0.0;
   for(k = 0; k<N*ncol; k++)
       dev += fabs(A[k] - B[k]); 
   return dev;
}

/**
 * Method to calculate devk
 *   @parameter: double *XMOD(1:N+1,1:N+1,1:KK)
 *   @parameter: double *XIJ(1:N+1, 1:N+1)
 *   @parameter: int nrow
 *   @parameter: int ncol
 *   ------------------------------------------
 *   @return     double devk
 */

double calcDevK(double *XMOD, double *XIJ, int nrow, int ncol)
{
   int i,j,k;
   int N = (nrow-1);
   double devk, soverk;

   devk = 0.0;
   for(i = 0 ; i<N ; i++)
   {
       for(j = 0; j<N; j++)
       {
           soverk = 0.0; 
           for(k=0; k<ncol; k++)
               soverk += XMOD[i*nrow*ncol + j*ncol + k];    
           devk += fabs(soverk - XIJ[i*nrow + j]);
       }
   }
   return devk;
}

