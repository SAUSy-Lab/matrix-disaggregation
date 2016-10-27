#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**  
 *  Function which creates the matrix SIK/DJK
 *    @parameter  char *filename: Name of the file
 *    @parameter  int nrow:       # Nrows (N+1)
 *    @parameter  int ncol:       # Ncols (KK)
 *    @return     matrix of size (NROW,NCOLS)        
 */
int * createSIK(char * filename, int nrow, int ncol) 
{
    int i,j,k;
    int jcol, irow, ival; 
    int *SIK;
    int numel=(nrow-1)*ncol;
    int numtotal=nrow*ncol;
    FILE *fp;

    // A. Allocate matrix and set all matrix el. to zero
    SIK=(int*)malloc(sizeof(int)*numtotal);
    if(SIK==NULL)
    {
       printf("createSIK_DJK:: The vector SIK/DJK can't be allocated\n");
       exit(-1);
    } 
 
    // B. Read ALL elements from file into SIK
    fp=fopen(filename,"r");
    if(fp==NULL)
    {
       printf("The file '%s' can not be opened\n",filename);
       exit(-1);
    }

    printf("  Start reading data from '%s'\n",filename);
    for(k=0; k<numel; k++)
    {
        fscanf(fp,"%d %d %d\n",&irow,&jcol,&ival); 
        i=irow-1;
        j=jcol-1; 
        SIK[i*ncol+j]=ival;  
    } 
    fclose(fp); 
    printf("  Finished reading from '%s'\n\n",filename);

    // C. For \all (i,k) where i \in [0,N[ AND k \in [0,KK[ 
    //    if SIK(i,k)==0 => 1
    for(k=0; k<numel; k++)
    {
       if(SIK[k]==0) SIK[k]=1;
    } 

   // D. All el in Last Row must be 1
   for(j=0; j<ncol; j++)
       SIK[numel+j]=1;

   return SIK;
} 


/**
 * Method which contracts a matrix (SIK/DJK) into a vector
 * over ALL BUT THE LAST ROW of the SIK/DJK matrix
 *   @parameter  int *SIK: SIK/DJK matrix
 *   @parameter  int nrow: #Rows (N+1)
 *   @parameter  int ncol: #Columns (KK)
 *   @return     Vector of length N (NOT N+1)!!!
 */
int * createSIKred(int *SIK, int nrow, int ncol)
{
   int irow, jcol;
   int N=nrow-1; 
   int shift;
   int *vec;

   // Allocate a vec of length (NROW-1)
   vec=(int*)malloc(sizeof(int)*N);
   if(vec==NULL)
   {
      printf("  createSIKred:: The vector can't be allocated\n");
      exit(-1);
   }

   // Contract the first N rows
   for(irow=0; irow<N ; irow++)
   {
       vec[irow]=0;
       shift = irow*ncol;
       for(jcol=0; jcol < ncol; jcol++)
           vec[irow] += SIK[shift + jcol]; 
   }
   return vec;
}

/**
 * Method which contracts a matrix (SIK/DJK) into a scalar
 * over ALL ELEMENTS BUT THE LAST ROW
 *   @parameter  int *SIKred: SIKred/DJKred vector
 *   @parameter  int nrow: #Rows (N+1)
 *   @return     int:      Scalar
 */
int createSIKscal(int *SIKred, int nrow)
{
    int k, res;
    int N = nrow-1;
    res = 0;
    for(k = 0; k < N; k++) 
        res += SIKred[k];
    return res;
} 


/**
 * Method to create the XIJ Matrix
 *   @parameter  char *filename   : Filename
 *   @parameter  int nrow         : (N+1) rows and (N+1) columns 
 *   @parameter  double THRESHOLD : If XIJ(i,j)==0 XIJ(i,j)=THRESHOLD 
 *                                     for all i,j \in (1,N)
 *   @parameter  int *SIK_S       : Vector of length N
 *   @param      int *DJK_S       : Vector of length N 
 *   @return     double *XIJ      : Matrix of (N+1,N+1)
 *
 *  NOTE THAT: XIJ(N+1,N+1) has NOT been initialized in the BASIC script 
 *             We EXPLICITLY initialized XIJ(N+1,N+1) = 0.0 
 */
 
double * createXIJ(char *filename, int nrow, double THRESHOLD, int *SIK_S, int * DJK_S)
{
    int N=nrow-1;
    int irow, jcol, ival;
    int i,j,k;
    double *XIJ, *TEMPI, *TEMPJ;
    FILE *fp;

    // Allocation of the matrices
    XIJ=  (double*)malloc(sizeof(double)*nrow*nrow);  // XIJ(1:N+1,1:N+1) = 0.0

    // Read ALL elements from file into the (1:N,1:N) SUB-matrix
    fp=fopen(filename,"r");
    if(fp==NULL)
    {
       printf("  createXIJ:: file '%s' can not be opened\n",filename);
       exit(-1);
    }

    printf("  Start reading from '%s'\n",filename);
    for(k=0; k< (N*N); k++)
    {
        fscanf(fp,"%d %d %d\n",&irow,&jcol,&ival);
        i=irow-1;
        j=jcol-1;

        if(ival==0)
           XIJ[i*nrow + j] = (double)THRESHOLD;
        else
           XIJ[i*nrow + j] = (double)ival;
    }
    fclose(fp);
    printf("  Finished reading from '%s'\n\n",filename);

    // Create TEMPI && TEMPJ from Submatrix XIJ(N,N)
    TEMPI=(double*)calloc(N,sizeof(double));  // TEMPI(1:N)       = 0.0
    TEMPJ=(double*)calloc(N,sizeof(double)); // TENPJ(1:N)       = 0.0
    for(irow=0; irow<N ; irow++)
    {
        for(jcol=0; jcol<N; jcol++) 
        {
            k=irow*nrow + jcol;
            TEMPI[irow] += XIJ[k];
            TEMPJ[jcol] += XIJ[k];
        }
    }  

    // Boundary conditions
    for(i=0;i<N; i++)
    {
        XIJ[i*nrow + N] = SIK_S[i] - TEMPI[i];
        XIJ[N*nrow + i] = DJK_S[i] - TEMPJ[i];
    } 
    XIJ[nrow*nrow-1]=0.00;  // NOT in ORIGINAL BASIC SCRIPT

    // Delete TEMPI & TEMPJ
    free(TEMPI);
    free(TEMPJ);
    return XIJ;
} 

/**
 * Method to create the DIS Matrix (N+1)x(N+1)
 *   @parameter   char *filename        : file to be read
 *   @parameter   int nrow              : #rows/#columns (N+1)
 *   @parameter   double BOUNDARY_VALUE : Value used to set boundary row and column
 *   @return      double *DIS           : Matrix (N+1,N+1)
 * NOTE: the DIS(N+1,N+1) element is NOT initialized in ORIGINAL BASIC code
 *       but here: DIS(N+1,N+1) = BOUNDARY_VALUE
 */ 
double * createDIS(char *filename, int nrow, double BOUNDARY_VALUE)
{
    int i,j,k;
    int irow, jcol;
    int N=nrow-1;

    double value;
    double *DIS;
    FILE *fp;
    DIS=(double*)malloc(sizeof(double)*nrow*nrow);


    // Read the (NxN) DIS submatrix from file 
    fp=fopen(filename,"r");
    if(fp==NULL)
    {
       printf("  createDIS:: file '%s' can not be opened\n",filename);
       exit(-1);
    }

    printf("  Start reading data from '%s'\n",filename);
    for(k=0; k< (N*N); k++)
    {
        fscanf(fp,"%d %d %lf\n",&irow,&jcol,&value);
        i=irow-1;
        j=jcol-1;
        DIS[i*nrow + j]=value;
    }
    fclose(fp);
    printf("  Finished reading from '%s'\n\n",filename);

    // DIS :: Boundary conditions
    for(i=0; i<N; i++)
    {
        DIS[i*nrow + N] = BOUNDARY_VALUE;
        DIS[N*nrow + i] = BOUNDARY_VALUE;
    }
    DIS[nrow*nrow-1]=BOUNDARY_VALUE;  // NOT in ORIGINAL BASIC SCRIPT
    return DIS;
}

/**
 * Method to create/initialize the BJK matrix
 *   @parameter int nrow: (N+1) rows
 *   @parameter int ncol: KK columns
 *   @oarameter double START_VALUE 
 *   @return    double* BJK(1:N+1,1:KK)=START_VALUE
 */ 
double * createBJK(int nrow, int ncol, double START_VALUE)
{
    int k;
    int nel=nrow*ncol;
    double *vec;

    vec=(double*)malloc(sizeof(double)*nel);
    if(vec==NULL)
    {
       printf("  createBJK:: The BJK vector can't be allocated\n");
       exit(-1);
    }
   
    for(k=0; k<nel;k++)
        vec[k]=START_VALUE;
    return vec; 
}

/**
 * Method to create the Q Matrix
 *   @parameter  int nrows  : #Rows (N+1)
 *   @parameter  double beta: Decay parameter 
 *   @parameter  double *DIS: DIS matrix (1:N+1,1:N+1)
 *   @return     double *Q  : Q matrix 
 *
 * NOTE: THE Q-MATRIX depends on DIS
 *       At DIS(N+1,N+1) the DIS is not defined (see above)! 
 *       The AIK, BJK and CIJ depend on the Q-MATRIX.
 */ 
double *createQ(int nrow, double BETA, double *DIS)
{
   int k;
   int nel=nrow*nrow;
   double *Q;

   Q=(double*)malloc(sizeof(double)*nrow*nrow);
   if(Q==NULL)
   {
      printf("createQ:: The matrix Q can't be allocated\n");
      exit(-1);
   }
    
   for(k=0; k<nel ; k++)
       Q[k]=exp(-BETA*DIS[k]);

   return Q;
}

/**
 * Method to create the CIJ Matrix
 *   @parameter double   *Q(1:N+1, 1:N+1)
 *   @parameter int    *SIK(1:N+1, 1:KK) 
 *   @parameter int    *DJK(1:N+1, 1:KK) 
 *   @parameter int    nrow(N+1)
 *   @parameter int    ncol(KK)
 *   @return    double *CIJ(1:N+1, 1:N+1)
 */
double *createCIJ(double *Q, int *SIK, int *DJK, int nrow, int ncol)
{
   int i,j,k,p;
   double *CIJ;

   CIJ=(double*)malloc(sizeof(double)*nrow*nrow);
   if(CIJ==NULL)
   {
      printf("createCIJ:: The matrix CIJ can't be allocated\n");
      exit(-1);
   }

   // Calculate CIJ
   for(i=0; i<nrow; i++)
   {
       for(j=0; j<nrow; j++)
       {
           p = i*nrow + j;
           CIJ[p]=0.0;
           for(k=0; k<ncol; k++) 
               CIJ[p] += SIK[i*ncol + k] * DJK[j*ncol + k]; 
           CIJ[p] *= Q[p];
  
           // Invert the element
           CIJ[p] = 1.0/CIJ[p];
       }
   } 
   return CIJ;
}

