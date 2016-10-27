#include <stdio.h>
#include <stdlib.h>

/**
 * Write the output into a file
 *   @parameter char *filename: Name of the file
 *   @parameter double *XIJ(1:N+1,1:N+1)
 *   @parameter double *XMOD(1:N+1,1:N+1,1:KK)
 *   @parameter int nrow
 *   @parameter int ncol 
 *   @parameter double THRESHOLD
 */
void dumpOutput(char *filename, double *XIJ, double *XMOD,
                                int nrow, int ncol, double THRESHOLD)
{
   int i,j,k;
   FILE *fp;

   printf("  Writing data to '%s'\n",filename);
   fp=fopen(filename,"w");
   if(fp==NULL)
   {
      printf("ERROR: File '%s' can not be opened\n",filename);
      exit(-1);
   }
  
   // Output Indices starting at 1:w
   for(i=0; i<nrow ; i++)
   {
       for(j=0; j<nrow ; j++)
       {
           if(XIJ[i*nrow + j] > THRESHOLD)
           {
              fprintf(fp,"%d %d %12.4lf  ",i+1,j+1,XIJ[i*nrow + j]);
              for(k=0; k<ncol; k++)
                  fprintf(fp,"  %16.10lf", XMOD[i*nrow*ncol + j*ncol + k]);
               fprintf(fp,"\n");
           }
       }
   }
   fclose(fp);
   printf("  Finished writing data to '%s'\n",filename);
   return;
}
