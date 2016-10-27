#include "geo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * Method to print all the command line options for the program
 *   @return: void
 **/
void printAllOptionsCmdLineInput(void)
{
   printf("\n");
   printf("  Input Options:\n");
   printf("    --version            shows program's version number and exit\n");
   printf("    -h/--help            shows the help message and exit\n");
   // MANDATORY parameters
   printf("    --oik  $oikfile      file used to generate the SIK vector (MNADATORY)\n");
   printf("    --djk  $djkfile      file used to generate the DJK vector (MANDATORY)\n");
   printf("    --od   $odfile       file used to generate the XIJ vector (MANDATORY)\n");
   printf("    --dis  $disfile      file used to generate the DIS vector (MANDATORY)\n");
   printf("    --out  $outputfile   name of the output file              (MANDATORY)\n");
   printf("    --nrow $x            #rows N (int)                        (MANDATORY)\n");
   printf("    --ncol $x            #columns KK (int)                    (MANDATORY)\n");
   printf("    --beta $x            damping factor (double)              (MANDATORY)\n");
   printf("    --verbose $VERBOSITY verbosity level [low] (choose between low/medium/high) --> NOT IMPLEMENTED YET\n\n");
   // OPTIONAL parameters
   printf("    --maxiter    $x      max. #iterations (int)                        [99]\n");
   printf("    --thresh_xij $x      threshold value the XIJ vector (double)       [0.001]\n");
   printf("    --bc_dis     $x      boundary values for the DIS vector (double)   [50.0]\n");
   printf("    --start_bjk  $x      starting value for the BJK vector (double)    [1.00]\n");
   printf("    --convi      $x      convergence threshold I (double)              [50.0]\n");
   printf("    --convj      $x      convergence threshold J (double)              [50.0]\n");
   printf("    --convk      $x      convergence threshold K (double)              [100.0]\n");
   return;
}

/**
 * Method to print the selected input options
 *   @parameter input   Input object
 *   @return    void    
 **/
void printInputParameters(Input inp)
{
    printf("\n");
    printf("  Mandatory parameters:\n");
    printf("     #Rows:%d\n", inp->N);
    printf("     #Colums:%d\n",inp->KK);
    printf("     Beta:%10.4lf\n",inp->Beta);
    printf("     File OIK '%s'\n", inp->fileOik);
    printf("     File DJK '%s'\n", inp->fileDjk);
    printf("     File OD  '%s'\n", inp->fileOd);
    printf("     File DIS '%s'\n", inp->fileDis);
    printf("     Output File '%s'\n", inp->fileOut);
    printf("  Optional parameters:\n");
    printf("     Max. #iterations:%d\n", inp->maxIter);
    printf("     Threshold XIJ vector:%10.4lf\n",inp->threshXij);   
    printf("     Boundary Condition DIS:%10.4lf\n", inp->bcDis);
    printf("     Init. value BJK vector:%10.4lf\n", inp->initBjk);
    printf("     ConvI:%10.4lf\n",inp->convI);
    printf("     ConvJ:%10.4lf\n",inp->convJ);
    printf("     ConvK:%10.4lf\n",inp->convK);
    printf("\n");
    return;   
}

/**
 * Method which parses the command line and generates an Input Object
 *   @parameter  int argc     #Input arguments
 *   @parameter  char **argv  Array list of arguments 
 *   @return     Input
 */

Input parseCmdLineInput(int argc,char **argv)
{
    int i;
    int ret, itemp;
    double dtemp;
    char name_method[MAXLENGTH]="parseCmdLineInput";
    Input inp;

    // Mandatory Input tests
    bool foundOik=false;
    bool foundDjk=false;
    bool foundOd=false;
    bool foundDis=false;
    bool foundOut=false;
    bool foundN=false;
    bool foundKK=false;
    bool foundBeta=false;

    // Allocation of an Input Object
    inp=(Input)malloc(sizeof(struct INPUT));

    // Default Values
    inp->maxIter=99;
    inp->threshXij=0.001;
    inp->bcDis=50.0;
    inp->initBjk=1.00;
    inp->convI=50.0;
    inp->convJ=50.0;
    inp->convK=100.0;
    inp->verbosity=0;
    
    printf("  Invoked Command Line:\n     ");
    for(i=0; i<argc; i++)
        printf("%s ",argv[i]);
    printf("\n");
    i=1;
    while(i<argc)
    {
       if((strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"--help")==0))
       {
           printAllOptionsCmdLineInput();
           exit(-1); 
       }
       else if(strcmp(argv[i],"--version")==0)
       {
           printf("  Version Code:%s\n", GEOVERSION);
           exit(-1);
       }
       // MANDATORY VARIABLES
       else if(strcmp(argv[i],"--oik")==0)
       {  
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --oik option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              strcpy(inp->fileOik,argv[i]); 
           foundOik=true;
           i++;
       }
       else if(strcmp(argv[i],"--djk")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --djk option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              strcpy(inp->fileDjk,argv[i]);
           foundDjk=true;
           i++;
       }
       else if(strcmp(argv[i],"--od")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --od option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              strcpy(inp->fileOd,argv[i]);
           foundOd=true;
           i++;
       }
       else if(strcmp(argv[i],"--dis")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --dis option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              strcpy(inp->fileDis,argv[i]);
           foundDis=true;
           i++;
       }
       else if(strcmp(argv[i],"--out")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --out option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              strcpy(inp->fileOut,argv[i]);
           foundOut=true;
           i++;
       }

       else if(strcmp(argv[i],"--nrow")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --nrow option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%d", &itemp);

           // Valid number?
           if(itemp<1)
           {
              printf("  ERROR in '%s'  Invalid argument (<1) for the --nrow option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           foundN=true;
           inp->N=itemp;
           i++;
       }

       else if(strcmp(argv[i],"--ncol")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --ncol option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%d", &itemp);

           // Valid number?
           if(itemp<1)
           {
              printf("  ERROR in '%s'  Invalid argument (<1) for the --ncol option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           foundKK=true;
           inp->KK=itemp;
           i++;
       }

       else if(strcmp(argv[i],"--beta")==0)
       {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --beta option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);

           // Valid Damping Factor
           if(dtemp<0.0)
           {
              printf("  ERROR in '%s'  Invalid argument (<0.0) for the --beta option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           foundBeta=true;
           inp->Beta=dtemp;
           i++;
      }

      /* OPTIONAL ARGUMENTS */
      else if(strcmp(argv[i],"--maxiter")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --maxiter option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%d", &itemp);

           /* Valid number? */
           if(itemp<1)
           {
              printf("  ERROR in '%s'  Invalid argument (<1) for the --maxiter option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           inp->maxIter=itemp;
           i++;
      }
      else if(strcmp(argv[i],"--thresh_xij")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --thresh_xij option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->threshXij=dtemp;
           i++;
      }

      else if(strcmp(argv[i],"--bc_dis")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --bc_dis option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->bcDis=dtemp;
           i++;
      }

      else if(strcmp(argv[i],"--start_bjk")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --start_bjk option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->initBjk=dtemp;
           i++;
      }
 
      else if(strcmp(argv[i],"--convi")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --convi option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->convI=dtemp;
           i++;
      }

      else if(strcmp(argv[i],"--convj")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --convj option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->convJ=dtemp;
           i++;
      }
 
      else if(strcmp(argv[i],"--convk")==0)
      {
           i++;
           if(argc==i)
           {
              printf("  ERROR in '%s'  Missing argument for the --convk option\n",name_method);
              printAllOptionsCmdLineInput();
              exit(-1);
           }
           else
              ret=sscanf(argv[i], "%lf", &dtemp);
           inp->convK=dtemp;
           i++;
      }

      else
      {
         printf("  ERROR in '%s'  Unknown option:'%s'\n",name_method,argv[i]);
         printAllOptionsCmdLineInput();
         exit(-1);
      }
    } // End of While loop
   
    /* Check Existence of Mandatory Inputs */
    if(!foundOik)
    {
       printf("  ERROR in '%s'  Oik file has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundDjk)
    {
       printf("  ERROR in '%s'  Djk file has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundOd)
    {
       printf("  ERROR in '%s'  Od file has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundDis)
    {
       printf("  ERROR in '%s'  Dis file has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundOut)
    {
       printf("  ERROR in '%s'  The output file has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundN)
    {
       printf("  ERROR in '%s'  The #rows has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundKK)
    {
       printf("  ERROR in '%s'  The #cols has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }

    if(!foundBeta)
    {
       printf("  ERROR in '%s'  Beta has not been specified\n",name_method);
       printAllOptionsCmdLineInput();
       exit(-1);
    }
    return inp;
}

