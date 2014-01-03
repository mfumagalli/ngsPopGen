
// ngsFST computes an approximation of the expected FST=a/(a+b) (see Reynolds et al 1983) from expected a and expected b and correct for the ratio of a and b; it works in blocks of sites so it is memory efficient; it can receive in input a file from ANGSD -realSFS 1 or ANGSD -realSFS1 + optimSFS + sfstools.

/// IMPORTANT NOTE

// if you give as input a prior 2D-SFS file, then -postfile are from -realSFS 1 (they are in log)
// if you give no priorfiles I assume you have run sfstools and -postfiles are not in log
// if you give 2 priorfiles the program will take the product of of -sfsfile (which must be from -realSFS 1 and in log) and -priorfiles

// Ouput: a text file, each row is a site and columsn are tab separated: a, a+b, FACT, theta, pvar

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>
#include <vector> 
#include <math.h>

#include "ngsFST.hpp"
 
// to compile: g++ ngsFST.cpp -Wall -o bin/ngsFST -lm -lz -O0 -g

int main (int argc, char *argv[]) {
  
   if (argc==1) {
    info();
   return 0;   
  }

  /// DECLARE AND INITIALIZE VARIABLES
  
  char *sfsfile1=NULL; // posterior probabilities or sample allele frequency likelihoods
  char *sfsfile2=NULL;
  char *priorfile1=NULL; // marginal priors
  char *priorfile2=NULL;
  char *priorfile12=NULL; // joint prior, it is 2D-SFS

  FILE *outpost;
  char *outfile=NULL;
  char *foutpost=NULL;
  
  int argPos = 1, increment = 0, nind = 0, nind1 = 0, nind2 = 0, nsites = 0, verbose = 0, nsums = 1, block_size = 0, isfold=0, firstbase=1, islog=0;

  /// READ AND ASSIGN INPUT PARAMETERS
  
   while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-postfiles")==0) {
      sfsfile1 = argv[argPos+1];
      sfsfile2 = argv[argPos+2];
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-priorfile")==0) {
      priorfile12 = argv[argPos+1];
    }
    else if(strcmp(argv[argPos],"-priorfiles")==0) {
      priorfile1 = argv[argPos+1];
      priorfile2 = argv[argPos+2];
      increment = 1;
    }     
    else if(strcmp(argv[argPos],"-outfile")==0) outfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-nind")==0) {
      nind1 = atoi(argv[argPos+1]);
      nind2 = atoi(argv[argPos+2]);
      nind = nind1 + nind2;
      increment = 1;
    }
    else if(strcmp(argv[argPos],"-nsites")==0) nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-verbose")==0) verbose = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-block_size")==0) block_size = atoi(argv[argPos+1]);    
    else if(strcmp(argv[argPos],"-firstbase")==0) firstbase = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-isfold")==0) isfold = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-islog")==0) islog = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
    }
    argPos = argPos + 2 + increment;
  }
  
  /// CHECK INPUT
  if((sfsfile1 == NULL) & (sfsfile2 == NULL) ) {
    fprintf(stderr,"\nMust supply -postfiles.\n");
    info();
     return 0;
  }
  if(outfile == NULL) {
    fprintf(stderr,"\nMust supply -outfile.\n");
    info();
    return 0;
  }
  if ((priorfile12!=NULL) & (islog==0)) {
    fprintf(stderr, "\nIf you provide a 2D-SFS as a prior then -postfiles should be in log froamt (from -realSFS 1, not from sfstools).\n");
    info();
    return 0;
  }

  if ((priorfile12==NULL) & (priorfile1==NULL)) {
    fprintf(stderr, "\nYou are not providing any prior information. Have you already run sfstools?\n");
  }

  /// OUTPUT
  foutpost = append(outfile, "");
  fprintf(stderr,"\t->Dumping file: %s\n", foutpost);
  outpost = getFILE(foutpost, "w");
  
  // print input arguments
  if(verbose==1) fprintf(stderr,"\t->Using some of these args: -nind %d -nind1 %d -nind2 %d -nsites %d -postfiles %s %s -priorfiles %s %s -priorfile %s -outfile %s -verbose %d -firstbase %d -isfold %d -block_size %d -islog %d\n", nind, nind1, nind2, nsites, sfsfile1, sfsfile2, priorfile1, priorfile2, priorfile12, foutpost, verbose, firstbase, isfold, block_size, islog);

  // READ PRIORS (if provided)
  // marginal spectra
  array<double> prior1;
  array<double> prior2;
  if (priorfile1 != NULL) {
    if (verbose==1) fprintf(stderr, "\nReading priors...");
    prior1 = readArray(priorfile1, nind1, isfold);
    prior2 = readArray(priorfile2, nind2, isfold);
  }
  // 2D-SFS
  matrix<double> prior12;
  if ((priorfile12==NULL)==0) {
    if (verbose==1) fprintf(stderr, "\nReading 2D prior...");
    prior12 = readPrior12(priorfile12, nind1*2+1, nind2*2+1); // same dimensions even if it is folded
    if (verbose==2) {
      fprintf(stderr, "\nPrior 2d:\n");
      writematrix(prior12, stderr);
    }
    //// the difference with this prior is that I don't add the prior, but I add the prior directly at computeFST step
  }

  /// GET POSITIONS OF BLOCKS
  if (block_size>(nsites-firstbase+1)) block_size=(nsites-firstbase+1);
  if (block_size==0) block_size=nsites-firstbase+1;
  array<int> start; array<int> end;
  start=getStart(nsites, firstbase, block_size);
  end=getEnd(nsites, firstbase, block_size);
  int nwin= (nsites-firstbase+1)/block_size;
  if ( ( (nsites-firstbase+1) % block_size)!=0) nwin++;

  if (verbose==1) fprintf(stderr, "\n num win %d win0 is %d %d\n", nwin, start.data[0], end.data[0]);

  /// ITERATE OVER EACH BLOCK
  for (int n=0; n<nwin; n++) {

    if (verbose==1) fprintf(stderr, "Block %d out of %d from %d to %d\n", n, (nwin-1), start.data[n], end.data[n]);
  
    // READ POSTERIOR PROBABILITIES FILES
    matrix<double> post1;
    matrix<double> post2;
    post1 = readFileSub(sfsfile1, nind1, start.data[n], end.data[n], isfold);
    post2 = readFileSub(sfsfile2, nind2, start.data[n], end.data[n], isfold);

   // print first post probs
    if (verbose==1) {
      fprintf(stderr, "initial post probs 1: %f %f %f\n", post1.data[0][0], post1.data[0][1], post1.data[0][nind1]);
      fprintf(stderr, "initial post probs 2: %f %f %f\n", post2.data[0][0], post2.data[0][1], post2.data[0][nind2]);
    }

    // NORM FROM LOG if necessary
    if (islog) {
      // if from -realSFS 1 they are -log
      normSFS(post1, 1); // 2nd argument is islog
      normSFS(post2, 1);
    }

   if (verbose==1) {
      fprintf(stderr, "mid post probs 1: %f %f %f\n", post1.data[0][0], post1.data[0][1], post1.data[0][nind1]);
      fprintf(stderr, "mid post probs 2: %f %f %f\n", post2.data[0][0], post2.data[0][1], post2.data[0][nind2]);
    }

    // ADD PRIOR is marginal priors
    if (priorfile1!=NULL) {
      addPrior(post1, prior1);
      addPrior(post2, prior2);
      normSFS(post1, 0); // 2nd argument is islog
      normSFS(post2, 0);
    }

    // print first post probs
    if (verbose==1) {
      fprintf(stderr, "final post probs 1: %f %f %f\n", post1.data[0][0], post1.data[0][1], post1.data[0][nind1]);
      fprintf(stderr, "final post probs 2: %f %f %f\n", post2.data[0][0], post2.data[0][1], post2.data[0][nind2]);
    }

    // CALCULATE FST
    if (verbose==1) fprintf(stderr,"Computing FST.\n");
    if (priorfile12==NULL) {
      if (isfold) {
        computeVarReyFold(post1, post2, verbose, outpost, nsums);
      } else {
        computeVarRey(post1, post2, verbose, outpost, nsums);
      }
    } else {
      if (verbose==1) fprintf(stderr,"Using 2D-SFS as prior.\n");
      if (isfold) {
        computeVarRey12Fold(post1, post2, verbose, outpost, nsums, prior12);
      } else {
        computeVarRey12(post1, post2, verbose, outpost, nsums, prior12);
      }


    }

    cleanup(post1);
    cleanup(post2);
    
  } // end blocks iterations
 
  delete [] start.data;
  delete [] end.data;

  if ((priorfile12==NULL)==0) cleanup(prior12);

  if ((priorfile1==NULL)==0) delete [] prior1.data;
  if ((priorfile2==NULL)==0) delete [] prior2.data;

  fclose(outpost);
  free(foutpost);
  
  return 0;

} // end main





