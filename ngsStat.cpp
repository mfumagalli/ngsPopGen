
// this program receives as input the .sfs files from realSFS and sfstools and outputs several statistics for each window: expected number of segregating sites, expected average heterozygosity, expected number of fixed differendes (if 2 pops are provided)

#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <math.h>
#include "ngsStat.hpp"

// to compile: g++ -Wall -O0 -g ngsStat.cpp -o ngsStat

int main (int argc, char *argv[]) {

  if (argc==1) {
    info();
    return 0;
  }

  /// DECLARE AND INITIALIZE VARIABLES

  char *sfsfile1=NULL; // posterior probabilities or sample allele frequency likelihoods
  char *sfsfile2=NULL;

  FILE *outpost;
  char *outfile=NULL;
  char *foutpost=NULL;

  int argPos = 1, increment = 0, npop = 1, nind1 = 0, nind2 = 0, nsites = 0, verbose = 0, block_size = 0, isfold=0, firstbase=1, iswin = 0, islog = 0;

  // if iswin==1 then block_size is the window size and print the valeus for each window
  // if iswin==0 then block_size is just for efficiency and print values for each site

  /// READ AND ASSIGN INPUT PARAMETERS

  while (argPos<argc) {
    increment = 0;
    if(strcmp(argv[argPos],"-postfiles")==0) {
      sfsfile1 = argv[argPos+1];
      if (npop==2) {
        sfsfile2 = argv[argPos+2];
        increment = 1;
      }
    }
    else if(strcmp(argv[argPos],"-outfile")==0) outfile = argv[argPos+1];
    else if(strcmp(argv[argPos],"-npop")==0) npop = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nind")==0) {
      nind1 = atoi(argv[argPos+1]);
      if (npop==2)  {
        nind2 = atoi(argv[argPos+2]);
        increment = 1;
      }
    }
    else if(strcmp(argv[argPos],"-nsites")==0) nsites = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-verbose")==0) verbose = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-block_size")==0) block_size = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-firstbase")==0) firstbase = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-isfold")==0) isfold = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-islog")==0) islog = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-iswin")==0) iswin = atoi(argv[argPos+1]);
    else {
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
    }
    argPos = argPos + 2 + increment;
  }

  /// CHECK INPUT
  if(sfsfile1 == NULL) {
    fprintf(stderr,"\nMust supply -postfiles.\n");
    info();
     return 0;
  }
  if(outfile == NULL) {
    fprintf(stderr,"\nMust supply -outfile.\n");
    info();
    return 0;
  }

  /// OUTPUT
  foutpost = append(outfile, "");
  if (verbose==1) fprintf(stderr,"\t->Dumping file: %s\n", foutpost);
  outpost = getFILE(foutpost, "w");

  /// GET POSITIONS OF BLOCKS
  // if block_size longer than nsites
  if (block_size>(nsites-firstbase+1)) block_size=(nsites-firstbase+1);
  if (block_size==0) block_size=nsites-firstbase+1;

  array<int> start; array<int> end;
  start=getStart(nsites, firstbase, block_size);
  end=getEnd(nsites, firstbase, block_size);
  int nwin= (nsites-firstbase+1)/block_size;
  if ( ( (nsites-firstbase+1) % block_size)>0) nwin=nwin+1;
  if (verbose==1) fprintf(stderr, "\nLen %d and win %d and rest is %d", (nsites-firstbase+1), block_size, ( (nsites-firstbase+1) % block_size));

  if (verbose==1) fprintf(stderr, "\nNum of win %d and win[0] is %d %d\n", nwin, start.data[0], end.data[0]);

  /// ITERATE OVER EACH BLOCK
  for (int n=0; n<nwin; n++) {

    if (verbose==1) fprintf(stderr, "Block %d out of %d from %d to %d\n", n, (nwin-1), start.data[n], end.data[n]);

    // READ POSTERIOR PROBABILITIES FILES
    matrix<double> post1;
    matrix<double> post2;
    post1 = readFileSub(sfsfile1, nind1, start.data[n], end.data[n], isfold);
    if (npop==2) {
      post2 = readFileSub(sfsfile2, nind2, start.data[n], end.data[n], isfold);
    }

    // NORM FROM LOG if necessary
    if (islog) {
      // if from -realSFS 1 they are -log
      normSFS(post1, 1); // 2nd argument is islog
      if(npop==2) normSFS(post2, 1);
    }

    if (npop==1) computeStats(post1, verbose, outpost, iswin, isfold, start.data[n]);
    if (npop==2) computeStats2Pops(post1, verbose, outpost, iswin, isfold, start.data[n], post2);

    cleanup(post1);
    if (npop==2) cleanup(post2);

  } // end for n

  delete [] start.data;
  delete [] end.data;

  fclose(outpost);
  free(foutpost);

  return 0;


} // end main

 


