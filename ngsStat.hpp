#include "shared.hpp"

// print help
void info() {
  fprintf(stdout, "\nInput:\n-npop: how many pops (1 or 2)\n-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population (with or without running sfstools)\n-outfile: name of the output file\n-nind: number of individuals for each population\n-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit\n-verbose: level of verbosity, if 0 suppress all messages\n-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk\n-firstbase: in case you want to analyze a subset of your sites this is the lower limit\n-iswin: if 1 then print the value computed for each non-overlapping window defined by block_size\n\n");
}


// compute summary stats for each block in case of 1 pop and print the results
void computeStats(matrix<double> &post1, int verbose, FILE *outpost, int iswin, int start) {

  int nind = (post1.y-1)/2; // sample size
  int nsites = post1.x;

  // init
  array<double> segsites;
  segsites.x=nsites;
  double *tmp1 = new double [nsites];
  for (int i=0; i<nsites; i++) { 
    tmp1[i]=0.0; 
  }
  segsites.data=tmp1;

  array<double> hetero;
  hetero.x=nsites;
  double *tmp2 = new double [nsites];
  for (int i=0; i<nsites; i++) { 
    tmp2[i]=0.0; 
  }
  hetero.data=tmp2;

  double temp=0.0;

  // to convert from 0 base to 1 base  
  int start0=start+1;

  double sum_segsites = 0.0, sum_hetero = 0.0;

  for (int s=0; s<nsites; s++) {

    start++;

    segsites.data[s] = 1.0 - post1.data[s][0];

    if (verbose==2) fprintf(stderr, "\n%f %f", post1.data[s][0], segsites.data[s]);
    segsites.data[s] = segsites.data[s] - post1.data[s][post1.y-1];
    if (verbose==2) fprintf(stderr, " %f %f", post1.data[s][post1.y-1],  segsites.data[s]);

    temp=0.0;
    for (int j=0; j<post1.y; j++) { temp = temp+ 2.0*(j/(nind*2.0))*((nind*2.0-j)/(nind*2.0))*post1.data[s][j];  }
    hetero.data[s]=temp;

    if (iswin==0) fprintf(outpost, "%d\t%d\t%f\t%f\n", start, start, segsites.data[s], hetero.data[s]);

    sum_segsites=sum_segsites+segsites.data[s];
    sum_hetero=sum_hetero+hetero.data[s];

  } // end cycle s

  // compute sum across all sites 
  if (iswin==1) fprintf(outpost, "%d\t%d\t%f\t%f\n", start0, start, sum_segsites, sum_hetero);

  delete [] segsites.data;
  delete [] hetero.data;

}


// compute summary stats for each block in case of 2 pops and print the results
void computeStats2Pops(matrix<double> &post1, int verbose, FILE *outpost, int iswin, int start, matrix<double> &post2) {
  // sample sizes
  int nind1 = (post1.y-1)/2;
  int nind2 = (post2.y-1)/2;
  int nsites=post1.x;

  // init
  array<double> segsites1;
  segsites1.x=nsites;
  array<double> hetero1;
  hetero1.x=nsites;
  array<double> segsites2;
  segsites2.x=nsites;
  array<double> hetero2;
  hetero2.x=nsites;
  array<double> fixed;
  fixed.x=nsites;
  array<double> dxy;
  dxy.x=nsites;

  double *tmp1= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp1[i]=0.0; }
  segsites1.data=tmp1;
  double *tmp2= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp2[i]=0.0; }
  hetero1.data=tmp2;
  double *tmp3= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp3[i]=0.0; }
  segsites2.data=tmp3;
  double *tmp4= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp4[i]=0.0; }
  hetero2.data=tmp4;
  double *tmp5= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp5[i]=0.0; }
  fixed.data=tmp5;
  double *tmp6= new double [nsites];
  for (int i=0; i<nsites; i++) { tmp6[i]=0.0; }
  dxy.data=tmp6;

  double temp=0.0;

  // to convert from 0 base to 1 base  
  int start0=start+1;

  double sum_segsites1=0.0, sum_hetero1=0.0;
  double sum_segsites2=0.0, sum_hetero2=0.0;
  double sum_fixed=0.0;
  double sum_dxy=0.0;

  for (int s=0; s<nsites; s++) {

    start++;

    segsites1.data[s]= 1.0 - post1.data[s][0];
    segsites1.data[s]=segsites1.data[s] - post1.data[s][post1.y-1];
    temp=0.0;
    for (int j=0; j<post1.y; j++) { temp=temp + 2.0*(j/(nind1*2.0))*((nind1*2.0-j)/(nind1*2.0))*post1.data[s][j]; }
    hetero1.data[s]=temp;
    sum_segsites1=sum_segsites1+segsites1.data[s];
    sum_hetero1=sum_hetero1+hetero1.data[s];

    segsites2.data[s]= 1.0 - post2.data[s][0];
    segsites2.data[s]=segsites2.data[s] - post2.data[s][post2.y-1];
    temp=0.0;
    for (int j=0; j<post2.y; j++) { temp=temp + 2.0*(j/(nind2*2.0))*((nind2*2.0-j)/(nind2*2.0))*post2.data[s][j]; }
    hetero2.data[s]=temp;
    sum_segsites2=sum_segsites2+segsites2.data[s];
    sum_hetero2=sum_hetero2+hetero2.data[s];

    fixed.data[s]=post1.data[s][0]*post2.data[s][post2.y-1] + post2.data[s][0]*post1.data[s][post1.y-1]; 
    sum_fixed = sum_fixed + fixed.data[s];

    dxy.data[s]=0.0;
    for (int i=0; i<post1.y; i++)
      for (int j=0; j<post2.y; j++)
	dxy.data[s] = dxy.data[s] + ( ( (i/(nind1*2.0))*(((nind2*2.0)-j)/(nind2*2.0)) + (((nind1*2.0)-i)/(nind1*2.0))*(j/(nind2*2.0)) ) * post1.data[s][i] * post2.data[s][j] ) ;
    sum_dxy = sum_dxy + dxy.data[s];


    if (iswin==0) fprintf(outpost, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", start, start, segsites1.data[s], hetero1.data[s], segsites2.data[s], hetero2.data[s], fixed.data[s], dxy.data[s]);

  } // end cycle s

  // compute sum across all sites 
  if (iswin==1) fprintf(outpost, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", start0, start, sum_segsites1, sum_hetero1, sum_segsites2, sum_hetero2, sum_fixed, sum_dxy);

  delete [] segsites1.data;
  delete [] hetero1.data;
  delete [] segsites2.data;
  delete [] hetero2.data;
  delete [] fixed.data;
  delete [] dxy.data;

}

