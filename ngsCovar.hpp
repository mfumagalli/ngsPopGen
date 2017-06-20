
#include "shared.hpp"

// print help
void info() {
  fprintf(stdout, "\nInput:\n-probfile: file with genotype posterior probabilities [required]\n-outfile: name of output file [required], currently it is a text file, tab separated with n*n cells\n-sfsfile: file with SFS posterior probabilities [required if you want to weight each site by its probability of being variable]\n-nind: nr of individuals [required]\n-nsites: nr of sites [required]\n-norm: if 0 no normalization, if 1 matrix is normalized by (p(1-p)) as in Patterson et al 2006, if 2 normalization is 2p(1-p) [0]\n-verbose: level of verbosity [0]\n-block_size: how many sites per block when reading the input file [0]\n-call: whether calling genotypes (1) or not (0) [0]\n-offset: starting position of subset analysis [1]\n-minmaf: filter out sites with estimated MAF less than minmaf or greater than 1-minmaf [0] (this filtering will be ignored when using the weighting approach\n-genoquality: text file with nsites lines; each line has a 0 and 1; if 0 the program will ignore this site [NULL]\n\n");
}

// compute estimated allele frequencies from genotype posterior probabilities
array<double> getAlleFreq (matrix<double> &m) {
  // m is dimensions: nsites * (nind*3)
  int nsites = m.x;
  int nind = m.y/3;
  double somma;
  double *tmp = new double[nsites];
  for (int i=0; i<nsites; i++) {
    somma = 0.0;
    for (int j=0; j<nind; j++) {
      somma = somma + m.data[i][(j*3)+1] + (2*m.data[i][(j*3)+2]);
    }
    tmp[i] = somma / (nind*2);    
  }  
  array<double> ret;
  ret.x = nsites;
  ret.data = tmp;
  return ret;
}

// get the covariance for a pair of individual; output is the effective number of sites (expected or passed the filter)
double calcCovarUp (matrix<double> &m, array<double> a, matrix<double> &covar, double minmaf, array<int> good, int start, int norm) {
  int nsites = m.x;
  int nind = m.y/3;
  if (nsites != a.x) {
    fprintf(stderr, "\n prob file and alle freq dimensions disagree. Terminate.");
    exit(-1);
  }
  double eff_nsites=0.0;
  for (int s=0; s<nsites; s++) {
    if ((a.data[s]>minmaf) & (a.data[s]<(1-minmaf))) {
      eff_nsites=eff_nsites+1.0;
    }
  }
  double somma = 0.0, subsomma = 0.0;
  double **data = new double*[nind];
  for (int i=0; i<nind; i++) {
    double *tmp = new double[nind];
    for (int j=0; j<(i+1); j++) {
      somma = 0.0;
      for (int s=0; s<nsites; s++) {
	subsomma = 0.0;
        if ((a.data[s]>minmaf) & (a.data[s]<(1-minmaf))) {
	  for (int C1=0; C1<3; C1++) {
            for (int C2=0; C2<3; C2++) {
             subsomma = subsomma + (C1-(2*a.data[s]))*(C2-(2*a.data[s]))*m.data[s][(i*3)+C1]*m.data[s][(j*3)+C2];
            }
	  }
          subsomma = subsomma * good.data[start+s];
          if (norm==1) subsomma = subsomma /  ((a.data[s]*(1-a.data[s])));
          if (norm==2) subsomma = subsomma /  (2*(a.data[s]*(1-a.data[s])));
        }
	if (isnan(subsomma)==0) somma = somma + subsomma;
      }
      tmp[j] = somma;
    }
    data[i] = tmp;
  }
  fprintf(stderr, "Message (not error/warning): for this window nsites is %d and effective is %f\n", nsites, eff_nsites);
  matrix<double> ret;
  ret.x = nind;
  ret.y = nind;
  ret.data = data;
  for (int i=0;i<covar.x;i++) {
    for (int j=0;j<(i+1);j++) {
      covar.data[i][j]=covar.data[i][j]+ret.data[i][j];
      covar.data[j][i]=covar.data[j][i]+ret.data[i][j];
    }
  }
  cleanup(ret);
  return eff_nsites;  
}


// get the covariance for a pair of individual by weighting by
void calcCovarUpProb (matrix<double> &m, array<double> a, matrix<double> &covar, array<double> pvar, array<int> good, int start, int norm) {
  int nsites = m.x;
  int nind = m.y/3;
  if (nsites != a.x) {
    fprintf(stderr, "\n prob file and alle freq dimensions disagree. Terminate.");
    exit(-1);
  }
  double somma = 0.0, subsomma = 0.0;
  double **data = new double*[nind];
  for (int i=0; i<nind; i++) {
    double *tmp = new double[nind];
    for (int j=0; j<(i+1); j++) {
      somma = 0.0;
      for (int s=0; s<nsites; s++) {
	    subsomma = 0.0;
	    for (int C1=0; C1<3; C1++) {
              for (int C2=0; C2<3; C2++) {
	        subsomma = subsomma + (C1-(2*a.data[s]))*(C2-(2*a.data[s]))*m.data[s][(i*3)+C1]*m.data[s][(j*3)+C2];
              }
	    }
            subsomma=subsomma * good.data[start+s] * pvar.data[s]; 
            if (norm==1) subsomma = subsomma /  ((a.data[s]*(1-a.data[s])));
            if (norm==2) subsomma = subsomma /  (2*(a.data[s]*(1-a.data[s])));
	    if (isnan(subsomma)==0) somma = somma + subsomma;
      }
      tmp[j] = somma;
    }
    data[i] = tmp;
  }
  matrix<double> ret;
  ret.x = nind;
  ret.y = nind;
  ret.data = data;
  for (int i=0;i<covar.x;i++) {
    for (int j=0;j<(i+1);j++) {
      covar.data[i][j]=covar.data[i][j]+(ret.data[i][j]);
      covar.data[j][i]=covar.data[j][i]+(ret.data[i][j]);
    }
  }
  cleanup(ret);
}


