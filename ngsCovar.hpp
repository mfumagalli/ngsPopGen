/// TEMPLATES

#include <string> //for str operations
using namespace std;

// a general matrix style structure
template <typename T>
struct matrix{
  int x;
  int y;
  T** data;
};

template <typename T>
struct array{
  int x;
  T* data;
};

template <typename T>
T *collapse(std::vector<T> &v){
  T *tmp = new T[v.size()];
  for(int i=0;i<v.size();i++)
    tmp[i] = v[i];
  return tmp;
}

//function to cleanup our generic matrix structure
template <typename T>
void cleanup(matrix<T> &m){//using a reference to avoid copying the data
  for(int i=0;i<m.x;i++)
    delete [] m.data[i];
  delete [] m.data;
}

/// FUNCTIONS

array<int> getStart(int nsites, int firstbase, int block_size) {

  // note that firstbase and nsites are 1-based
  int len = nsites-firstbase+1;

  int nwin = len/block_size;
  if ( (len % block_size)!=0) nwin=nwin+1;

  array<int> start;
  start.x=nwin;
  
  int *tStart= new int [nwin];
  for (int i=0; i<nwin; i++) {
    tStart[i]=(i)*block_size;
  }

  // if you dont start from beginning
  if (firstbase>1) {
    for (int i=0; i<nwin; i++) {
      tStart[i]=tStart[i]+firstbase-1; // -1 because firstbase is 1-based
    }
  }
  start.data=tStart;
  return start;
}

array<int> getEnd(int nsites, int firstbase, int block_size) {

  // note that firstbase and nsites are 1-based
  int len = nsites-firstbase+1;

  int nwin = len/block_size;
  if ( (len % block_size)!=0) nwin=nwin+1;
  
  array<int> end;
  end.x=nwin;

  int *tEnd= new int [nwin];
  for (int i=0; i<nwin; i++) {
    tEnd[i]=(i)*block_size+block_size-1;
  }
  tEnd[nwin-1]=nsites-1; // nsites is 1 based

 
  // if you dont start from beginning
  if (firstbase>0) {
    for (int i=0; i<nwin; i++) {
      tEnd[i]=tEnd[i]+firstbase-1;
    }
  }

  end.data=tEnd;
  return end;
}


// get the filesize of afile
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

// find out if a file exists
int fexists(const char* str) {
  struct stat buffer ;
  return (stat(str, &buffer )==0 );
}

// a nice wrapper for getting files
FILE *getFILE(const char*fname,const char* mode) {
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  if(writeFile&&fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

// read a file into a matrix but only for a specific subsets of positions (0-based notation)
matrix<double> readFileSub(char *fname, int nInd, int start, int end, int isfold) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if (isfold==0) {
    if((filesize %(sizeof(double)*(2*nInd+1)) )) {
      fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
      exit(-1);
    }
  } else {
    if((filesize %(sizeof(double)*(nInd+1)) )) {
      fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
      exit(-1);
    }	  
  }
  int nsites = end-start+1;  
  double **data = new double*[nsites];
  if (isfold) {
	  fseek(fp, sizeof(double)*(nInd+1)*start, SEEK_SET);
  } else {
	  fseek(fp, sizeof(double)*(2*nInd+1)*start, SEEK_SET);
  }
  if (isfold) {
    for(int i=0; i<nsites; i++) {
      double *tmp = new double[nInd+1];
      fread(tmp,sizeof(double),nInd+1,fp);
      data[i]= tmp;
    }
  } else {
    for(int i=0; i<nsites; i++) {
      double *tmp = new double[2*nInd+1];
      fread(tmp,sizeof(double),2*nInd+1,fp);
      data[i]= tmp;
    }
  }
  fclose(fp);
  matrix<double> ret;
  ret.x = nsites;
  if (isfold) {
    ret.y = nInd+1;
  } else {
    ret.y = 2*nInd+1;
  }
  ret.data = data;
  return ret;
}

// read genotype posterior probabilities from angsd (-dogeno 64), but only for a specific subsets of positions (0-based notation)
matrix<double> readEstiSub(char *fname, int nInd, int start, int end) {
  FILE *fp = getFILE(fname,"rb");
  size_t filesize =fsize(fname);
  if((filesize %(sizeof(double)*(3*nInd)) )) {
    fprintf(stderr,"\n\t-> Possible error,binaryfiles might be broken\n");
    exit(-1);
  }
  int nsites = end-start+1;
  double **data = new double*[nsites];
  fseek(fp, sizeof(double)*(3*nInd)*start, SEEK_SET);
  for(int i=0; i<nsites; i++) {
    double *tmp = new double[3*nInd];
    fread(tmp,sizeof(double),3*nInd,fp);
    data[i]= tmp;
  }
  fclose(fp);
  matrix<double> ret;
  ret.x = nsites;
  ret.y = 3*nInd;
  ret.data = data;
  return ret;
}

// read genotype quality (boolean), analysis only on sites to be kept (1) and discard the rest (0)
array<int> readGenoQuality(const char *fname, int nsites) {
  // nsites is how many sites you want
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  int *tmp = new int[nsites];
  char *buf = new char[filesize];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atoi(strtok(buf,"\t \n"));
  for(int i=1;i<(nsites);i++)
    tmp[i] = atoi(strtok(NULL,"\t \n"));
  fclose(fp);
  array<int> allvalues;
  allvalues.x = nsites;
  allvalues.data = tmp;
  return allvalues;
}

// return max value position of a row from a matrix of geno likes or post probs for a specific individual
int maxposarr(matrix<double> &m, int row, int ind) {
  int i=0, res = ind;
  double val;
  val = m.data[row][ind];
  for (i = ind; i < (ind+3); i++) {
    if (m.data[row][i] > val) {
      res = i;
      val = m.data[row][i];
    }
  }
  return res;
}

// normalize SFS and exp it if log (if from -realSFS 1)
void normSFS(matrix<double> &sfs, int islog) {
  int nsites = sfs.x;
  int ncol = sfs.y;
  double somma = 0.0;
  for (int j=0; j<nsites; j++) {
    // get the sum of values (do exp if they are in log scale)
    somma = 0;
    for(int i=0; i<ncol; i++) {
      if (islog) {
        somma = somma + exp(sfs.data[j][i]);
        } else {
        somma = somma + sfs.data[j][i];
      }
    }
    // divide each value for the sum
    for(int i=0; i<ncol; i++) {
      if (islog) {
        sfs.data[j][i] = exp(sfs.data[j][i]) / somma;
      } else {
        sfs.data[j][i] = sfs.data[j][i] / somma;
      }
    }
   }
}


// get probability of being variable from posterior probabilities of sample allele frequencies
void getPvar(matrix<double> &sfs, array<double> pvar, int isfold) {
  if (isfold==0) {
    for (int i=0; i<sfs.x; i++)
      pvar.data[i]=1-sfs.data[i][0]-sfs.data[i][2*(sfs.y-1)];
  } else {
    for (int i=0; i<sfs.x; i++)
      pvar.data[i]=1-sfs.data[i][0];
  }
}

// write an array into a file
void writearray(array<double> &m,FILE *fp){
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%f\t",m.data[i]);
  fprintf(fp,"\n");
  
}

// write an array of Int into a file
void writearrayInt(array<int> &m,FILE *fp){
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%d\t",m.data[i]);
  fprintf(fp,"\n");
}

// write a matrix into a file
void writematrix(matrix<double> &m,FILE *fp){
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%f\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// write a matrix of Int into a file
void writematrixInt(matrix<int> &m,FILE *fp){
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%d\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}


// print help
void info() {
  fprintf(stdout, "\nInput:\n-probfile: file with genotype posterior probabilities [required]\n-outfile: name of output file [required], currently it is a text file, tab separated with n*n cells\n-sfsfile: file with SFS posterior probabilities [required if you want to weight each site by its probability of being variable]\n-nind: nr of individuals [required]\n-nsites: nr of sites [required]\n-norm: if 0 no normalization, if 1 matrix is normalized by (p(1-p)) as in Patterson et al 2006, if 2 normalization is 2p(1-p) [0]\n-verbose: level of verbosity [0]\n-block_size: how many sites per block when reading the input file [0]\n-call: whether calling genotypes (1) or not (0) [0]\n-offset: starting position of subset analysis [1]\n-minmaf: filter out sites with estimated MAF less than minmaf or greater than 1-minmaf [0] (this filtering will be ignored when using the weighting approach\n-genoquality: text file with nsites lines; each line has a 0 and 1; if 0 the program will ignore this site [NULL]\n-isfold: whether data in -sfsfile is folded [0]\n-islog: whether data in -sfsfile is in log values [0]\n\n");
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
  fprintf(stderr, "\n nsites is %d but effective is %f", nsites, eff_nsites);
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


