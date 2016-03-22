/// TEMPLATES

#include <string> //for str operations
#include <unistd.h>

using namespace std;

// a general matrix style structure
template <typename T>
struct matrix{
  int x;
  int y;
  T** data;
};

// a general array style structure
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

// get the filesize of a file
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

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// a nice wrapper for getting files
FILE *getFILE(const char* fname,const char* mode) {
  FILE *fp;
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  if(writeFile && fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  
  if(strcmp(fname, "-") == 0){
    if(writeFile)
      fp = stdout;
    else
      fp = stdin;
  } else {
    if(NULL==(fp=fopen(fname,mode))){
      fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
      exit(0);
    }
  }

  if(isatty(fileno(fp))) {
    fprintf(stderr, "Your stdin/stdout is not a pipe, this might not be what you want!\n");
    exit(0);
  }
  
  return fp;
}

// read a file of prior into an array
array<double> readArray(const char *fname, int nInd) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  int len = 2*nInd+1;
  char *buf = new char[filesize];
  double *tmp = new double[len];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atof(strtok(buf,"\t \n"));
  for(int i=1;i<(len);i++)
      tmp[i] = atof(strtok(NULL,"\t \n"));
  fclose(fp);
  delete [] buf;
  array<double> ret;
  ret.x = len;
  ret.data = tmp;
  return ret;
}

// read 2d sfs
matrix<double> readPrior12(const char *fname, int nrow, int ncol) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  double *tmp = new double[nrow*ncol];
  char *buf = new char[filesize];
  fread(buf,sizeof(char),filesize,fp);
  tmp[0] = atof(strtok(buf,"\t \n"));
  for(int i=1;i<(nrow*ncol);i++)
    tmp[i] = atof(strtok(NULL,"\t \n"));
  fclose(fp);
  delete [] buf;
  array<double> allvalues;
  allvalues.x = nrow*ncol;
  allvalues.data = tmp;
  int index=0;
  double **data = new double*[nrow];
  for(int i=0;i<nrow;i++){
    double *tmp2 = new double[ncol];
    for(int k=0;k<ncol;k++) {
      tmp2[k]=allvalues.data[index];
      index=index+1;
    }
    data[i]= tmp2;
  }
  matrix<double> ret;
  ret.x=nrow;
  ret.y=ncol;
  ret.data=data;
  delete [] allvalues.data;
  return ret;
}


// read a file of posterior probabilities into a matrix but only for a specific subsets of positions (0-based notation)
matrix<double> readFileSub(char *fname, int nInd, int start, int end) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  int n_categ = 2*nInd+1;

  if( strcmp(fname,"-")!=0 ) {
    if( filesize % (n_categ*sizeof(float)) != 2*sizeof(float) ) {
      fprintf(stderr,"\n\t-> Possible error reading SFS, binary file might be broken...\n");
      exit(-1);
    }
  }

  int nsites = end-start+1;
  double **data = new double*[nsites];
  float float_tmp = 0;

  // Locate data to read
  fseek(fp, (2+n_categ*start)*sizeof(float), SEEK_SET);

  // Read data
  for(int i=0; i<nsites; i++) {
    double *tmp = new double[n_categ];
    for(int c=0; c<n_categ; c++){
      fread(&float_tmp, sizeof(float), 1, fp);
      tmp[c] = (double) float_tmp;
    }
    data[i] = tmp;
  }

  fclose(fp);
  matrix<double> ret;
  ret.x = nsites;
  ret.y = n_categ;
  ret.data = data;
  return ret;
}

// read genotype posterior probabilities from angsd (-dogeno 64), but only for a specific subsets of positions (0-based notation)
matrix<double> readEstiSub(char *fname, int nInd, int start, int end) {
  FILE *fp = getFILE(fname,"rb");
  size_t filesize =fsize(fname);
  if( strcmp(fname,"-")!=0 && (filesize %(sizeof(double)*3*nInd)) ) {
    fprintf(stderr,"\n\t-> Possible error read GENO, binary file might be broken...\n");
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
  if( strcmp(fname,"-")!=0 && filesize==0){
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
int maxposarr(matrix<double> &m, int row, int start, int end) {
  double val = m.data[row][start];

  for (int i = start; i < end; i++) {
    if (m.data[row][i] > val) {
      start = i;
      val = m.data[row][i];
    }
  }

  return start;
}

// normalize SFS and exp it if log
void normSFS(matrix<double> &sfs, bool islog) {
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
void getPvar(matrix<double> &sfs, array<double> pvar) {
  for (int i=0; i<sfs.x; i++)
    pvar.data[i]=1-sfs.data[i][0]-sfs.data[i][2*(sfs.y-1)];
}

// write an array of doubles into a file
void writearray(array<double> &m,FILE *fp) {
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%f\t",m.data[i]);
  fprintf(fp,"\n");
}

// write an array of ints into a file
void writearrayInt(array<int> &m,FILE *fp) {
  for(int i=0;i<m.x;i++)
    fprintf(fp,"%d\t",m.data[i]);
  fprintf(fp,"\n");
}

// write a matrix of doubles into a file
void writematrix(matrix<double> &m,FILE *fp) {
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%f\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// write a matrix of ints into a file
void writematrixInt(matrix<int> &m,FILE *fp) {
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%d\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}

// add prior
void addPrior(matrix<double> &sfs, array<double> prior) {
  int nsites = sfs.x;
  int ncol = sfs.y;
  for (int j=0; j<nsites; j++) {
    for(int i=0; i<ncol; i++) {
        sfs.data[j][i] = sfs.data[j][i]*prior.data[i];
    }
  }
}
