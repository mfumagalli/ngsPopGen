
/// TEMPLATES

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

array<int> getStart(int nsites, int firstbase, int block_size) {
  // note that firstbase and nsites are 1-based
  int len = nsites-firstbase+1;
  int nwin = len/block_size;
  if ( (len % block_size)!=0) { nwin=nwin+1; }
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
  if ( (len % block_size)!=0) { nwin=nwin+1; }
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

// a nice wrapper for getting files
FILE *getFILE(const char*fname,const char* mode) {
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  if(writeFile&&fexists(fname)){//DRAGON
    fprintf(stderr,"\t-> File exists: %s exiting...\n",fname);
    exit(0);
  }
  FILE *fp = fopen(fname, mode);
  if(fp==NULL){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    fclose(fp);
    exit(0);
  }
  return fp;
}

// read a file of posterior probabilities into a matrix but only for a specific subsets of positions (0-based notation)
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

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// print help
void info() {
  fprintf(stdout, "\nInput:\n-npop: how many pops (1 or 2)\n-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population (with or without running sfstools)\n-outfile: name of the output file\n-nind: number of individuals for each population\n-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit\n-verbose: level of verbosity, if 0 suppress all messages\n-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk\n-firstbase: in case you want to analyze a subset of your sites this is the lower limit\n-isfold: boolean, is your data folded or not?\n-islog: boolean, are postfiles in log (from -realSFS 1 only, required if 2D-SFS is given)? If you use sfstools then set -islog 1\n-iswin: if 1 then print the value computed for each non-overlapping window defined by block_size\n\n");
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

// compute summary stats for each block in case of 1 pop and print the results
void computeStats(matrix<double> &post1, int verbose, FILE *outpost, int iswin, int isfold, int start) {

  int nind=0; // sample size
  if (isfold) {
    nind=(post1.y-1);
  } else {
    nind=(post1.y-1)/2;
  }
  int nsites=post1.x;

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
    if (isfold==0) segsites.data[s] = segsites.data[s] - post1.data[s][post1.y-1];
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
void computeStats2Pops(matrix<double> &post1, int verbose, FILE *outpost, int iswin, int isfold, int start, matrix<double> &post2) {

  int nind1=0,nind2=0; // sample sizes
  if (isfold) {
    nind1=(post1.y-1);
    nind2=(post2.y-1);
  } else {
    nind1=(post1.y-1)/2;
    nind2=(post2.y-1)/2;
  }
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

  double temp=0.0;

  // to convert from 0 base to 1 base  
  int start0=start+1;

  double sum_segsites1=0.0, sum_hetero1=0.0;
  double sum_segsites2=0.0, sum_hetero2=0.0;
  double sum_fixed=0.0;

  for (int s=0; s<nsites; s++) {

    start++;

    segsites1.data[s]= 1.0 - post1.data[s][0];
    if (isfold==0) segsites1.data[s]=segsites1.data[s] - post1.data[s][post1.y-1];
    temp=0.0;
    for (int j=0; j<post1.y; j++) { temp=temp + 2.0*(j/(nind1*2.0))*((nind1*2.0-j)/(nind1*2.0))*post1.data[s][j]; }
    hetero1.data[s]=temp;
    sum_segsites1=sum_segsites1+segsites1.data[s];
    sum_hetero1=sum_hetero1+hetero1.data[s];

    segsites2.data[s]= 1.0 - post2.data[s][0];
    if (isfold==0) segsites2.data[s]=segsites2.data[s] - post2.data[s][post2.y-1];
    temp=0.0;
    for (int j=0; j<post2.y; j++) { temp=temp + 2.0*(j/(nind2*2.0))*((nind2*2.0-j)/(nind2*2.0))*post2.data[s][j]; }
    hetero2.data[s]=temp;
    sum_segsites2=sum_segsites2+segsites2.data[s];
    sum_hetero2=sum_hetero2+hetero2.data[s];

    if (isfold) {
      fixed.data[s]=-999.9;
    } else {
      fixed.data[s]=post1.data[s][0]*post2.data[s][post2.y-1] + post2.data[s][0]*post1.data[s][post1.y-1];
    }
    sum_fixed = sum_fixed + fixed.data[s];

    if (iswin==0) fprintf(outpost, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", start, start, segsites1.data[s], hetero1.data[s], segsites2.data[s], hetero2.data[s], fixed.data[s]);

  } // end cycle s

  // compute sum across all sites 
  if (iswin==1) fprintf(outpost, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", start0, start, sum_segsites1, sum_hetero1, sum_segsites2, sum_hetero2, sum_fixed);

  delete [] segsites1.data;
  delete [] hetero1.data;
  delete [] segsites2.data;
  delete [] hetero2.data;
  delete [] fixed.data;

}





