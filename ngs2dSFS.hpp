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

/// FUNCTIONS

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

// to append names
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

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

// return max value position of a row from a matrix
int maxposarr(matrix<double> &m, int row) {
  int i=0, res = 0;
  double val;
  val = m.data[row][0];
  for (i = 1; i < m.y; i++) {
    if (m.data[row][i] > val) {
      res = i;
      val = m.data[row][i];
    }
  }
  return res;
}

// comput 2d-sfs, adding each site the joint max post prob of allele frequencies
void sumSpectrum(matrix<double> &spec, matrix<double> &m1, matrix<double> &m2, int maxlike) {
  int nsites = 0;
  nsites = m1.x;
  if (maxlike==1) { // get the max like joint alle freq
    for (int s=0; s<nsites; s++) spec.data[maxposarr(m1, s)][maxposarr(m2, s)]+=1.0;
  } else { // sum the products
   for (int s=0; s<nsites; s++) {
    for (int i=0; i<m1.y; i++) {
      for (int j=0; j<m2.y; j++) {
        spec.data[i][j]+=(m1.data[s][i]*m2.data[s][j]); // in case of folded, the "derived" stated will have 0 frequency
      }
    }
   }
  }
}

// write a matrix into a file
void writematrix(matrix<double> &m,FILE *fp){ 
  for(int i=0;i<m.x;i++){
    for(int j=0;j<m.y;j++)
      fprintf(fp,"%f\t",m.data[i][j]);
    fprintf(fp,"\n");
  }
}






