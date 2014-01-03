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

// read a file of prior into an array
array<double> readArray(const char *fname, int nInd, int isfold) {
  FILE *fp = getFILE(fname,"r");
  size_t filesize =fsize(fname);
  if(filesize==0){
    fprintf(stderr,"file:%s looks empty\n",fname);
    exit(0);
  }
  int len = 2*nInd+1;
  if (isfold) len = nInd+1;
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
  fprintf(stdout, "\nInput:\n-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population [required]\n-priorfile: 2D-SFS to be used as a prior; you can use ngs2DSFS with parameter -relative set to 1 [NULL]\n-priorfiles: marginal spectra to be used as a prior; you can use optimSFS in ANGSD [NULL]\n-outfile: name of the output file [required]\n-nind: number of individuals for each population [required]\n-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit [required]\n-verbose: level of verbosity [0]\n-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk [0]\n-firstbase: in case you want to analyze a subset of your sites this is the lower limit [1]\n-isfold: boolean, is your data folded? [0]\n-islog: boolean, are postfiles in -log (from -realSFS 1 only, required if 2D-SFS is given)? If you use sfstools then set to 1 [0]\n\n");
}

/// 

// compute a and ab (using short-cut formulas) given allele frequencies and sample sizes
array<double> calcAB(int s1, int s2, int n1, int n2, int debug) {
  if (debug) fprintf(stderr, "\nInside calcAB");
  double alfa1 = 0.0, alfa2 = 0.0, p1 = 0.0, p2 = 0.0, a = 0.0, ab = 0.0;
  array<double> res; res.x = 2;
  double *values= new double [2];
  if (debug) fprintf(stderr, "%f\t%f\t%f\t%f\t%f\t%f\n", alfa1, alfa2, p1,p2,a,ab);
  p1 = static_cast<double>(s1) / (2*n1);
  p2 = static_cast<double>(s2) / (2*n2);
  if (debug) fprintf(stderr, "%f\t%f\n", p1, p2);
  alfa1 = 2*p1*(1-p1);  
  alfa2 = 2*p2*(1-p2);  
  if (debug) fprintf(stderr, "%f\t%f\n", alfa1, alfa2);
  a = ((p1-p2)*(p1-p2)) - (((n1+n2) * (n1 * alfa1 + n2 * alfa2)) / (4*n1*n2*(n1+n2-1)));
  ab = ((p1-p2)*(p1-p2)) + ( ( (4*n1*n2 - n1 - n2) * (n1 * alfa1 + n2 * alfa2) ) / ((4*n1*n2)*(n1+n2-1)) );
  values[0]=a;
  values[1]=ab;
  if (debug) fprintf(stderr, "%f\t%f\n", a, ab);
  res.data = values;
  return res;
}

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, but no correction from a first guess of fst
void computeVarRey(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==2) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1)/2; // nind1
  n2 = (m2.y-1)/2; // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  array<double> temp;
  for (int s=0; s<nsites; s++) {
    if (verbose==2) fprintf(stderr, "\t s %d", s);
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1*2)+1;
    A.y=AB.y=(n2*2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1*2)+1];
    double **dataAB = new double*[(n1*2)+1];
    // get also the probability of site being variable
    double pvar = 0.0;
    pvar = 1 - m1.data[s][0]*m2.data[s][0] - m1.data[s][2*n1]*m2.data[s][2*n2];
    if (verbose==2) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(2*n1+1); i++) {
      if (verbose==2) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2*2)+1];
      double *bufAB = new double[(n2*2)+1];
      for (int j=0; j<(2*n2+1); j++) {
       if (verbose==2) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    if (EA<0.0) EA=0.0; // protect against -0.000 cases
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(2*n1+1); i++) {
        for (int j=0; j<(2*n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
  } // end for s in nsites
} // end


// compute a and ab estimate for all possible combinations of sample size, then weight by their prob, , but no correction from a first guess of fst, here if folded!
void computeVarReyFold(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums) {
  // m1 and m2 are post probs
  /// DEFINITION (DECLARATION AND INITIALIZATION)
  if (verbose==2) fprintf(stderr, "\nInside computeVarRey");
  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites
  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1); // nind1
  n2 = (m2.y-1); // nind2
  nsites = m1.x; // nsites from file
  double VAR = 0.0, COVAR = 0.0, FACT = 0.0;
  matrix<double> A;
  matrix<double> AB;
  array<double> temp;
  double pvar = 0.0;
  for (int s=0; s<nsites; s++) {
    if (verbose==2) fprintf(stderr, "\t s %d", s);
    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    A.x=AB.x=(n1)+1;
    A.y=AB.y=(n2)+1;
    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1)+1];
    double **dataAB = new double*[(n1)+1];
    // get also the probability of site being variable
    pvar = 1 - (m1.data[s][0]*m2.data[s][0]);
    if (verbose==2) fprintf(stderr, "\t first cycle");
    double EA = 0.0, EAB = 0.0;
    for (int i=0; i<(n1+1); i++) {
      if (verbose==2) fprintf(stderr, "\ti%d",i);
      double *bufA = new double[(n2)+1];
      double *bufAB = new double[(n2)+1];
      for (int j=0; j<(n2+1); j++) {
       if (verbose==2) fprintf(stderr, "\tj%d",j);
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];
        EA = EA + temp.data[0]*m1.data[s][i]*m2.data[s][j];
        EAB = EAB + temp.data[1]*m1.data[s][i]*m2.data[s][j];
        delete [] temp.data;
      } // end for in j
      dataA[i]=bufA;
      dataAB[i]=bufAB;
    } // end for in i
    A.data=dataA;
    AB.data=dataAB;
    if (EA<0.0) EA=0.0; // protect against -0.000 cases
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;
    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);
    for (int q=1; q<=nsums; q++) {
      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(n1+1); i++) {
        for (int j=0; j<(n2+1); j++) {
          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m1.data[s][i]*m2.data[s][j];
        } // end for in j
      } // end for in i
      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));
      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);
    } // end for in nsums
    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);
    cleanup(A);
    cleanup(AB);
  } // end for s in nsites
} // end

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

// compute a and ab estimate for all possible combinations of sample size, then weight by their prob12 computed from post1, pos12 and prior12 (normalize it)
void computeVarRey12(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums, matrix<double> &p12) {

  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites

  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1)/2; // nind1
  n2 = (m2.y-1)/2; // nind2
  nsites = m1.x; // nsites from file

  double VAR = 0.0, COVAR = 0.0, FACT = 0.0, pvar = 0.0;

  for (int s=0; s<nsites; s++) {

    // m1 and m2 are post probs, p12 is the 2d sfs
    matrix<double> m12;
    m12.x=m1.y;
    m12.y=m2.y;
    double **ddata = new double*[m12.x];
      for(int i=0;i<m12.x;i++){
        double *dtmp = new double[m12.y];
        for(int k=0;k<m12.y;k++) {
          dtmp[k]=0.0;
        }
        ddata[i]= dtmp;
      }
    m12.data=ddata;

    // multiply
    for (int j=0; j<m12.x; j++) {
      for (int i=0;i<m12.y; i++) {
        m12.data[j][i] = m1.data[s][j]*m2.data[s][i]*p12.data[j][i];
      }
    }

    // get the sum
    double somma=0.0;
      for (int j=0; j<m12.x; j++) {
        for (int i=0;i<m12.y; i++) {
          somma=somma+m12.data[j][i];
      }
    }

    // divide
    for (int j=0; j<m12.x; j++) {
      for (int i=0;i<m12.y; i++) {
        m12.data[j][i] = m12.data[j][i]/somma;
      }
    }

    //fprintf(stderr, "final2d %f %f %f %f\n", m12.data[0][0], m12.data[0][1], m12.data[1][0], m12.data[1][1]);  

    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    matrix<double> A;
    matrix<double> AB;

    A.x=AB.x=(n1*2)+1;
    A.y=AB.y=(n2*2)+1;

    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[(n1*2)+1];
    double **dataAB = new double*[(n1*2)+1];

    // get also the probability of site being variable
    pvar = 1 - m12.data[0][0] - m12.data[2*n1][2*n2];

    if (verbose==2) fprintf(stderr, "\t first cycle");

    double EA = 0.0, EAB = 0.0;

    for (int i=0; i<(2*n1+1); i++) {

      if (verbose==2) fprintf(stderr, "\ti%d",i);

      double *bufA = new double[(n2*2)+1];
      double *bufAB = new double[(n2*2)+1];

      for (int j=0; j<(2*n2+1); j++) {

        array<double> temp;
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];

        EA = EA + temp.data[0]*m12.data[i][j];
        EAB = EAB + temp.data[1]*m12.data[i][j];

        delete [] temp.data;

      } // end for in j

      dataA[i]=bufA;
      dataAB[i]=bufAB;

    } // end for in i

    A.data=dataA;
    AB.data=dataAB;

    if (EA<0.0) EA=0.0; // protect against -0.000 cases
    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;

    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);

    for (int q=1; q<=nsums; q++) {

      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(2*n1+1); i++) {
        for (int j=0; j<(2*n2+1); j++) {

          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) )*m12.data[i][j];

          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m12.data[i][j];

        } // end for in j
      } // end for in i

      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));

      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);

    } // end for in nsums

    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);

    cleanup(A);
    cleanup(AB);
    cleanup(m12);

  } // end for s in nsites



} // end


// compute a and ab estimate for all possible combinations of sample size, then weight by their prob12 computed from post1, pos12 and prior12 (normalize it)
void computeVarRey12Fold(matrix<double> &m1, matrix<double> &m2, int verbose, FILE *fname, int nsums, matrix<double> &p12) {

  int n1 = 0, n2 = 0; // sample sizes
  int nsites = 0; // nsites

  // estimates nsites and pop sizes from matrices
  n1 = (m1.y-1); // nind1
  n2 = (m2.y-1); // nind2
  nsites = m1.x; // nsites from file

  double VAR = 0.0, COVAR = 0.0, FACT = 0.0, pvar = 0.0;

  for (int s=0; s<nsites; s++) {

    // m1 and m2 are post probs, p12 is the 2d sfs
    matrix<double> m12;
    m12.x=m1.y;
    m12.y=m2.y;
    double **ddata = new double*[m12.x];
      for(int i=0;i<m12.x;i++){
        double *dtmp = new double[m12.y];
        for(int k=0;k<m12.y;k++) {
          dtmp[k]=0.0;
        }
        ddata[i]= dtmp;
      }
    m12.data=ddata;

    // multiply
    for (int j=0; j<m12.x; j++) {
      for (int i=0;i<m12.y; i++) {
        m12.data[j][i] = m1.data[s][j]*m2.data[s][i]*p12.data[j][i];
      }
    }

    // get the sum
    double somma=0.0;
      for (int j=0; j<m12.x; j++) {
        for (int i=0;i<m12.y; i++) {
          somma=somma+m12.data[j][i];
      }
    }

    // divide
    for (int j=0; j<m12.x; j++) {
      for (int i=0;i<m12.y; i++) {
        m12.data[j][i] = m12.data[j][i]/somma;
      }
    }

    //fprintf(stderr, "final2d %f %f %f %f\n", m12.data[0][0], m12.data[0][1], m12.data[1][0], m12.data[1][1]);

    // for each possible value of freq 1 and freq 2 compute the FST, so compute A, AB, VAR, COVAR (see Price paper for its meaning)
    matrix<double> A;
    matrix<double> AB;

    A.x=AB.x=n1+1;
    A.y=AB.y=n2+1;

    // FIRST CYCLE: get expected A and AB and retain matrices of A and AB
    double **dataA = new double*[n1+1];
    double **dataAB = new double*[n1+1];

    // get also the probability of site being variable
    pvar = 1 - m12.data[0][0];

    if (verbose==2) fprintf(stderr, "\t first cycle");

    double EA = 0.0, EAB = 0.0;

    for (int i=0; i<(n1+1); i++) {

      if (verbose==2) fprintf(stderr, "\ti%d",i);

      double *bufA = new double[(n2)+1];
      double *bufAB = new double[(n2)+1];

      for (int j=0; j<(n2+1); j++) {

        array<double> temp;
        temp = calcAB(i, j, n1, n2, 0);
        bufA[j]=temp.data[0];
        bufAB[j]=temp.data[1];

        EA = EA + temp.data[0]*m12.data[i][j];
        EAB = EAB + temp.data[1]*m12.data[i][j];

        delete [] temp.data;

      } // end for in j

      dataA[i]=bufA;
      dataAB[i]=bufAB;

    } // end for in i

    A.data=dataA;
    AB.data=dataAB;

    if (EA<0.0) EA=0.0; // protect against -0.000 cases

    // SECOND CYCLE: get VAR and COVAR, and then the correcting FACTor, according to number of sums to retain
    VAR = 0.0, COVAR = 0.0, FACT = 0.0;

    if ((verbose==4) & (s==0)) fprintf(stderr, "\t second cycle %d", nsums);

    for (int q=1; q<=nsums; q++) {

      VAR = 0.0; COVAR = 0.0;
      for (int i=0; i<(n1+1); i++) {
        for (int j=0; j<(n2+1); j++) {

          VAR = VAR + pow((AB.data[i][j]-EAB), static_cast <double> (q) )*m12.data[i][j];

          COVAR = COVAR + pow((AB.data[i][j]-EAB)*(A.data[i][j]-EA), static_cast <double> (q) ) *m12.data[i][j];

        } // end for in j
      } // end for in i

      FACT = FACT + pow( (-1.0), static_cast <double> (q)) *( (EA*VAR + COVAR) / pow(EAB, static_cast <double> (q+1)));

      if ((verbose==4) & (s==0)) fprintf(stderr, "\n q %d v %f c %f f %f",q,VAR,COVAR,FACT);

    } // end for in nsums

    // print results
    fprintf(fname, "%f\t%f\t%f\t%f\t%f\n", EA, EAB, FACT, (EA/EAB)+FACT, pvar);

    cleanup(A);
    cleanup(AB);
    cleanup(m12);

  } // end for s in nsites



} // end

