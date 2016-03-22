
#include "shared.hpp"

// print help
void info() {
  fprintf(stdout, "\nInput:\n-postfiles: .sfs files with posterior probabilities of sample allele frequencies for each population [required]\n-priorfile: 2D-SFS to be used as a prior; you can use ngs2DSFS with parameter -relative set to 1 [NULL]\n-priorfiles: marginal spectra to be used as a prior; you can use optimSFS in ANGSD [NULL]\n-outfile: name of the output file [required]\n-nind: number of individuals for each population [required]\n-nsites: total number of sites; in case you want to analyze a subset of sites this is the upper limit [required]\n-verbose: level of verbosity [0]\n-block_size: to be memory efficient, set this number as the number of sites you want to analyze at each chunk [0]\n-firstbase: in case you want to analyze a subset of your sites this is the lower limit [1]\n\n");
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
