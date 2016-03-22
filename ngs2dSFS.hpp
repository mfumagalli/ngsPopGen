#include "shared.hpp"

// comput 2d-sfs, adding each site the joint max post prob of allele frequencies
void sumSpectrum(matrix<double> &spec, matrix<double> &m1, matrix<double> &m2, int maxlike) {
  int nsites = 0;
  nsites = m1.x;
  if (maxlike==1) { // get the max like joint alle freq
    for (int s=0; s<nsites; s++) spec.data[maxposarr(m1,s,0,m1.y)][maxposarr(m2,s,0,m2.y)]+=1.0;
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
