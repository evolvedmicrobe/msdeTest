////////////////////////////////////////////////////////////////

#ifndef mvnUtils_h
#define mvnUtils_h 1


#include <Rcpp.h>
#include "swappableDefinitions.h"
using namespace Rcpp;


// a few utility functions for multivariate normals

////////////////////////////////////////////////////

// class for proposal mean and variances (or rather sd's)
// also create a dummy space of size mean.
class propMV {
public:
  int nDims;
  double *mean, *sd, *z;
  propMV(int);
  ~propMV();
};

////////////////////////////////////////////////////

/*
bool inRange(double x[], double lBound[], double uBound[], int n);
bool isNaN(double x[], int n);
*/

// cholesky decomposition of a symmetric, positive definite matrix.
// returns a vector of the *Upper Triangular* cholesy factor, leaving the other elements of the array unchanged.
// in other words, U' %*% U \neq A, but lowTri(U') %*% upTri(U) = A.
// both U and A are stacked by COLUMN
// can be performed IN-PLACE, i.e., with U and A refering to same memory location.
void chol_decomp(double *U, double *A, int n);

// x = sd * z + mean.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void xmvn(double *x, double *z, double *mean, double *cholSd, int n);

// z = sd^{-1} * (x - mean).  only calculates first nMax values of z.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void zmvn(double *z, double *x, double *mean, double *cholSd, int n, int nMax);

// log-normal density evaluation.  z[] is required as temporary storage of residuals.
// i.e., z = sd^{-1} * (x - mean)
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
// TODO: include pi factor
double lmvn(double *x, double *z, double *mean, double *cholSd, int n);

//////////////////////////////////////////////

// log-density of a Gaussian copula
// TODO: include pi factor
double lgcop(double *x, double *qnorm, int *nBreaks, double *range,
	     double *dx, double *pdf, double *lpdf, double *cdf,
	     double *mean, double *sd, double *RhoCholSd, int n);

//////////////////////////////////////////////

// scalar multiplication

// multiply vector by scalar a
inline void v_mult(double *v, double a, int n) {
  for(int ii=0; ii<n; ii++) {
    v[ii] *= a;
  }
  return;
}
// multiply an upper triangular matrix by scalar a
inline void U_mult(double *U, double a, int n) {
  int ii,jj,colI;
  for(ii=0; ii<n; ii++) {
    colI = ii*n;
    for(jj=0; jj<=ii; jj++) {
      U[colI+jj] *= a;
    }
  }
  return;
}

#endif
