////////////////////////////////////////////////////////////////

// a few utility functions for multivariate normals

#include "mvnUtils.h"

////////////////////////////////////////////////////////////////

// proMV class: mean, variance, and temporary z storage
// constructor and destructor for propMV class

propMV::propMV(int d) {
  nDims = d;
  mean = new double[nDims];
  sd = new double[nDims*nDims];
  z = new double[nDims];
  // initialize
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    mean[ii] = 0.0;
    z[ii] = 0.0;
  }
  for(ii = 0; ii < nDims*nDims; ii++) {
    sd[ii] = 0.0;
  }
}

propMV::~propMV() {
  delete [] mean;
  delete [] sd;
  delete [] z;
}

////////////////////////////////////////////////////////////////

/*
bool inRange(double x[], double lBound[], double uBound[], int n) {
  int ii = 0;
  while(x[ii] >= lBound[ii] && x[ii] <= uBound[ii] && ii < n) ii++;
  return(ii == n);
}

bool isNaN(double x[], int n) {
  int ii = 0;
  while(x[ii] == x[ii] &&
	((2.0 * x[ii] != x[ii] && -2.0 * x[ii] != x[ii]) || x[ii] == 0.0)
	&& ii < n) ii++;
  return(ii < n);
}
*/

// cholesky decomposition of a symmetric, positive definite matrix.
// returns a vector of the *Upper Triangular* cholesy factor, leaving the other elements of the array unchanged.
// in other words, U' %*% U \neq A, but lowTri(U') %*% upTri(U) = A.
// both U and A are stacked by COLUMN
// can be performed IN-PLACE, i.e., with U and A refering to same memory location.
void chol_decomp(double *U, double *A, int n) {
  int ii, jj, kk, colI, colJ;
  double tmpSum, tmpInv;
  for(ii = 0; ii < n; ii++) {
    colI = ii*n;
    tmpSum = 0.0;
    for(kk = 0; kk < ii; kk++) {
      tmpSum += U[colI + kk] * U[colI + kk];
    }
    tmpInv = sqrt(A[colI + ii] - tmpSum);
    U[colI + ii] = tmpInv;
    tmpInv = 1.0/tmpInv;
    for(jj = ii+1; jj < n; jj++) {
      colJ = jj*n;
      tmpSum = 0.0;
      for(kk = 0; kk < ii; kk++) tmpSum += U[colJ + kk] * U[colI + kk];
      U[colJ + ii] = tmpInv * (A[colJ + ii] - tmpSum);
    }
  }
  return;
}

// x = sd * z + mean.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void xmvn(double *x, double *z, double *mean, double *cholSd, int n) {
  int ii, jj, colI;
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    x[ii] = 0;
    for(jj = 0; jj <= ii; jj++) x[ii] += cholSd[colI + jj] * z[jj];
    x[ii] += mean[ii];
  }
  return;
}

// z = sd^{-1} * (x - mean).  only calculates first nMax values of z.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
void zmvn(double *z, double *x, double *mean, double *cholSd, int n, int nMax) {
  int ii, jj, colI;
  double tmpSum;
  for(ii = 0; ii < nMax; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii = 0; ii < nMax; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  return;
}

// log-normal density evaluation.  z[] is required as temporary storage of residuals.
// i.e., z = sd^{-1} * (x - mean)
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
// TODO: include pi factor
double lmvn(double *x, double *z, double *mean, double *cholSd, int n) {
  int ii, jj, colI;
  double tmpSum;
  for(ii = 0; ii < n; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) tmpSum += z[ii] * z[ii];
  tmpSum *= 0.5;
  for(ii = 0; ii < n; ii++) tmpSum += log(cholSd[n*ii + ii]);
  return(-tmpSum);
}

// log-density of a Gaussian copula
// TODO: include pi factor
double lgcop(double *x, double *qNorm, int *nBreaks, double *range,
	     double *dx, double *pdf, double *lpdf, double *cdf,
	     double *mean, double *sd, double *RhoCholSd, int n) {
  int ii, jj, colI, densElt, start;
  double lp = 0.0;
  double tmpSum = 0.0;
  // normal quantiles and marginal components
  start = 0;
  for(ii = 0; ii < n; ii ++) {
    densElt = (int) floor((x[ii]-range[2*ii])/dx[ii]);
    if((densElt >= 0) & (densElt < nBreaks[ii])) {
      lp += lpdf[densElt + start];
      qNorm[ii] = x[ii] - (range[2*ii] + densElt * dx[ii]);
      qNorm[ii] *= pdf[densElt + start];
      qNorm[ii] = Rf_qnorm5(cdf[densElt + start + ii] +  qNorm[ii], 0.0, 1.0, 1, 0);
    }
    else {
      lp += Rf_dnorm4(x[ii], mean[ii], sd[ii], 1);
      qNorm[ii] = (x[ii] - mean[ii])/sd[ii];
    }
    start += nBreaks[ii];
  }
  // copula components
  // iid standard normal densities
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  lp += 0.5 * tmpSum;
  // multivariate normal density
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) {
      tmpSum += RhoCholSd[colI + jj] * qNorm[jj];
    }
    qNorm[ii] = (qNorm[ii] - tmpSum)/RhoCholSd[colI + ii];
  }
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  tmpSum *= 0.5;
  for(ii = 0; ii < n; ii++) {
    tmpSum += log(RhoCholSd[n*ii + ii]);
  }
  lp -= tmpSum;
  return(lp);
}
