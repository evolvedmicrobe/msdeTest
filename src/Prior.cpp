///////////////////////////////////////

#include "Prior.h"

// prior classes
//
// there are four types of priors: flat, normal, gaussian copula, and custom
// each prior calls two arguments: x and theta;
// tuning parameters are set by the constructor.
// in addition to tuning parameters, the normal and gcop constructors
// take the following dimension arguments: nD = nDims, nTh = nTheta,
// and nRv <= nDims + nTh, the number of "active variables".
// the inputs are ligned up as tmpX = (x, theta), and the active variables are the _last_
// nRV of these.

///////////////////////////////////////

// constructor, destructor, and logPrior for Flat Prior

FlatPrior::FlatPrior() {
}
FlatPrior::~FlatPrior() {
}
double FlatPrior::logPrior(double *theta, double *x) {
  return 0.0;
}


///////////////////////////////////////

// constructor, destructor, and logPrior for Normal Prior

NormalPrior::NormalPrior(List priorParams, int nRv) {
  nTotalRV = nDims + nParams;
  nActiveRV = nRv;
  nSkippedRV = nTotalRV - nActiveRV;
  mean = REAL(priorParams["Mu"]);
  cholSd = REAL(priorParams["V"]);
  tmpX = new double[nTotalRV];
  tmpZ = new double[nTotalRV];
}

NormalPrior::~NormalPrior() {
  delete [] tmpX;
  delete [] tmpZ;
}

double NormalPrior::logPrior(double *theta, double *x) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nParams; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lmvn(&tmpX[nSkippedRV], &tmpZ[nSkippedRV], mean, cholSd, nActiveRV);
  return(lp);
}


///////////////////////////////////////

// constructor, destructor, and logPrior for Gaussian Copula Prior

GCopPrior::GCopPrior(List priorParams, int nRv) {
  nTotalRV = nDims + nParams;
  nActiveRV = nRv;
  nSkippedRV = nTotalRV - nActiveRV;
  nBreaks = INTEGER(priorParams["nbreaks"]);
  range = REAL(priorParams["rx"]);
  dx = REAL(priorParams["dx"]);
  pdf = REAL(priorParams["dens.y"]);
  lpdf = REAL(priorParams["ldens.y"]);
  cdf = REAL(priorParams["Dens.y"]);
  mean = REAL(priorParams["mean"]);
  sd = REAL(priorParams["sd"]);
  RhoCholSd = REAL(priorParams["Rho"]);
  tmpX = new double[nTotalRV];
  tmpQ = new double[nTotalRV];
}

GCopPrior::~GCopPrior() {
  delete [] tmpX;
  delete [] tmpQ;
}

double GCopPrior::logPrior(double *theta, double *x) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nParams; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lgcop(&tmpX[nSkippedRV], &tmpQ[nSkippedRV],
	     nBreaks, range, dx, pdf, lpdf, cdf,
	     mean, sd, RhoCholSd, nActiveRV);
  return(lp);
}
