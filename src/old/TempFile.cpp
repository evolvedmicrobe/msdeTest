#include <Rcpp.h>
using namespace Rcpp;
/////////////////////////////////////////



// sdeModel class: constains functions unique to a given sde

class sdeModel {
  // put private storage variables here
  // end private storage
public:
  static const int nParams = 5;
  static const int nDims = 2;
  friend void sdeDr(double *dr, double *x, double *theta, sdeModel *sde);
  friend void sdeDf(double *df, double *x, double *theta, sdeModel *sde);
  static bool isValidData(double *x);
  static bool isValidParams(double *theta);
  sdeModel();
  ~sdeModel();
};

////////////////////////////////////////////////////////////////



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

////////////////////////////////////////////////////////

// utility functions for sdes



// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double *mean, double *sd,
	     double *x0, double t, double t2, double *theta,
	     sdeModel *sde);

// loglikelihood function (actually log-density, minus factor of pi)
double loglik(double *theta, double *x, double *dT, double *sqrtDT,
	      propMV **mvX, sdeModel *sde, int nComp);

/////////////////////////////////////////////////////////////////////////////



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

//base prior class
class Prior {
 protected:
  static const int nDims = sdeModel::nDims;
  static const int nParams = sdeModel::nParams;
 public:
  // the type of prior
  enum Type {
    Flat = 1,
    Normal = 2,
    GCop = 3,
    Custom = 4
  };

  //pure virtual function must be defined in each derived class
  virtual double logPrior(double *theta, double *x)=0;
};

// flat prior
class FlatPrior : public Prior {
 public:
  FlatPrior();
  ~FlatPrior();
  double logPrior(double *theta, double *x);
};

// normal prior
class NormalPrior : public Prior {
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  NormalPrior(List priorParams, int nRv);
  ~NormalPrior();
  double logPrior(double *theta, double *x);
};

// gaussian copula prior
class GCopPrior : public Prior {
  // gaussian copula parameters
  int *nBreaks;
  double *range, *dx, *pdf, *lpdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  double *tmpX, *tmpQ;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  GCopPrior(List priorParams, int nRv);
  ~GCopPrior();
  // parameter + unobserved first states gaussian copula prior.
  double logPrior(double *theta, double *x);
};

/////////////////////////////////////////

// Utilities for MCMC



// x is either:
// T/F: in which case it's a yes or no
// integer, in which case it's every so many iterations (i.e., iter is a multiple of x)
// fraction between 0 an 1, in which case it's randomly updated that fraction of the time
bool updateComponent(double x, int iter);

////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////

// class and functions for MCMC object.
// class contains all the storage etc. required to run the update methods.

// in terms of storage in {curr/prop}Full, let's put theta before x.

class sdeMCMC {
  int *missInd;
  int nMiss, nMiss0;
  Prior *prior;
public:
  int nComp, nDims, nParams;
  double *currFull, *propFull, *propAccept;
  double *currX, *propX, *currTheta, *propTheta;
  double *dT, *sqrtDT, *B, *sqrtB;
  int *nObsComp;
  bool *fixedTheta;
  propMV **mvX;
  propMV *mvTheta;
  sdeModel *sde;
  void missGibbsUpdate(double *jumpSd, int *gibbsAccept, int *paramAccept);
  void paramVanillaUpdate(double *jumpSd, int *paramAccept);
  // internal loglikelihood
  double loglik(double *theta, double *x) {
    return ::loglik(theta, x, dT, sqrtDT, mvX, sde, nComp);
  };
  sdeMCMC(int N, double *dt, double *xInit, double *thetaInit,
	  int *xIndex, bool *thetaIndex, Prior *priorIn);
  // thin initialization, for loglikelihoods only
  sdeMCMC(int N, double *dt);
  ~sdeMCMC();
};

////////////////////////////////////////////////////

// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
void mvEraker(double *mean, double *sd,
	      double *x0, double *x2, double b, double b2, double *theta,
	      sdeModel *sde);



///////////////////////////////////////////////////////////


// custom functions for heston's model

// constructor
sdeModel::sdeModel() {
}

// destructor
sdeModel::~sdeModel() {} // do nothing

// todo: drop t and sqrtT in these functions
void sdeDr(double *dr, double *x, double *theta, sdeModel *sde) {
  dr[0] = (theta[0] - .125 * x[1]*x[1]); // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x[1]); // z
  return;
}

void sdeDf(double *df, double *x, double *theta, sdeModel *sde) {
  df[0] = .5 * x[1];
  df[2] = theta[3];
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

bool sdeModel::isValidData(double *x) {
  return(x[1] > 0.0);
}

bool sdeModel::isValidParams(double *theta) {
  int isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}
//////////////////////////////////////////////////////

// utility functions for sdes


// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double *mean, double *sd,
	     double *x0, double t, double t2, double *theta,
	     sdeModel *sde) {
  sdeDr(mean, x0, theta, sde);
  v_mult(mean, t, sdeModel::nDims);
  for(int jj = 0; jj < sdeModel::nDims; jj++) {
    mean[jj] += x0[jj];
  }
  sdeDf(sd, x0, theta, sde);
  U_mult(sd, t2, sdeModel::nDims);
  return;
}

////////////////////////////////////////////////////////////////////////////////

// loglikelihood function

double loglik(double *theta, double *x, double *dT, double *sqrtDT,
	      propMV **mvX, sdeModel *sde, int nComp) {
  double ll = 0;
  // *** PARALLELIZABLE FOR-LOOP ***
  for(int ii = 0; ii < nComp-1; ii++) {
    mvEuler(mvX[ii]->mean, mvX[ii]->sd, &x[ii*sdeModel::nDims],
	    dT[ii], sqrtDT[ii], theta, &sde[ii]);
    ll += lmvn(&x[(ii+1)*sdeModel::nDims], mvX[ii]->z, mvX[ii]->mean, mvX[ii]->sd,
	       sdeModel::nDims);
  }
  return(ll);
}
///////////////////////////////////////


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
////////////////////////////////////////////////////////////////

// a few utility functions for multivariate normals


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
/////////////////////////////////////////

// Utilities for MCMC


// x is either:
// T/F: in which case it's a yes or no
// integer, in which case it's every so many iterations (i.e., iter is a multiple of x)
// fraction between 0 an 1, in which case it's randomly updated that fraction of the time
bool updateComponent(double x, int iter) {
  bool update = false;
  if(x > 0.0) {
    if(x < 1.0) {
      update = unif_rand() <= x;
    }
    else {
      update = (iter % (int) x) == 0;
    }
  }
  return update;
}
////////////////////////////////////////////


// functions for class sdeMCMC and associated functions

/*
class sdeMCMC {
  int *missInd;
  int nMiss, nMiss0;
  Prior *prior;
public:
  int nComp, nDims, nParams;
  double *currFull, *propFull;
  double *currX, *propX, *currTheta, *propTheta;
  double *dT, *sqrtDT, *B, *sqrtB;
  int *nObsComp;
  bool *fixedTheta;
  propMV *mvX;
  propMV *mvTheta;
  sdeModel *sde;
  void missGibbsUpdate(double *jumpSd, int *gibbsAccept, int *paramAccept);
  void paramVanillaUpdate(double *jumpSd, int *paramsAccept);
  double loglik(double *theta, double *x);
  sdeMCMC(int N, double *dt, double *xInit, double *thetaInit,
	  int *xIndex, int *thetaIndex);
  ~sdeMCMC();
};
*/

///////////////////////////////////////////////////////////////

// class constructor & destructor

// small initialization for loglikelihoods
sdeMCMC::sdeMCMC(int N, double *dt) {
  int ii;
  // problem dimensions
  nComp = N;
  nDims = sdeModel::nDims;
  nParams = sdeModel::nParams;
  //times
  dT = dt;
  sqrtDT = new double[nComp];
  B = new double[1];
  sqrtB = new double[1];
  for(ii=0; ii<nComp-1; ii++) {
    sqrtDT[ii] = sqrt(dT[ii]);
  }
  // data
  currFull = new double[1];
  propFull = new double[1];
  currX = currFull + nParams;
  propX = propFull + nParams;
  propAccept = new double[1];
  mvX = new propMV*[nComp];
  for(ii=0; ii<nComp; ii++) {
    mvX[ii] = new propMV(nDims);
  }
  missInd = new int[1];
  // parameters
  currTheta = currFull;
  propTheta = propFull;
  // TODO: prior
  // sde
  sde = new sdeModel[nComp];
}

sdeMCMC::sdeMCMC(int N, double *dt, double *xInit, double *thetaInit,
		 int *xIndex, bool *thetaIndex, Prior *priorIn) {
  int ii, jj;
  // problem dimensions
  nComp = N;
  nDims = sdeModel::nDims;
  nParams = sdeModel::nParams;
  // times
  dT = dt;
  sqrtDT = new double[nComp];
  B = new double[nComp];
  sqrtB = new double[nComp];
  for(ii=0; ii<nComp-1; ii++) {
    sqrtDT[ii] = sqrt(dT[ii]);
    if(ii > 0) {
      B[ii] = dT[ii]/(dT[ii-1] + dT[ii]);
      sqrtB[ii] = sqrt((1-B[ii]) * dT[ii]);
    }
  }
  // data
  currFull = new double[nComp*nDims + nParams];
  propFull = new double[nComp*nDims + nParams];
  propAccept = new double[nComp];
  currX = currFull + nParams;
  propX = propFull + nParams;
  mvX = new propMV*[nComp];
  for(ii=0; ii<nComp; ii++) {
    propAccept[ii] = 0.0;
    mvX[ii] = new propMV(nDims);
    for(jj=0; jj<nDims; jj++) {
      currX[ii*nDims + jj] = xInit[ii*nDims + jj];
      propX[ii*nDims + jj] = currX[ii*nDims + jj];
    }
  }
  // missing data
  nObsComp = xIndex;
  int nMiss0 = nDims-nObsComp[0]; // unobserved states in first observation
  int nMissN = nDims-nObsComp[nComp-1]; // unobserved states in last observation
  // identify missing data indices, i.e. at least one component to update
  nMiss = 0;
  for(ii = 0; ii < nComp; ii++) {
    nMiss += (nObsComp[ii] < nDims);
  }
  missInd = new int[nMiss + (nMiss == 0)];
  jj = 0;
  for(ii = 0; ii < nComp; ii++) {
    if(nObsComp[ii] < nDims) {
      missInd[jj++] = ii;
    }
  }
  // parameters
  fixedTheta = thetaIndex;
  currTheta = currFull;
  propTheta = propFull;
  for(ii=0; ii<nParams; ii++) {
    currTheta[ii] = thetaInit[ii];
    propTheta[ii] = currTheta[ii];
  }
  // prior
  prior = priorIn;
  // sde
  sde = new sdeModel[nComp];
}

sdeMCMC::~sdeMCMC() {
  delete [] sqrtDT;
  delete [] B;
  delete [] sqrtB;
  delete [] currFull;
  delete [] propFull;
  delete [] propAccept;
  delete [] missInd;
  for(int ii=nComp-1; ii>=0; ii--) {
    delete mvX[ii];
  }
  delete [] mvX;
  delete [] sde;
}

///////////////////////////////////////////////////////////////

// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
void mvEraker(double *mean, double *sd,
	      double *x0, double *x2, double b, double b2, double *theta,
	      sdeModel *sde) {
  for(int jj = 0; jj < sdeModel::nDims; jj++) {
    mean[jj] = x0[jj] * b + x2[jj] * (1-b);
  }
  sdeDf(sd, x0, theta, sde);
  U_mult(sd, b2, sdeModel::nDims);
  return;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////

// componentwise data updates.


void sdeMCMC::missGibbsUpdate(double *jumpSd, int *gibbsAccept, int *paramAccept) {
  int ii, II, jj, JJ;
  int startII, endII;
  //int foundNaN = 0;
  // only elements in missInd are updated.
  // first and last components are handled separately.
  startII = (missInd[0] == 0) ? 1 : 0;
  endII = (missInd[nMiss-1] == nComp-1) ? nMiss-1 : nMiss;
  // Markov chain elements are conditionally independent,
  // so every other can be updated in parallel
  for(JJ = 0; JJ < 2; JJ++) {
    // *** PARALLELIZABLE FOR-LOOP ***
    for(II = startII+JJ; II < endII; II = II+2) {
      ii = missInd[II];
      // intermediate data points
      if((0 < ii) && (ii < nComp-1)) {
	mvEraker(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
		 &currX[(ii+1)*nDims], B[ii], sqrtB[ii], currTheta, &sde[ii]);
	// partial observations
	if(nObsComp[ii] == 0) {
	  for(jj = 0; jj < nDims; jj++) mvX[ii]->z[jj] = norm_rand();
	}
	else {
	  zmvn(mvX[ii]->z, &currX[ii*nDims], mvX[ii]->mean,
	       mvX[ii]->sd, nDims, nObsComp[ii]);
	  for(jj = nObsComp[ii]; jj < nDims; jj++) {
	    mvX[ii]->z[jj] = norm_rand();
	  }
	}
	// proposals
	xmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
	     mvX[ii]->sd, nDims);
	//if(isNaN(&propX[ii*nDims], nDims)) {
	//  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
	//  foundNaN = 1;
	//  goto stop;
	//}
	// only calculate acceptance rate if proposal is valid
	if(sdeModel::isValidData(&propX[ii*nDims])) {
	  // acceptance rate
	  // proposal
	  propAccept[ii] = lmvn(&currX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				mvX[ii]->sd, nDims);
	  propAccept[ii] -= lmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // target 1
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
		  dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  propAccept[ii] -= lmvn(&currX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // target 2
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &propX[ii*nDims], dT[ii],
		  sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[ii*nDims], dT[ii],
		  sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // evaluate mh ratio
	  if(exp(propAccept[ii]) >= unif_rand()) {
	    for(jj = 0; jj < nDims; jj++) {
	      currX[ii*nDims + jj] = propX[ii*nDims + jj];
	    }
	    gibbsAccept[ii]++;
	  }
	}
      }
    }
  }
  // special treatment for last datapoint
  if(missInd[nMiss-1] == nComp-1) {
    ii = nComp-1;
    mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
    // partial observations
    if(nObsComp[ii] == 0) {
      for(jj = 0; jj < nDims; jj++) {
	mvX[ii]->z[jj] = norm_rand();
      }
    }
    else {
      zmvn(mvX[ii]->z, &currX[ii*nDims], mvX[ii]->mean,
	   mvX[ii]->sd, nDims, nObsComp[ii]);
      for(jj = nObsComp[ii]; jj < nDims; jj++) {
	mvX[ii]->z[jj] = norm_rand();
      }
    }
    // proposals
    xmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
	 mvX[ii]->sd, nDims);
    //if(isNaN(&propX[ii*nDims], nDims)) {
    //  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
    //  foundNaN = 1;
    //  goto stop;
    //}
    // acceptance is 100% as long as the proposal is valid
    if(sdeModel::isValidData(&propX[ii*nDims])) {
      for(jj = 0; jj < nDims; jj++) {
	currX[ii*nDims + jj] = propX[ii*nDims + jj];
      }
      gibbsAccept[ii]++;
    }
  }
  // special treatment for first datapoint
  if(missInd[0] == 0) {
    ii = 0;
    // initialize
    for(jj = 0; jj < nDims; jj++) {
      propX[jj] = currX[jj];
    }
    // random walk metropolis
    for(jj = 0; jj < nMiss0; jj++) {
      // proposal
      propX[nObsComp[0]+jj] = currX[nObsComp[0]+jj] + jumpSd[nParams+jj] * norm_rand();
      if(sdeModel::isValidData(&propX[ii*nDims])) {
	// acceptance rate.
	// target 1
	propAccept[ii] = prior->logPrior(currTheta, propX);
	propAccept[ii] -= prior->logPrior(currTheta, currX);
	// target 2
	mvEuler(mvX[ii]->mean, mvX[ii]->sd, &propX[ii*nDims], dT[ii],
		sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] += lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
			       mvX[ii]->sd, nDims);
	mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
			       mvX[ii]->sd, nDims);
	// evaluate mh ratio
	if(exp(propAccept[ii]) >= unif_rand()) {
	  currX[ii*nDims + nObsComp[0]+jj] = propX[ii*nDims + nObsComp[0]+jj];
	  paramAccept[nParams + jj]++;
	}
	else {
	  propX[ii*nDims + nObsComp[0]+jj] = currX[ii*nDims + nObsComp[0]+jj];
	}
      }
    }
  }
  //stop:
  //return foundNaN;
}
////////////////////////////////////////////


// componentwise vanilla MH parameter updates
void sdeMCMC::paramVanillaUpdate(double jumpSd[], int paramAccept[]) {
  double acc, currLoglik, propLoglik;
  int ii;
  // initialize
  for(ii = 0; ii < nParams; ii++) {
    propTheta[ii] = currTheta[ii];
  }
  currLoglik = loglik(currTheta, currX);
  // random walk metropolis updates
  for(ii = 0; ii < nParams; ii++) {
    if(!fixedTheta[ii]) {
    // proposal
    propTheta[ii] = currTheta[ii] + jumpSd[ii] * norm_rand();
    // only calculate acceptance if valid
    if(sdeModel::isValidParams(propTheta)) {
      // likelihood
      propLoglik = loglik(propTheta, currX);
      acc = propLoglik - currLoglik;
      // prior
      acc += prior->logPrior(propTheta, currX);
      acc -= prior->logPrior(currTheta, currX);
      // acceptance rate
      if(exp(acc) >= unif_rand()) {
	currTheta[ii] = propTheta[ii];
	currLoglik = propLoglik;
	paramAccept[ii]++;
      }
    }
    // propTheta and currTheta should only differ by one element
    propTheta[ii] = currTheta[ii];
    }
  }
  return;
}
// debugging functions (i.e. check the code)
// also provides loglikelihood


//[[Rcpp::export("sde.model$drift")]]
NumericVector sdeDrift(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sdeModel *sde = new sdeModel[nReps];
  // *** parallelizable for-loop ***
  for(int ii = 0; ii < nReps; ii++) {
    sdeDr(&dr[ii*nDims], &x[ii*nDims], &theta[ii*nParams], &sde[ii]);
  }
  delete [] sde;
  return drOut;
}

//[[Rcpp::export("sde.model$diff")]]
NumericVector sdeDiff(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims*nDims);
  double *df = REAL(dfOut);
  sdeModel *sde = new sdeModel[nReps];
  // *** parallelizable for-loop ***
  for(int ii = 0; ii < nReps; ii++) {
    sdeDf(&df[ii*nDims*nDims], &x[ii*nDims], &theta[ii*nParams], &sde[ii]);
  }
  delete [] sde;
  return dfOut;
}

// SDE log-likelihood evaluation.
//[[Rcpp::export("sde.model$loglik")]]
NumericVector sdeLogLik(NumericVector xIn, NumericVector thetaIn, NumericVector dT,
			int nComp, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  sdeMCMC sdeLL(nComp, REAL(dT));
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*nParams], &x[ii*nDims*nComp]);
  }
  return llOut;
}

// Prior
//[[Rcpp::export("sde.model$logprior")]]
NumericVector sdePrior(NumericVector xIn, NumericVector thetaIn,
		       List priorParams, int priorType, int nRv, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector lpiOut(nReps);
  double *lpi = REAL(lpiOut);
  Prior *prior = NULL;
  if((Prior::Type)priorType == Prior::Flat) {
    prior = new FlatPrior();
  }
  else if((Prior::Type)priorType == Prior::Normal) {
    prior = new NormalPrior(priorParams, nRv);
  }
  else if((Prior::Type)priorType == Prior::GCop) {
    prior = new GCopPrior(priorParams, nRv);
  }
  //else if((Prior::Type)priorType == Prior::Custom) {
  //  prior = new CustomPrior(priorParams, nRv);
  //}
  else {
    throw Rcpp::exception("ERROR: Unrecognized prior type\n");
  }
  //
  for(int ii=0; ii<nReps; ii++) {
    lpi[ii] = prior->logPrior(&theta[nParams*ii], &x[nDims*ii]);
  }
  delete prior;
  return lpiOut;
}
/////////////////////////////////////////////////////

// Euler Simulation of Multivariate SDEs

/////////////////////////////////////////////////////


//[[Rcpp::export("sde.model$sim")]]
List sdeEulerSim(int nDataOut,
		 int N, int reps, int r, double delta,
		 int MAXBAD, NumericVector initData, NumericVector params) {
  RNGScope scope;

  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double sqrtDelta = sqrt(delta);
  int ii, jj, kk;
  int bad = 0;
  // output
  NumericVector dataOut(nDataOut);
  int nBadDraws;

  double *X = new double[nDims];
  double *tmpX = new double[nDims];
  propMV **mvX = new propMV*[reps];
  for(ii=0; ii<reps; ii++) {
    mvX[ii] = new propMV(nDims);
  }
  sdeModel *sde = new sdeModel[reps];
  // initialize
  for(ii = 0; ii < nDims; ii++) {
    X[ii] = 0.0;
    tmpX[ii] = 0.0;
  }

  // *** PARALLELIZABLE FOR-LOOP ***
  for(ii = 0; ii < reps; ii++) {
    for(jj = 0; jj < N*r; jj++) {
      // initialize chains
      if(jj == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = initData[ii*nDims + kk];
	}
      }
      else {
	mvEuler(mvX[ii]->mean, mvX[ii]->sd, X, delta, sqrtDelta,
		&params[ii*nParams], &sde[ii]);
	for(kk = 0; kk < nDims; kk++) {
	  mvX[ii]->z[kk] = norm_rand();
	}
	xmvn(tmpX, mvX[ii]->z, mvX[ii]->mean, mvX[ii]->sd, nDims);
	// validate draw
	while(!sdeModel::isValidData(tmpX) && bad < MAXBAD) {
	  for(kk = 0; kk < nDims; kk++) {
	    mvX[ii]->z[kk] = norm_rand();
	  }
	  xmvn(tmpX, mvX[ii]->z, mvX[ii]->mean, mvX[ii]->sd, nDims);
	  bad++;
	}
	if (bad == MAXBAD) {
	  goto stop;
	}
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = tmpX[kk];
	}
      }
      // store
      if(jj % r == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  dataOut[ii*N*nDims + (jj/r)*nDims + kk] = X[kk];
	}
      }
    }
  }

 stop:
  nBadDraws = bad;

  delete [] X;
  delete [] tmpX;
  for(ii = reps-1; ii>=0; ii--) {
    delete mvX[ii];
  }
  delete [] mvX;
  delete [] sde;

  return List::create(_["dataOut"] = dataOut, _["nBadDraws"] = nBadDraws);
}
////////////////////////////////////////////////////

// Basic MCMC for SDE models

// Parameters and initial missing data are Vanilla Metropolis-within-Gibbs
// Inner missing data are Independence Metropolis-within-Gibbs using the proposal of
// Eraker (2001).

////////////////////////////////////////////////////


//[[Rcpp::export("sde.model$post")]]
List sdeEulerMCMC(NumericVector initParams, NumericVector initData,
		  NumericVector dT, IntegerVector nDimsPerObs, LogicalVector fixedParams,
		  int nSamples, int burn,
		  int nParamsOut, int nDataOut,
		  IntegerVector dataOutRow, IntegerVector dataOutCol,
		  double updateParams, double updateData,
		  int priorType, List priorParams, NumericVector rwJumpSd,
		  int updateLogLik, int nLogLikOut,
		  int updateLastMiss, int nLastMissOut) {
  RNGScope scope;
  int ii, jj, kk;

  // problem dimensions
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  int nComp = initData.length()/nDims;
  int nCompOut = dataOutCol.length();
  int nMiss0 = nDims-nDimsPerObs[0]; // unobserved states in first observation
  int nMissN = nDims-nDimsPerObs[nComp-1]; // unobserved states in last observation

  // output variables
  NumericVector paramsOut(nParamsOut);
  NumericVector dataOut(nDataOut);
  IntegerVector paramAcceptOut(nParams + nMiss0);
  IntegerVector gibbsAcceptOut(nComp);
  NumericVector logLikOut(nLogLikOut);
  NumericVector lastMissOut(nLastMissOut);
  NumericVector lastIter(nParams + nComp*nDims);
  // pointers to acceptance rate counters for internal use
  int *paramAccept = INTEGER(paramAcceptOut);
  int *gibbsAccept = INTEGER(gibbsAcceptOut);

  // initialize prior
  double *jumpSd = REAL(rwJumpSd); // random walk jump sizes
  Prior *prior = NULL;
  if((Prior::Type)priorType == Prior::Flat) {
    prior = new FlatPrior();
  }
  else if((Prior::Type)priorType == Prior::Normal) {
    prior = new NormalPrior(priorParams, nParams + nMiss0);
  }
  else if((Prior::Type)priorType == Prior::GCop) {
    prior = new GCopPrior(priorParams, nParams + nMiss0);
  }
  //else if((Prior::Type)priorType == Prior::Custom) {
  //  prior = new CustomPrior(priorParams, nParams + nMiss0);
  //}
  else {
    throw Rcpp::exception("ERROR: Unrecognized prior type\n");
  }

  // initialize MCMC
  sdeMCMC mcmc(nComp, REAL(dT), REAL(initData), REAL(initParams),
	       INTEGER(nDimsPerObs), (bool*)LOGICAL(fixedParams), prior);

  // main MCMC loop
  jj = 0;
  for(int smp = -burn; smp < nSamples; smp++) {
    // user interrupt
    if(smp % (int) 1e4) {
      Rcpp::checkUserInterrupt();
    }
    // missing data update
    if(updateComponent(updateData, smp)) {
      mcmc.missGibbsUpdate(jumpSd, gibbsAccept, paramAccept);
    }
    // parameter update
    if(updateComponent(updateParams, smp)) {
      mcmc.paramVanillaUpdate(jumpSd, paramAccept);
    }
    // log-likelihood
    // TODO: keep track of this interally after every MCMC step
    if(smp >= 0) {
      if(updateLogLik) logLikOut[smp] = mcmc.loglik(mcmc.currTheta, mcmc.currX);
    }
    // storage
    if(smp == dataOutRow[jj]) {
      if(updateData > 0.0) {
	for(ii=0; ii<nCompOut; ii++) {
	  for(kk=0; kk<nDims; kk++) {
	    dataOut[jj*nDims*nCompOut+ii*nDims+kk] = mcmc.currX[dataOutCol[ii]*nDims+kk];
	  }
	}
      }
      jj++;
    }
    if((updateParams > 0.0) && (smp >= 0)) {
      for(ii=0; ii<nParams; ii++) paramsOut[smp*nParams + ii] = mcmc.currTheta[ii];
    }
    if(updateLastMiss && smp >= 0) {
      for(ii=0; ii<nMissN; ii++) {
	lastMissOut[smp*nMissN+ii] = mcmc.currX[(nComp-1)*nDims+nDimsPerObs[nComp-1]+ii];
      }
    }
    // NaN breakpoint
    //if(foundNaN) {
    //  Rprintf("smp = %i\n", smp);
    //  goto stop;
    //}
  }

  // store last iteration, to resume MCMC later if needed
  for(ii=0; ii<nParams; ii++) {
    lastIter[ii] = mcmc.currTheta[ii];
  }
  for(ii=0; ii<nDims*nComp; ii++) {
    lastIter[nParams+ii] = mcmc.currX[ii];
  }

  //stop:
  // delete dynamic variables
  delete prior;
  //delete [] sqrtDT;
  //delete [] B;
  //delete [] sqrtB;
  //delete [] currData;
  //delete [] propData;
  //delete [] propAccept;
  //delete [] currParams;
  //delete [] propParams;
  //delete [] missInd;
  //delete [] pastDT;
  //delete [] pastSqrtDT;
  //delete [] logMultiAcc;
  //delete prior;
  //delete pastPrior;

  return List::create(_["paramsOut"] = paramsOut, _["dataOut"] = dataOut,
		      _["paramAccept"] = paramAcceptOut,
		      _["gibbsAccept"] = gibbsAcceptOut,
		      _["logLikOut"] = logLikOut,
		      _["lastMissOut"] = lastMissOut, _["lastIter"] = lastIter);
}



