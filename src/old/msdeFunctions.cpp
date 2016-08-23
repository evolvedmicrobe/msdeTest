////////////////////////////////////////

// some functions and classes for msde inference

// design notes:
// Don't use Eigen at top level, as its generally slower for small matrices
// sdeModel class will have two different examples: one without eigen and one with

/////////////////////////////////////////

// sdeModel class: constains functions unique to a given sde

class sdeModel {
  // put private storage variables here
  // end private storage
public:
  static int nParams = 5;
  static int nDims = 2;
  friend void sdeDr(double *dr, double *x, double t, double *theta);
  friend void sdeDf(double *df, double *x, double sqrtT, double *theta);
  static bool isValidData(double *x);
  static bool isValidParams(double *theta);
  sdeModel();
  ~sdeModel();
};


///////////////////////////////////////////////////////////

// examples for heston's model

// constructor
sdeModel::sdeModel() {
}

// destructor
sdeModel::~sdeModel() {} // do nothing

// todo: drop t and sqrtT in these functions
void sdeDr(double *dr, double *x, double t, double *theta) {
  dr[0] = (theta[0] - .125 * x[1]*x[1]) * t; // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x(1)) * t; // z
  return;
}

void sdeDf(double *df, double *x, double sqrtT, double *theta) {
  df[0] = .5 * x[1] * sqrtT;
  df[2] = theta[3] * sqrtT;
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

static bool sdeModel::isValidData(double *x) {
  return(x[1] > 0.0);
}

static bool sdeModel::isValidParams(double *theta) {
  int isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}

////////////////////////////////////////////////////////////////

// a few utility functions for multivariate normals

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
double lgcop(double *x, double *qnorm, int *nBreaks, double *range,
	     double *dx, double *pdf, double *lpdf, double *cdf,
	     double *mean, double *RhoCholSd, int n) {
  int ii, jj, colI, densElt, start;
  double lp = 0.0;
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
  tmpSum = 0.0;
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


//////////////////////////////////////////////////////////////////////////

// class for proposal mean and variances (or rather sd's)
// also create a dummy space of size mean.

class propMV {
public:
  int nDims;
  double *mean, *sd, *z;
  propMV(int);
  ~propMV();
};

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
 public:
  // the type of prior
  enum Type {
    Flat = 1,
    Normal = 2,
    GCop = 3,
    Custom = 4
  };

  //pure virtual function must be defined in each derived class
  virtual double logPrior(double *x, double *theta)=0;
};

// flat prior
class FlatPrior : public Prior
{
 public:
  FlatPrior();
  ~FlatPrior();
  double logPrior(double *x, double *theta);
};

// normal prior
class NormalPrior : public Prior
{
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
  int nDims, nTheta;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  NormalPrior(double **priorParams, int nD, int nTh, int nRv);
  ~NormalPrior();
  double logPrior(double *x, double *theta);
};

NormalPrior::NormalPrior(double **priorParams, int nD, int nTh, int nRv) {
  nDims = nD;
  nTheta = nTh;
  nTotalRV = nDims + nTheta;
  nActiveRV = nRV;
  nSkippedRV = nTotalRV - nActiveRV;
  mean = priorParams[0];
  sd = priorParams[1];
  tmpX = new double[nTotalRV];
  tmpZ = new double[nTotalRV];
}
NormalPrior::~NormalPrior() {
  delete [] tmpX;
  delete [] tmpZ;
}

double NormalPrior::logPrior(double *x, double *theta) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nTheta; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lmvn(&(tmpX+nSkippedRV), &(tmpZ+nSkippedRV), mean, cholSd, nActiveRV);
  return(lp);
}


// gaussian copula prior
class GCopPrior : public Prior
{
  // gaussian copula parameters
  int *nBreaks;
  double *range, *dx, *pdf, *logPdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  int nActiveRV, nSkippedRV, nTotalRV;
  int nDims, nTheta;
 public:
  GCopPrior(double **priorParamsD, int *priorParamsI, int nD, int nTh, int nRv);
  ~GCopPrior();
  // parameter + unobserved first states gaussian copula prior.
  double logPrior(double x[], double theta[]);
};


GCopPrior::GCopPrior(double **priorParamsD, int *priorParamsI,
		     int nD, int nTh, int nRv) {
  nDims = nD;
  nTheta = nTh;
  nTotalRV = nDims + nTheta;
  nActiveRV = nRV;
  nSkippedRV = nTotalRV - nActiveRV;
  nBreaks = priorParamsI;
  range = priorParamsD[0];
  dx = priorParamsD[0];
  pdf = priorParamsD[0];
  logPdf = priorParamsD[0];
  cdf = priorParamsD[0];
  mean = priorParamsD[0];
  sd = priorParamsD[0];
  RhoCholSd = priorParamsD[0];
  tmpX = new double[nTotalRV];
  tmpQ = new double[nTotalRV];
}
GCopPrior::~GCopPrior() {
  delete [] tmpX;
  delete [] tmpQ;
}

GCopPrior::logPrior(double *x, double *theta) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nTheta; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lgcop(&(tmpX+nSkippedRV), &(tmpZ+nSkippedRV),
	     nBreaks, range, dx, pdf, lpdf, cdf,
	     mean, RhoCholSd, nActiveRV);
  return(lp);
}

// custom prior
class CustomPrior : public Prior {
  // private custom prior parameters
  // end custom prior parameters
 public:
  CustomPrior(double **priorParamsD, int **priorParamsI, int nDims, int nTheta);
  ~CustomPrior();
  double logPrior(double *x, double *theta);
};

////////////////////////////////////////////////////////////////////////////

// ok now an MCMC object.
// contains all the storage etc. required to run the update methods.

// in terms of storage in {curr/prop}Full, let's put theta before x.


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
  sdeMCMC(int N, int *xIndex, int *thetaIndex);
  ~sdeMCMC();
};


///////////////////////////////////////////////////////////////

// a few utility functions for SDEs,
// namely mean and cholesky factor for Eraker and Euler normals

// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
void mvEraker(double *mean, double *sd,
	      double *x0, double *x2, double b, double b2, double *theta,
	      sdeModel *sde) {
  for(int jj = 0; jj < sdeModel::nDims; jj++) {
    mean[jj] = x0[jj] * b + x2[jj] * (1-b);
  }
  sde->sdeDf(sd, x0, b2, theta);
  return;
}

// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double *mean, double *sd,
	     double *x0, double t, double t2, double theta*,
	     sdeModel *sde) {
  sde->sdeDr(mean, x0, t, theta);
  for(int jj = 0; jj < sdeModel::nDims; jj++) {
    mean[jj] += x0[jj];
  }
  sde->sdeDf(sd, x0, t2, theta);
  return;
}

////////////////////////////////////////////////////////////////////////////////

// loglikelihood function

double sdeMCMC:loglik(double *theta, double *x) {
  double ll = 0;
  // *** PARALLELIZABLE FOR-LOOP ***
  for(int ii = 0; ii < nComp-1; ii++) {
    mvEuler(mvX[ii]->mean, mvX[ii]->sd, &x[ii*nDims],
	    dT[ii], sqrtDT[ii], theta, &sde[ii]);
    ll += lmvn(&x[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean, mvX[ii]->sd, nDims);
  }
  return(ll);
}

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

// componentwise vanilla MH parameter updates
void sdeMCMC::paramVanillaUpdate(double jumpSd[], int paramsAccept[]) {
  double acc, currLoglik, propLoglik;
  int ii;
  // initialize
  for(ii = 0; ii < nParams; ii++) {
    propTheta[ii] = currTheta[ii];
  }
  currLoglik = loglik(currTheta, currX);
  // random walk metropolis updates
  for(ii = 0; ii < nParams; ii++) {
    if(!fixedParams[ii]) {
    // proposal
    propTheta[ii] = currTheta[ii] + jumpSD[ii] * norm_rand();
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
