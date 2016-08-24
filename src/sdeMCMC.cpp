////////////////////////////////////////////

#include "sdeMCMC.h"

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
  sde->sdeDf(sd, x0, theta);
  U_mult(sd, b2, sdeModel::nDims);
  return;
}


////////////////////////////////////////////////////////////////////////////////
