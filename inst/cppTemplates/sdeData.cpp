#include "sdeData.h"

// functions for class sdeData

///////////////////////////////////////////////////////////////

// constructor
sdeData::sdeData(int N, double *dt) {
  int ii;
  // problem dimensions
  nComp = N;
  nDims = sdeModel::nDims;
  nParams = sdeModel::nParams;
  //times
  dT = dt;
  sqrtDT = new double[nComp];
  for(ii=0; ii<nComp-1; ii++) {
    sqrtDT[ii] = sqrt(dT[ii]);
  }
  mvX = new propMV*[nComp];
  for(ii=0; ii<nComp; ii++) {
    mvX[ii] = new propMV(nDims);
  }
  // sde
  sde = new sdeModel[nComp];
}

// destructor
sdeData::~sdeData() {
  delete [] sqrtDT;
  for(int ii=nComp-1; ii>=0; ii--) {
    delete mvX[ii];
  }
  delete [] mvX;
  delete [] sde;
}

///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
double sdeData::loglik(double *theta, double *x) {
  double ll = 0;
  // *** PARALLELIZABLE FOR-LOOP ***
  for(int ii = 0; ii < nComp-1; ii++) {
    mvEuler(mvX[ii]->mean, mvX[ii]->sd, &x[ii*nDims], theta, ii);
    ll += lmvn(&x[(ii+1)*nDims], mvX[ii]->z,
	       mvX[ii]->mean, mvX[ii]->sd, nDims);
  }
  return(ll);
}
