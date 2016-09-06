/////////////////////////////////////////////////////

// Euler Simulation of Multivariate SDEs

/////////////////////////////////////////////////////

#include "mvnUtils.h"
#include "sdeModel.h"
#include "sdeUtils.h"

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
	  mvX[ii]->z[kk] = NORM_RAND();
	}
	xmvn(tmpX, mvX[ii]->z, mvX[ii]->mean, mvX[ii]->sd, nDims);
	// validate draw
	while(!sdeModel::isValidData(tmpX) && bad < MAXBAD) {
	  for(kk = 0; kk < nDims; kk++) {
	    mvX[ii]->z[kk] = NORM_RAND();
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
