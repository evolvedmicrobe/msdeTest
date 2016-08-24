//////////////////////////////////////////////////////

// utility functions for sdes

#include "sdeUtils.h"

// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double *mean, double *sd,
	     double *x0, double t, double t2, double *theta,
	     sdeModel *sde) {
  sde->sdeDr(mean, x0, theta);
  v_mult(mean, t, sdeModel::nDims);
  for(int jj = 0; jj < sdeModel::nDims; jj++) {
    mean[jj] += x0[jj];
  }
  sde->sdeDf(sd, x0, theta);
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
