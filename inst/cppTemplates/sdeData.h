/////////////////////////////////////////

#ifndef sdeData_h
#define sdeData_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"
#include "sdeModel.h"

// Class to store data and loglikelihood

class sdeData {
  // mean and variance of a forward euler step from a given observation
//  void mvEuler(double *mean, double *sd, double *x0, double *theta,
//	       int iObs0);
 public:
  int nComp, nDims, nParams;
  double *dT, *sqrtDT;
  propMV **mvX;
  sdeModel *sde;
  sdeData(int N, double *dt);
  ~sdeData();
  // log-density
  double loglik(double *theta, double *x);
private:
  // euler approximation mean and standard deviation
  // NOTE: sde = upper triangular cholesky factor
  inline void mvEuler(double *mean, double *sd,
                               double *x0, double *theta, int iObs0) {
    sde[iObs0].sdeDr(mean, x0, theta);
    v_mult(mean, dT[iObs0], nDims);
    for(int jj = 0; jj < nDims; jj++) {
      mean[jj] += x0[jj];
    }
    sde[iObs0].sdeDf(sd, x0, theta);
    U_mult(sd, sqrtDT[iObs0], nDims);
    return;
  }
};

#endif
