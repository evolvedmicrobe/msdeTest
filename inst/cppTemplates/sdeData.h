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
  void mvEuler(double *mean, double *sd, double *x0, double *theta,
	       int iObs0);
 public:
  int nComp, nDims, nParams;
  double *dT, *sqrtDT;
  propMV **mvX;
  sdeModel *sde;
  sdeData(int N, double *dt);
  ~sdeData();
  // log-density
  double loglik(double *theta, double *x);
};

#endif
