////////////////////////////////////////////////////////

// utility functions for sdes

#ifndef sdeUtils_h
#define sdeUtils_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeModel.h"

// euler approximation mean and standard deviation
// NOTE: sde = upper triangular cholesky factor
void mvEuler(double *mean, double *sd,
	     double *x0, double t, double t2, double *theta,
	     sdeModel *sde);

// loglikelihood function (actually log-density, minus factor of pi)
double loglik(double *theta, double *x, double *dT, double *sqrtDT,
	      propMV **mvX, sdeModel *sde, int nComp);

#endif
