////////////////////////////////////////////////////////////////////////////

#ifndef sdeMCMC_h
#define sdeMCMC_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"
#include "sdeModel.h"
#include "sdeUtils.h"
#include "Prior.h"

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


#endif

