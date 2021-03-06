////////////////////////////////////////////

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"
#include "sdeMCMC.h"

// componentwise vanilla MH parameter updates
void sdeMCMC::paramVanillaUpdate(double *jumpSd, int *paramAccept) {
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
