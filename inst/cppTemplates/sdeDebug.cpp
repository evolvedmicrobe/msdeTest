// debugging functions (i.e. check the code)
// also provides loglikelihood

#include "sdeModel.h"
#include "Prior.h"
#include "sdeData.h"
#include "sdeMCMC.h"

//[[Rcpp::export("sde.model$drift")]]
NumericVector sdeDrift(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sdeModel *sde = new sdeModel[nReps];
  // *** parallelizable for-loop ***
  for(int ii = 0; ii < nReps; ii++) {
    sde[ii].sdeDr(&dr[ii*nDims], &x[ii*nDims], &theta[ii*nParams]);
  }
  delete [] sde;
  return drOut;
}

//[[Rcpp::export("sde.model$diff")]]
NumericVector sdeDiff(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims*nDims);
  double *df = REAL(dfOut);
  sdeModel *sde = new sdeModel[nReps];
  // *** parallelizable for-loop ***
  for(int ii = 0; ii < nReps; ii++) {
    sde[ii].sdeDf(&df[ii*nDims*nDims], &x[ii*nDims], &theta[ii*nParams]);
  }
  delete [] sde;
  return dfOut;
}

// SDE log-likelihood evaluation.
//[[Rcpp::export("sde.model$loglik")]]
NumericVector sdeLogLik(NumericVector xIn, NumericVector thetaIn,
			NumericVector dT,
			int nComp, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  sdeData sdeLL(nComp, REAL(dT));
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*nParams], &x[ii*nDims*nComp]);
  }
  return llOut;
}

// Prior
//[[Rcpp::export("sde.model$logprior")]]
NumericVector sdePrior(NumericVector xIn, NumericVector thetaIn,
		       List priorParams, int priorType, int nRv, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector lpiOut(nReps);
  double *lpi = REAL(lpiOut);
  Prior *prior = NULL;
  if((Prior::Type)priorType == Prior::Flat) {
    prior = new FlatPrior();
  }
  else if((Prior::Type)priorType == Prior::Normal) {
    prior = new NormalPrior(priorParams, nRv);
  }
  else if((Prior::Type)priorType == Prior::GCop) {
    prior = new GCopPrior(priorParams, nRv);
  }
  //else if((Prior::Type)priorType == Prior::Custom) {
  //  prior = new CustomPrior(priorParams, nRv);
  //}
  else {
    throw Rcpp::exception("ERROR: Unrecognized prior type\n");
  }
  //
  for(int ii=0; ii<nReps; ii++) {
    lpi[ii] = prior->logPrior(&theta[nParams*ii], &x[nDims*ii]);
  }
  delete prior;
  return lpiOut;
}
