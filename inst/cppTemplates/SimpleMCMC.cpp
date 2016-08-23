////////////////////////////////////////////////////

// Basic MCMC for SDE models

// Parameters and initial missing data are Vanilla Metropolis-within-Gibbs
// Inner missing data are Independence Metropolis-within-Gibbs using the proposal of
// Eraker (2001).

////////////////////////////////////////////////////

#include "mcmcUtils.h"
#include "sdeModel.h"
#include "sdeMCMC.h"

//[[Rcpp::export("sde.model$post")]]
List sdeEulerMCMC(NumericVector initParams, NumericVector initData,
		  NumericVector dT, IntegerVector nDimsPerObs, LogicalVector fixedParams,
		  int nSamples, int burn,
		  int nParamsOut, int nDataOut,
		  IntegerVector dataOutRow, IntegerVector dataOutCol,
		  double updateParams, double updateData,
		  int priorType, List priorParams, NumericVector rwJumpSd,
		  int updateLogLik, int nLogLikOut,
		  int updateLastMiss, int nLastMissOut) {
  RNGScope scope;
  int ii, jj, kk;

  // problem dimensions
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  int nComp = initData.length()/nDims;
  int nCompOut = dataOutCol.length();
  int nMiss0 = nDims-nDimsPerObs[0]; // unobserved states in first observation
  int nMissN = nDims-nDimsPerObs[nComp-1]; // unobserved states in last observation

  // output variables
  NumericVector paramsOut(nParamsOut);
  NumericVector dataOut(nDataOut);
  IntegerVector paramAcceptOut(nParams + nMiss0);
  IntegerVector gibbsAcceptOut(nComp);
  NumericVector logLikOut(nLogLikOut);
  NumericVector lastMissOut(nLastMissOut);
  NumericVector lastIter(nParams + nComp*nDims);
  // pointers to acceptance rate counters for internal use
  int *paramAccept = INTEGER(paramAcceptOut);
  int *gibbsAccept = INTEGER(gibbsAcceptOut);

  // initialize prior
  double *jumpSd = REAL(rwJumpSd); // random walk jump sizes
  Prior *prior = NULL;
  if((Prior::Type)priorType == Prior::Flat) {
    prior = new FlatPrior();
  }
  else if((Prior::Type)priorType == Prior::Normal) {
    prior = new NormalPrior(priorParams, nParams + nMiss0);
  }
  else if((Prior::Type)priorType == Prior::GCop) {
    prior = new GCopPrior(priorParams, nParams + nMiss0);
  }
  //else if((Prior::Type)priorType == Prior::Custom) {
  //  prior = new CustomPrior(priorParams, nParams + nMiss0);
  //}
  else {
    throw Rcpp::exception("ERROR: Unrecognized prior type\n");
  }

  // initialize MCMC
  sdeMCMC mcmc(nComp, REAL(dT), REAL(initData), REAL(initParams),
	       INTEGER(nDimsPerObs), (bool*)LOGICAL(fixedParams), prior);

  // main MCMC loop
  jj = 0;
  for(int smp = -burn; smp < nSamples; smp++) {
    // user interrupt
    if(smp % (int) 1e4) {
      Rcpp::checkUserInterrupt();
    }
    // missing data update
    if(updateComponent(updateData, smp)) {
      mcmc.missGibbsUpdate(jumpSd, gibbsAccept, paramAccept);
    }
    // parameter update
    if(updateComponent(updateParams, smp)) {
      mcmc.paramVanillaUpdate(jumpSd, paramAccept);
    }
    // log-likelihood
    // TODO: keep track of this interally after every MCMC step
    if(smp >= 0) {
      if(updateLogLik) logLikOut[smp] = mcmc.loglik(mcmc.currTheta, mcmc.currX);
    }
    // storage
    if(smp == dataOutRow[jj]) {
      if(updateData > 0.0) {
	for(ii=0; ii<nCompOut; ii++) {
	  for(kk=0; kk<nDims; kk++) {
	    dataOut[jj*nDims*nCompOut+ii*nDims+kk] = mcmc.currX[dataOutCol[ii]*nDims+kk];
	  }
	}
      }
      jj++;
    }
    if((updateParams > 0.0) && (smp >= 0)) {
      for(ii=0; ii<nParams; ii++) paramsOut[smp*nParams + ii] = mcmc.currTheta[ii];
    }
    if(updateLastMiss && smp >= 0) {
      for(ii=0; ii<nMissN; ii++) {
	lastMissOut[smp*nMissN+ii] = mcmc.currX[(nComp-1)*nDims+nDimsPerObs[nComp-1]+ii];
      }
    }
    // NaN breakpoint
    //if(foundNaN) {
    //  Rprintf("smp = %i\n", smp);
    //  goto stop;
    //}
  }

  // store last iteration, to resume MCMC later if needed
  for(ii=0; ii<nParams; ii++) {
    lastIter[ii] = mcmc.currTheta[ii];
  }
  for(ii=0; ii<nDims*nComp; ii++) {
    lastIter[nParams+ii] = mcmc.currX[ii];
  }

  //stop:
  // delete dynamic variables
  delete prior;
  //delete [] sqrtDT;
  //delete [] B;
  //delete [] sqrtB;
  //delete [] currData;
  //delete [] propData;
  //delete [] propAccept;
  //delete [] currParams;
  //delete [] propParams;
  //delete [] missInd;
  //delete [] pastDT;
  //delete [] pastSqrtDT;
  //delete [] logMultiAcc;
  //delete prior;
  //delete pastPrior;

  return List::create(_["paramsOut"] = paramsOut, _["dataOut"] = dataOut,
		      _["paramAccept"] = paramAcceptOut,
		      _["gibbsAccept"] = gibbsAcceptOut,
		      _["logLikOut"] = logLikOut,
		      _["lastMissOut"] = lastMissOut, _["lastIter"] = lastIter);
}

