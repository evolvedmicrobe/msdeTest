/////////////////////////////////////////////////////////////////////////////

#ifndef Prior_h
#define Prior_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"
#include "sdeModel.h"

// prior classes
//
// there are four types of priors: flat, normal, gaussian copula, and custom
// each prior calls two arguments: x and theta;
// tuning parameters are set by the constructor.
// in addition to tuning parameters, the normal and gcop constructors
// take the following dimension arguments: nD = nDims, nTh = nTheta,
// and nRv <= nDims + nTh, the number of "active variables".
// the inputs are ligned up as tmpX = (x, theta), and the active variables are the _last_
// nRV of these.

//base prior class
class Prior {
 protected:
  static const int nDims = sdeModel::nDims;
  static const int nParams = sdeModel::nParams;
 public:
  // the type of prior
  enum Type {
    Flat = 1,
    Normal = 2,
    GCop = 3,
    Custom = 4
  };

  //pure virtual function must be defined in each derived class
  virtual double logPrior(double *theta, double *x)=0;
};

// flat prior
class FlatPrior : public Prior {
 public:
  FlatPrior();
  ~FlatPrior();
  double logPrior(double *theta, double *x);
};

// normal prior
class NormalPrior : public Prior {
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  NormalPrior(List priorParams, int nRv);
  ~NormalPrior();
  double logPrior(double *theta, double *x);
};

// gaussian copula prior
class GCopPrior : public Prior {
  // gaussian copula parameters
  int *nBreaks;
  double *range, *dx, *pdf, *lpdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  double *tmpX, *tmpQ;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  GCopPrior(List priorParams, int nRv);
  ~GCopPrior();
  // parameter + unobserved first states gaussian copula prior.
  double logPrior(double *theta, double *x);
};

#endif
