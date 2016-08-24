/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"

// sdeModel class: constains functions unique to a given sde

class sdeModel {
  // put private storage variables here
  // end private storage
public:
  static const int nParams = 5;
  static const int nDims = 2;
  void sdeDr(double *dr, double *x, double *theta);
  void sdeDf(double *df, double *x, double *theta);
  static bool isValidData(double *x);
  static bool isValidParams(double *theta);
  sdeModel();
  ~sdeModel();
};

#endif
