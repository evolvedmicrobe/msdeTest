///////////////////////////////////////////////////////////

#include "sdeModel.h"

// custom functions for heston's model

// constructor
sdeModel::sdeModel() {
}

// destructor
sdeModel::~sdeModel() {} // do nothing

// todo: drop t and sqrtT in these functions
void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = (theta[0] - .125 * x[1]*x[1]); // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x[1]); // z
  return;
}

void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = .5 * x[1];
  df[2] = theta[3];
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

bool sdeModel::isValidData(double *x) {
  return(x[1] > 0.0);
}

bool sdeModel::isValidParams(double *theta) {
  int isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}
