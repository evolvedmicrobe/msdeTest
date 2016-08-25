#--- test coefficients of SDE model ---------------------------------------------

require(msdeTest)
require(msde)

hest.model <- sde.make.model(model = hestList)

# heston drift and diffusion functions

# to save time encode these with minimal efficiency
hest.dr <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .125 * Z^2, theta[,3]/Z - .5*theta[,2]*Z)
}
hest.df <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cv <- .5*theta[,5]*theta[,4]*Z
  ans <- cbind(.25 * Z^2, cv, cv, theta[,4]^2)
  t(apply(ans, 1, function(x) chol(matrix(x,2,2))))
}

nReps <- 1e3
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

dr1 <- hest.dr(X = X0[,1], Z = X0[,2], theta = Theta)
df1 <- hest.df(X = X0[,1], Z = X0[,2], theta = Theta)
dr2 <- hest.drift(model = hest.model, x = X0, theta = Theta)
df2 <- hest.diff(model = hest.model, x = X0, theta = Theta)

range(dr1-dr2)
range(df1-df2)
