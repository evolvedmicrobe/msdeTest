#' SDE Loglikelihood Function
#'
#' @export
sde.loglik2 <- function(model, x, dt, theta, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(debug) browser()
  if(is.matrix(x)) {
    ncomp <- nrow(x)
    x <- array(t(x), dim = c(ndims,ncomp,1))
  } else {
    ncomp <- dim(x)[2]
    x <- aperm(x, perm = 3:1)
  }
  if(!is.matrix(theta)) theta <- matrix(theta, ncol = 1) else theta <- t(theta)
  nreps <- max(dim(x)[3], ncol(theta))
  if(dim(x)[3] == 1) x <- array(x, dim = c(ndims,ncomp,nreps))
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(!dim(x)[3] == ncol(theta)) {
    stop("x and theta must have the same number of samples.")
  }
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  # compute
  ans <- model$loglik(xIn = as.double(x),
                      thetaIn = as.double(theta), dT = as.double(dt),
                      nComp = as.integer(ncomp),
                      nReps = as.integer(nreps))
  ans
}
