#' SDE Drift Function
#'
#' @export
<<<<<<< HEAD
sde.drift2 <- function(model, x, theta) {
=======
sde.drift <- function(model, x, theta) {
>>>>>>> rename
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(!is.matrix(x)) x <- matrix(x, ncol = 1) else x <- t(x)
  if(!is.matrix(theta)) theta <- matrix(theta, ncol = 1) else theta <- t(theta)
  nreps <- max(ncol(x), ncol(theta))
  if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(ncol(x) != ncol(theta)) {
    stop("x and theta must have the same number of rows.")
  }
  # compute
  ans <- model$drift(x = as.double(x), theta = as.double(theta),
                     nReps = as.integer(nreps))
  dr <- matrix(ans, nreps, ndims, byrow = TRUE)
  dr
}
