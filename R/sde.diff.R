#' SDE Diffusion Function
#'
#' @export
sde.diff2 <- function(model, x, theta) {
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
  ans <- model$diff(x = as.double(x), theta = as.double(theta),
                    nReps = as.integer(nreps))
  df <- matrix(ans, nreps, ndims^2, byrow = TRUE)
  # put zeros into the off-triangular elements
  df[,lower.tri(diag(ndims))] <- 0
  df
}
