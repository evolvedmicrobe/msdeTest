#'@name sde.init
#'@title Posterior Sampler Initialization
#'@description  initialize posterior sampler for a given resolution level k, \cr
#'              or m missing data points per interval, i.e. m = 2^k-1. \cr
#'              currently only supports constant interobservation times.
#'@param data An array of data
#'@param dt A vector of Delta-t. Must have \code{length(dt) == 1 || length(dt) == nrow(data)-1}
#'@param k An integer representing a resolution level (only used if \code{m} is missing)
#'@param m An integerrepresenting the number of missing data points per interval
#'@param par.index An integer
#'@param params A vector of parameters that, if provided, will be included in the output of this function.
#'@param debug a boolean (FALSE by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return a list containing: data, dt, par.index, params (only if \code{params} was passed to sde.init)
#'@examples
#'# Create the model
#'hest.model <- sde.make.model(list = hestList, model.name = "hest")
#'
#'theta <- c(alpha = .1, gamma = 5, beta = .8, sigma = .6, rho = -.7)
#'Y0 <- c(X = log(100), Z = .1)
#'
#'# simulate data
#'N <- 10
#'burn <- 10
#'dT <- 1/252
#'
#'hsim <- sde.sim(model = hest.model, init.data = Y0, params = theta, dt = dT, dt.sim = dT/100,
#'                N = N, burn = burn, nreps = 1)
#'
#'k <- 1
#'par.index <- 1
#'
#'init <- sde.init(data = hsim$data, dt = dT, k = k, par.index = par.index, params = theta)
#'@export
sde.init <- function(data, dt, k, m, par.index, params, debug = FALSE) {
  nobs <- nrow(data)
  ndims <- ncol(data)
  if(missing(m)) m <- 2^k-1
  ncomp <- (nobs-1)*(m+1)+1
  init.data <- matrix(NA, ncomp, ndims)
  colnames(init.data) <- colnames(data)
  # interpolation to create missing data
  if(debug) browser()
  dt1 <- length(dt) == 1
  if(dt1) dt <- rep(dt, nobs-1)
  if(length(dt) != nobs-1) stop("Incorrect specification of dt.")
  dtnew <- rep(dt/(m+1), each = m+1)
  told <- cumsum(c(0, dt))
  tnew <- cumsum(c(0, dtnew))
  for(ii in 1:ndims) {
    init.data[,ii] <- approx(x = told, y = data[,ii],
                             xout = tnew)$y
  }
  if(dt1) dtnew <- dtnew[1]
  if(missing(par.index)) par.index <- ndims
  par.ind <- rep(0, ncomp)
  par.ind[seq(1, ncomp, len = nobs)] <- par.index
  ans <- list(data = init.data, dt = dtnew, par.index = par.ind)
  if(!missing(params)) ans$params <- params
  ans
}
