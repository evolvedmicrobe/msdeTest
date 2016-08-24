#' Basic Forward Simulation
#'
#' @export
sde.sim <- function(model,
                    init.data, params, dt, dt.sim,
                    N, burn = 0, nreps = 1,
                    max.bad.draws = 5e3, verbose = TRUE,
                    debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(!is.matrix(init.data)) init.data <- t(as.matrix(init.data))
  if(!is.matrix(params)) params <- t(as.matrix(params))
  if(ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(nparams != ncol(params))
    stop("params does not have the right length/number of columns.")
  # data
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != data.names))
      stop("Incorrect data.names.")
  }
  if(nrow(init.data) == 1) {
    init.data <- matrix(init.data, nrow = nreps, ncol = ndims, byrow = TRUE)
  }
  if(nrow(init.data) != nreps) stop("init.data does not have the right number of rows.")
  colnames(init.data) <- data.names
  # params
  init.params <- params
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != param.names))
      stop("Incorrect param.names.")
  }
  if(nrow(params) == 1) {
    params <- matrix(params, nrow = nreps, ncol = nparams, byrow = TRUE)
  }
  if(nrow(params) != nreps) stop("params does not have the right number of rows.")
  colnames(init.params) <- param.names
  # time
  if(dt.sim <= dt) {
    r <- ceiling(dt/dt.sim)
    t <- dt/r
  } else {
    r <- 1
    t <- dt
  }
  if(burn < 1) burn <- N*burn
  burn <- floor(burn)
  if(verbose) {
    message("Normal draws required: ", round((N+burn)*r*nreps, 2))
    if(verbose > 1) {
      ans <- readline("Press Q to quit, any other key to proceed: ")
      if(substr(ans, 1, 1) == "Q") {
        message("Ended by user.")
        return()
      }
    }
    message("Running simulation...")
  }
  if(debug) browser()
  tm <- chrono()
  # forward simulation
  ans <- model$sim(nDataOut = as.integer((N+burn)*ndims*nreps),
                   N = as.integer(N+burn),
                   reps = as.integer(nreps),
                   r = as.integer(r),
                   delta = as.double(t),
                   MAXBAD = as.integer(max.bad.draws),
                   initData = as.double(t(init.data)),
                   params = as.double(t(params)))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("dataOut", "nBadDraws")
  if(verbose) message("Bad Draws = ", ans$nBadDraws)
  data <- aperm(array(ans$dataOut, dim = c(ndims, N+burn, nreps)),
                perm = 3:1)
  dimnames(data) <- list(NULL, NULL, data.names)
  out <- list(data = data[,burn+1:N,], params = init.params[,], dt = dt,
              dt.sim = t, nbad = ans$nBadDraws)
  out
}
