#' Basic MCMC sampler for Multivariate SDEs
#'
#' @export
hest.post <- function(model, init,
                      init.data, init.params, par.index, dt,
                      nsamples, burn,
                      data.out.ind, prior, fixed.params,
                      rw.jump.sd = NULL,
                      update.data = TRUE, update.params = TRUE,
                      loglik.out = FALSE, last.miss.out = FALSE,
                      verbose = TRUE, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initial parameters
  if(!missing(init)) {
    if(!is.null(init$data)) init.data <- init$data
    if(!is.null(init$dt)) dt <- init$dt
    if(!is.null(init$par.index)) par.index <- init$par.index
    if(!is.null(init$params)) init.params <- init$params
  }
  # parse inputs
  ncomp <- nrow(init.data)
  nmiss0 <- ndims - par.index[1]
  nparams2 <- nparams+nmiss0
  if(ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != data.names))
      stop("Incorrect data.names.")
  }
  if(nparams != length(init.params)) {
    stop("init.params does not have the right length.")
  }
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != param.names))
      stop("Incorrect param.names.")
  }
  # time
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # storage
  # data can be subset both by time and by iteration
  if(missing(data.out.ind)) data.out.ind <- 2e3
  if(!is.list(data.out.ind)) {
    data.out.row <- data.out.ind
    data.out.col <- 1:ncomp
  } else {
    data.out.row <- data.out.ind$row
    data.out.col <- data.out.ind$col
  }
  if(length(data.out.row) == 1)
    data.out.row <- unique(floor(seq(1, nsamples, len = data.out.row)))
  if(is.logical(data.out.row)) data.out.row <- which(data.out.row)
  nsamples.out <- length(data.out.row)
  if(is.logical(data.out.col)) data.out.col <- which(data.out.col)
  ncomp.out <- length(data.out.col)
  data.out.ind <- list(row = data.out.row, col = data.out.col)
  if(all(par.index == ndims)) update.data <- FALSE
  if(update.data) {
    ndata.out <- nsamples.out*ndims*ncomp.out
  } else {
    ndata.out <- 1
    data.out.ind$row <- 1:nsamples
  }
  nparams.out <- ifelse(update.params, nsamples*nparams, 1)
  # prior specification
  prior2 <- prior
  if(debug) browser()
  prior <- format.prior(model = model, prior = prior)
  # random walk jump size
  if(is.null(rw.jump.sd)) {
    rw.jump.sd <- abs(init.params)/4
    if(nmiss0 > 0) rw.jump.sd <- c(rw.jump.sd, abs(init.data[par.index[1]+(1:nmiss0)])/4)
  }
  if(length(rw.jump.sd) < nparams + nmiss0)
    stop("Incorrectly specified random walk jump sd's (need at least nparams + nmiss0).")
  if(missing(fixed.params)) fixed.params <- rep(FALSE, nparams)
  # last missing output
  nmissN <- ndims - par.index[ncomp]
  if(nmissN == 0) last.miss.out <- FALSE
  if(verbose) {
    message("Output size:")
    if(update.params) message("params = ", round(nparams.out, 2))
    if(update.data) message("data = ", round(ndata.out, 2))
    #if(update.multi & log.multi.acc)
    #  message("multi.acc = ", round(nmulti.out, 2))
    if(verbose > 1) {
      ans <- readline("Press Q to quit, any other key to proceed: ")
      if(substr(ans, 1, 1) == "Q") {
        message("Ended by user.")
        return()
      }
    }
    message("Running posterior sampler...")
  }
  if(debug) browser()
  tm <- chrono()
  # compute
  ans <- .hestPost(initParams = as.double(init.params),
                   initData = as.double(t(init.data)),
                   dT = as.double(dt),
                   nDimsPerObs = as.integer(par.index),
                   fixedParams = as.logical(fixed.params),
                   nSamples = as.integer(nsamples),
                   burn = as.integer(burn),
                   nParamsOut = as.integer(nparams.out),
                   nDataOut = as.integer(ndata.out),
                   dataOutRow = as.integer(data.out.row-1),
                   dataOutCol = as.integer(data.out.col-1),
                   updateParams = as.double(update.params),
                   updateData = as.double(update.data),
                   priorType = which(prior$type == .PriorTypes),
                   priorParams = prior$args,
                   rwJumpSd = as.double(rw.jump.sd),
                   updateLogLik = as.integer(loglik.out),
                   nLogLikOut = as.integer(ifelse(loglik.out, nsamples, 1)),
                   updateLastMiss = as.integer(last.miss.out),
                   nLastMissOut = as.integer(ifelse(last.miss.out, nsamples*nmissN, 1)))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("paramsOut", "dataOut", "paramAccept", "gibbsAccept",
                  "logLikOut", "lastMissOut", "lastIter")
  # acceptance rates
  if(debug) browser()
  accept <- NULL
  if(update.data) {
    accept <- c(accept, list(data = ans$gibbsAccept/(nsamples+burn)))
    if(verbose) {
      message("Gibbs accept: ", signif(mean(accept$data[par.index < ndims])*100,3), "%")
    }
  }
  if(update.params) {
    accept <- c(accept, list(params = ans$paramAccept[1:nparams]/(nsamples+burn)))
    if(verbose) {
      for(ii in 1:nparams)
        message(param.names[ii], " accept: ", signif(accept$params[ii]*100,3), "%")
    }
  }
  if(update.data && (nmiss0 > 0)) {
    accept <- c(accept,
                list(miss0 = ans$paramAccept[nparams+(1:nmiss0)]/(nsamples+burn)))
    for(ii in 1:nmiss0) {
      message(data.names[par.index[1] + ii], "0 accept: ",
              signif(accept$miss0[ii]*100,3), "%")
    }
  }
  out <- list()
  if(update.params) {
    out <- c(out, list(params = matrix(ans$paramsOut,
                         ncol = nparams, byrow = TRUE,
                         dimnames = list(NULL, param.names))))
  } else out <- c(out, list(params = init.params))
  if(update.data) {
    out <- c(out, list(data = aperm(array(ans$dataOut,
                         dim = c(ndims, ncomp.out, nsamples.out),
                         dimnames = list(data.names, NULL, NULL)), perm = 3:1)))
  } else out <- c(out, list(data = init.data))
  if(loglik.out) out <- c(out, list(loglik = ans$logLikOut))
  out <- c(out, list(dt = dt, par.index = par.index, data.out.ind = data.out.ind,
                     init.data = init.data, init.params = init.params,
                     rw.jump.sd = rw.jump.sd, prior = prior2))
  last.iter <- list(params = ans$lastIter[1:nparams],
                    data = matrix(ans$lastIter[nparams + 1:(ncomp*ndims)],
                      ncomp, ndims, byrow = TRUE))
  names(last.iter$params) <- param.names
  colnames(last.iter$data) <- data.names
  out <- c(out, list(last.iter = last.iter))
  if(last.miss.out) {
    last.miss <- matrix(ans$lastMissOut, nsamples, nmissN, byrow = TRUE)
    colnames(last.miss) <- data.names[par.index[ncomp]+(1:nmissN)]
    out <- c(out, list(last.miss = last.miss))
  }
  out <- c(out, list(accept = accept))
  out
}
