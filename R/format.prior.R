#--- prior ----------------------------------------------------------------------

# names of provided priors
.PriorTypes <- c("flat", "normal")
.NormalNames <- sort(c("Mu", "V"))
.GCopNames <- sort(c("dens.x", "dx", "rx", "dens.y", "ldens.y", "Dens.y",
                     "Rho", "mean", "sd"))
# format the prior
format.prior <- function(model, prior, debug = FALSE) {
  nparams <- model$nparams
  ndims <- model$ndims
  nmax.rv <- nparams + ndims
  # check if there's a custom prior
  if(!is.null(model$logCustomPrior)) {
    has.custom.prior <- TRUE
    custom.names <- model$custom.names
    if(!is.null(custom.names)) custom.names <- sort(custom.names)
  } else {
    has.custom.prior <- FALSE
    custom.names <- NULL
  }
  # make sure custom.names doesn't conflict with normal or gcop
  # although this shoudl probably be done during sde.make.model
  if(!is.null(custom.names)) {
    if(identical(custom.names, .NormalNames) || identical(custom.names, .GCopNames)) {
      stop("custom.prior parameter names conflict with those of Normal or GCop prior.")
    }
  }
  # prior specification
  prior.names <- names(prior)
  if(!is.null(prior.names)) prior.names <- sort(prior.names)
  #prior2 <- prior
  prior.type <- NULL
  if(debug) browser()
    if(is.null(prior)) {
    # flat prior
    prior.type <- "flat"
    prior <- list()
  } else if (identical(prior.names, .NormalNames)) {
    # mv prior
    prior.type <- "normal"
    nrv <- length(prior$Mu)
    mv.dim <- c(Mu = nrv, V = nrv^2)
    if(nrv > nmax.rv) {
      stop("Too many random variables in prior specification.")
    }
    if(any(sapply(prior, length)[prior.names] != mv.dim[prior.names])) {
      stop("Prior elements have incompatible dimensions.")
    }
    prior$V <- chol(prior$V)
  } else if(identical(prior.names, .GCopNames)) {
    # gcop prior
    prior.type <- "gcop"
    nrv <- length(prior$dens.x)
    gcop.dim <- c(dens.x = nrv, dx = nrv, rx = 2*nrv,
                  dens.y = nrv, ldens.y = nrv, Dens.y = nrv,
                  Rho = nrv^2, mean = nrv, sd = nrv)
    if(nrv > nmax.rv) {
      stop("Too many random variables in prior specification.")
    }
    if(any(sapply(prior, length)[prior.names] != gcop.dim[prior.names])) {
      stop("Prior elements have incompatible dimensions.")
    }
    prior$Rho <- chol(prior$Rho)
    prior$nbreaks <- sapply(prior$dens.x, length)
    prior <- sapply(prior, unlist)
  } else {
    # custom prior
    prior.type <- "custom"
    if(!has.custom.prior) {
      stop("prior is not recognized default and no custom prior supplied.")
    }
    if(!is.null(custom.names) && !identical(prior.names, custom.names))
      stop("Supplied custom prior names don't match those defined by sde.model.")
    prior <- sapply(prior, unlist, recursive = TRUE)
  }
  # debug feature
  if(!prior.type %in% .PriorTypes) {
    stop("Specified prior type not supported.")
  }
  list(type = prior.type, args = prior)
}
