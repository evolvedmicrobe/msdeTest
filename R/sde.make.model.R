#' Create sde.model object
#'
#' @export
<<<<<<< HEAD
sde.make.model2 <- function(code, ndims, nparams,
                            data.names, param.names, custom.names,
                            cpp.out = TRUE, ..., debug = FALSE) {
=======
sde.make.model <- function(code, ndims, nparams,
                           data.names, param.names, custom.names,
                           cpp.out = TRUE, ..., debug = FALSE) {
>>>>>>> rename
  # default data and parameter names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) stop("Incorrect data.names.")
  if(length(param.names) != nparams) stop("Incorrect param.names.")
  if(missing(custom.names)) custom.names <- NULL
  sde.model <- list(ndims = ndims, nparams = nparams,
                    data.names = data.names, param.names = param.names,
                    custom.names = custom.names)
  # check for custom prior
  #if(missing(logCustomPrior)) logCustomPrior <- NULL
  #if(missing(custom.names)) custom.names <- NULL
  #sde.model <- list(ndims = ndims, nparams = nparams,
  #                  data.names = data.names, param.names = param.names,
  #                  sdeDr = sdeDr, sdeDf = sdeDf, isValidData = isValidData, isValidParams = isValidParams)
  #if(!is.null(logCustomPrior)) {
  #  sde.model <- c(sde.model, list(logCustomPrior = logCustomPrior))
  #  if(!is.null(custom.names)) sde.model <- c(sde.model, list(custom.names = custom.names))
  #} else {
  #  logCustomPrior <- "double CustomPrior::logPrior(double params[], double x[]) {
  #return(0.0);
<<<<<<< HEAD
#}"
=======
  #}"
>>>>>>> rename
  #}
  # create c++ code
  header.names <- c("sdeModel.h", "mvnUtils.h", "sdeUtils.h",
                    "Prior.h", "mcmcUtils.h", "sdeMCMC.h")
  cpp.names <- c("sdeUtils.cpp", "Prior.cpp",
                 "mvnUtils.cpp", "mcmcUtils.cpp", "sdeMCMC.cpp",
                 "missGibbsUpdate.cpp", "paramVanillaUpdate.cpp",
                 "sdeDebug.cpp", "EulerSim.cpp", "SimpleMCMC.cpp")
  #.msdeCppPath <- "C:/Users/Jerome/Documents/R/msdeTest/inst/cppTemplate"
  hLines <- sapply(header.names, function(nm) {
    readLines(file.path(.msdeCppPath, nm))
  })
  cppLines <- sapply(cpp.names, function(nm) {
    readLines(file.path(.msdeCppPath, nm))
  })
  cpp.code <- c(unlist(hLines), readLines(code), unlist(cppLines),"\n")
  # strip redundant include statements
  cpp.code <- cpp.code[!sapply(cpp.code, grepl, pattern = "^#.*")]
  cpp.code <- cpp.code[!sapply(cpp.code, grepl,
                               pattern = "using namespace Rcpp;")]
  # reappend #include Rcpp
  cpp.code <- c("#include <Rcpp.h>", "using namespace Rcpp;", cpp.code)
  if(cpp.out) sde.model <- c(sde.model, list(cpp.code = cpp.code))

  # compile c++ code
  if(debug) browser()
  sourceCpp(code = paste(cpp.code, collapse = "\n"),
            env = environment(), ...)
  environment(sde.model$sim) <- globalenv()
  environment(sde.model$post) <- globalenv()
  environment(sde.model$drift) <- globalenv()
  environment(sde.model$diff) <- globalenv()
  environment(sde.model$loglik) <- globalenv()

  class(sde.model) <- "sde.model"
  sde.model
}
