#' Simulation and inference for multivariate stochastic differential equation models
#'
#' @useDynLib msdeTest
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name msdeTest
NULL

# tell package where to find templates folder
.onLoad <- function(libname, pkgname) {
  assign(".msdeCppPath", file.path(libname, pkgname, "cppTemplates"), envir = parent.env(environment()))
}
