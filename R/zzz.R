# tell package where to find templates folder
.onLoad <- function(libname, pkgname) {
  assign(".msdeCppPath", file.path(libname, pkgname, "cppTemplates"), envir = parent.env(environment()))
}
