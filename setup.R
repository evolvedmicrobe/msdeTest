###########################################

# the MSDE package, final version

# mlysy and tkitt 2014

# create the package skeleton etc for the package.

require(devtools)
require(Rcpp)

pkg <- "C:/Users/Jerome/Documents/R/msdeTest"

#--- install package ------------------------------------------------------------

compileAttributes(pkg)
document(pkg)
install(pkg)

#--- build ----------------------------------------------------------------------

build()

#--- create package skeleton (do this only once) --------------------------------

# Rcpp.package.skeleton(name = pkg.name, cpp_files = "msdeFunctions.cpp",
#                       example_code = FALSE)

