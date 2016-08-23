###########################################

# the MSDE package, final version

# mlysy and tkitt 2014

# create the package skeleton etc for the package.

require(devtools)
require(Rcpp)

<<<<<<< HEAD
pkg.name <- "msdeTest"
pkg.dir <- getwd()

#--- create package skeleton (do this only once) --------------------------------

# Rcpp.package.skeleton(name = pkg.name, cpp_files = "msdeFunctions.cpp",
#                       example_code = FALSE)


#--- install package ------------------------------------------------------------

compileAttributes()
document()
install()
=======
pkg <- "C:/Users/Jerome/Documents/R/msdeTest"

#--- install package ------------------------------------------------------------

compileAttributes(pkg)
document(pkg)
install(pkg)
>>>>>>> rename


#--- build ----------------------------------------------------------------------

build()
<<<<<<< HEAD
=======

#--- create package skeleton (do this only once) --------------------------------

# Rcpp.package.skeleton(name = pkg.name, cpp_files = "msdeFunctions.cpp",
#                       example_code = FALSE)


>>>>>>> rename
