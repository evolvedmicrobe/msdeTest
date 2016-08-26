#--- speed comparisons with old version of package ------------------------------

# the old msde package was written without OOP and (*) merged all cpp/h
# files into 1 before compile on-the-fly with Rcpp
# the new version is written with OOP and checks the speed effect of (*)

require(msde) # old version
require(msdeTest) # new version
pkg1 <- function(fun) {
  eval(parse(text = paste0("msde::", fun)))
}
pkg2 <- function(fun) {
  eval(parse(text = paste0("msdeTest::", fun)))
}

# the test model is called heston's model, which is an sde with
# two components y = (x,z), and drift (2x1 vector) and diffusion
# (2x2 upper tri chol matrix) coefficients given by the following
# functions.  for more info check vignette("msde-vignette")
hest.dr <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .125 * Z^2, theta[,3]/Z - .5*theta[,2]*Z)
}
hest.df <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cv <- .5*theta[,5]*theta[,4]*Z
  ans <- cbind(.25 * Z^2, cv, cv, theta[,4]^2)
  # output as Nx4 matrix, shitty coding to save time
  t(apply(ans, 1, function(x) chol(matrix(x,2,2))))
}

# create model objects (i.e., compile c++ code)

# old version (msde)
hmod <- pkg1("sde.make.model")(list = hestList,
                               showOutput = TRUE, rebuild = TRUE)

# new version (msdeTest -- interface still in progress)
ndims <- 2
nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod2 <- file.path(pkg2(":.msdeCppPath"), "hestModel.cpp")
hmod2 <- pkg2("sde.make.model")(code = hmod2,
                                ndims = ndims, nparams = nparams,
                                param.names = param.names,
                                data.names = data.names,
                                showOutput = TRUE, rebuild = TRUE)


#--- check drift/diffusion in the c++ code --------------------------------------
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# R code
dr.R <- hest.dr(X = X0[,1], Z = X0[,2], theta = Theta)
df.R <- hest.df(X = X0[,1], Z = X0[,2], theta = Theta)
# msde
dr.m <- pkg1("sde.drift")(model = hmod, x = X0, theta = Theta)
df.m <- pkg1("sde.diff")(model = hmod, x = X0, theta = Theta)
# msdeTest (separate .o files)
dr.mto <- hest.drift(model = hmod, x = X0, theta = Theta)
df.mto <- hest.diff(model = hmod, x = X0, theta = Theta)
# msdeTest (single .o file)
dr.mt <- pkg1("sde.drift")(model = hmod2, x = X0, theta = Theta)
df.mt <- pkg2("sde.diff")(model = hmod2, x = X0, theta = Theta)

# check calcs
sapply(list(m = dr.m, mt = dr.mt, mto = dr.mto),
       function(x) range(dr.R - x))
sapply(list(m = df.m, mt = df.mt, mto = df.mto),
       function(x) range(df.R - x))

# ok.  speed comparisons
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 5e6
Theta <- matrix(theta, nReps, nparams, byrow = TRUE)
colnames(Theta) <- param.names
X0 <- matrix(x0, nReps, ndims, byrow = TRUE)
colnames(X0) <- data.names

time.m <- system.time({
  dr.m <- pkg1("sde.drift")(model = hmod, x = X0, theta = Theta)
  df.m <- pkg1("sde.diff")(model = hmod, x = X0, theta = Theta)
})
time.mto <- system.time({
  dr.mto <- hest.drift(model = hmod, x = X0, theta = Theta)
  df.mto <- hest.diff(model = hmod, x = X0, theta = Theta)
})
time.mt <- system.time({
  dr.mt <- pkg2("sde.drift")(model = hmod2, x = X0, theta = Theta)
  df.mt <- pkg2("sde.diff")(model = hmod2, x = X0, theta = Theta)
})

# nice boost
c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
##  m.elapsed  mt.elapsed mto.elapsed
##  1.0000000   0.8846154   1.0000000

#--- loglikelihood evaluations --------------------------------------------------

nReps <- 1e4
nObs <- 1e3
dT <- 1/252

# initialize
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# generate data
hsim <- pkg1("sde.sim")(model = hmod, init.data = X0, params = Theta,
                        dt = dT, dt.sim = dT, N = nObs, nreps = nReps)

# loglikelihood calcs
time.m <- system.time({
  ll.m <- pkg1("sde.loglik")(model = hmod, x = X0, dt = dT,
                             theta = Theta)
})
time.mto <- system.time({
  ll.mto <- hest.loglik(model = hmod, x = X0, dt = dT,
                        theta = Theta)
})
time.mt <- system.time({
  ll.mt <- pkg2("sde.loglik")(model = hmod2, x = X0, dt = dT,
                              theta = Theta)
})

range(ll.m-ll.mt)
range(ll.m-ll.mto)

# still getting overhead. now with loglik as a method of sdeData
# way worse performance
c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
##  m.elapsed  mt.elapsed mto.elapsed
##   1.000000    1.346966    1.245383


#--- forward simulation ---------------------------------------------------------

# I've checked correctness of code in test-sde.sim.  it's a bit sloppy but
# it's there.
# if you want to do speed comparisons with exactly the same output
# (otherwise output is stochastic) then change this to true
same.rnd <- TRUE
SEED <- 3843

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# each of these takes about 15s on my lappy
nObs <- 1e3
nReps <- 1e3
dT <- 1/252
if(same.rnd) set.seed(SEED)
time.m <- system.time({
  hsim.m <- pkg1("sde.sim")(model = hmod, init.data = x0, params = theta,
                            dt = dT, dt.sim = dT/100,
                            N = nObs, nreps = nReps)
})
if(same.rnd) set.seed(SEED)
time.mt <- system.time({
  hsim.mt <- pkg2("sde.sim")(model = hmod2, init.data = x0, params = theta,
                             dt = dT, dt.sim = dT/100,
                             N = nObs, nreps = nReps)
})
if(same.rnd) set.seed(SEED)
time.mto <- system.time({
  hsim.mto <- hest.sim(model = hmod2, init.data = x0, params = theta,
                       dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
})

if(same.rnd) {
  rbind(mt = sapply(names(hsim.mt),
          function(nm) max(abs(hsim.m[[nm]] - hsim.mt[[nm]]))),
        mto = sapply(names(hsim.mt),
          function(nm) max(abs(hsim.m[[nm]] - hsim.mto[[nm]]))))
}


# irritating loss of performance
c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
##  m.elapsed  mt.elapsed mto.elapsed
##   1.000000    1.028704    1.085337


#--- MCMC draws -----------------------------------------------------------------

same.rnd <- TRUE
SEED <- 2531

# simulate data
nObs <- 250
dT <- 1/252
hsim <- pkg1("sde.sim")(model = hmod, init.data = x0, params = theta,
                        dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

# initialize MCMC rv's
init <- pkg2("sde.init")(data = hsim$data, dt = dT, m = 2,
                         par.index = c(2, rep(nObs-1, 1)),
                         params = theta)

# prior
prior <- list(Mu = c(.1, .35, 1.0, .5, -.81),
              V = crossprod(matrix(rnorm(25),5)))
prior$V <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$V)
prior$V <- prior$V %*% sqrt(diag(c(.1, 8, .15, .002, .002)))

# mcmc specs
rw.jump.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
update.params <- TRUE
update.data <- TRUE
nsamples <- ifelse(update.data, 2e4, 4e4)
burn <- 1e3

if(same.rnd) set.seed(SEED)
time.m <- system.time({
  hpost.m <- pkg1("sde.post")(model = hmod, init = init,
                              nsamples = nsamples, burn = burn,
                              prior = prior,
                              rw.jump.sd = rw.jump.sd,
                              update.params = update.params,
                              update.data = update.data)
})

if(same.rnd) set.seed(SEED)
time.mt <- system.time({
  hpost.mt <- pkg2("sde.post")(model = hmod2, init = init,
                               nsamples = nsamples, burn = burn,
                               prior = prior,
                               rw.jump.sd = rw.jump.sd,
                               update.params = update.params,
                               update.data = update.data)
})

if(same.rnd) set.seed(SEED)
time.mto <- system.time({
  hpost.mto <- hest.post(model = hmod, init = init,
                         nsamples = nsamples, burn = burn,
                         prior = prior,
                         rw.jump.sd = rw.jump.sd,
                         update.params = update.params,
                         update.data = update.data)
})

if(same.rnd) {
  rbind(mt = sapply(c("params", "data"),
          function(nm) max(abs(hpost.m[[nm]] - hpost.mt[[nm]]))),
        mto = sapply(c("params", "data"),
          function(nm) max(abs(hpost.m[[nm]] - hpost.mto[[nm]]))))
}


# gross
c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
##  m.elapsed  mt.elapsed mto.elapsed
##   1.000000    1.087065    1.210614
