#--- speed comparisons with old version of package ------------------------------

# the old msde package was written without OOP and (*) merged all cpp/h
# files into 1 before compile on-the-fly with Rcpp
# the new version is written with OOP and checks the speed effect of (*)

require(msde) # old version
require(msdeTest) # new version

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
hest.model <- sde.make.model(list = hestList)

# new version (msdeTest -- interface still in progress)
ndims <- 2
nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hest.model2 <- sde.make.model2(code = file.path(msdeTest:::.msdeCppPath,
                                 "hestModel.cpp"),
                               ndims = ndims, nparams = nparams,
                               param.names = param.names,
                               data.names = data.names)


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
dr.m <- sde.drift(model = hest.model, x = X0, theta = Theta)
df.m <- sde.diff(model = hest.model, x = X0, theta = Theta)
# msdeTest (separate .o files)
dr.mto <- hest.drift(model = hest.model, x = X0, theta = Theta)
df.mto <- hest.diff(model = hest.model, x = X0, theta = Theta)
# msdeTest (single .o file)
dr.mt <- sde.drift2(model = hest.model2, x = X0, theta = Theta)
df.mt <- sde.diff2(model = hest.model2, x = X0, theta = Theta)

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
  dr.m <- sde.drift(model = hest.model, x = X0, theta = Theta)
  df.m <- sde.diff(model = hest.model, x = X0, theta = Theta)
})
time.mto <- system.time({
  dr.mto <- hest.drift(model = hest.model, x = X0, theta = Theta)
  df.mto <- hest.diff(model = hest.model, x = X0, theta = Theta)
})
time.mt <- system.time({
  dr.mt <- sde.drift2(model = hest.model2, x = X0, theta = Theta)
  df.mt <- sde.diff2(model = hest.model2, x = X0, theta = Theta)
})

# should be imperceptible
c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]

#--- forward simulation ---------------------------------------------------------

# I've checked correctness of code in test-sde.sim.  it's a bit sloppy but
# it's there.
# if you want to do speed comparisons with exactly the same output
# (otherwise output is stochastic) then change this to true
same.rnd <- FALSE

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# each of these takes about 15s on my lappy
nObs <- 1e3
nReps <- 1e3
dT <- 1/252
if(!same.rnd) set.seed(3843)
time.m <- system.time({
  hsim.m <- sde.sim(model = hest.model, init.data = x0, params = theta,
                    dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
})
if(!same.rnd) set.seed(3843)
time.mt <- system.time({
  hsim.mt <- sde.sim2(model = hest.model2, init.data = x0, params = theta,
                      dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
})
if(!same.rnd) set.seed(3843)
time.mto <- system.time({
  hsim.mto <- sde.sim2(model = hest.model2, init.data = x0, params = theta,
                       dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
})

c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
# why the fuck is the OOP version 5% slower??

#--- MCMC draws -----------------------------------------------------------------

same.rnd <- FALSE

# simulate data
nObs <- 250
dT <- 1/252
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

# initialize MCMC rv's
init <- sde.init(data = hsim$data, dt = dT, m = 2,
                 par.index = c(2, rep(nObs-1, 1)),
                 params = theta)

# prior
prior <- list(Mu = c(.1, .35, 1.0, .5, -.81),
              V = crossprod(matrix(rnorm(25),5)))
prior$V <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$V)
prior$V <- prior$V %*% sqrt(diag(c(.1, 8, .15, .002, .002)))

# mcmc specs
nsamples <- 2e4
burn <- 1e3
rw.jump.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
update.params <- TRUE
update.data <- TRUE

if(!same.rnd) set.seed(2531)
time.m <- system.time({
  hpost.m <- sde.post(model = hest.model, init = init,
                      nsamples = nsamples, burn = burn, prior = prior,
                      rw.jump.sd = rw.jump.sd,
                      update.params = update.params,
                      update.data = update.data)
})

if(!same.rnd) set.seed(2531)
time.mt <- system.time({
  hpost.mt <- sde.post2(model = hest.model2, init = init,
                        nsamples = nsamples, burn = burn, prior = prior,
                        rw.jump.sd = rw.jump.sd,
                        update.params = update.params,
                        update.data = update.data)
})

if(!same.rnd) set.seed(2531)
time.mto <- system.time({
  hpost.mto <- hest.post(model = hest.model, init = init,
                         nsamples = nsamples, burn = burn, prior = prior,
                         rw.jump.sd = rw.jump.sd,
                         update.params = update.params,
                         update.data = update.data)
})

c(m = time.m[3], mt = time.mt[3], mto = time.mto[3])/time.m[3]
# why is this lil bitch 7% slower??? that's 5mins slower per hour!
# why the FUCK do i learn how to program better if it just makes shit
# slower!!!
