#--- check that two versions of same package produce identical results -----------

pkg1 <- "msdeTest"
pkg2 <- "msdeTest2"

parsePkg <- function(pkg, fun) {
  eval(parse(text = paste0(pkg, "::", fun)))
}

require(pkg1, character.only = TRUE)
require(pkg2, character.only = TRUE)

# compile code
ndims <- 2
nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")

hmod <- parsePkg(pkg1, "sde.make.model2")
hmod <- hmod(code = file.path(msdeTest:::.msdeCppPath,
               "hestModel.cpp"),
             ndims = ndims, nparams = nparams,
             param.names = param.names,
             data.names = data.names)

hmod2 <- parsePkg(pkg2, "sde.make.model")
hmod2 <- hmod2(code = file.path(msdeTest:::.msdeCppPath,
                 "hestModel.cpp"),
               ndims = ndims, nparams = nparams,
               param.names = param.names,
               data.names = data.names)

#--- drift/diffusion ------------------------------------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

dr <- parsePkg(pkg1, "sde.drift2")(model = hmod, x = X0, theta = Theta)
df <- parsePkg(pkg1, "sde.diff2")(model = hmod, x = X0, theta = Theta)

dr2 <- parsePkg(pkg2, "sde.drift")(model = hmod2, x = X0, theta = Theta)
df2 <- parsePkg(pkg2, "sde.diff")(model = hmod2, x = X0, theta = Theta)

range(dr - dr2)
range(df - df2)

#--- forward simulation ---------------------------------------------------------

nReps <- 1e3

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# each of these takes about 15s on my lappy
nObs <- 1e3
dT <- 1/252

set.seed(3843)
hsim <- parsePkg(pkg1, "sde.sim2")
hsim <- hsim(model = hmod, init.data = X0, params = Theta,
             dt = dT, dt.sim = dT/10, N = nObs, nreps = nReps)

set.seed(3843)
hsim2 <- parsePkg(pkg2, "sde.sim")
hsim2 <- hsim2(model = hmod, init.data = X0, params = Theta,
               dt = dT, dt.sim = dT/10, N = nObs, nreps = nReps)

identical(hsim, hsim2)

#--- loglikelihood --------------------------------------------------------------

ll1 <- parsePkg(pkg1, "sde.loglik2")
ll1 <- ll1(model = hmod, x = hsim$data, dt = dT,
           theta = hsim$params)

ll2 <- parsePkg(pkg2, "sde.loglik")
ll2 <- ll2(model = hmod2, x = hsim$data, dt = dT,
           theta = hsim$params)

range(ll1-ll2)

#--- posterior sampling ---------------------------------------------------------

nObs <- 250
dT <- 1/252
hsim <- msdeTest::sde.sim2(model = hmod, init.data = x0,
                           params = theta,
                           dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

# initialize MCMC rv's
init <- msdeTest2::sde.init(data = hsim$data, dt = dT, m = 2,
                            par.index = c(2, rep(nObs-1, 1)),
                            params = theta)

# prior
prior <- list(Mu = c(.1, .35, 1.0, .5, -.81),
              V = crossprod(matrix(rnorm(25),5)))
prior$V <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$V)
prior$V <- prior$V %*% sqrt(diag(c(.1, 8, .15, .002, .002)))

# mcmc specs
nsamples <- 1e3
burn <- 1e2
rw.jump.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
update.params <- TRUE
update.data <- TRUE

set.seed(531)
hpost <- parsePkg(pkg1, "sde.post2")
hpost <- hpost(model = hmod, init = init,
               nsamples = nsamples, burn = burn, prior = prior,
               rw.jump.sd = rw.jump.sd,
               update.params = update.params,
               update.data = update.data)

set.seed(531)
hpost2 <- parsePkg(pkg2, "sde.post")
hpost2 <- hpost2(model = hmod2, init = init,
                 nsamples = nsamples, burn = burn, prior = prior,
                 rw.jump.sd = rw.jump.sd,
                 update.params = update.params,
                 update.data = update.data)

identical(hpost, hpost2)
