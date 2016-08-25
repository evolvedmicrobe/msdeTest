#--- on-the-fly compile with sourceCpp ------------------------------------------

require(msdeTest)
require(msde)
require(multiplot)

hest.model <- sde.make.model(list = hestList)

ndims <- 2
nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")

hest.model2 <- sde.make.model2(ndims = ndims, nparams = nparams,
                               param.names = param.names,
                               data.names = data.names, debug = FALSE)


#--- ok speed test --------------------------------------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nObs <- 20
dT <- 1/252
set.seed(382)
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)
set.seed(382)
hsim2 <- sde.sim2(model = hest.model2, init.data = x0, params = theta,
                  dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

range(hsim$data - hsim2$data)
range(hsim$params - hsim2$params)

# missing data test

init <- sde.init(data = hsim$data, dt = dT, k = 0,
                 par.index = c(2,sample(0:2, nObs-1, replace = TRUE)),
                 params = theta)

set.seed(979)
nsamples <- 1e3
burn <- 0
hpost1 <- sde.post(model = hest.model, init = init,
                   nsamples = nsamples, burn = burn, prior = NULL,
                   update.params = FALSE, data.out.ind = nsamples)
set.seed(979)
hpost2 <- sde.post2(model = hest.model2, init = init,
                    nsamples = nsamples, burn = burn, prior = NULL,
                    update.params = FALSE, data.out.ind = nsamples)

range(hpost1$data - hpost2$data)
range(hpost1$params - hpost2$params)

# ok now with parameters

nObs <- 250
dT <- 1/252
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

init <- sde.init(data = hsim$data, dt = dT, m = 2,
                 par.index = c(2, rep(nObs-1, 1)),
                 params = theta)

prior <- list(Mu = c(.1, .35, 1.0, .5, -.81),
              V = crossprod(matrix(rnorm(25),5)))
prior$V <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$V)
prior$V <- prior$V %*% sqrt(diag(c(.1, 8, .15, .002, .002)))

nsamples <- 5e4
burn <- 1e3
rw.jump.sd <- c(.1, 1, .1, .01, .01)
update.params <- TRUE
update.data <- TRUE

set.seed(2531)
hpost1 <- sde.post(model = hest.model, init = init,
                   nsamples = nsamples, burn = burn, prior = prior,
                   rw.jump.sd = rw.jump.sd,
                   update.params = update.params,
                   update.data = update.data)
set.seed(2531)
hpost2 <- sde.post2(model = hest.model2, init = init,
                    nsamples = nsamples, burn = burn, prior = prior,
                    rw.jump.sd = rw.jump.sd,
                    update.params = update.params,
                    update.data = update.data)

range(hpost1$data - hpost2$data)
range(hpost1$params - hpost2$params)

mplot.dens(list(msde = hpost1$params, msdeTest = hpost2$params))
