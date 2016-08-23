# some unit tests

require(msdeTest)
require(msde)
require(multiplot)

hest.model <- sde.make.model(list = hestList)

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

#--- drift and diffusion functions ----------------------------------------------

nReps <- 1e3
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

system.time({
  dr1 <- sde.drift(model = hest.model, x = X0, theta = Theta)
  df1 <- sde.diff(model = hest.model, x = X0, theta = Theta)
})

system.time({
  dr2 <- hest.drift(model = hest.model, x = X0, theta = Theta)
  df2 <- hest.diff(model = hest.model, x = X0, theta = Theta)
})

# this is just for speed comparisons
system.time({
  dr3 <- matrix(hest.model$drift(xIn = t(X0), thetaIn = t(Theta), nReps = nReps),
                nrow = nReps, byrow = TRUE)
  df3 <- matrix(hest.model$diff(xIn = t(X0), thetaIn = t(Theta), nReps = nReps),
                nrow = nReps, byrow = TRUE)
})

range(dr1-dr2)
range(df1-df2)

#--- loglikelihood --------------------------------------------------------------

nReps <- 1e3
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

dT <- 1/252
N <- 1e3

hsim <- sde.sim(model = hest.model, init.data = X0, params = Theta,
                dt = dT, dt.sim = dT, N = N, nreps = nReps)

system.time({
  ll1 <- sde.loglik(model = hest.model, x = hsim$data, dt = dT, theta = hsim$params)
})
system.time({
  ll2 <- sde.loglik(model = hest.model, x = hsim$data, dt = dT, theta = hsim$params)
  #ll2 <- hest.loglik(xIn = aperm(hsim$data, 3:1), thetaIn = t(hsim$params),
  #                   dT = rep(dT, N-1), nComp = N, nReps = nReps)
})

range(ll1-ll2)

#--- prior ----------------------------------------------------------------------

hest.model <- sde.make.model(list = hestList)

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 1e5
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)
dT <- 1/252
N <- 1e1
hsim <- sde.sim(model = hest.model, init.data = X0, params = Theta,
                dt = dT*N, dt.sim = dT, N = 2, nreps = nReps)


X <- hsim$data[,2,]
RV <- cbind(X, Theta)
Mu <- colMeans(RV)
V <- var(RV)

nRV <- sample(ncol(RV), 1)
ind <- ncol(RV)-nRV+1:nRV
lpi1 <- dmvnorm(x = RV[,ind,drop=FALSE], mean = Mu[ind,drop=FALSE],
                sigma = V[ind,ind,drop=FALSE], log = TRUE)
lpi2 <- hest.logprior(xIn = t(X), thetaIn = t(Theta),
                      priorParams = list(Mu = Mu[ind], V = chol(V[ind,ind])),
                      priorType = 2, nRv = nRV, nReps = nReps)
range(diff(lpi1-lpi2))

#--- ok mcmc algorithm ----------------------------------------------------------

require(msdeTest)
require(msde)
require(multiplot)

hest.model <- sde.make.model(list = hestList)

# ok just missing data, no first observation

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nObs <- 20
dT <- 1/252
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

init <- sde.init(data = hsim$data, dt = dT, k = 0,
                 par.index = c(2,sample(0:2, nObs-1, replace = TRUE)),
                 params = theta)

nsamples <- 1e5
burn <- 100
hpost1 <- sde.post(model = hest.model, init = init,
                   nsamples = nsamples, burn = burn, prior = NULL,
                   update.params = FALSE, data.out.ind = nsamples)

hpost2 <- hest.post(model = hest.model, init = init,
                    nsamples = nsamples, burn = burn, prior = NULL,
                    update.params = FALSE, data.out.ind = nsamples)


ind <- sample(nObs, 1)
mplot.dens(list(msde = hpost1$data[,ind,], msdeTest = hpost2$data[,ind,]))

# ok now with parameters

nObs <- 50
dT <- 1/252
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

init <- sde.init(data = hsim$data, dt = dT, k = 0,
                 par.index = 2,
                 params = theta)

prior <- list(Mu = c(.1, .35, 1.0, .5, -.81),
              V = crossprod(matrix(rnorm(25),5)))
prior$V <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$V)
prior$V <- prior$V %*% sqrt(diag(c(.1, 8, .15, .002, .002)))

nsamples <- 1e5
burn <- 1e3
rw.jump.sd <- c(.1, 1, .1, .01, .01)
hpost1 <- sde.post(model = hest.model, init = init,
                   nsamples = nsamples, burn = burn, prior = prior,
                   rw.jump.sd = rw.jump.sd,
                   update.params = TRUE, update.data = FALSE)

hpost2 <- hest.post(model = hest.model, init = init,
                    nsamples = nsamples, burn = burn, prior = prior,
                    rw.jump.sd = rw.jump.sd,
                    update.params = TRUE, update.data = FALSE)

mplot.dens(list(msde = hpost1$params, msdeTest = hpost2$params))
