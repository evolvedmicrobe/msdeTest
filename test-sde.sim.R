#--- test simulation ------------------------------------------------------------

require(msdeTest)
require(msde)
require(multiplot)

hest.model <- sde.make.model(list = hestList)
ndims <- 2
nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hest.model2 <- sde.make.model2(code = file.path(msdeTest:::.msdeCppPath,
                                 "hestModel.cpp"),
                               ndims = ndims, nparams = nparams,
                               param.names = param.names,
                               data.names = data.names)

# ok just missing data, no first observation

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

nObs <- 20
dT <- 1/252

# msde
set.seed(777)
hsim <- sde.sim(model = hest.model, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
# msdeTest - multi-o
set.seed(777)
hsim2 <- hest.sim(model = hest.model, init.data = x0, params = theta,
                  dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
# msdeTest - single-o
set.seed(777)
hsim3 <- sde.sim2(model = hest.model2, init.data = x0, params = theta,
                  dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)

range(hsim$data-hsim2$data)
range(hsim$params-hsim2$params)
range(hsim$data-hsim3$data)
range(hsim$params-hsim3$params)
