#context("test-hidden")

# Data simulation
#= Number of sites
nsite <- 50
#= Number of species
nsp <- 5
#= Set seed for repeatability
seed <- 1234
#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1
#= Number of latent variables
n_latent <- 2

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
#= Latent variables W
W <- matrix(rnorm(nsite*n_latent,0,1), nsite, n_latent)
#= Factor loading lambda
lambda.target <- matrix(0, n_latent, nsp)
mat <- t(matrix(runif(nsp*n_latent, -2, 2), byrow=TRUE, nrow=nsp))
lambda.target[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]

log.theta <- X %*% t(beta.target) + W%*%lambda.target
theta <- exp(log.theta)
set.seed(seed)
Y.pois <- apply(theta, 2, rpois, n=nsite) 


#= Number of visits associated to each site
visits <- rpois(nsite,3)
visits[visits==0] <- 1
logit.theta <- X %*% t(beta.target)
theta <- jSDM::inv_logit(logit.theta)
set.seed(seed)
Y.bin <- apply(theta, 2, rbinom, n=nsite, size=visits)


test_that("checks works", {
  expect_equal(check.mcmc.parameters(1000, 1000, 1), 0)
  expect_equal(check.verbose(1), 0)
  expect_equal(check.Y.poisson(Y.pois), 0)
  expect_equal(check.Y.binomial(Y.bin, T=visits), 0)
  expect_equal(form.beta.start(0,np), rep(0,np))
  expect_equal(form.beta.start(rep(0,np),np), rep(0,np))
  expect_equal(form.beta.start.sp(0,np, nsp), matrix(0, np, nsp))
  expect_equal(form.beta.start.sp(0,np, nsp), matrix(0, np, nsp))
  expect_equal(check.Vbeta.mat(1,np), diag(rep(1,np)))
  expect_equal(check.Vbeta.mat(rep(1,np),np), diag(rep(1,np)))
  expect_equal(check.Vbeta.mat(diag(rep(1,np)),np), diag(rep(1,np)))
  expect_equal(check.Vlambda.mat(1, n_latent), diag(rep(1,n_latent)))
  expect_equal(check.Vlambda.mat(rep(1, n_latent), n_latent), diag(rep(1,n_latent)))
  expect_equal(check.Vlambda.mat(diag(rep(1,n_latent)), n_latent), diag(rep(1,n_latent)))
})