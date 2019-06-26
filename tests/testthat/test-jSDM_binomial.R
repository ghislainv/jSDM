context("test-jSDM_binomial")

#==================
#== Data simulation

#= Number of sites
nsite <- 200

#= Set seed for repeatability
seed <- 1234

#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
set.seed(seed)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1)
logit.theta <- X %*% beta.target
theta <- inv_logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Data-sets
data.obs <- data.frame(Y,visits,x1,x2)

#==================
#== Fit the model 
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial(presences=data.obs$Y,
                      trials=data.obs$visits,
                      suitability=~x1+x2,
                      data=data.obs,
                      burnin=burnin, mcmc=mcmc, thin=thin,
                      beta_start=0,
                      mubeta=0, Vbeta=1.0E6,
                      seed=1234, ropt=0.44, verbose=1)

test_that("jSDM_binomial works", {
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,1))
  expect_equal(dim(mod$mcmc),c(nsamp,ncol(X)+1))
  expect_equal(sum(is.na(mod$mcmc)),0)
})