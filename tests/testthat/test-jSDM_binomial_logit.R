context("test-jSDM_binomial_logit")

#==================
#== Data simulation

#= Number of sites
nsite <- 100
#= Number of species
nsp <- 50
#= Set seed for repeatability
seed <- 1234

#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
logit.theta <- X %*% t(beta.target)
theta <- inv_logit(logit.theta)
set.seed(seed)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)

#= Site-occupancy model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_logit(# Chains
                          burnin, mcmc, thin,
                          # Response variable
                          presence_site_sp=Y, trials=visits,
                          # Explanatory variables
                          site_suitability=~x1+x2, site_data=X,
                          # Starting values
                          beta_start=0,
                          # Priors
                          mu_beta=0, V_beta=1.0E6,
                          # Various
                          seed=1234, ropt=0.44, verbose=1)

test_that("jSDM_binomial_logit works", {
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(unique(lapply(mod$mcmc.betas,dim))[[1]],c(nsamp,np))
  expect_equal(length(mod$mcmc.betas),nsp)
  expect_equal(sum(is.na(mod$mcmc.betas)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
})