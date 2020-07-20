context("test-jSDM_binomial_logit_rand_site")

#==================
#== Data simulation

#= Number of species
nsp <- 50
#= Number of sites
nsite <- 100
#= Number of latent variables
n_latent <- 0
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of visits associated to each site
visits<- rpois(nsite,3)
visits[visits==0] <- 1

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
logit.theta <- X %*% t(beta.target) + alpha.target
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


#= Site-occupancy model
burnin <- 100
mcmc <- 100
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_logit_rand_site(# Chains
  burnin, mcmc, thin,
  # Response variable
  presence_site_sp=Y, trials=visits,
  # Explanatory variables
  site_suitability=~x1+x2, site_data=X,
  # Starting values
  alpha_start=0,
  beta_start=0,
  V_alpha_start=1,
  # Priors
  shape=0.5, rate=0.0005,
  mu_beta=0, V_beta=1.0E6,
  # Various
  seed=1234, ropt=0.44, verbose=1)

test_that("jSDM_binomial_logit_rand_site works", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,ncol(X)+n_latent))
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
})