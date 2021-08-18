context("test-jSDM_binomial_logit")

#================ Single species distribution model (SDM) =======================

# Data simulation
#= Number of sites
nsite <- 50
#= Number of species
nsp <- 1
#= Set seed for repeatability
seed <- 1234
#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

# Ecological process (suitability)
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
mod <- jSDM_binomial_logit(burnin, mcmc, thin,# Chains
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           # Starting values
                           beta_start=0,
                           # Priors
                           mu_beta=0, V_beta=1.0E6,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)

test_that("jSDM_binomial_logit works with one species", {
  expect_equal(sum(is.na(mod$theta_latent)),0)
  #' theta_latent \tab Predictive posterior mean of the probability associated to the suitability process for each observation. \cr
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,np))
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
})

#=============== Joint species distribution model (JSDM) =====================

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

# Ecological process (suitability)
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

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_logit(burnin, mcmc, thin,# Chains
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           # Starting values
                           beta_start=0,
                           # Priors
                           mu_beta=0, V_beta=1.0E6,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)
# Tests
test_that("jSDM_binomial_logit works ", {
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,np))
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
})

#============= JSDM with fixed site effect =================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0
logit.theta <- X %*% t(beta.target) + alpha.target
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model 
mod <- jSDM_binomial_logit(burnin, mcmc, thin, # Chains
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           site_effect="fixed", 
                           # Starting values
                           alpha_start=0,
                           beta_start=0,
                           # Priors
                           V_alpha=10,
                           mu_beta=0, V_beta=1.0E6,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)

# Tests
test_that("jSDM_binomial_logit works with fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
})
#========== JSDM with random site effect ====================

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


# Fit the model 
mod <- jSDM_binomial_logit(burnin, mcmc, thin, # Chains
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           site_effect="random", 
                           # Starting values
                           alpha_start=0,
                           beta_start=0,
                           V_alpha=1,
                           # Priors
                           shape=0.5, rate=0.0005,
                           mu_beta=0, V_beta=1.0E6,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)

# Tests
test_that("jSDM_binomial_logit works with random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
})

#============ JSDM with latent variables ==================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
set.seed(2*seed)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <-  ncol(W)
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*2-3,-2,2)
lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
logit.theta <- X %*% t(beta.target) + W %*% t(lambda.target)
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model
mod <- jSDM_binomial_logit(burnin, mcmc, thin,# Chains
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           n_latent=n_latent,
                           # Starting values
                           beta_start=0, lambda_start = 0,
                           W_start=0,
                           # Priors
                           mu_beta=0, V_beta=1.0E6,
                           mu_lambda=0, V_lambda=10,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)
# Tests
test_that("jSDM_binomial_logit works with latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,ncol(X)+n_latent))
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
})

#========= JSDM with latent variables and fixed site effect ===================================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
set.seed(2*seed)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <- ncol(W)
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*2-3,-2,2)
lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0
logit.theta <- X %*% t(beta.target) + W %*% t(lambda.target) + alpha.target
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model 
mod <- jSDM_binomial_logit(burnin, mcmc, thin, # Chains 
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           n_latent= n_latent,
                           site_effect ="fixed",
                           # Starting values
                           alpha_start=0,
                           beta_start=0,
                           lambda_start = 0,
                           W_start=0,
                           # Priors
                           V_alpha=10,
                           mu_beta=0, V_beta=1.0E6,
                           mu_lambda=0, V_lambda=10,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)
# Tests 
test_that("jSDM_binomial_logit works with fixed site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
  expect_equal(unique(lapply(mod$mcmc.sp,dim))[[1]],c(nsamp,ncol(X)+n_latent))
  expect_equal(sum(is.na(mod$mcmc.sp)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
})

#=========== JSDM with latent variables and random site effect ========================================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
np <- ncol(X)
set.seed(2*seed)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <- ncol(W)
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*2-3,-2,2)
lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
logit.theta <- X %*% t(beta.target) + W %*% t(lambda.target) + alpha.target
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model 
mod <- jSDM_binomial_logit(burnin, mcmc, thin, # Chains 
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~x1+x2, site_data=X,
                           n_latent= n_latent,
                           site_effect ="random",
                           # Starting values
                           alpha_start=0,
                           beta_start=0,
                           lambda_start = 0,
                           W_start=0,
                           V_alpha=1,
                           # Priors
                           shape=0.5, rate=0.0005,
                           mu_beta=0, V_beta=1.0E6,
                           mu_lambda=0, V_lambda=10,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)
# Tests 
test_that("jSDM_binomial_logit works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
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

#=== JSDM with intercept only, latent variables and random site effect =================================

# Ecological process (suitability)
X <- matrix(1,nsite,1)
colnames(X) <- "Int"
np <- ncol(X)
set.seed(2*seed)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <- ncol(W)
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*2-3,-2,2)
lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp)
beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
logit.theta <- X %*% t(beta.target) + W %*% t(lambda.target) + alpha.target
theta <- inv_logit(logit.theta)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model 
mod <- jSDM_binomial_logit(burnin, mcmc, thin, # Chains 
                           # Response variable
                           presence_data=Y, trials=visits,
                           # Explanatory variables
                           site_formula=~Int-1, site_data=X,
                           n_latent= n_latent,
                           site_effect ="random",
                           # Starting values
                           alpha_start=0,
                           beta_start=0,
                           lambda_start = 0,
                           W_start=0,
                           V_alpha=1,
                           # Priors
                           shape=0.5, rate=0.0005,
                           mu_beta=0, V_beta=1.0E6,
                           mu_lambda=0, V_lambda=10,
                           # Various
                           seed=1234, ropt=0.44, verbose=1)
# Tests 
test_that("jSDM_binomial_logit works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(dim(mod$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$logit_theta_latent)),0)
  expect_equal(dim(mod$logit_theta_latent),c(nsite,nsp))
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