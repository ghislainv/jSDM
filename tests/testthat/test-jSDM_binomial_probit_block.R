context("test-jSDM_binomial_probit_block")

#=========================================
#= Simple species distribution model (SDM) 

# Data simulation
#= Number of sites
nsite <- 50
#= Set seed for repeatability
seed <- 1234
set.seed(seed)
#= Number of species
nsp <- 1

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
V <- 1
probit_theta <- as.matrix(X) %*% beta.target 
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}
# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                        presence_site_sp = Y ,
                                        site_suitability = ~ x1 + x2,
                                        site_data = X, 
                                        beta_start=0,
                                        mu_beta=0, V_beta=1.0E6,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with one species", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=========================================
#= joint species distribution model (jSDM)

# Data simulation
#= Number of sites
nsite <- 50
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 5

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
V <- 1
probit_theta <- as.matrix(X) %*% beta.target 
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                        presence_site_sp = Y ,
                                        site_suitability = ~ x1 + x2,
                                        site_data = X, 
                                        beta_start=0,
                                        mu_beta=0, V_beta=1.0E6,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============================
#= jSDM with latent variables

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <- ncol(W)
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
param.target <- rbind(beta.target,lambda.target)
V <- 1
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target 
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                  presence_site_sp = Y,
                                  site_suitability = ~ x1 + x2,
                                  site_data = X, site_effect="none",
                                  n_latent=n_latent,
                                  beta_start=0,
                                  lambda_start=0, W_start=0,
                                  mu_beta=0, V_beta=1.0E6,
                                  mu_lambda=0, V_lambda=10,
                                  seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#===============================
#= jSDM with fixed site effect 

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
V <- 1
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- as.matrix(X) %*% beta.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                        presence_site_sp=Y, n_latent=0,
                                        site_suitability = ~ x1 + x2,
                                        site_data = X, site_effect="fixed",
                                        alpha_start=0, beta_start=0,
                                        V_alpha=10,
                                        shape=0.5, rate=0.0005,
                                        mu_beta=0, V_beta=1.0E6,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#===============================
#= jSDM with random site effect 

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
V <- 1
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                        presence_site_sp=Y, n_latent=0,
                                        site_suitability = ~ x1 + x2,
                                        site_data = X, site_effect="random",
                                        alpha_start=0, beta_start=0,
                                        V_alpha=1,
                                        shape=0.5, rate=0.0005,
                                        mu_beta=0, V_beta=1.0E6,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#====================================================
#= jSDM with fixed site effect and latent variables 

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
param.target <- rbind(beta.target,lambda.target)
V <- 1
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(presence_site_sp=Y,
                                        site_suitability=~x1+x2,
                                        site_data=X, n_latent=2, 
                                        site_effect = "fixed", 
                                        burnin=burnin, mcmc=mcmc, thin=thin,
                                        alpha_start=0, beta_start=0,
                                        lambda_start=0, W_start=0,
                                        V_alpha=10,
                                        shape=0.5, rate=0.0005,
                                        mu_beta=0, V_beta=1.0E6,
                                        mu_lambda=0, V_lambda=10,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with fixed site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#====================================================
#= jSDM with random site effect and latent variables 

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
param.target <- rbind(beta.target,lambda.target)
Valpha.target <- 0.5
V <- 1
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(presence_site_sp=Y,
                                        site_suitability=~x1+x2,
                                        site_data=X, n_latent=2, 
                                        site_effect = "random", 
                                        burnin=burnin, mcmc=mcmc, thin=thin,
                                        alpha_start=0, beta_start=0,
                                        lambda_start=0, W_start=0,
                                        V_alpha=1,
                                        shape=0.5, rate=0.0005,
                                        mu_beta=0, V_beta=1.0E6,
                                        mu_lambda=0, V_lambda=10,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_block works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(dim(mod$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod$probit_theta_pred)),0)
  expect_equal(dim(mod$probit_theta_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})