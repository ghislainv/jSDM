context("test-jSDM_binomial_probit_long_format")

#=== Without traits =======

##============== Simple species distribution model (SDM) =======


# Data simulation
#= Number of sites
nsite <- 50
#= Set seed for repeatability
seed <- 1234
set.seed(seed)
#= Number of species
nsp <- 1

#= Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Int","x1","x2")
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
probit_theta <- X %*% beta.target 
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
colnames(X_supObs) <- c("Int","x1","x2")
probit_theta_supObs <- X_supObs%*%beta.target 
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp) 
Id_sp <- rep(1:nsp,each=nsite) 
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data,
                                              site_formula = ~ (x1 + x2):species,
                                              beta_start=0,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with one species", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

##======joint species distribution model (jSDM)=====

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
X <- cbind(rep(1,nsite), x1, x2)
colnames(X) <- c("Int","x1","x2")
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))

probit_theta <- X %*% beta.target 
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
colnames(X_supObs) <- c("Int","x1","x2")
probit_theta_supObs <- X_supObs%*%beta.target 
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1 
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data = data,
                                              site_formula = ~ (x1 + x2):species,
                                              beta_start=0,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

##===== JSDM with latent variables ======


# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite), x1, x2)
colnames(X) <- c("Int","x1","x2")
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


probit_theta <- X %*% beta.target + W %*% lambda.target 
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
probit_theta_supObs <- X_supObs%*%beta.target + W%*%lambda.target 
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                        data=data,
                                        site_formula = ~ x1:species + x2:species + species,
                                        site_effect="none",
                                        n_latent=n_latent,
                                        beta_start=0,
                                        lambda_start=0, W_start=0,
                                        mu_beta=0, V_beta=100,
                                        mu_lambda=0, V_lambda=10,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=========JSDM with fixed site effect=============

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite), x1, x2)
colnames(X) <- c("Int","x1","x2")
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))

alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- X %*% beta.target + alpha.target
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
probit_theta_supObs <- X_supObs%*%beta.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data, n_latent=0,
                                              site_formula = ~(x1 + x2):species,
                                              site_effect="fixed",
                                              alpha_start=0, beta_start=0,
                                              V_alpha=10,
                                              shape=0.5, rate=0.0005,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============JSDM with random site effect==================

#= Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Int","x1","x2")
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
V_alpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
probit_theta <- X %*% beta.target + alpha.target
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
probit_theta_supObs <- X_supObs%*%beta.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp) 
Id_sp <- rep(1:nsp,each=nsite) 
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data, n_latent=0,
                                              site_formula = ~ (x1 + x2):species,
                                              site_effect="random",
                                              alpha_start=0, beta_start=0,
                                              V_alpha=1,
                                              shape=0.5, rate=0.0005,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============JSDM with fixed site effect and latent variables===============

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Int","x1","x2")
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- X %*% beta.target + W %*% lambda.target + alpha.target
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
probit_theta_supObs <- X_supObs%*%beta.target + W%*%lambda.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ (x1 + x2):species,
                                              n_latent=2, 
                                              site_effect = "fixed", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=10,
                                              shape=0.5, rate=0.0005,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with fixed site effect and latent variables", {
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
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with random site effect and latent variables =============

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Int","x1","x2")
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))

Valpha.target <- 0.5

alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- X %*% beta.target + W %*% lambda.target + alpha.target
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
set.seed(2*seed)
X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
probit_theta_supObs <- X_supObs%*%beta.target + W%*%lambda.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
                   x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ (x1 + x2):species,
                                              n_latent=2, 
                                              site_effect = "random", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=1,
                                              shape=0.5, rate=0.0005,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect and latent variables", {
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
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with random intercept only, site effect and latent variables =============

# Ecological process (suitability)
X <- matrix(1,nsite,1)
colnames(X) <- c("Int")
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
data <- cbind (X,W)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))

Valpha.target <- 0.5

alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- X %*% beta.target + W %*% lambda.target + alpha.target
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=Id_site, species=Id_sp, Y=Y, Int=1)
# missing observation
data <- data[-1,]
nobs <- nobs -1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ Int:species - species ,
                                              n_latent=2, 
                                              site_effect = "random", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=1,
                                              shape=0.5, rate=0.0005,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect and latent variables", {
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
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

context("test-jSDM_binomial_probit_long_format")


#=== With traits =======

##============== Simple species distribution model (SDM) =======


# Data simulation
#= Number of sites
nsite <- 50
#= Set seed for repeatability
seed <- 1234
set.seed(seed)
#= Number of species
nsp <- 1

# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
gamma.target <-runif(nd,-1,1)
# constraint of identifiability
beta.target[,1] <- 0.0
#= Variance of random site effect
V_alpha.target <- 0.5
#= Random site effect
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
## probit_theta
probit_theta <- c(X %*% beta.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data,
                                              site_formula = ~ x1.SLA + x1:species + x1.2:species,
                                              gamma_start=0,
                                              mu_gamma=0, V_gamma=100,
                                              beta_start=0,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with one species", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

##======joint species distribution model (jSDM)=====

# Data simulation
#= Number of sites
nsite <- 50
#= Number of species
nsp <- 5
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
gamma.target <-runif(nd,-1,1)
## probit_theta
probit_theta <- c(X %*% beta.target) + as.matrix(D) %*% gamma.target 
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + as.matrix(D_supObs) %*% gamma.target 
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1 
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data = data,
                                              site_formula = ~ (x1 + x1.2):species + x1.SLA,
                                              beta_start=0,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

##===== JSDM with latent variables ======
#= Number of latent variables
n_latent <- 2
#'
# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## W
W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame( x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
gamma.target <-runif(nd,-1,1)
# constraints of identifiability
beta.target[,1] <- 0.0
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
## probit_theta
probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) + as.matrix(D) %*% gamma.target 
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + as.matrix(D_supObs) %*% gamma.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                        data=data,
                                        site_formula = ~ (x1 + x1.2):species + x1.SLA,
                                        site_effect="none",
                                        n_latent=n_latent,
                                        gamma_start=0,
                                        mu_gamma=0, V_gamma=100,
                                        beta_start=0,
                                        lambda_start=0, W_start=0,
                                        mu_beta=0, V_beta=100,
                                        mu_lambda=0, V_lambda=10,
                                        seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=========JSDM with fixed site effect=============

# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
gamma.target <-runif(nd,-1,1)
#= Fixed site effect
alpha.target <- runif(nsite,-2,2)
# constraints of identifiability
beta.target[,1] <- 0.0
alpha.target[1] <- 0 
## probit_theta
probit_theta <- c(X %*% beta.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data, n_latent=0,
                                              site_formula = ~(x1 + x1.2):species + x1.SLA,
                                              site_effect="fixed",
                                              alpha_start=0, gamma_start=0,
                                              beta_start=0,
                                              V_alpha=10,
                                              shape=0.5, rate=0.0005,
                                              mu_gamma=0, V_gamma=100,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============JSDM with random site effect==================

# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
gamma.target <-runif(nd,-1,1)
#= Variance of random site effect
V_alpha.target <- 0.5
#= Random site effect
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
# constraints of identifiability
beta.target[,1] <- 0.0
## probit_theta
probit_theta <- c(X %*% beta.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target)+ as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(burnin=burnin, mcmc=mcmc, thin=thin,
                                              data=data, n_latent=0,
                                              site_formula = ~ (x1 + x1.2):species + x1.SLA,
                                              site_effect="random",
                                              alpha_start=0, gamma_start=0,
                                              V_alpha=1, beta_start=0,
                                              shape=0.5, rate=0.0005,
                                              mu_gamma=0, V_gamma=100,
                                              mu_beta=0, V_beta=100,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============JSDM with fixed site effect and latent variables===============
#= Number of latent variables
n_latent <- 2
#'
# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## W
W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
gamma.target <-runif(nd,-1,1)
#= Fixed site effect
alpha.target <- runif(nsite,-2,2)
# constraints of identifiability
beta.target[,1] <- 0.0
alpha.target[1] <- 0 
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
## probit_theta
probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ (x1 + x1.2):species + x1.SLA,
                                              n_latent=2, 
                                              site_effect = "fixed", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, gamma_start=0,
                                              beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=10,
                                              shape=0.5, rate=0.0005,
                                              mu_gamma=0, V_gamma=100,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with fixed site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with random site effect and latent variables =============

#= Number of latent variables
n_latent <- 2
#'
# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## W
W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
gamma.target <-runif(nd,-1,1)
# constraints of identifiability
beta.target[,1] <- 0.0
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
#= Variance of random site effect
V_alpha.target <- 0.5
#= Random site effect
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
## probit_theta
probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ (x1 + x1.2):species + x1.SLA,
                                              n_latent=2, 
                                              site_effect = "random", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, gamma_start=0,
                                              beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=1,
                                              shape=0.5, rate=0.0005,
                                              mu_gamma=0, V_gamma=100,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with intercept only, random site effect and latent variables =============

#= Number of latent variables
n_latent <- 2
#'
# Ecological process (suitability)
## X
X <- matrix(1,nsite,1)
colnames(X) <- c("Int")
np <- ncol(X)
## W
W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
## D
SLA <- runif(nsp,-1,1)
x1 <- rnorm(nsite,0,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
gamma.target <-runif(nd,-1,1)
# constraints of identifiability
beta.target[,1] <- 0.0
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
#= Variance of random site effect
V_alpha.target <- 0.5
#= Random site effect
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
## probit_theta
probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
X_supObs <- X
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   Int=1, x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1

# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_long_format(data=data,
                                              site_formula= ~ Int:species + x1.SLA -species,
                                              n_latent=2, 
                                              site_effect = "random", 
                                              burnin=burnin, mcmc=mcmc, thin=thin,
                                              alpha_start=0, gamma_start=0,
                                              beta_start=0,
                                              lambda_start=0, W_start=0,
                                              V_alpha=1,
                                              shape=0.5, rate=0.0005,
                                              mu_gamma=0, V_gamma=100,
                                              mu_beta=0, V_beta=100,
                                              mu_lambda=0, V_lambda=10,
                                              seed=1234, verbose=1)
# Tests
test_that("jSDM_binomial_probit_long_format works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.gamma),c(nsamp,ncol(D)))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Z_latent)),0)
  expect_equal(sum(is.infinite(mod$Z_latent)),0)
  expect_equal(length(mod$Z_latent),nobs)
  expect_equal(sum(is.na(mod$theta_latent)),0)
  expect_equal(length(mod$theta_latent),nobs)
  expect_equal(sum(is.na(mod$probit_theta_latent)),0)
  expect_equal(length(mod$probit_theta_latent),nobs)
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})
