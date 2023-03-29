#context("test-jSDM_gaussian")
#== Without traits ======
#===== Simple species distribution model (SDM) ==========

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
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
mu.target <- as.matrix(X) %*% beta.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data = Y ,
                           site_formula = ~ x1 + x2,
                           site_data = X, 
                           beta_start=0,
                           mu_beta=0, V_beta=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with one species", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#======== Joint species distribution model (JSDM) ====================

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
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
mu.target <- as.matrix(X) %*% beta.target 
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data = Y ,
                           site_formula = ~ x1 + x2,
                           site_data = X, 
                           beta_start=0,
                           mu_beta=0, V_beta=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============= JSDM with latent variables ===============
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target 
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data = Y,
                           site_formula = ~ x1 + x2,
                           site_data = X, site_effect="none",
                           n_latent=n_latent,
                           beta_start=0,
                           lambda_start=0, W_start=0,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============== JSDM with fixed site effect =================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
alpha.target <- runif(nsite,-1,1)
alpha.target[1] <- 0 
mu.target <- as.matrix(X) %*% beta.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data=Y, n_latent=0,
                           site_formula = ~ x1 + x2,
                           site_data = X, site_effect="fixed",
                           alpha_start=0, beta_start=0,
                           V_alpha=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           mu_beta=0, V_beta=10,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with random site effect ================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data=Y, n_latent=0,
                           site_formula = ~ x1 + x2,
                           site_data = X, site_effect="random",
                           alpha_start=0, beta_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#======= JSDM with fixed site effect and latent variables ==============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Fixed site effect
alpha.target <- runif(nsite,-1,1)
alpha.target[1] <- 0 
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)



# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=~x1+x2,
                           site_data=X, n_latent=2, 
                           site_effect = "fixed", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with fixed site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============ JSDM with random site effect and latent variables ==================================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=~x1+x2,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#== JSDM with intercept only, random site effect and latent variables ===============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=~Int-1,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#== With traits ===========
#======== Joint species distribution model (JSDM) ====================

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
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
form.Tr <- function(trait_formula, trait_data,X){
  data <- trait_data
  # add column of 1 with names of covariables in site_data 
  data[,colnames(X)] <- 1
  mf.suit.tr <- model.frame(formula=trait_formula, data=data)
  # full design matrix corresponding to formula
  mod.mat <- model.matrix(attr(mf.suit.tr,"terms"), data=mf.suit.tr)
  # Remove duplicated columns to get design matrix for traits 
  Tr <- as.matrix(mod.mat[,!duplicated(mod.mat,MARGIN=2)])
  colnames(Tr) <- colnames(mod.mat)[!duplicated(mod.mat,MARGIN=2)]
  # Rename columns according to considered trait 
  for(p in 1:np){
    if(sum(colnames(Tr)==colnames(X)[p])==0){
      colnames(Tr) <- gsub(pattern=paste0(":",colnames(X)[p]), replacement="",
                           x=colnames(Tr), fixed=TRUE)
      colnames(Tr) <- gsub(pattern=paste0(colnames(X)[p],":"), replacement="",
                           x=colnames(Tr), fixed=TRUE)
    }
  }
  nt <- ncol(Tr)
  n_Tint <- sum(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
  col_Tint <- which(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
  gamma_zeros <- matrix(0,nt,np)
  rownames(gamma_zeros) <- colnames(Tr)
  colnames(gamma_zeros) <- colnames(X)
  for(t in 1:nt){
    for(p in 1:np){
      term <- c(grep(paste0(colnames(X)[p],":"), colnames(mod.mat), value=TRUE, fixed=TRUE),grep(paste0(":",colnames(X)[p]), colnames(mod.mat), value=TRUE, fixed=TRUE))
      if(length(term)==0) next
      # fixed=TRUE pattern is a string to be matched as is 
      # not a regular expression because of special characters in formula (^, /, [, ...) 
      gamma_zeros[t,p] <- length(c(grep(paste0(":",colnames(Tr)[t]), term, fixed=TRUE),grep(paste0(colnames(Tr)[t],":"), term, fixed=TRUE)))
    }
    gamma_zeros[t,1] <- length(which(colnames(mod.mat)==colnames(Tr)[t]))  
  }
  gamma_zeros[col_Tint,] <- 1
  return(list(gamma_zeros=gamma_zeros,Tr=Tr))
}
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect 
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
mu.target <- as.matrix(X) %*% beta.target 
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data = Y ,
                           site_formula = site_formula,
                           site_data = X,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           beta_start=0,
                           mu_beta=0, V_beta=10,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============= JSDM with latent variables ===============
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target 
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data = Y,
                           site_formula = site_formula,
                           site_data = X, site_effect="none",
                           n_latent=n_latent,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           beta_start=0,
                           lambda_start=0, W_start=0,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits, latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============== JSDM with fixed site effect =================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
alpha.target <- runif(nsite,-1,1)
alpha.target[1] <- 0 
mu.target <- as.matrix(X) %*% beta.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data=Y, n_latent=0,
                           site_formula = site_formula,
                           site_data = X, site_effect="fixed",
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           alpha_start=0, beta_start=0,
                           V_alpha=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           mu_beta=0, V_beta=10,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits, fixed site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#=============== JSDM with random site effect ================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data=Y, n_latent=0,
                           site_formula =site_formula,
                           site_data = X, site_effect="random",
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           alpha_start=0, beta_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits, random site effect", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]), c(nsamp,ncol(X)))
  expect_equal(length(mod$mcmc.gamma), ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#======= JSDM with fixed site effect and latent variables ==============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Fixed site effect
alpha.target <- runif(nsite,-1,1)
alpha.target[1] <- 0 
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(burnin=burnin, mcmc=mcmc, thin=thin,
                           response_data=Y,
                           site_formula=site_formula,
                           site_data=X, n_latent=2, 
                           site_effect = "fixed", 
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=10,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits, fixed site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#============ JSDM with random site effect and latent variables ==================================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=site_formula,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with traits, random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#== JSDM with intercept only in X, random site effect and latent variables ===============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + I(WSD^2) + I(SLA^2) 
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=~Int-1,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=10,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=10,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with  intercept only in X, random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#== JSDM with intercept only in Tr, random site effect and latent variables ===============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(Int=rep(1,nsp))
trait_formula <- ~. -1
# trait_formula <- ~ Int + x1:Int + x2:Int + I(x1^2):Int + I(x2^2):Int -1
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=site_formula,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=1,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=1,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with  intercept only in Tr, traits, random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(as.matrix(sapply(mod$mcmc.gamma,colMeans))!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})

#== JSDM with intercept only in Tr and X, random site effect and latent variables ===============================
#= Number of latent variables
n_latent <- 2
# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
np <- ncol(X)
# Species traits dataset
trait_data <- data.frame(Int=rep(1,nsp))
trait_formula <- ~. -1
# trait_formula <- ~ Int + x1:Int + x2:Int + I(x1^2):Int + I(x2^2):Int -1
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
# Interaction traits environment parameters
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# Species effect
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
# Latent variables 
W <- matrix(rnorm(nsite*n_latent,0,1), byrow=TRUE, nrow=nsite)
#= Species effect 
beta.target <- t(matrix(runif(nsp*ncol(X),-1,1), byrow=TRUE, nrow=nsp))
# Factor loadings  
mat <- t(matrix(runif(nsp*n_latent,-1,1), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# Random site effect
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
mu.target <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
# Variance of residuals
V.target <- 0.2
Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_gaussian(response_data=Y,
                           site_formula=~Int-1,
                           site_data=X, n_latent=2, 
                           site_effect = "random", 
                           burnin=burnin, mcmc=mcmc, thin=thin,
                           trait_formula = trait_formula,
                           trait_data = trait_data,
                           gamma_start=0,
                           mu_gamma=0, V_gamma=1,
                           alpha_start=0, beta_start=0,
                           lambda_start=0, W_start=0,
                           V_alpha=1,
                           V_start=1,
                           shape_V=0.5, rate_V=0.0005,
                           shape_Valpha=0.5, rate_Valpha=0.0005,
                           mu_beta=0, V_beta=1,
                           mu_lambda=0, V_lambda=1,
                           seed=1234, verbose=0)
# Tests
test_that("jSDM_gaussian works with intercept only in Tr and X, random site effect and latent variables", {
  expect_equal(length(mod$mcmc.sp),nsp)
  expect_equal(dim(mod$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod$mcmc.gamma),ncol(X))
  expect_equal(dim(mod$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(as.matrix(sapply(mod$mcmc.gamma,colMeans))!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(dim(mod$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod$Y_pred)),0)
  expect_equal(sum(is.infinite(mod$Y_pred)),0)
  expect_equal(dim(mod$Y_pred),c(nsite,nsp))
  expect_equal(sum(is.na(mod$mcmc.V)),0)
  expect_equal(dim(mod$mcmc.V),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.V_alpha)),0)
  expect_equal(dim(mod$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod$mcmc.Deviance)),0)
  expect_equal(dim(mod$mcmc.Deviance),c(nsamp,1))
})