#context("test-jSDM_binomial_probit_sp_constrained")
#== Without traits ======
#======== Joint species distribution model (JSDM) ====================
#============= JSDM with latent variables ===============
# Data simulation
#= Number of sites
nsite <- 10
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 5 

#= Number of latent variables
n_latent <- 2

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- matrix(rnorm(nsite*n_latent,0,1),nsite)
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target 
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(burnin=burnin, mcmc=mcmc, thin=thin,
                                                 presence_data = Y,
                                                 site_formula = ~ x1 + x2,
                                                 site_data = X, site_effect="none",
                                                 n_latent=n_latent,
                                                 beta_start=0,
                                                 lambda_start=0, W_start=0,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#======= JSDM with fixed site effect and latent variables ==============================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
                                                 site_formula=~x1+x2,
                                                 site_data=X, n_latent=2, 
                                                 site_effect = "fixed", 
                                                 burnin=burnin, mcmc=mcmc, thin=thin,
                                                 alpha_start=0, beta_start=0,
                                                 lambda_start=0, W_start=0,
                                                 V_alpha=10,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with fixed site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#============ JSDM with random site effect and latent variables ==================================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
                                                 site_formula=~x1+x2,
                                                 site_data=X, n_latent=2, 
                                                 site_effect = "random", 
                                                 burnin=burnin, mcmc=mcmc, thin=thin,
                                                 alpha_start=0, beta_start=0,
                                                 lambda_start=0, W_start=0,
                                                 V_alpha=1,
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#== JSDM with intercept only, random site effect and latent variables ===============================

# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
beta.target <- t(matrix(runif(nsp*ncol(X),-2,2), byrow=TRUE, nrow=nsp))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
                                                 site_formula=~Int-1,
                                                 site_data=X, n_latent=2, 
                                                 site_effect = "random", 
                                                 burnin=burnin, mcmc=mcmc, thin=thin,
                                                 alpha_start=0, beta_start=0,
                                                 lambda_start=0, W_start=0,
                                                 V_alpha=1,
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#== With traits ===========
#======== Joint species distribution model (JSDM) ====================
form.Tr <- function(trait_formula, trait_data,X){
  data <- trait_data
  # add column of 1 with names of covariates in site_data 
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

#============= JSDM with latent variables ===============
# Data simulation
#= Number of sites
nsite <- 10
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 5 

#= Number of latent variables
n_latent <- 2

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- matrix(rnorm(nsite*n_latent,0,1),nsite)
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target 
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(burnin=burnin, mcmc=mcmc, thin=thin,
                                                 presence_data = Y,
                                                 site_formula = site_formula,
                                                 site_data = X, site_effect="none",
                                                 n_latent=n_latent,
                                                 trait_formula = trait_formula,
                                                 trait_data = trait_data,
                                                 gamma_start=0,
                                                 mu_gamma=0, V_gamma=10,
                                                 beta_start=0,
                                                 lambda_start=0, W_start=0,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with traits, latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod[[1]]$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#============== JSDM with fixed site effect =================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- as.matrix(X) %*% beta.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}



#======= JSDM with fixed site effect and latent variables ==============================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
alpha.target <- runif(nsite,-2,2)
alpha.target[1] <- 0 
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(burnin=burnin, mcmc=mcmc, thin=thin,
                                                 presence_data=Y,
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
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with traits, fixed site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod[[1]]$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#============ JSDM with random site effect and latent variables ==================================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
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
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with traits, random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod[[1]]$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#== JSDM with intercept only in X, random site effect and latent variables ===============================

# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + I(WSD^2) + I(SLA^2) 
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
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
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=10,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with  intercept only in X, random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(sapply(mod[[1]]$mcmc.gamma,colMeans)!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#== JSDM with intercept only in Tr, random site effect and latent variables ===============================

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
site_data <- data.frame(x1=x1,x2=x2)
site_formula <- ~ x1 + x2 + I(x1^2) + I(x2^2)
X <- model.matrix(site_formula, site_data)
np <- ncol(X)
trait_data <- data.frame(Int=rep(1,nsp))
trait_formula <- ~. -1
# trait_formula <- ~ Int + x1:Int + x2:Int + I(x1^2):Int + I(x2^2):Int -1
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
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
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=1,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with  intercept only in Tr, traits, random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(as.matrix(sapply(mod[[1]]$mcmc.gamma,colMeans))!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})

#== JSDM with intercept only in Tr and X, random site effect and latent variables ===============================
# Ecological process (suitability)
X <- data.frame(Int=rep(1,nsite))
np <- ncol(X)
trait_data <- data.frame(Int=rep(1,nsp))
trait_formula <- ~. -1
# trait_formula <- ~ Int + x1:Int + x2:Int + I(x1^2):Int + I(x2^2):Int -1
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
l.zero <- 0
l.diag <- runif(2,0,2)
l.other <- runif(nsp*n_latent-3,-2,2)
lambda.target <- t(matrix(c(l.diag[1],l.zero,
                            l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))
probit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}


# Fit the model
burnin <- 500
mcmc <- 500
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_sp_constrained(presence_data=Y,
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
                                                 shape_Valpha=0.5, rate_Valpha=0.0005,
                                                 mu_beta=0, V_beta=1,
                                                 mu_lambda=0, V_lambda=1,
                                                 verbose=0)
# Tests
test_that("jSDM_binomial_probit_sp_constrained works with intercept only in Tr and X, random site effect and latent variables", {
  expect_equal(length(mod[[1]]$mcmc.sp),nsp)
  expect_equal(dim(mod[[1]]$mcmc.sp[["sp_1"]]),c(nsamp,ncol(X)+n_latent))
  expect_equal(length(mod[[1]]$mcmc.gamma),ncol(X))
  expect_equal(dim(mod[[1]]$mcmc.gamma[[1]]),c(nsamp,ncol(Tr)))
  expect_equal(which(as.matrix(sapply(mod[[1]]$mcmc.gamma,colMeans))!=0), which(gamma_zeros!=0))
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_1),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_1)),0)
  expect_equal(dim(mod[[1]]$mcmc.latent$lv_2),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.latent$lv_2)),0)
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.alpha),c(nsamp,nsite))
  expect_equal(sum(is.na(mod[[1]]$mcmc.alpha)),0)
  expect_equal(sum(is.na(mod[[1]]$Z_latent)),0)
  expect_equal(sum(is.infinite(mod[[1]]$Z_latent)),0)
  expect_equal(dim(mod[[1]]$Z_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$probit_theta_latent)),0)
  expect_equal(dim(mod[[1]]$probit_theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$theta_latent)),0)
  expect_equal(dim(mod[[1]]$theta_latent),c(nsite,nsp))
  expect_equal(sum(is.na(mod[[1]]$mcmc.V_alpha)),0)
  expect_equal(dim(mod[[1]]$mcmc.V_alpha),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$mcmc.Deviance)),0)
  expect_equal(dim(mod[[1]]$mcmc.Deviance),c(nsamp,1))
  expect_equal(sum(is.na(mod[[1]]$sp_constrained)),0)
  expect_equal(length(mod[[1]]$sp_constrained),n_latent)
})
