#context("test-predict.jSDM")


# Data simulation
#= Number of sites
nsite <- 50
#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 5

#= Number of latent variables 
n_latent <- 2

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
#### Binomial logit ###########

#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

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
logit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
theta <- jSDM::inv_logit(logit_theta)
set.seed(seed)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)


# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_logit(presence_data=Y, trials=visits,
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
                                 shape=0.5, rate=0.0005,
                                 mu_beta=0, V_beta=10,
                                 mu_lambda=0, V_lambda=1,
                                 seed=1234, verbose=1)

# Sites and species concerned by predictions :
## half of sites 
npred <- round(nrow(Y)/2)
Id_sites <- sample.int(nrow(Y), npred)
## All species 
Id_species <- colnames(mod$model_spec$presence_data)
# Simulate new observations of covariates on those sites 
np <- ncol(mod$model_spec$site_data)
simdata <- matrix(1, nrow=npred, ncol = np)
colnames(simdata) <- colnames(mod$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
# Predictions 
theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                      Id_sites=Id_sites, type="mean")
# Tests
test_that("predict.jSDM works with jSDM_binomial_logit results considering intercept only in Tr and X, random site effect and latent variables.", {
  expect_equal(dim(theta_pred),c(npred,nsp))
  expect_equal(sum(is.na(theta_pred)),0)
  expect_equal(sum(theta_pred<0),0)
  expect_equal(sum(theta_pred>1),0)
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
logit_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
theta <- jSDM::inv_logit(logit_theta)
set.seed(seed)
Y <- apply(theta, 2, rbinom, n=nsite, size=visits)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_logit(presence_data=Y, trials=visits,
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
                                 shape=0.5, rate=0.0005,
                                 mu_beta=0, V_beta=1,
                                 mu_lambda=0, V_lambda=1,
                                 seed=1234, verbose=1)

# Sites and species concerned by predictions :
## half of sites 
npred <- round(nrow(Y)/2)
Id_sites <- sample.int(nrow(Y), npred)
## All species 
Id_species <- colnames(mod$model_spec$presence_data)
# Simulate new observations of covariates on those sites 
np <- ncol(mod$model_spec$site_data)
simdata <- matrix(rnorm(npred*np), nrow=npred, ncol = np)
colnames(simdata) <- colnames(mod$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
# Predictions 
theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                      Id_sites=Id_sites, type="mean")
# Tests
test_that("predict.jSDM works with jSDM_binomial_logit results considering traits, random site effect and latent variables", {
  expect_equal(dim(theta_pred),c(npred,nsp))
  expect_equal(sum(is.na(theta_pred)),0)
  expect_equal(sum(theta_pred<0),0)
  expect_equal(sum(theta_pred>1),0)
})
#### Binomial probit ##########

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
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit(presence_data=Y,
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
                                  shape=0.5, rate=0.0005,
                                  mu_beta=0, V_beta=10,
                                  mu_lambda=0, V_lambda=1,
                                  seed=1234, verbose=1)

# Sites and species concerned by predictions :
## half of sites 
npred <- round(nrow(Y)/2)
Id_sites <- sample.int(nrow(Y), npred)
## All species 
Id_species <- colnames(mod$model_spec$presence_data)
# Simulate new observations of covariates on those sites 
np <- ncol(mod$model_spec$site_data)
simdata <- matrix(rnorm(npred*np), nrow=npred, ncol = np)
colnames(simdata) <- colnames(mod$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
# Predictions 
theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                      Id_sites=Id_sites, type="mean")
# Tests
test_that("predict.jSDM works with jSDM_binomial_probit results considering traits, random site effect and latent variables", {
  expect_equal(dim(theta_pred),c(npred,nsp))
  expect_equal(sum(is.na(theta_pred)),0)
  expect_equal(sum(theta_pred<0),0)
  expect_equal(sum(theta_pred>1),0)
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
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit(presence_data=Y,
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
                                  shape=0.5, rate=0.0005,
                                  mu_beta=0, V_beta=1,
                                  mu_lambda=0, V_lambda=1,
                                  seed=1234, verbose=1)
# Tests
  # Sites and species concerned by predictions :
  ## half of sites 
  npred <- round(nrow(Y)/2)
  Id_sites <- sample.int(nrow(Y), npred)
  ## All species 
  Id_species <- colnames(mod$model_spec$presence_data)
  # Simulate new observations of covariates on those sites 
  np <- ncol(mod$model_spec$site_data)
  simdata <- matrix(1, nrow=npred, ncol = np)
  colnames(simdata) <- colnames(mod$model_spec$site_data)
  rownames(simdata) <- Id_sites
  simdata <- as.data.frame(simdata)
  # Predictions 
  theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                        Id_sites=Id_sites, type="mean")
  # Tests
  test_that("predict.jSDM works with jSDM_binomial_probit results considering intercept only in Tr and X, random site effect and latent variables", {
    expect_equal(dim(theta_pred),c(npred,nsp))
    expect_equal(sum(is.na(theta_pred)),0)
    expect_equal(sum(theta_pred<0),0)
    expect_equal(sum(theta_pred>1),0)
  })
  
#### Poisson log ##############
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
log_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
theta <- exp(log_theta)
set.seed(seed)
Y <- apply(theta, 2, rpois, n=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_poisson_log(count_data=Y,
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
                              shape=0.5, rate=0.0005,
                              mu_beta=0, V_beta=1,
                              mu_lambda=0, V_lambda=1,
                              seed=1234, verbose=1)
# Sites and species concerned by predictions :
## half of sites 
npred <- round(nrow(Y)/2)
Id_sites <- sample.int(nrow(Y), npred)
## All species 
Id_species <- colnames(mod$model_spec$count_data)
# Simulate new observations of covariates on those sites 
np <- ncol(mod$model_spec$site_data)
simdata <- matrix(rnorm(npred*np), nrow=npred, ncol = np)
colnames(simdata) <- colnames(mod$model_spec$site_data)
rownames(simdata) <- Id_sites
simdata <- as.data.frame(simdata)
# Predictions 
theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                      Id_sites=Id_sites, type="mean")
# Tests
test_that("predict.jSDM works with jSDM_poisson_log results considering traits, random site effect and latent variables", {
  expect_equal(dim(theta_pred),c(npred,nsp))
  expect_equal(sum(is.na(theta_pred)),0)
  expect_equal(sum(theta_pred<0),0)
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
log_theta <- as.matrix(X) %*% beta.target + W %*% lambda.target + alpha.target
theta <- exp(log_theta)
set.seed(seed)
Y <- apply(theta, 2, rpois, n=nsite)

# Fit the model
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_poisson_log(count_data=Y,
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
                              shape=0.5, rate=0.0005,
                              mu_beta=0, V_beta=1,
                              mu_lambda=0, V_lambda=1,
                              seed=1234, verbose=1)
# Tests
  # Sites and species concerned by predictions :
  ## half of sites 
  npred <- round(nrow(Y)/2)
  Id_sites <- sample.int(nrow(Y), npred)
  ## All species 
  Id_species <- colnames(mod$model_spec$count_data)
  # Simulate new observations of covariates on those sites 
  np <- ncol(mod$model_spec$site_data)
  simdata <- matrix(1, nrow=npred, ncol = np)
  colnames(simdata) <- colnames(mod$model_spec$site_data)
  rownames(simdata) <- Id_sites
  simdata <- as.data.frame(simdata)
  # Predictions 
  theta_pred <- predict(mod, newdata=simdata, Id_species=Id_species,
                        Id_sites=Id_sites, type="mean")
  # Tests
  test_that("predict.jSDM works with jSDM_poisson_log results considering intercept only in Tr and X, random site effect and latent variables.", {
    expect_equal(dim(theta_pred),c(npred,nsp))
    expect_equal(sum(is.na(theta_pred)),0)
    expect_equal(sum(theta_pred<0),0)
  })

########## Binomial probit long format ###############

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

# Sites and species concerned by predictions :
## half of sites 
nsite_pred <- round(nsite/2)
## All species 
sites <- sample(as.character(unique(data$site)), nsite_pred)
Id_sites <- rep(sites, nsp)
species <- unique(data$species)
Id_species <- rep(species,each=nsite_pred)
nobs <- length(Id_species)
# Simulate new observations of covariates on those sites 
simdata <- data.frame(site=Id_sites, species=Id_species)
x1 <- rnorm(nsite_pred)
p2 <- poly(x1, 2)
simdata$x1 <- p2[,1]
simdata$x1.2 <- p2[,2]
simdata$x1.SLA <- scale(c(p2[,1] %*% t(SLA)))
# Predictions 
theta_pred <- predict(mod, newdata=simdata,
                      Id_species=Id_species,   Id_sites=Id_sites, type="mean")
  # Tests
  test_that("predict.jSDM works with jSDM_binomial_probit_long_format results considering traits, random site effect and latent variables", {
    expect_equal(length(theta_pred), nsite_pred*nsp)
    expect_equal(sum(is.na(theta_pred)),0)
    expect_equal(sum(theta_pred<0),0)
    expect_equal(sum(theta_pred>1),0)
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
  # Sites and species concerned by predictions :
  ## half of sites 
  nsite_pred <- round(nsite/2)
  ## All species 
  sites <- sample(as.character(unique(data$site)), nsite_pred)
  Id_sites <- rep(sites, nsp)
  species <- unique(data$species)
  Id_species <- rep(species,each=nsite_pred)
  nobs <- length(Id_species)
  # Simulate new observations of covariates on those sites 
  simdata <- data.frame(site=Id_sites, species=Id_species)
  x1 <- rnorm(nsite_pred)
  p2 <- poly(x1, 2)
  simdata$Int <- 1
  simdata$x1.SLA <- scale(c(p2[,1] %*% t(SLA)))
  # Predictions 
  theta_pred <- predict(mod, newdata=simdata,
                        Id_species=Id_species,   Id_sites=Id_sites, type="mean")
  # Tests
  test_that("predict.jSDM works with jSDM_binomial_probit_long_format results considering intercept only, random site effect and latent variables", {
    expect_equal(length(theta_pred), nsite_pred*nsp)
    expect_equal(sum(is.na(theta_pred)),0)
    expect_equal(sum(theta_pred<0),0)
    expect_equal(sum(theta_pred>1),0)
  })

  