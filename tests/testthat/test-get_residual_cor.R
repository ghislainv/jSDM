#context("test-get_residual_cor")

#============= JSDM with latent variables ===============

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

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
x2 <- rnorm(nsite,0,1)
X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#= Number of latent variables
n_latent <- ncol(W)
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
burnin <- 1000
mcmc <- 1000
thin <- 1
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit(burnin=burnin, mcmc=mcmc, thin=thin,
                                  presence_data = Y,
                                  site_formula = ~ x1 + x2,
                                  site_data = X, site_effect="none",
                                  n_latent=n_latent,
                                  beta_start=0,
                                  lambda_start=0, W_start=0,
                                  mu_beta=0, V_beta=10,
                                  mu_lambda=0, V_lambda=1,
                                  seed=1234, verbose=0)
########## get average residual correlation  ################
results <- get_residual_cor(mod, prob=0.95, type="mean")
# Tests
test_that("get_residual_cor works", {
  expect_equal(length(results),8)
  expect_equal(dim(results$cov.mean),c(nsp,nsp))
  expect_equal(dim(results$cor.mean),c(nsp,nsp))
  expect_equal(dim(results$cor.lower),c(nsp,nsp))
  expect_equal(dim(results$cor.upper),c(nsp,nsp))
  expect_equal(dim(results$cor.sig),c(nsp,nsp))
  expect_equal(dim(results$cov.lower),c(nsp,nsp))
  expect_equal(dim(results$cov.upper),c(nsp,nsp))
  expect_equal(dim(results$cov.sig),c(nsp,nsp))
  expect_equal(sum(is.na(results$cor.mean)),0)
  expect_equal(sum(is.na(results$cov.mean)),0)
  expect_equal(sum(is.na(results$cor.lower)),0)
  expect_equal(sum(is.na(results$cor.upper)),0)
  expect_equal(sum(is.na(results$cor.sig)),0)
  expect_equal(sum(is.na(results$cov.lower)),0)
  expect_equal(sum(is.na(results$cov.upper)),0)
  expect_equal(sum(is.na(results$cov.sig)),0)
  expect_equal(sum(diag(results$cov.mean)<=0),0)
  expect_equal(sum(diag(results$cor.mean)!=1),0)
  expect_equal(sum( (results$cor.mean < -1) | (results$cor.mean >1)),0)
})
########## get median residual correlation  ################
results <- get_residual_cor(mod, prob=0.95, type="median")
# Tests
test_that("get_residual_cor works", {
  expect_equal(length(results),8)
  expect_equal(dim(results$cov.median),c(nsp,nsp))
  expect_equal(dim(results$cor.median),c(nsp,nsp))
  expect_equal(dim(results$cor.lower),c(nsp,nsp))
  expect_equal(dim(results$cor.upper),c(nsp,nsp))
  expect_equal(dim(results$cor.sig),c(nsp,nsp))
  expect_equal(dim(results$cov.lower),c(nsp,nsp))
  expect_equal(dim(results$cov.upper),c(nsp,nsp))
  expect_equal(dim(results$cov.sig),c(nsp,nsp))
  expect_equal(sum(is.na(results$cor.median)),0)
  expect_equal(sum(is.na(results$cov.median)),0)
  expect_equal(sum(is.na(results$cor.lower)),0)
  expect_equal(sum(is.na(results$cor.upper)),0)
  expect_equal(sum(is.na(results$cor.sig)),0)
  expect_equal(sum(is.na(results$cov.lower)),0)
  expect_equal(sum(is.na(results$cov.upper)),0)
  expect_equal(sum(is.na(results$cov.sig)),0)
  expect_equal(sum(diag(results$cov.median)<=0),0)
  expect_equal(sum(diag(results$cor.median)!=1),0)
  expect_equal(sum( (results$cor.median < -1) | (results$cor.median >1)),0)
})