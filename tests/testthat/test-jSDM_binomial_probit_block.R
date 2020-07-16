context("test-jSDM_binomial_probit_block")

#==================
#== Data simulation

#= Number of sites
nsite <- 50

#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp <- 5

#= Ecological process (suitability)
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

#==================
#== Fit the model
burnin <- 3000
mcmc <- 3000
thin <- 3
nsamp <- mcmc/thin
mod <- jSDM::jSDM_binomial_probit_block(burnin=burnin, mcmc=mcmc, thin=thin,
                                        presence_site_sp = Y ,
                                        site_suitability = ~ x1 + x2,
                                        site_data = X, 
                                        beta_start=0,
                                        mu_beta=0, V_beta=1.0E6,
                                        seed=1234, verbose=1)

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
