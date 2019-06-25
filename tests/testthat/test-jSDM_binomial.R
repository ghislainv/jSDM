context("test-jSDM_binomial")

#==================
#== Data simulation

#= Number of sites
nsite <- 200

#= Set seed for repeatability
seed <- 1234

#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1

#= Ecological process (suitability)
set.seed(seed)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
beta.target <- c(-1,1,-1)
logit.theta <- X %*% beta.target
theta <- inv_logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Data-sets
data.obs <- data.frame(Y,visits,x1,x2)

# Fit the model 
mod <- jSDM::jSDM_binomial(presences=data.obs$Y,
                      trials=data.obs$visits,
                      suitability=~x1+x2,
                      data=data.obs,
                      burnin=1000, mcmc=1000, thin=1,
                      beta_start=0,
                      mubeta=0, Vbeta=1.0E6,
                      seed=1234, ropt=0.44, verbose=1)

test_that("jSDM_binomial works", {
  expect_equal(mean(mod$theta_latent),0.3227354614)
  expect_equal(as.numeric(mod$mcmc[1,"beta_(Intercept)"]),-1.104863817)
  expect_equal(as.numeric(mod$mcmc[1,"Deviance"]),335.7975026)
})