\name{jSDM_probit_block}
\alias{jSDM_probit_block}

\title{Binomial probit regression model}

\description{The \code{jSDM_probit_block} function performs a Binomial probit regression in a Bayesian framework. The function calls a Gibbs sampler written in C++ code which uses conjugate priors
to estimate the conditional posterior distribution of
model's parameters.}

\usage{jSDM_probit_block(presence_site_sp, site_suitability,
site_data, n_latent=2, burnin=5000, mcmc=15000, thin=10,
alpha_start=0, beta_start=0, lambda_start=0, W_start=0,
V_alpha_start=1, shape=0.5, rate=0.0005, 
mu_beta=0, V_beta=1.0E6, mu_lambda=0, V_lambda=10,
seed=1234, verbose=1)}

\arguments{
  
  \item{presence_site_sp}{A matrix \eqn{n_{site} \times n_{species}}{n_site x n_species} indicating the presence by a 1 (or the absence by a 0) of each species on each site.}

  \item{n_latent}{An integer indicating the number of latent variables.}

  \item{site_suitability}{A one-sided formula of the form '~x1+...+xp' with p terms specifying the explicative variables for the suitability process of the model.}
  
  \item{site_data}{A data frame containing the model's explicative
    variables by site.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to \code{burnin+mcmc}.     \code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.} 

  \item{beta_start}{Starting values for beta parameters of the
  suitability process for each species must be either a scalar or a \eqn{p \times n_{species}}{p x n_species} matrix. If \code{beta_start} takes a scalar value, then that value will serve for all of the betas.}
  
    \item{lambda_start}{Starting values for lambda parameters corresponding to the latent variables for each species must be either a scalar or a \eqn{n_{latent} \times n_{species}}{n_latent x n_species} upper triangular matrix with strictly positive values on the diagonal. If \code{lambda_start} takes a scalar value, then that value will serve for all of the lambdas except those concerned by the constraints explained above.}
  
    \item{alpha_start}{Starting values for random site effect parameters must be either a scalar or a nsite-length vector. If \code{alpha_start} takes a scalar value, then that value will serve for all of the alphas.}

    \item{V_alpha_start}{Starting value for variance of random site effect must be a stricly positive scalar.}

    \item{W_start}{Starting values for latent variables must be either a scalar or a \eqn{nsite \times n_latent}{n_site x n_latent} matrix. If \code{W_start} takes a scalar value, then that value will serve for all of the Ws.}  
    
    \item{shape}{Shape parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}. Must be a stricly positive scalar. Default to 0.5 for weak informative prior.}

    \item{rate}{Rate parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}. Must be a stricly positive scalar. Default to 0.0005 for weak informative prior.}

  \item{mu_beta}{Means of the Normal priors for the \eqn{\beta}{beta} parameters of the suitability process. \code{mubeta} must be either a scalar or a p-length vector. If \code{mubeta} takes a scalar value, then that value will serve as the prior mean for all of the betas. The default value is set to 0 for an uninformative prior.}

  \item{V_beta}{Variances of the Normal priors for the \eqn{\beta}{beta} parameters of the suitability process. \code{Vbeta} must be either a scalar or a \eqn{p \times p}{p x p} symmetric positive semi-definite square matrix. If \code{Vbeta} takes a scalar value,
then that value will serve as the prior variance for all of the
betas, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. The default variance is large and set to 1.0E6 for an uninformative flat prior.}

  \item{mu_lambda}{Means of the Normal priors for the \eqn{\lambda}{lambda} parameters corresponding to the latent variables. \code{mulambda} must be either a scalar or a n_latent-length vector. If \code{mulambda} takes a scalar value, then that value will serve as the prior mean for all of the lambdas. The default value is set to 0 for an uninformative prior.}

  \item{V_lambda}{Variances of the Normal priors for the \eqn{\lambda}{lambda} parameters corresponding to the latent variables. \code{Vlambda} must be either a scalar or a \eqn{n_{latent} \times n_{latent}}{n_latent x n_latent} symmetric positive semi-definite square matrix. If \code{Vlambda} takes a scalar value, then that value will serve as the prior variance for all of the lambdas, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. The default variance is large and set to 10 for an uninformative flat prior.}
 
  \item{seed}{The seed for the random number generator. Default to 1234.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.}
    
}

\value{
 
  \item{mcmc.alpha}{An mcmc object that contains the posterior samples for alphas. This object can be summarized by functions provided by the coda package.}
  
  \item{mcmc.Valpha}{An mcmc object that contains the posterior samples for variance of random site effect.}
  
  \item{mcmc.latent}{A list by latent variable of mcmc objects that contains the posterior samples for latent variables Ws.}
    
  \item{mcmc.sp}{A list by species of mcmc objects that contains the posterior samples for betas and lambdas.}

   \item{mcmc.Deviance}{The posterior sample of the deviance \eqn{D}{D}, with \eqn{D=-2\log(\prod_{i,j} P(y_{i,j}|\beta_j,\lambda_j, \alpha_i, W_i))}{D=-2log(\Pi_{i,j} P(y_{i,j}|\beta_j,\lambda_j, \alpha_i, W_i))}, is also provided.} 
   
  \item{Z_latent}{Predictive posterior mean of the latent variable Z.}
    
  \item{probit_theta_pred}{Predictive posterior mean of the probability to each species to be present on each site, transformed by probit link function.}
    
  \item{model_spec}{Model's specifications}
      
}

\details{

  We model an ecological process where the presence or absence
  of the species is explained by habitat suitability.

  \bold{Ecological process:}
  \deqn{y_{i,j} \sim \mathcal{B}ernoulli(\theta_{i,j})}{y_{i,j} ~ Bernoulli(\theta_{i,j})}
  \deqn{probit(\theta_{i,j}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j + \alpha_i }{probit(\theta_i) = \beta_{0j} + X_i * \beta_j +  W_i * \lambda_j + \alpha_i}
  
}

\references{

  Chib, S. and Greenberg, E. (1998) Analysis of multivariate probit models. \emph{Biometrika}, 85, 347-361.
  
  Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.

}

\author{
  Jeanne Clement \email{jeanne.clement16@laposte.net}
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}
}

\seealso{
  \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
}

\examples{

#==============================================
# jSDM_probit_block()
# Example with simulated data
#==============================================

#=================
#== Load libraries
library(jSDM)

#==================
#== Data simulation

#= Number of sites
nsite <- 50

#= Set seed for repeatability
seed <- 1234
set.seed(seed)

#= Number of species
nsp<- 5

#= Number of latent variables
n_latent <- 2

#= Ecological process (suitability)
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
probit_theta <- as.matrix(X) \%*\% beta.target + W \%*\% lambda.target + alpha.target
e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
Z_true <- probit_theta + e
Y <- matrix (NA, nsite,nsp)
for (i in 1:nsite){
  for (j in 1:nsp){
    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
    else {Y[i,j] <- 0}
  }
}

#==================================
#== Site-occupancy model

# Increase number of iterations (burnin and mcmc) to get convergence 
mod_jSDM_probit_block <- jSDM::jSDM_probit_block ( presence_site_sp = Y ,
                                             site_suitability = ~ x1 + x2,
                                             site_data = X[,-1], n_latent=2,
                                             burnin=3000, mcmc=3000, thin=3,
                                             alpha_start=0, beta_start=0,
                                             lambda_start=0, W_start=0,
                                             V_alpha_start=1,
                                             shape=0.5, rate=0.0005,
                                             mu_beta=0, V_beta=1.0E6,
                                             mu_lambda=0, V_lambda=10,
                                             seed=1234, verbose=1)

# ===================================================
# Result analysis
# ===================================================

#==========
#== Outputs

#= Parameter estimates

## alpha
summary(mod_jSDM_probit_block$mcmc.alpha)
pdf(file=file.path(tempdir(), "Posteriors_alpha_jSDM_probit_block.pdf"))
plot(alpha.target,
summary(mod_jSDM_probit_block$mcmc.alpha)[[1]][,"Mean"],
xlab ="alphas target", ylab ="alphas estimated")
abline(a=0,b=1,col='red')
dev.off()

## Valpha
summary(mod_jSDM_probit_block$mcmc.Valpha)
pdf(file=file.path(tempdir(), "Posteriors_Valpha_jSDM_probit_block.pdf"))
par(mfrow=c(1,2))
coda::traceplot(mod_jSDM_probit_block$mcmc.Valpha)
coda::densplot(mod_jSDM_probit_block$mcmc.Valpha)
abline(v=Valpha.target,col='red')
dev.off()


## beta_j
summary(mod_jSDM_probit_block$mcmc.sp$sp_1[,1:ncol(X)])
pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_probit_block.pdf"))
par(mfrow=c(ncol(X),2))
for (j in 1:nsp) {
  for (p in 1:ncol(X)) {
    coda::traceplot(coda::as.mcmc(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]][,p]))
    coda::densplot(coda::as.mcmc(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]][,p]), 
    main = paste(colnames(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]])[p],", species : ",j))
    abline(v=beta.target[p,j],col='red')
  }
}
dev.off()

## lambda_j
summary(mod_jSDM_probit_block$mcmc.sp$sp_1[,(ncol(X)+1):(ncol(X)+n_latent)])
summary(mod_jSDM_probit_block$mcmc.sp$sp_2[,(ncol(X)+1):(ncol(X)+n_latent)])
pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_probit_block.pdf"))
par(mfrow=c(n_latent*2,2))
for (j in 1:nsp) {
  for (l in 1:n_latent) {
    coda::traceplot(coda::as.mcmc(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]][,ncol(X)+l]))
    coda::densplot(coda::as.mcmc(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]][,ncol(X)+l]), 
    main=paste(colnames(mod_jSDM_probit_block$mcmc.sp[[paste0("sp_",j)]])[ncol(X)+l],",
    species : ",j))
    abline(v=lambda.target[l,j],col='red')
  }
}
dev.off()

## W latent variables
pdf(file=file.path(tempdir(), "Posteriors_lv_jSDM_probit_block.pdf"))
par(mfrow=c(1,2))
for (l in 1:n_latent) {
plot(W[,l],
summary(mod_jSDM_probit_block$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
main = paste0("Latent variable W_", l),
xlab =paste0("W_", l, " target"), ylab =paste0("W_", l, " estimated"))
abline(a=0,b=1,col='red')
}
dev.off()

## Deviance
summary(mod_jSDM_probit_block$mcmc.Deviance)
plot(mod_jSDM_probit_block$mcmc.Deviance)

#= Predictions

pdf(file=file.path(tempdir(), "Pred-Init.pdf"))
## probit_theta
summary(mod_jSDM_probit_block$probit_theta_pred)
par(mfrow=c(1,1))
plot(probit_theta,mod_jSDM_probit_block$probit_theta_pred)
abline(a=0,b=1,col='red')

## Z
summary(mod_jSDM_probit_block$Z_latent)
plot(Z_true,mod_jSDM_probit_block$Z_latent)
abline(a=0,b=1,col='red')
dev.off()

}

\keyword{Binomial probit regression}
\keyword{biodiversity}
\keyword{joint species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Gibbs Sampling}
