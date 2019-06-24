\name{jSDM_binomial}
\alias{jSDM_binomial}

\title{Binomial logistic regression model}

\description{The \code{jSDM_binomial} function performs a Binomial
  logistic regression in a Bayesian framework. The function calls
  a Gibbs sampler written in C++ code which uses an adaptive Metropolis
  algorithm to estimate the conditional posterior distribution of
  model's parameters.}

\usage{jSDM_binomial(presences, trials, suitability, data, 
burnin = 5000, mcmc = 10000, thin = 10,
beta_start, mubeta = 0, Vbeta = 1e+06, seed = 1234, ropt = 0.44, verbose = 1)}

\arguments{
  
  \item{presences}{A vector indicating the number of successes (or
  presences) for each observation.}

  \item{trials}{A vector indicating the number of trials for each
  observation. \eqn{t_n} should be superior or equal to \eqn{y_n}, the
  number of successes for observation \eqn{n}. If \eqn{t_n=0}, then
  \eqn{y_n=0}.}

  \item{suitability}{A one-sided formula of the form '~x1+...+xp' with p
  terms specifying the explicative variables for the suitability process
  of the model.}
  
  \item{data}{A data frame containing the model's explicative
    variables.}

  \item{burnin}{The number of burnin iterations for the sampler.}
    
  \item{mcmc}{The number of Gibbs iterations for the sampler. Total
    number of Gibbs iterations is equal to
    \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and
    superior or equal to 100 so that the progress bar can be displayed.}
    
  \item{thin}{The thinning interval used in the simulation. The number
    of mcmc iterations must be divisible by this value.} 

  \item{beta_start}{Starting values for beta parameters of the
  suitability process. If \code{beta_start} takes a scalar value, then
  that value will serve for all of the betas.}

  \item{mubeta}{Means of the priors for the \eqn{\beta}{beta} parameters
  of the suitability process. \code{mubeta} must be either a scalar or a
  p-length vector. If \code{mubeta} takes a scalar value, then that value will
  serve as the prior mean for all of the betas. The default value is set
  to 0 for an uninformative prior.}

  \item{Vbeta}{Variances of the Normal priors for the \eqn{\beta}{beta}
  parameters of the suitability process. \code{Vbeta} must be either a
  scalar or a p-length vector. If \code{Vbeta} takes a scalar value,
  then that value will serve as the prior variance for all of the
  betas. The default variance is large and set to 1.0E6 for an
  uninformative flat prior.}
 
  \item{seed}{The seed for the random number generator. Default to 1234.}
  
  \item{ropt}{Target acceptance rate for the adaptive Metropolis algorithm. Default to 0.44.}
  
  \item{verbose}{A switch (0,1) which determines whether or not the
    progress of the sampler is printed to the screen. Default is 1: a
    progress bar is printed, indicating the step (in \%) reached by the
    Gibbs sampler.}
    
}

\value{
  
  \item{mcmc}{An mcmc object that contains the posterior sample. This
    object can be summarized by functions provided by the coda
    package. The posterior sample of the deviance \eqn{D}{D}, with
    \eqn{D=-2\log(\prod_i P(y_i|\beta,t_i))}{%
      D=-2log(prod_i P(y_i|beta,t_i))}, is also provided.}

  \item{theta_latent}{Predictive posterior mean of the probability
    associated to the suitability process for each observation.}
    
  \item{spec}{Model's specifications}
      
}

\details{

  We model an ecological process where the presence or absence
  of the species is explained by habitat suitability.

  \bold{Ecological process:}
  \deqn{y_i \sim \mathcal{B}inomial(\theta_i,t_i)}{y_i ~ Binomial(theta_i,t_i)}
  \deqn{logit(\theta_i) = X_i \beta}{logit(theta_i) = X_i * beta}
  
}

\references{

  Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and
  Rebelo, A. G. (2005) Modelling species diversity through species level
  hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
  
  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building
  statistical models to analyze species distributions. \emph{Ecological
  Applications}, 16, 33-50.

}

\author{
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}
}

\seealso{
  \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
}

\examples{

#==============================================
# jSDM_binomial()
# Example with simulated data
#==============================================

#=================
#== Load libraries
library(jSDM)

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
logit.theta <- X \%*\% beta.target
theta <- inv_logit(logit.theta)
set.seed(seed)
Y <- rbinom(nsite,visits,theta)

#= Data-sets
data.obs <- data.frame(Y,visits,x1,x2)

#==================================
#== Site-occupancy model

mod_jSDM_binomial <- jSDM_binomial(presences=data.obs$Y,
                                   trials=data.obs$visits,
                                   suitability=~x1+x2,
                                   data=data.obs,
                                   burnin=1000, mcmc=1000, thin=1,
                                   beta_start=0,
                                   mubeta=0, Vbeta=1.0E6,
                                   seed=1234, ropt=0.44, verbose=1)

#==========
#== Outputs

#= Parameter estimates
summary(mod_jSDM_binomial$mcmc)
pdf(file=file.path(tempdir(), "Posteriors_jSDM_binomial.pdf"))
plot(mod_jSDM_binomial$mcmc)
dev.off()

#== glm resolution to compare
mod_glm <- glm(cbind(Y,visits-Y)~x1+x2,family="binomial",data=data.obs)
summary(mod_glm)

#= Predictions
summary(mod_jSDM_binomial$theta_latent)
pdf(file=file.path(tempdir(), "Pred-Init.pdf"))
plot(theta, mod_jSDM_binomial$theta_latent)
abline(a=0 ,b=1, col="red")
dev.off()

}

\keyword{Binomial logistic regression}
\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
