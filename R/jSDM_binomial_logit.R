## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name jSDM_binomial_logit
#' @aliases jSDM_binomial_logit
#' @title Binomial logistic regression to fit joint species distribution model
#' @description The \code{jSDM_binomial_logit} function performs a Binomial logistic regression in a Bayesian framework. The function calls a Gibbs sampler written in C++ code which uses an adaptive Metropolis algorithm to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param presence_site_sp A vector indicating the number of successes (or presences) and the absence by a zero for each species at studied sites.
#' @param site_suitability A one-sided formula of the form '~x1+...+xp' with p terms specifying the explicative variables for the suitability process of the model.
#' @param site_data data frame containing the model's explicative variables.
#' @param trials A vector indicating the number of trials for each site. \eqn{t_i} should be superior or equal to \eqn{y_{ij}}{y_ij}, the number of successes for observation \eqn{n}. If \eqn{t_i
#' =0}, then \eqn{y_{ij}=0}{y_ij=0}.
#' @param beta_start Starting values for beta parameters of the suitability process must be either a scalar or a \eqn{p \times n_{species}}{p x n_species} matrix. If \code{beta_start} takes a scalar value, then that value will serve for all of the betas.
#' @param mu_beta Means of the priors for the \eqn{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a p-length vector. If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the betas. The default value is set for an uninformative prior.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta} parameters of the suitability process. \code{V_beta} must be either a scalar or a p-length vector. If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the betas. The default variance is large and set to 1.0E6 for an uninformative flat prior.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param ropt Target acceptance rate for the adaptive Metropolis algorithm. Default to 0.44.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including : \tabular{ll}{
#' mcmc.sp \tab An mcmc object that contains the posterior sample of estimated species effects. This object can be summarized by functions provided by the coda package. \cr
#' mcmc.Deviance \tab The posterior sample of the deviance \eqn{D}, with \eqn{D=-2\log(\prod_{ij} P(y_{ij}|\beta_j,t_i))}{D=-2\log(\prod_ij P(y_ij|\beta_j,t_i))}, is also provided. \cr
#' theta_latent \tab Predictive posterior mean of the probability associated to the suitability process for each observation. \cr
#' model_spec \tab Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function and family, trial sizes, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin. \cr
#'}
#' @details   We model an ecological process where the presence or absence of species \eqn{j} on site \eqn{i} is explained by habitat suitability.
#'
#' \bold{Ecological process:}
#' \deqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij},t_i)}{y_ij ~ Binomial(theta_ij,t_i)}
#' \deqn{logit(\theta_{ij}) = \beta_{0j} + X_i \beta_j}{logit(\theta_ij) = \beta_{0j} + X_i \beta_j}
#' @examples #==============================================
#'# jSDM_binomial_logit()
#'# Example with simulated data
#'#==============================================
#'
#'#=================
#'#== Load libraries
#'library(jSDM)
#'
#'#==================
#'#== Data simulation
#'
#'#= Number of sites
#' nsite <- 100
#'#= Number of species
#' nsp <- 50
#'#= Set seed for repeatability
#'seed <- 1234
#'
#'#= Number of visits associated to each site
#'set.seed(seed)
#'visits <- rpois(nsite,3)
#'visits[visits==0] <- 1
#'
#'#= Ecological process (suitability)
#'x1 <- rnorm(nsite,0,1)
#'set.seed(2*seed)
#'x2 <- rnorm(nsite,0,1)
#'X <- cbind(rep(1,nsite),x1,x2)
#' np <- ncol(X)
#'beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
#'logit.theta <- X %*% t(beta.target)
#'theta <- inv_logit(logit.theta)
#'set.seed(seed)
#'Y <- apply(theta, 2, rbinom, n=nsite, size=visits)
#'
#'#= Site-occupancy model
#'
#' mod <- jSDM_binomial_logit(# Chains
#'                            burnin=100,
#'                            mcmc=100,
#'                            thin=1,
#'                            # Response variable 
#'                            presence_site_sp=Y,
#'                            trials=visits,
#'                            # Explanatory variables
#'                            site_suitability=~x1+x2,
#'                            site_data=X,
#'                            # Starting values 
#'                            beta_start=0,
#'                            # Priors 
#'                            mu_beta=0, V_beta=1.0E6,
#'                            # Various
#'                            seed=1234, ropt=0.44, verbose=1)
#'
#' #== Outputs
#' # Parameter estimates
#' # summary(mod$mcmc)
#'pdf(file=file.path(tempdir(), "Posteriors_jSDM_binomial.pdf"))
#'plot(mod$mcmc.Deviance)
#'for(j in 1:nsp) plot(mod$mcmc.sp[[j]])
#'dev.off()
#' 
#'#= Predictions
#' # summary(mod$theta_latent)
#'pdf(file=file.path(tempdir(), "Pred-Init.pdf"))
#'plot(theta, mod$theta_latent,
#'     main="theta",xlab="obs", ylab="fitted")
#' abline(a=0 ,b=1, col="red")
#' dev.off()
#'@references \tabular{l}{
#' Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and Rebelo, A. G. (2005) Modelling species diversity through species level hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.\cr
#'Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building statistical models to analyze species distributions. \emph{Ecological Applications}, 16, 33-50.\cr
#'}
#' @author \tabular{l}{
#'  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr }
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#' @keywords multivariate logistic regression model binomial biodiversity MCMC, Metropolis algorithm 
#' @export
#' 

jSDM_binomial_logit <- function(# Iteration
                               burnin=5000, mcmc=10000, thin=5,
                               # Data and suitability process
                               presence_site_sp, site_suitability, site_data, trials,
                               # Prior and starting values
                               beta_start=0, mu_beta=0, V_beta=1.0E6,
                               # Various 
                               ropt=0.44, seed=1234, verbose=1)
  
{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  
  #======== 
  # Form response, covariate matrices and model parameters
  #========
  
  #= Response
  Y <- as.matrix(presence_site_sp)
  nsp <- ncol(Y)
  nsite <- nrow(Y)
  nobs <- nsite*nsp
  if(!is.null(trials)){
    T <- as.vector(trials)
  } else {
    T <- rep(1,nobs)
  }  
  #= Suitability
  mf.suit <- model.frame(formula=site_suitability, data=as.data.frame(site_data))
  X <- model.matrix(attr(mf.suit,"terms"), data=mf.suit)
  np <- ncol(X)
  
  #= Iterations
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin
  
  #========== 
  # Check data
  #==========
  check.T.binomial(T, nsite)
  check.Y.binomial(c(as.matrix(Y)), replicate(nsp,T))
  check.X(as.matrix(X), nsite)
  
  #========
  # Initial starting values for M-H
  #========
  beta_start <- form.beta.start.sp(beta_start, np, nsp)
  
  #========
  # Form and check priors
  #========
  mu_beta <- check.mubeta(mu_beta,np)
  V_beta <- check.Vbeta(V_beta,np)
  
  #========
  # call Rcpp function
  #========
  mod <- Rcpp_jSDM_binomial_logit(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                  Y=Y,T=T, X=as.matrix(X),
                                  beta_start=beta_start, mu_beta = mu_beta, V_beta=V_beta,
                                  ropt=ropt, seed=seed, verbose=verbose)
  
  
  
  
  #= Transform Sample list in an MCMC object
  MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
  colnames(MCMC.Deviance) <- "Deviance"
  MCMC.sp <- list()
  for (j in 1:nsp) {
    ## beta_j
    MCMC.beta_j <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
    colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
    MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_j, start=nburn+1, end=ngibbs, thin=nthin)
  }
  #= Model specification, site_suitability,
  model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                     presences=presence_site_sp,
                     site_suitability=site_suitability,
                     site_data=site_data,
                     beta_start=beta_start, mu_beta=mu_beta, V_beta=V_beta,
                     family="binomial", link="logit",
                     ropt=ropt, seed=seed, verbose=verbose)
  
  #= Output
  output <- list(mcmc.sp= MCMC.sp, mcmc.Deviance=MCMC.Deviance,
                 theta_latent=mod$theta_latent,
                 model_spec=model_spec)
  
  class(output) <- "jSDM"
  # return S3 object output belonging to class jSDM
  # acting like list
  return(output)
  
}

# End
