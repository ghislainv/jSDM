## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' Binomial logistic regression joint species distribution model
#' @description The \code{jSDM_binomial_logit} function performs a Binomial logistic regression in a Bayesian framework. The function calls a Gibbs sampler written in C++ code which uses an adaptive Metropolis algorithm to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}. \code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param presence_site_sp A vector indicating the number of successes (or presences) for each observation..
#' @param site_suitability A one-sided formula of the form '~x1+...+xp' with p terms specifying the explicative variables for the suitability process of the model.
#' @param site_data data frame containing the model's explicative variables.
#' @param trials A vector indicating the number of trials for each observation. \eqn{t_n} should be superior or equal to \eqn{y_n}, the number of successes for observation \eqn{n}. If \eqn{t_n
#' =0}, then \eqn{y_n=0}.
#' @param beta_start Starting values for beta parameters of the suitability process. If \code{beta_start} takes a scalar value, then that value will serve for all of the betas.
#' @param mu_beta Means of the priors for the \eqn{\beta}{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a p-length vector. If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the betas. The default value is set for an uninformative prior.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta}{\beta} parameters of the suitability process. \code{V_beta} must be either a scalar or a p-length vector. If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the betas. The default variance is large and set to 1.0E6 for an uninformative flat prior.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param ropt Target acceptance rate for the adaptive Metropolis algorithm. Default to 0.44.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including : \tabular{ll}{
#' mcmc \tab An mcmc object that contains the posterior sample of estimated parameters betas and the posterior sample of the deviance \eqn{D}, with \eqn{D=-2\log(\prod_n P(y_n|\beta,t_n))}. This object can be summarized by functions provided by the coda package. \cr
#' theta_latent \tab Predictive posterior mean of the probability associated to the suitability process for each observation. \cr
#' model_spec \tab Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function and family, trial sizes, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin. \cr
#' }
#' @details   We model an ecological process where the presence or absence of species \eqn{j} on site \eqn{i} is explained by habitat suitability.
#'
#' \bold{Ecological process:}
#' \deqn{y_n \sim \mathcal{B}inomial(\theta_n,t_n)}{y_n ~ Binomial(\theta_n,t_n)}
#' \deqn{logit(\theta_n) = X_n \beta}{logit(\theta_n) = X_n \beta}
#' 
#' @examples 
#' #==============================================
#'# jSDM_binomial_logit_one_species()
#'# Example with simulated data
#'#==============================================
#'
#'#=================
#'#== Load libraries
#'library(jSDM)
#' library(coda)
#'#==================
#'#== Data simulation
#'
#'#= Number of sites
#' nsite <- 200
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
#'beta.target <- beta.target <- c(-1,1,-1)
#'logit.theta <- X %*% beta.target
#'theta <- inv_logit(logit.theta)
#'set.seed(seed)
#'Y <- rbinom(nsite,visits,theta)
#'
#'#= Site-occupancy model
#'
#' mod <- jSDM_binomial_logit_one_species(# Chains 
#'                                        burnin=100,
#'                                        mcmc=100,
#'                                        thin=1,
#'                                        # Response variable 
#'                                        presence_site_sp=Y,
#'                                        trials=visits,
#'                                        # Explanatory variables
#'                                        site_suitability=~x1+x2,
#'                                        site_data=X,
#'                                        # Starting values 
#'                                        beta_start=0,
#'                                        # Prior hyperparameters
#'                                        mu_beta=0, V_beta=1.0E6,
#'                                        # Various
#'                                        seed=1234, ropt=0.44, verbose=1)
#'
#' #== Outputs
#' # Parameter estimates
#'summary(mod$mcmc)
#'pdf(file=file.path(tempdir(), "Posteriors_jSDM_binomial.pdf"))
#'plot(mod$mcmc[,"Deviance"])
#'plot(mod$mcmc[,grepl("beta", colnames(mod$mcmc))])
#'dev.off()
#' 
#'#= Predictions
#'summary(mod$theta_latent)
#'pdf(file=file.path(tempdir(), "Pred-Init.pdf"))
#'plot(theta, mod$theta_latent,
#'     main="theta",xlab="obs", ylab="fitted")
#' abline(a=0 ,b=1, col="red")
#' dev.off()
#' @references \tabular{l}{
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
jSDM_binomial_logit_one_species <- function (# Chains
                                             burnin=5000, mcmc=10000, thin=10,
                                             # Response variable
                                             presence_site_sp, trials,
                                             # Explanatory variables 
                                             site_suitability, site_data,
                                             # Starting values
                                             beta_start,
                                             # Priors
                                             mu_beta=0, V_beta=1.0E6,
                                             # Various
                                             seed=1234, ropt=0.44, verbose=1)
  
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
  Y <- presence_site_sp
  nobs <- length(Y)
  T <- trials
  #= Suitability
  mf.suit <- model.frame(formula=site_suitability,data=as.data.frame(site_data))
  X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  np <- ncol(X)
  #= Iterations
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin
  
  #========== 
  # Check data
  #==========
  check.T.binomial(T,nobs)
  check.Y.binomial(as.vector(Y),T)
  check.X(X,nobs)
  
  #========
  # Initial starting values for M-H
  #========
  beta_start <- form.beta.start(beta_start,np)
  
  #========
  # Form and check priors
  #========
  mu_beta <- check.mubeta(mu_beta,np)
  V_beta <- check.Vbeta(V_beta,np)
  
  #========
  # call Rcpp function
  #========
  mod <- Rcpp_jSDM_binomial_logit_one_species(ngibbs, nthin, nburn,
                                              Y=as.vector(Y), T=T, X=as.matrix(X),
                                              beta_start=beta_start, mu_beta=mu_beta, V_beta=V_beta,
                                              seed=seed, ropt=ropt, verbose=verbose)
  
  #= Matrix of MCMC samples
  Matrix <- cbind(mod$beta, mod$Deviance)
  names.fixed <- paste("beta_",colnames(X),sep="")
  colnames(Matrix) <- c(names.fixed,"Deviance")
  
  #= Transform Sample list in an MCMC object
  MCMC <- coda::mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
  
  #= Model specification
  model_spec <- list(presences=presence_site_sp, trials=trials,
                     site_suitability=site_suitability,
                     site_data=site_data,
                     burnin=burnin, mcmc=mcmc, thin=thin,
                     beta_start=beta_start, mu_beta=mu_beta, V_beta=V_beta,
                     family="binomial", link="logit",
                     seed=seed, ropt=ropt, verbose=verbose)
  
  #= Output
  output <- list(mcmc=MCMC, theta_latent=mod$theta_latent,
                 model_spec=model_spec)
  
  class(output) <- "jSDM"
  # return S3 object output belonging to class jSDM
  # acting like list
  return(output)
  
}

# End
