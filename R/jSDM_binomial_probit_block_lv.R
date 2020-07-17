## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name jSDM_binomial_probit_block_lv
#' @aliases jSDM_binomial_probit_block_lv
#' @title Binomial probit regression model to fit joint species distribution model with latent variables
#' @description The \code{jSDM_binomial_probit_block_lv} function performs a Binomial probit regression in a Bayesian framework. The function calls a Gibbs sampler written in C++ code which uses conjugate priors to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}.\code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param presence_site_sp A matrix \eqn{n_{site} \times n_{species}}{n_site x n_species} indicating the presence by a 1 (or the absence by a 0) of each species on each site.
#' @param n_latent An integer indicating the number of latent variables.
#' @param site_suitability A one-sided formula of the form '~x1+...+xp' with p terms specifying the explicative variables for the suitability process of the model.
#' @param site_data A data frame containing the model's explicative variables by site.
#' @param beta_start Starting values for beta parameters of the suitability process for each species must be either a scalar or a \eqn{p \times n_{species}}{p x n_species} matrix. If \code{beta_start} takes a scalar value, then that value will serve for all of the betas.
#' @param lambda_start Starting values for lambda parameters corresponding to the latent variables for each species must be either a scalar or a \eqn{n_{latent} \times n_{species}}{n_latent x n_species} upper triangular matrix with strictly positive values on the diagonal. If \code{lambda_start} takes a scalar value, then that value will serve for all of the lambdas except those concerned by the constraints explained above.
#' @param W_start Starting values for latent variables must be either a scalar or a \eqn{n_{site} \times n_{latent}}{n_site x n_latent} matrix. If \code{W_start} takes a scalar value, then that value will serve for all of the Ws.
#' @param mu_beta Means of the Normal priors for the \eqn{\beta}{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a p-length vector. If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the betas. The default value is set to 0 for an uninformative prior.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta}{\beta} parameters of the suitability process. \code{V_beta} must be either a scalar or a \eqn{p \times p}{p x p} symmetric positive semi-definite square matrix. If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the betas, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. The default variance is large and set to 1.0E6 for an uninformative flat prior.
#' @param mu_lambda Means of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. \code{mu_lambda} must be either a scalar or a n_latent-length vector. If \code{mu_lambda} takes a scalar value, then that value will serve as the prior mean for all of the lambdas. The default value is set to 0 for an uninformative prior.
#' @param V_lambda Variances of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. \code{V_lambda} must be either a scalar or a \eqn{n_{latent} \times n_{latent}}{n_latent x n_latent} symmetric positive semi-definite square matrix. If \code{V_lambda} takes a scalar value, then that value will serve as the prior variance for all of the lambdas, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. The default variance is large and set to 10 for an uninformative flat prior.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including : \tabular{ll}{
#' mcmc.latent \tab A list by latent variable of mcmc objects that contains the posterior samples for latent variables Ws.\cr
#' mcmc.sp \tab A list by species of mcmc objects that contains the posterior samples for betas and lambdas.\cr
#' mcmc.Deviance \tab The posterior sample of the deviance \eqn{D}{D}, with \eqn{D=-2\log(\prod_{ij} P(y_{ij}|\beta_j,\lambda_j, W_i))}{D=-2log(\prod_ij P(y_ij | \beta_j, \lambda_j, W_i))}, is also provided.\cr 
#' Z_latent \tab Predictive posterior mean of the latent variable Z. \cr
#' probit_theta_pred \tab Predictive posterior mean of the probability to each species to be present on each site, transformed by probit link function.\cr
#' model_spec \tab Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function, family and number of latent variables, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin.\cr}
#' @details We model an ecological process where the presence or absence of the species is explained by habitat suitability.
#' \bold{Ecological process:}
#' \deqn{y_{ij} \sim \mathcal{B}ernoulli(\theta_{ij})}{y_ij ~ Bernoulli(\theta_ij)}
#' \deqn{probit(\theta_{ij}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j}{probit(\theta_ij) = \beta_0j + X_i x \beta_j +  W_i x \lambda_j }
#' #' @references \tabular{l}{
#' Chib, S. and Greenberg, E. (1998) Analysis of multivariate probit models. \emph{Biometrika}, 85, 347-361. \cr
#' Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.\cr}
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr> \cr
#' Jeanne Clément <jeanne.clement16@laposte.net> \cr }
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}
#' @examples
#'#==============================================
#'# jSDM_binomial_probit_block_lv()
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
#'nsite <- 50
#'
#'#= Set seed for repeatability
#'seed <- 1234
#'set.seed(seed)
#'
#'#= Number of species
#'nsp<- 5
#'
#'#= Number of latent variables
#'n_latent <- 2
#'
#'#= Ecological process (suitability)
#'x1 <- rnorm(nsite,0,1)
#'x2 <- rnorm(nsite,0,1)
#'X <- data.frame(Int=rep(1,nsite),x1=x1,x2=x2)
#'W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#'data <- cbind (X,W)
#'beta.target <- t(matrix(runif(nsp*ncol(X),-2,2),
#'                 byrow=TRUE, nrow=nsp))
#'l.zero <- 0
#'l.diag <- runif(2,0,2)
#'l.other <- runif(nsp*n_latent-3,-2,2)
#'lambda.target <- t(matrix(c(l.diag[1],l.zero,
#'l.other[1],l.diag[2],l.other[-1]), byrow=TRUE, nrow=nsp))
#'param.target <- rbind(beta.target,lambda.target)
#'V <- 1
#'probit_theta <- as.matrix(X) %*% beta.target 
#' + W %*% lambda.target 
#'e <- matrix(rnorm(nsp*nsite,0,sqrt(V)),nsite,nsp)
#'Z_true <- probit_theta + e
#' Y <- matrix (NA, nsite,nsp)
#' for (i in 1:nsite){
#'  for (j in 1:nsp){
#'    if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
#'    else {Y[i,j] <- 0}
#'  }
#'}
#'
#'#==================================
#'#== Site-occupancy model
#'
#'# Increase number of iterations (burnin and mcmc) to get convergence 
#' mod<-jSDM_binomial_probit_block_lv(# Chains
#'                                    burnin=100,
#'                                    mcmc=100,
#'                                    thin=1,
#'                                    # Response variable
#'                                    presence_site_sp = Y ,
#'                                    # Explanatory variables 
#'                                    site_suitability = ~x1+x2,
#'                                    site_data = X,
#'                                    n_latent=2,
#'                                    # Starting values 
#'                                    beta_start=0,
#'                                    lambda_start=0,
#'                                    W_start=0,
#'                                    # Priors 
#'                                    mu_beta=0, V_beta=1.0E6,
#'                                    mu_lambda=0, V_lambda=10,
#'                                    seed=1234, verbose=1)
#' # ===================================================
#' # Result analysis
#' # ===================================================
#' 
#' #==========
#' #== Outputs
#' 
#' #= Parameter estimates
#' ## beta_j
#' # summary(mod$mcmc.sp$sp_1[,1:ncol(X)])
#' pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_probit_block.pdf"))
#' par(mfrow=c(ncol(X),2))
#' for (j in 1:nsp) {
#'   for (p in 1:ncol(X)) {
#'     coda::traceplot(coda::as.mcmc(
#'     mod$mcmc.sp[[paste0("sp_",j)]][,p]))
#'     coda::densplot(coda::as.mcmc(
#'     mod$mcmc.sp[[paste0("sp_",j)]][,p]), 
#'     main = paste(colnames(
#'     mod$mcmc.sp[[paste0("sp_",j)]])[p],
#'     ", species : ",j))
#'     abline(v=beta.target[p,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#'## lambda_j
#' # summary(mod$mcmc.sp$sp_1[,(ncol(X)+1):(ncol(X)+n_latent)])
#' # summary(mod$mcmc.sp$sp_2[,(ncol(X)+1):(ncol(X)+n_latent)])
#' pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_probit_block.pdf"))
#' par(mfrow=c(n_latent*2,2))
#' for (j in 1:nsp) {
#'   for (l in 1:n_latent) {
#'   coda::traceplot(coda::as.mcmc(mod$mcmc.sp
#'                                 [[paste0("sp_",j)]][,ncol(X)+l]))
#'   coda::densplot(coda::as.mcmc(mod$mcmc.sp
#'                               [[paste0("sp_",j)]][,ncol(X)+l]), 
#'                  main=paste(colnames(mod$mcmc.sp[[paste0("sp_",j)]])
#'                  [ncol(X)+l],", species : ",j))
#'  abline(v=lambda.target[l,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#' ## W latent variables
#' pdf(file=file.path(tempdir(), "Posteriors_lv_jSDM_probit_block.pdf"))
#' par(mfrow=c(1,2))
#' for (l in 1:n_latent) {
#' plot(W[,l],
#' summary(mod$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
#' main = paste0("Latent variable W_", l),
#' xlab ="obs", ylab="fitted")
#' abline(a=0,b=1,col='red')
#' }
#' dev.off()
#' 
#' ## Deviance
#' summary(mod$mcmc.Deviance)
#' plot(mod$mcmc.Deviance)
#' 
#' #= Predictions
#' 
#' pdf(file=file.path(tempdir(), "Pred-Init.pdf"))
#' ## probit_theta
#' # summary(mod$probit_theta_pred)
#' par(mfrow=c(1,1))
#' plot(probit_theta,mod$probit_theta_pred)
#' abline(a=0,b=1,col='red')
#' 
#' ## Z
#' # summary(mod$Z_latent)
#' plot(Z_true,mod$Z_latent)
#' abline(a=0,b=1,col='red')
#' dev.off()
#' @keywords Binomial probit regression biodiversity JSDM hierarchical Bayesian models MCMC Markov Chains Monte Carlo Gibbs Sampling
#' @export 

jSDM_binomial_probit_block_lv <- function (burnin=5000, mcmc=15000, thin=10,
                                           presence_site_sp, site_suitability,
                                           site_data, n_latent=2,
                                           lambda_start=0, W_start=0,
                                           beta_start=0, 
                                           mu_beta=0, V_beta=1.0E6,
                                           mu_lambda=0, V_lambda=10,
                                           seed=1234, verbose=1)
  
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
  nsp <- ncol(Y)
  nsite <- nrow(Y)
  nobs <- nsite*nsp
  T <- matrix(1, nsite, nsp)
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
  check.T.binomial(c(T), nobs)
  check.Y.binomial(c(as.matrix(Y)), c(T))
  check.X(X, nsite)
  
  #========
  # Initial starting values for M-H
  #========
  beta_start <- form.beta.start.sp(beta_start, np, nsp)
  lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
  W_start <-form.W.start.sp(W_start, nsite, n_latent)
  param_start = rbind(beta_start,lambda_start)
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mu_beta,np)
  Vbeta <- check.Vbeta.mat(V_beta,np)
  mulambda <- check.mubeta(mu_lambda,n_latent)
  Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
  Vparam <- diag(c(diag(Vbeta),diag(Vlambda)))
  muparam <- c(mubeta,mulambda)
  V_W <- diag(rep(1,n_latent))
  
  #========
  # call Rcpp function
  #========
  mod <- Rcpp_jSDM_binomial_probit_block_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                            Y=as.matrix(Y), X=as.matrix(X),
                                            param_start= param_start, V_param=Vparam, mu_param = muparam,
                                            W_start=W_start, V_W=V_W,
                                            seed=seed, verbose=verbose)
  
  #= Transform Sample list in an MCMC object
  MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
  MCMC.sp <- list()
  for (j in 1:nsp) {
    ## beta_j
    MCMC.beta_j <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
    colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
    ## lambda_j
    MCMC.lambda_j <- coda::mcmc(mod$param[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
    colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
    
    MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
  }
  ## W latent variables 
  MCMC.latent <- list()
  for (l in 1:n_latent) {
    MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
    MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
  }
  
  #= Model specification, site_suitability,
  model_spec <- list(presences=presence_site_sp,
                     site_suitability=site_suitability,
                     site_data=site_data, n_latent=n_latent,
                     burnin=burnin, mcmc=mcmc, thin=thin,
                     beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                     lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                     W_start=W_start, V_W=V_W,
                     family="binomial", link="probit",
                     seed=seed, verbose=verbose)
  
  #= Output
  output <- list(mcmc.Deviance=MCMC.Deviance,
                 mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                 Z_latent=mod$Z_latent, 
                 probit_theta_pred=mod$probit_theta_pred,
                 model_spec=model_spec)
  
  class(output) <- "jSDM"
  # return S3 object output belonging to class jSDM
  # acting like list
  return(output)
}

# End
