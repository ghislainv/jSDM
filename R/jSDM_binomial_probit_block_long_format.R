## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name jSDM_binomial_probit_block_long_format
#' @aliases jSDM_binomial_probit_block_long_format
#' @title Binomial probit regression on long format data 
#' @description The \code{jSDM_binomial_probit_block_long_format} function performs a Binomial probit regression in a Bayesian framework. 
#' The function calls a Gibbs sampler written in C++ code which uses conjugate priors to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}.\code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param data A \code{data.frame} with at least the following columns :  
#' \tabular{ll}{
#' \code{Y} \tab \eqn{n_{obs}}{n_obs}-length vector indicating the presence by a 1 (or absence by a 0) of the species observed during each visit of the sites. \cr
#' \code{site} \tab numerc or character \eqn{n_{obs}}{n_obs}-length vector indicating the visited site (sites can be visited several times). \cr
#' \code{species} \tab numerc or character eqn{n_{obs}}{n_obs}-length vector indicating the species observed (species may not have been recorded at all sites) \cr 
#' \code{x1,...,xp} \tab columns of explicative variables for the suitability process of the model \cr
#' }
#' @param site_suitability A one-sided formula with the form '~x1 + ... + xd + species:x1 + ... + species:xp' with \eqn{d} terms related to \eqn{\beta} parameters and \eqn{p} terms related to species effects \eqn{\beta_sp} parameters,
#'  specifying the explicative variables for the suitability process of the model as formula used by the \code{lm} function.
#' @param n_latent An integer which specifies the number of latent variables to generate. Defaults to \code{0}.
#' @param site_effect A string indicating whether row effects are included as fixed effects (\code{"fixed"}), as random effects (\code{"random"}), or not included (\code{"none"}) in the model. 
#'  If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used.
#'  If random effects, they are drawn from a normal distribution with mean zero and unknown variance, analogous to a random intercept in mixed models. Defaults to \code{"none"}.
#' @param beta_start Starting values for beta parameters of the suitability process must be either a scalar or a \eqn{d}-length vector. If \code{beta_start} takes a scalar value, then that value will serve for all of the \eqn{\beta} parameters.
#' @param beta_sp_start Starting values for beta_sp parameters of the suitability process for each species must be either a scalar or a \eqn{p \times n_{species}}{p x n_species} matrix. If \code{beta_sp_start} takes a scalar value, then that value will serve for all of the \eqn{\beta_{sp}}{\beta_sp} parameters.
#' @param lambda_start Starting values for lambda parameters corresponding to the latent variables for each species must be either a scalar or a \eqn{n_{latent} \times n_{species}}{n_latent x n_species} upper triangular matrix with strictly positive values on the diagonal, ignored if \code{n_latent=0}.
#'  If \code{lambda_start} takes a scalar value, then that value will serve for all of the \eqn{\lambda} parameters except those concerned by the constraints explained above.
#' @param W_start Starting values for latent variables must be either a scalar or a \eqn{nsite \times n_latent}{n_site x n_latent} matrix, ignored if \code{n_latent=0}.
#'  If \code{W_start} takes a scalar value, then that value will serve for all of the \eqn{W_{il}}{W_il} with \eqn{i=1,\ldots,n_{site}}{l=1,...,n_site} and \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}.
#' @param alpha_start Starting values for random site effect parameters must be either a scalar or a \eqn{n_{site}}{n_site}-length vector, ignored if \code{site_effect="none"}.
#'  If \code{alpha_start} takes a scalar value, then that value will serve for all of the \eqn{\alpha} parameters.
#' @param V_alpha Starting value for variance of random site effect if \code{site_effect="random"} or constant variance of the Normal prior for the fixed site effect if \code{site_effect="fixed"}.
#' Must be a stricly positive scalar, ignored if \code{site_effect="none"}.
#' @param shape Shape parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}. 
#' Must be a stricly positive scalar. Default to 0.5 for weak informative prior.
#' @param rate Rate parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}
#' Must be a stricly positive scalar. Default to 0.0005 for weak informative prior.
#' @param mu_beta Means of the Normal priors for the \eqn{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a \eqn{d}-length vector.
#'  If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\beta} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta} parameters of the suitability process. \code{V_beta} must be either a scalar or a \eqn{d \times d}{d x d} symmetric positive semi-definite square matrix.
#'  If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\beta} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. 
#'  The default variance is large and set to \code{1e+06} for an uninformative flat prior.
#' @param mu_beta_sp Means of the Normal priors for the \eqn{\beta_{sp}}{\beta_sp} parameters of the suitability process. \code{mu_beta_sp} must be either a scalar or a \eqn{p}-length vector.
#'  If \code{mu_beta_sp} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\beta} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_beta_sp Variances of the Normal priors for the \eqn{\beta_{sp}}{\beta_sp} parameters of the suitability process. \code{V_beta_sp} must be either a scalar or a \eqn{p \times p}{p x p} symmetric positive semi-definite square matrix.
#'  If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\beta_{sp}}{\beta_sp} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. 
#' The default variance is large and set to \code{1e+06} for an uninformative flat prior.
#' @param mu_lambda Means of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. 
#' \code{mu_lambda} must be either a scalar or a \eqn{n_{latent}}{n_latent}-length vector. 
#' If \code{mu_lambda} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\lambda}{\lambda} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_lambda Variances of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. \code{V_lambda} must be either a scalar or a \eqn{n_{latent} \times n_{latent}}{n_latent x n_latent} symmetric positive semi-definite square matrix.
#' If \code{V_lambda} takes a scalar value, then that value will serve as the prior variance for all of \eqn{\lambda}{\lambda} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal.
#' The default variance is large and set to 10 for an uninformative flat prior.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including : \tabular{ll}{
#' mcmc.alpha \tab An mcmc object that contains the posterior samples for for site effects \eqn{\alpha_i}, not returned if \code{site_effect="none"}.\cr
#' mcmc.V_alpha \tab An mcmc object that contains the posterior samples for variance of random site effect, not returned if \code{site_effect="none"} or \code{site_effect="fixed"}.\cr
#' mcmc.latent \tab A list by latent variable of mcmc objects that contains the posterior samples for latent variables  \eqn{W_l} with \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}, not returned if \code{n_latent=0}.\cr
#' mcmc.sp \tab A list by species of mcmc objects that contains the posterior samples for species effects \eqn{\beta_{sp}}{\beta_sp} and the loading factors \eqn{\lambda} if \code{n_latent>0}.\cr
#' mcmc.Deviance \tab The posterior sample of the deviance \eqn{D}{D}, with \eqn{D=-2\log(\prod_{n} P(y_{n}|\beta_j,\lambda_j, \alpha_i, W_i))}{D=-2log(\prod_n P(y_n|\beta_j,\lambda_j, \alpha_i, W_i))}, is also provided.\cr 
#' Z_latent \tab Predictive posterior mean of the latent variable Z. \cr
#' probit_theta_pred \tab Predictive posterior mean of the probability to each species to be present on each site, transformed by probit link function.\cr
#' model_spec \tab Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function, family and number of latent variables, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin.\cr}
#' The \code{mcmc.} objects can be summarized by functions provided by the \code{coda} package. 
#' @details We model an ecological process where the presence or absence of species \eqn{j} on site \eqn{i} is explained by habitat suitability.
#'
#' \bold{Ecological process:}
#' \deqn{y_{n} \sim \mathcal{B}ernoulli(\theta_{n}) \text{ such as } species_n=j \text{ and } site_n=i}{y_n ~ Bernoulli(\theta_n) such as species_n=j and site_n=i}
#'  where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{n}) = D_n \beta + X_n \beta_{sp.j}}{(\theta_n) = D_n \beta + X_n \beta_sp.j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{n}) = D_n \beta + X_n \beta_{sp.j} + W_i \lambda_j}{(\theta_n) =  D_n \beta  + X_n \beta_sp.j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_{n}) = D_n \beta + X_n \beta_{sp.j}  + \alpha_i}{(\theta_n) = D_n \beta + X_n \beta_sp.j + \alpha_i}  and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_{n}) = D_n \beta + X_n \beta_{sp.j} + W_i \lambda_j + \alpha_i}{(\theta_n) = D_n \beta  + X_n  \beta_sp.j +  W_i \lambda_j + \alpha_i} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{n}) = D_n \beta  + X_n \beta_{sp.j}  + \alpha_i}{(\theta_n) = D_n \beta + X_n \beta_sp.j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{n}) = D_n \beta + X_n \beta_{sp.j} + W_i \lambda_j + \alpha_i}{(\theta_n) = D_n \beta  + X_n  \beta_sp.j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }
#' @references \tabular{l}{
#' Chib, S. and Greenberg, E. (1998) Analysis of multivariate probit models. \emph{Biometrika}, 85, 347-361. \cr
#' Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.\cr}
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr> \cr
#' Jeanne Clément <jeanne.clement16@laposte.net> \cr }
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}} \code{\link{jSDM_binomial_probit_block}} \code{\link{jSDM_binomial_logit}} \code{\link{jSDM_poisson_log}} 
#' @examples
#' #==============================================
#' # jSDM_binomial_probit_block_long_format()
#' # Example with simulated data
#' #==============================================
#' 
#' #=================
#' #== Load libraries
#' library(jSDM)
#' 
#' #==================
#' #== Data simulation
#' 
#' #= Number of sites
#' nsite <- 50
#' 
#' #= Set seed for repeatability
#' seed <- 1234
#' set.seed(seed)
#' 
#' #' #= Number of species
#' nsp <- 25
#' 
#' #= Number of latent variables
#' n_latent <- 2
#' #'
#' # Ecological process (suitability)
#' ## X
#' x1 <- rnorm(nsite,0,1)
#' x1.2 <- scale(x1^2)
#' X <- cbind(rep(1,nsite),x1,x1.2)
#' colnames(X) <- c("Int","x1","x1.2")
#' np <- ncol(X)
#' ## W
#' W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
#' ## D
#' SLA <- runif(nsp,-1,1)
#' D <- data.frame(Int=1, x1=x1, x1.2=x1.2, x1.SLA= scale(c(x1 %*% t(SLA))))
#' nd <- ncol(D)
#' ## parameters
#' beta_sp.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
#' mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
#' diag(mat) <- runif(n_latent,0,2)
#' lambda.target <- matrix(0,n_latent,nsp)
#' beta.target <-runif(nd,-1,1)
#' # constraintq of identifiability
#' beta_sp.target[,1] <- 0.0
#' lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
#' #= Variance of random site effect 
#' V_alpha.target <- 0.5
#' #= Random site effect 
#' alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
#' ## probit_theta
#' probit_theta <- c(X %*% beta_sp.target) + c(W %*% lambda.target)
#'                 + as.matrix(D) %*% beta.target + rep(alpha.target,nsp)
#' # Supplementary observation (each site have been visited twice)
#' # Environmental variables at the time of the second visit
#' x1_supObs <- rnorm(nsite,0,1)
#' x1.2_supObs <- scale(x1^2)
#' X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
#' D_supObs <- data.frame(Int=1, x1=x1_supObs, x1.2=x1.2_supObs, x1.SLA=scale(c(x1_supObs %*% t(SLA))))
#' probit_theta_supObs <- c(X_supObs%*%beta_sp.target) + c(W%*%lambda.target)
#'                        + as.matrix(D_supObs) %*% beta.target + alpha.target
#' probit_theta <- c(probit_theta, probit_theta_supObs)
#' nobs <- length(probit_theta)
#' e <- rnorm(nobs,0,1)
#' Z_true <- probit_theta + e
#' Y<-rep(0,nobs)
#' for (n in 1:nobs){
#'   if ( Z_true[n] > 0) {Y[n] <- 1}
#' }
#' Id_site <- rep(1:nsite,nsp)
#' Id_sp <- rep(1:nsp,each=nsite)
#' data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
#'                    x1=c(D$x1,D_supObs$x1), x1.2=c(D$x1.2,D_supObs$x1.2),
#'                    x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
#' # missing observation
#' data <- data[-1,]
#' 
#' #==================================
#' #== Site-occupancy model
#' 
#' # Increase number of iterations (burnin and mcmc) to get convergence
#' mod<-jSDM_binomial_probit_block_long_format( # Iteration
#'   burnin=500,
#'   mcmc=500,
#'   thin=1,
#'   # Response variable
#'   data=data,
#'   # Explanatory variables
#'   site_suitability=~ (1 + x1 + x1.2)*species + x1.SLA + site,
#'   n_latent=2,
#'   site_effect="random",
#'   # Starting values
#'   alpha_start=0,
#'   beta_start=0,
#'   beta_sp_start=0,
#'   lambda_start=0,
#'   W_start=0,
#'   V_alpha=1,
#'   # Priors
#'   shape=0.5, rate=0.0005,
#'   mu_beta_sp=0, V_beta_sp=10,
#'   mu_beta=0, V_beta=1,
#'   mu_lambda=0, V_lambda=10,
#'   seed=1234, verbose=1)
#' 
#' #= Parameter estimates
#' 
#' # beta 
#' par(mfrow=c(nd,2))
#' for(d in 1:nd){
#'  coda::traceplot(mod$mcmc.beta[,d])
#'  coda::densplot(mod$mcmc.beta[,d],
#'                 main = colnames(mod$mcmc.beta)[d])
#'  abline(v=beta.target[d],col='red')
#' }
#' ## beta_j
#' # summary(mod$mcmc.sp$sp_1[,1:ncol(X)])
#' mean_beta_sp <- matrix(0,nsp,ncol(X))
#' pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_probit_block.pdf"))
#' par(mfrow=c(ncol(X),2))
#' for (j in 1:nsp) {
#'   mean_beta_sp[j,] <- apply(mod$mcmc.sp[[paste0("sp_",j)]]
#'                          [,1:ncol(X)], 2, mean)
#'   for (p in 1:ncol(X)){
#'     coda::traceplot(coda::as.mcmc(
#'       mod$mcmc.sp[[paste0("sp_",j)]][,p]))
#'     coda::densplot(coda::as.mcmc(
#'       mod$mcmc.sp[[paste0("sp_",j)]][,p]),
#'       main = paste(colnames(
#'         mod$mcmc.sp[[paste0("sp_",j)]])[p],
#'         ", species : ",j))
#'     abline(v=beta_sp.target[p,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#' ## lambda_j
#' # summary(mod$mcmc.sp$sp_1[,(ncol(X)+1):(ncol(X)+n_latent)])
#' # summary(mod$mcmc.sp$sp_2[,(ncol(X)+1):(ncol(X)+n_latent)])
#' mean_lambda <- matrix(0,nsp,n_latent)
#' pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_probit_block.pdf"))
#' par(mfrow=c(n_latent*2,2))
#' for (j in 1:nsp) {
#'   mean_lambda[j,] <- apply(mod$mcmc.sp[[paste0("sp_",j)]]
#'                            [,(ncol(X)+1):(ncol(X)+n_latent)], 2, mean)
#'   for (l in 1:n_latent) {
#'     coda::traceplot(coda::as.mcmc(mod$mcmc.sp[[paste0("sp_",j)]]
#'                                   [,ncol(X)+l]))
#'     coda::densplot(coda::as.mcmc(mod$mcmc.sp[[paste0("sp_",j)]]
#'                                  [,ncol(X)+l]),
#'                    main=paste(colnames(mod$mcmc.sp[[paste0("sp_",j)]])
#'                               [ncol(X)+l],", species : ",j))
#'     abline(v=lambda.target[l,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#' # Species effects beta and factor loadings lambda
#' par(mfrow=c(1,2))
#' plot(t(beta_sp.target), mean_beta_sp,
#'      main="species effect beta",
#'      xlab ="obs", ylab ="fitted")
#' abline(a=0,b=1,col='red')
#' plot(t(lambda.target), mean_lambda,
#'      main="factor loadings lambda",
#'      xlab ="obs", ylab ="fitted")
#' abline(a=0,b=1,col='red')
#' 
#' ## W latent variables
#' par(mfrow=c(1,2))
#' for (l in 1:n_latent) {
#'   plot(W[,l],
#'        summary(mod$mcmc.latent[[paste0("lv_",l)]])[[1]][,"Mean"],
#'        main = paste0("Latent variable W_", l),
#'        xlab ="obs", ylab ="fitted")
#'   abline(a=0,b=1,col='red')
#' }
#' 
#' ## alpha
#' # summary(mod$mcmc.alpha)
#' par(mfrow=c(1,3))
#' plot(alpha.target, summary(mod$mcmc.alpha)[[1]][,"Mean"],
#'      xlab ="obs", ylab ="fitted", main="site effect alpha")
#' abline(a=0,b=1,col='red')
#' ## Valpha
#' coda::traceplot(mod$mcmc.V_alpha, main="Trace V_alpha")
#' coda::densplot(mod$mcmc.V_alpha,main="Density V_alpha")
#' abline(v=V_alpha.target,col='red')
#' 
#' ## Deviance
#' #summary(mod$mcmc.Deviance)
#' plot(mod$mcmc.Deviance)
#' 
#' #= Predictions
#' 
#' ## probit_theta
#' # summary(mod$probit_theta_pred)
#' par(mfrow=c(1,2))
#' plot(probit_theta[-1],mod$probit_theta_pred,
#'      main="probit(theta)",xlab="obs",ylab="fitted")
#' abline(a=0,b=1,col='red')
#' 
#' ## Z
#' # summary(mod$Z_latent)
#' plot(Z_true[-1],mod$Z_latent,
#'      main="Z_latent", xlab="obs", ylab="fitted")
#' abline(a=0,b=1,col='red')
#' 
#' @keywords Binomial probit regression biodiversity JSDM hierarchical Bayesian models MCMC Markov Chains Monte Carlo Gibbs Sampling
#' @export 

jSDM_binomial_probit_block_long_format <- function(data, site_suitability, n_latent=0,
                                                          site_effect="none",
                                                          burnin=5000, mcmc=10000, thin=10,
                                                          alpha_start=0, beta_start=0,
                                                          beta_sp_start=0,
                                                          lambda_start=0, W_start=0,
                                                          V_alpha=1,
                                                          shape=0.5, rate=0.0005,
                                                          mu_beta=0, V_beta=1.0E6,
                                                          mu_beta_sp=0, V_beta_sp=1.0E6,
                                                          mu_lambda=0, V_lambda=10,
                                                          seed=1234, verbose=1)
  
{   
  #========
  # Basic checks
  #========
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  
  #===================
  # Defining constants
  #===================
  nsp <- length(unique(data$species))
  nsite <- length(unique(data$site))
  nobs <- nrow(data)
  
  #======== 
  # Form response, covariate matrices and model parameters
  #========
  #= Response
  Y <- as.vector(data$Y)
  Id_sp <- rep(0,nrow(data))
  if(is.numeric(data$species)){
    if(min(data$species)==1 && max(data$species)==nsp){
      Id_sp <- data$species - 1
    }
    else if(min(data$species)==0 && max(data$species)==(nsp-1)){
      Id_sp <- data$species
    }
  }else{
    data$species <- as.character(data$species)
    for (j in 1:nsp) {
      Id_sp[grepl(unique(data$species)[j],data$species)] <- j-1
    }
  }
  Id_site <- rep(0,nrow(data))
  if(is.numeric(data$site)) {
    if(min(data$site)==1 && max(data$site)==nsite){
      Id_site <- data$site - 1
    }
    if(min(data$site)==0 && max(data$site)==(nsite-1)){
      Id_site <- data$site
    }
  }else{
    data$site <- as.character(data$site)
    for (i in 1:nsite) {
      Id_site[grepl(unique(data$site)[i],data$site)] <- i-1
    }
  }
  
  #= Suitability
  suitability <- site_suitability 
  if(site_suitability==~.) suitability <- ~. - site - Y
  mf.suit <- model.frame(formula=suitability, data=data)
  # design matrix X for species effects beta_sp
  Xterms <- stringi::stri_remove_empty(gsub(":?species:?", "", 
                                            grep("species", attr(attr(mf.suit,"terms"),"term.labels"), value=T)))
  Xformula <- paste0("~",paste0(Xterms, collapse="+"))
  mf.suit.X <- model.frame(formula=Xformula, data=data)
  attr(attr(mf.suit.X,"terms"),"intercept") <- ifelse(grepl("- *species", suitability[2]),0,1)
  X <- model.matrix(attr(mf.suit.X,"terms"), data=mf.suit.X)
  np <- ncol(X)
  # design matrix D for parameters beta 
  Dterms <- grep("species", grep("site", attr(attr(mf.suit,"terms"),"term.labels"), value=T, invert=T), value=T, invert=T)
  if(length(Dterms)!=0){
    Dformula <- paste0("~", paste0(Dterms, collapse="+"))
    mf.suit.D <- model.frame(formula=Dformula, data=data)
    attr(attr(mf.suit.D,"terms"),"intercept") <- attr(attr(mf.suit,"terms"),"intercept")
    D <- model.matrix(attr(mf.suit.D,"terms"), data=mf.suit.D)
    nd <- ncol(D)
    # Common covariables between X and D
    Id_common_var <- rep(0,ncol(X))
    for(p in 1:ncol(X)){
      if(length(which(colnames(D)==colnames(X)[p]))!=0){
      Id_common_var[p] <- which(colnames(D)==colnames(X)[p])
      }
    }
    if(sum(Id_common_var)!=0){
    Id_common_var <- Id_common_var[Id_common_var>0] -1
    } else {
    Id_common_var <- -1
    }
  }
  if(site_suitability==~. - site - Y) site_suitability <- ~.
  
  if(nrow(unique(X))==nsite){
    n_visit <- rep(1,nobs)
  } else {
    n_visit <- rep(0,nobs)
    for (n in 1:nobs){
      Id_obs <- which(Id_sp[n]==Id_sp & Id_site[n]==Id_site)
      n_visit[Id_obs] <- 1:length(Id_obs)
    }
  }
  
  #= Iterations
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin
  
  #========== 
  # Check data
  #==========
  check.T.binomial(n_visit, nobs)
  check.Y.binomial(Y, n_visit)
  check.X(X, nobs)
  if(length(Dterms)!=0){
    check.X(D, nobs)
  }
  
  if(n_latent==0 && site_effect=="none"){
    
    if(length(Dterms)==0){
      #========
      # Initial starting values for M-H
      #========
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      
      #========
      # Form and check priors
      #========
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                         Y=Y, X=as.matrix(X),
                                                         Id_sp=Id_sp, Id_site=Id_site,
                                                         beta_start=beta_sp_start,
                                                         V_beta=Vbeta_sp, mu_beta = mubeta_sp,
                                                         seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         data=data, n_latent=n_latent, 
                         site_suitability=site_suitability,
                         beta_sp_start=beta_sp_start, 
                         mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                         site_effect=site_effect,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.sp = MCMC.sp, 
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
    } else{
      #========
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                beta_start=beta_start,
                                                                V_beta=Vbeta, mu_beta=mubeta,
                                                                beta_sp_start=beta_sp_start,
                                                                V_beta_sp=Vbeta_sp, mu_beta_sp=mubeta_sp,
                                                                seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_spj
        MCMC.beta_spj <- coda::mcmc(mod$beta_sp[,j,], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp.",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         data=data, n_latent=n_latent, 
                         site_suitability=site_suitability,
                         beta_start=beta_start, 
                         beta_sp_start=beta_sp_start, 
                         mu_beta=mubeta, V_beta=Vbeta,
                         mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                         site_effect=site_effect,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta = MCMC.beta, 
                     mcmc.sp = MCMC.sp, 
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
      
    }
  }
  
  if(n_latent>0 && site_effect=="none"){
    
    if (nsp==1) {
      cat("Error: Unable to adjust latent variables from data about only one species.\n n_latent must be equal to 0 with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    
    if(length(Dterms)==0){
      #========
      # Initial starting values for M-H
      #========
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      param_start = rbind(beta_sp_start,lambda_start)
      
      #========
      # Form and check priors
      #========
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      Vparam <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
      muparam <- c(mubeta_sp,mulambda)
      V_W <- diag(rep(1,n_latent))
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                            Y=Y, X=as.matrix(X),
                                                            Id_sp=Id_sp, Id_site=Id_site,
                                                            param_start= param_start,
                                                            V_param=Vparam, mu_param=muparam,
                                                            W_start=W_start, V_W=V_W,
                                                            seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$param[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
        
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         site_effect=site_effect,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
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
    } else {
      #========
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      param_sp_start = rbind(beta_sp_start,lambda_start)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      Vparam_sp <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
      muparam_sp <- c(mubeta_sp,mulambda)
      V_W <- diag(rep(1,n_latent))
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                   Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                   Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                   beta_start=beta_start,
                                                                   V_beta=Vbeta, mu_beta=mubeta,
                                                                   param_sp_start= param_sp_start,
                                                                   V_param_sp=Vparam_sp, mu_param_sp=muparam_sp,
                                                                   W_start=W_start, V_W=V_W,
                                                                   seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_sp
        MCMC.beta_spj <- coda::mcmc(mod$param_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp.",j)
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$param_sp[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
        
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         site_effect=site_effect,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta=MCMC.beta,
                     mcmc.sp=MCMC.sp, mcmc.latent=MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
    }
  }
  
  if(n_latent==0 && site_effect=="fixed"){
    
    if (nsp==1) {
      cat("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to none with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    
    if(length(Dterms)==0){
      #========
      # Initial starting values for M-H
      #========
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_fixed_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                    Y=Y, X=as.matrix(X),
                                                                    Id_sp=Id_sp, Id_site=Id_site,
                                                                    beta_start=beta_sp_start, V_beta=Vbeta_sp, mu_beta=mubeta_sp,
                                                                    alpha_start=alpha_start, V_alpha=V_alpha,
                                                                    seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_sp_start=beta_sp_start, alpha_start=alpha_start,
                         V_alpha=V_alpha, site_effect=site_effect,
                         mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
    } else {
      #========
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_fixed_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                           Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                           Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                           beta_start=beta_start,
                                                                           V_beta=Vbeta, mu_beta=mubeta,
                                                                           beta_sp_start=beta_sp_start,
                                                                           V_beta_sp=Vbeta_sp, mu_beta_sp=mubeta_sp,
                                                                           alpha_start=alpha_start, V_alpha=V_alpha,
                                                                           seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",1:nsite)
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        MCMC.beta_spj <- coda::mcmc(mod$beta_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp.",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start,
                         beta_sp_start=beta_sp_start, alpha_start=alpha_start,
                         V_alpha=V_alpha, site_effect=site_effect,
                         mu_beta=mubeta, V_beta=Vbeta, 
                         mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta=MCMC.beta,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
      
    }
  }
  
  if(n_latent==0 && site_effect=="random"){
    
    if (nsp==1) {
      cat("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to none with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    
    if(length(Dterms)==0){
      #========
      # Initial starting values for M-H
      #========
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_rand_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                   Y=Y, X=as.matrix(X),
                                                                   Id_sp=Id_sp, Id_site=Id_site,
                                                                   beta_start=beta_sp_start, V_beta=Vbeta_sp, mu_beta=mubeta_sp,
                                                                   alpha_start=alpha_start, V_alpha_start=V_alpha,
                                                                   shape = shape, rate = rate,
                                                                   seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_sp_start=beta_sp_start, alpha_start=alpha_start,
                         V_alpha_start=V_alpha, shape=shape, rate=rate,
                         site_effect=site_effect, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
    } else {
      #========
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_rand_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                          Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                          Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                          beta_start=beta_start,
                                                                          V_beta=Vbeta, mu_beta=mubeta,
                                                                          beta_sp_start=beta_sp_start, 
                                                                          V_beta_sp=Vbeta_sp, mu_beta_sp=mubeta_sp,
                                                                          alpha_start=alpha_start, V_alpha_start=V_alpha,
                                                                          shape = shape, rate = rate,
                                                                          seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        MCMC.beta_spj <- coda::mcmc(mod$beta_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp.",j)
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(MCMC.beta_spj,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start,
                         beta_sp_start=beta_sp_start, alpha_start=alpha_start,
                         V_alpha_start=V_alpha, shape=shape, rate=rate,
                         site_effect=site_effect,
                         mu_beta=mubeta, V_beta=Vbeta, 
                         mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta=MCMC.beta, 
                     mcmc.alpha=MCMC.alpha,
                     mcmc.V_alpha=MCMC.V_alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
      
    }
  }
  
  
  if(n_latent>0 && site_effect=="random"){
    
    if (nsp==1) {
      cat("Error: Unable to adjust site effect and latent variables from data about only one species.\n site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    
    if(length(Dterms)==0){
    #========
    # Initial starting values for M-H
    #========
    beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
    lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
    alpha_start <- form.alpha.start.sp(alpha_start, nsite)
    W_start <-form.W.start.sp(W_start, nsite, n_latent)
    param_start = rbind(beta_sp_start,lambda_start)
    
    #========
    # Form and check priors
    #========
    mubeta_sp <- check.mubeta(mu_beta_sp,np)
    Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
    mulambda <- check.mubeta(mu_lambda,n_latent)
    Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
    Vparam <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
    muparam <- c(mubeta_sp,mulambda)
    V_W <- diag(rep(1,n_latent))
    V_alpha_start <- check.Valpha(V_alpha)
    
    #========
    # call Rcpp function
    #========
    mod <- Rcpp_jSDM_binomial_probit_block_rand_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                    Y=Y, X=as.matrix(X),
                                                                    Id_sp=Id_sp, Id_site=Id_site,
                                                                    param_start= param_start,
                                                                    V_param=Vparam, mu_param = muparam,
                                                                    W_start=W_start, V_W=V_W,
                                                                    alpha_start=alpha_start,
                                                                    V_alpha_start=V_alpha_start,
                                                                    shape = shape, rate = rate,
                                                                    seed=seed, verbose=verbose)
    
    #= Transform Sample list in an MCMC object
    MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
    colnames(MCMC.Deviance) <- "Deviance"
    MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
    colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
    MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
    colnames(MCMC.V_alpha) <- "V_alpha"
    MCMC.sp <- list()
    for (j in 1:nsp) {
      ## beta_j
      MCMC.beta_spj <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
      colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
      ## lambda_j
      MCMC.lambda_j <- coda::mcmc(mod$param[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
      colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
      
      MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
    }
    ## W latent variables 
    MCMC.latent <- list()
    for (l in 1:n_latent) {
      MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
      MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
    }
    
    #= Model specification, site_suitability,
    model_spec <- list(data=data,
                       site_suitability=site_suitability, n_latent=n_latent,
                       burnin=burnin, mcmc=mcmc, thin=thin,
                       beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                       lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                       alpha_start=alpha_start, V_alpha_start=V_alpha_start, 
                       shape=shape, rate=rate, site_effect=site_effect, W_start=W_start, V_W=V_W,
                       family="binomial", link="probit",
                       seed=seed)
    
    #= Output
    output <- list(mcmc.Deviance=MCMC.Deviance,
                   mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                   mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                   Z_latent=mod$Z_latent, 
                   probit_theta_pred=mod$probit_theta_pred,
                   model_spec=model_spec)
    
    } else {
      #========
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_sp_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      param_sp_start = rbind(beta_sp_start,lambda_start)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      Vparam_sp <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
      muparam_sp <- c(mubeta_sp,mulambda)
      V_W <- diag(rep(1,n_latent))
      V_alpha_start <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_rand_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                      Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                      Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                      beta_start=beta_start,
                                                                      V_beta=Vbeta, mu_beta=mubeta,
                                                                      param_sp_start= param_sp_start,
                                                                      V_param_sp=Vparam_sp, mu_param_sp = muparam_sp,
                                                                      W_start=W_start, V_W=V_W,
                                                                      alpha_start=alpha_start,
                                                                      V_alpha_start=V_alpha_start,
                                                                      shape = shape, rate = rate,
                                                                      seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$param_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$param_sp[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
        
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha_start=V_alpha_start, 
                         shape=shape, rate=rate, site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta=MCMC.beta,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec)
      
    }
  }
  
  if(n_latent>0 && site_effect=="fixed"){
    
    if (nsp==1) {
      cat("Error: Unable to adjust site effect and latent variables from data about only one species.\n
        site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    
    if(length(Dterms)==0){
    #=======
    # Initial starting values for M-H
    #========
    beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
    lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
    alpha_start <- form.alpha.start.sp(alpha_start, nsite)
    W_start <-form.W.start.sp(W_start, nsite, n_latent)
    param_start = rbind(beta_sp_start,lambda_start)
    
    #========
    # Form and check priors
    #========
    mubeta_sp <- check.mubeta(mu_beta_sp,np)
    Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
    mulambda <- check.mulambda(mu_lambda,n_latent)
    Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
    Vparam <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
    muparam <- c(mubeta_sp,mulambda)
    V_W <- diag(rep(1,n_latent))
    V_alpha <- check.Valpha(V_alpha)
    
    
    #========
    # call Rcpp function
    #========
    mod <- Rcpp_jSDM_binomial_probit_block_fixed_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                     Y=Y, X=as.matrix(X),
                                                                     Id_sp=Id_sp, Id_site=Id_site,
                                                                     param_start= param_start,
                                                                     V_param=Vparam, mu_param = muparam,
                                                                     W_start=W_start, V_W=V_W,
                                                                     alpha_start=alpha_start, V_alpha=V_alpha,
                                                                     seed=seed, verbose=verbose)
    
    #= Transform Sample list in an MCMC object
    MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
    colnames(MCMC.Deviance) <- "Deviance"
    MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
    colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
    MCMC.sp <- list()
    for (j in 1:nsp) {
      ## beta_j
      MCMC.beta_spj <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
      colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
      ## lambda_j
      MCMC.lambda_j <- coda::mcmc(mod$param[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
      colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
      
      MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
    }
    ## W latent variables 
    MCMC.latent <- list()
    for (l in 1:n_latent) {
      MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
      MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
    }
    
    #= Model specification, site_suitability,
    model_spec <- list(data=data,
                       site_suitability=site_suitability, n_latent=n_latent,
                       burnin=burnin, mcmc=mcmc, thin=thin,
                       beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                       lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                       alpha_start=alpha_start, V_alpha=V_alpha, 
                       shape=shape, rate=rate, site_effect=site_effect, W_start=W_start, V_W=V_W,
                       family="binomial", link="probit",
                       seed=seed)
    
    #= Output
    output <- list(mcmc.Deviance=MCMC.Deviance,
                   mcmc.alpha = MCMC.alpha,
                   mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                   Z_latent=mod$Z_latent, 
                   probit_theta_pred=mod$probit_theta_pred,
                   model_spec=model_spec)
    } else {
      #=======
      # Initial starting values for M-H
      #========
      beta_start <- form.beta.start(beta_start, nd)
      beta_sp_start <- form.beta.start.sp(beta_sp_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      param_sp_start = rbind(beta_sp_start,lambda_start)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,nd)
      Vbeta <- check.Vbeta.mat(V_beta,nd)
      mubeta_sp <- check.mubeta(mu_beta_sp,np)
      Vbeta_sp <- check.Vbeta.mat(V_beta_sp,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      Vparam_sp <- diag(c(diag(Vbeta_sp),diag(Vlambda)))
      muparam_sp <- c(mubeta_sp,mulambda)
      V_W <- diag(rep(1,n_latent))
      V_alpha <- check.Valpha(V_alpha)
      
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_block_traits_fixed_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                       Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                       Id_sp=Id_sp, Id_site=Id_site, Id_common_var,
                                                                       beta_start=beta_start,
                                                                       V_beta=Vbeta, mu_beta=mubeta,
                                                                       param_sp_start= param_sp_start,
                                                                       V_param_sp=Vparam_sp, mu_param_sp = muparam_sp,
                                                                       W_start=W_start, V_W=V_W,
                                                                       alpha_start=alpha_start, V_alpha=V_alpha,
                                                                       seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",unique(data$site))
      MCMC.beta <- coda::mcmc(mod$beta,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.beta) <- paste0("beta_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_spj <- coda::mcmc(mod$param_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_spj) <- paste0("beta_",colnames(X),"_sp",j)
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$param_sp[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent,"_sp",j)
        
        MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_spj, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_suitability,
      model_spec <- list(data=data,
                         site_suitability=site_suitability, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         beta_sp_start=beta_sp_start, mu_beta_sp=mubeta_sp, V_beta_sp=Vbeta_sp,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha=V_alpha, 
                         shape=shape, rate=rate, site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.beta=MCMC.beta,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_pred=mod$probit_theta_pred,
                     model_spec=model_spec) 
    }
  }
  
  class(output) <- "jSDM"
  # return S3 object output belonging to class jSDM
  # acting like list
  return(output) 
}

# End
