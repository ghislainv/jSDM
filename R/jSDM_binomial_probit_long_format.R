## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name jSDM_binomial_probit_long_format
#' @aliases jSDM_binomial_probit_long_format
#' @title Binomial probit regression on long format data 
#' @description The \code{jSDM_binomial_probit_long_format} function performs a Binomial probit regression in a Bayesian framework. 
#' The function calls a Gibbs sampler written in C++ code which uses conjugate priors to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}.\code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param data A \code{data.frame} with at least the following columns : \tabular{ll}{
#'  \code{Y} \tab \eqn{n_{obs}}{n_obs}-length vector indicating the presence by a 1 (or absence by a 0), \cr
#'  \tab of the species observed during each visit of the sites.\cr 
#'  \code{site} \tab numeric or character \eqn{n_{obs}}{n_obs}-length vector indicating the visited site, \cr
#'  \tab  (sites can be visited several times).\cr
#'  \code{species} \tab numeric or character eqn{n_{obs}}{n_obs}-length vector indicating the species observed, \cr 
#'  \tab (species may not have been recorded at all sites).\cr 
#'  \code{x1,...,xp} \tab columns of explanatory variables for the suitability process of the model.\cr}
#' @param site_formula A one-sided formula, as the formulas used by the \code{\link[stats]{lm}} function, of the form: '~ x1 + ... + xd + species:x1 + ... + species:xp' with \eqn{p} terms related to species effects \eqn{\beta},
#' specifying the explanatory variables for the suitability process of the model, including the intercept, different from the \eqn{d} terms related to \eqn{\gamma} parameters.
#' @param n_latent An integer which specifies the number of latent variables to generate. Defaults to \eqn{0}.
#' @param site_effect A string indicating whether row effects are included as fixed effects (\code{"fixed"}), as random effects (\code{"random"}), or not included (\code{"none"}) in the model. 
#'  If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used.
#'  If random effects, they are drawn from a normal distribution with mean zero and unknown variance, analogous to a random intercept in mixed models. Defaults to \code{"none"}.
#' @param gamma_start Starting values for gamma parameters of the suitability process must be either a scalar or a \eqn{d}-length vector. If \code{gamma_start} takes a scalar value, then that value will serve for all of the \eqn{\gamma} parameters.
#' @param beta_start Starting values for beta parameters of the suitability process for each species must be either a scalar or a \eqn{p \times n_{species}}{p x n_species} matrix. If \code{beta_start} takes a scalar value, then that value will serve for all of the \eqn{\beta} parameters.
#' @param lambda_start Starting values for lambda parameters corresponding to the latent variables for each species must be either a scalar or a \eqn{n_{latent} \times n_{species}}{n_latent x n_species} upper triangular matrix with strictly positive values on the diagonal, ignored if \code{n_latent=0}.
#'  If \code{lambda_start} takes a scalar value, then that value will serve for all of the \eqn{\lambda} parameters except those concerned by the constraints explained above.
#' @param W_start Starting values for latent variables must be either a scalar or a \eqn{nsite \times n_latent}{n_site x n_latent} matrix, ignored if \code{n_latent=0}.
#'  If \code{W_start} takes a scalar value, then that value will serve for all of the \eqn{W_{il}}{W_il} with \eqn{i=1,\ldots,n_{site}}{l=1,...,n_site} and \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}.
#' @param alpha_start Starting values for random site effect parameters must be either a scalar or a \eqn{n_{site}}{n_site}-length vector, ignored if \code{site_effect="none"}.
#'  If \code{alpha_start} takes a scalar value, then that value will serve for all of the \eqn{\alpha} parameters.
#' @param V_alpha Starting value for variance of random site effect if \code{site_effect="random"}  or constant variance of the Gaussian prior distribution for the fixed site effect if
#'  \code{site_effect="fixed"}. Must be a strictly positive scalar, ignored if \code{site_effect="none"}.
#' @param shape_Valpha Shape parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}. 
#' Must be a strictly positive scalar. Default to 0.5 for weak informative prior.
#' @param rate_Valpha Rate parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}
#' Must be a strictly positive scalar. Default to 0.0005 for weak informative prior.
#' @param mu_gamma Means of the Normal priors for the \eqn{\gamma} parameters of the suitability process. \code{mu_gamma} must be either a scalar or a \eqn{d}-length vector.
#'  If \code{mu_gamma} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\gamma} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_gamma Variances of the Normal priors for the \eqn{\gamma} parameters of the suitability process. \code{V_gamma} must be either a scalar or a \eqn{d \times d}{d x d} symmetric positive semi-definite square matrix.
#'  If \code{V_gamma} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\gamma} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. 
#'  The default variance is large and set to \code{1e+06} for an uninformative flat prior.
#' @param mu_beta Means of the Normal priors for the \eqn{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a \eqn{p}-length vector.
#'  If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\beta} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta} parameters of the suitability process. \code{V_beta} must be either a scalar or a \eqn{p \times p}{p x p} symmetric positive semi-definite square matrix.
#'  If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\beta} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. 
#' The default variance is large and set to \code{1e+06} for an uninformative flat prior.
#' @param mu_lambda Means of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. 
#' \code{mu_lambda} must be either a scalar or a \eqn{n_{latent}}{n_latent}-length vector. 
#' If \code{mu_lambda} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\lambda}{\lambda} parameters. The default value is set to 0 for an uninformative prior.
#' @param V_lambda Variances of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. \code{V_lambda} must be either a scalar or a \eqn{n_{latent} \times n_{latent}}{n_latent x n_latent} symmetric positive semi-definite square matrix.
#' If \code{V_lambda} takes a scalar value, then that value will serve as the prior variance for all of \eqn{\lambda}{\lambda} parameters, so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal.
#' The default variance is large and set to 10 for an uninformative flat prior.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including :
#'  \item{mcmc.alpha}{An mcmc object that contains the posterior samples for site effects \eqn{\alpha_i}, not returned if \code{site_effect="none"}.}
#'  \item{mcmc.V_alpha}{An mcmc object that contains the posterior samples for variance of random site effect, not returned if \code{site_effect="none"} or \code{site_effect="fixed"}.}
#'  \item{mcmc.latent}{A list by latent variable of mcmc objects that contains the posterior samples for latent variables  \eqn{W_l} with \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}, not returned if \code{n_latent=0}.}
#'  \item{mcmc.sp}{A list by species of mcmc objects that contains the posterior samples for species effects \eqn{\beta} and the loading factors \eqn{\lambda} if \code{n_latent>0}.}
#'  \item{mcmc.gamma}{An mcmc objects that contains the posterior samples for parameters \eqn{\gamma} not returned if \code{d=0}.}
#'  \item{mcmc.Deviance}{The posterior sample of the deviance \eqn{D} is also provided, with \eqn{D} defined as:\eqn{D=-2\log(\prod_{n} P(y_{n}|\beta_j,\lambda_j, \alpha_i, W_i))}{D=-2log(\prod_n P(y_n|\beta_j,\lambda_j, \alpha_i, W_i))}.} 
#'  \item{Z_latent}{Predictive posterior mean of the latent variable Z. }
#'  \item{probit_theta_latent}{Predictive posterior mean of the probability to each species to be present on each site, transformed by probit link function.}
#'  \item{theta_latent}{Predictive posterior mean of the probability to each species to be present on each site.}
#'  \item{model_spec}{Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function, family and number of latent variables, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin.}
#' The \code{mcmc.} objects can be summarized by functions provided by the \code{coda} package. 
#' @details We model an ecological process where the presence or absence of species \eqn{j} on site \eqn{i} is explained by habitat suitability.
#'
#' \bold{Ecological process:}
#' \deqn{y_n \sim \mathcal{B}ernoulli(\theta_n)}{y_n ~ Bernoulli(\theta_n)} such as \eqn{species_n=j} and \eqn{site_n=i},
#'  where :
#'   \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab probit\eqn{(\theta_n) = D_n \gamma+ X_n \beta_j + W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j  + \alpha_i}  and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j + W_i \lambda_j + \alpha_i} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab probit\eqn{(\theta_n) = D_n \gamma  + X_n \beta_j  + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j + W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} 
#' }
#' 
#' @references 
#' Chib, S. and Greenberg, E. (1998) Analysis of multivariate probit models. \emph{Biometrika}, 85, 347-361. 
#'
#' Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.
#' 
#' @author 
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr> 
#' 
#' Jeanne Clément <jeanne.clement16@laposte.net> 
#' 
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}} \code{\link{jSDM_binomial_probit}}
#'  \code{\link{jSDM_binomial_logit}} \code{\link{jSDM_poisson_log}} 
#' @importFrom stringi stri_remove_empty
#' @examples
#' #==============================================
#' # jSDM_binomial_probit_long_format()
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
#' D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
#' nd <- ncol(D)
#' ## parameters
#' beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
#' mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
#' diag(mat) <- runif(n_latent,0,2)
#' lambda.target <- matrix(0,n_latent,nsp)
#' gamma.target <-runif(nd,-1,1)
#' lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat,
#'                                                          diag=TRUE)]
#' #= Variance of random site effect 
#' V_alpha.target <- 0.5
#' #= Random site effect 
#' alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
#' ## probit_theta
#' probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) +
#'                 as.matrix(D) %*% gamma.target + rep(alpha.target, nsp)
#' # Supplementary observation (each site have been visited twice)
#' # Environmental variables at the time of the second visit
#' x1_supObs <- rnorm(nsite,0,1)
#' x1.2_supObs <- scale(x1^2)
#' X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
#' D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
#' probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + 
#'                        as.matrix(D_supObs) %*% gamma.target + alpha.target
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
#'                    x1=c(rep(x1,nsp),rep(x1_supObs,nsp)),
#'                    x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
#'                    x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
#' # missing observation
#' data <- data[-1,]
#' 
#' #==================================
#' #== Site-occupancy model
#' 
#' # Increase number of iterations (burnin and mcmc) to get convergence
#' mod<-jSDM_binomial_probit_long_format( # Iteration
#'   burnin=500,
#'   mcmc=500,
#'   thin=1,
#'   # Response variable
#'   data=data,
#'   # Explanatory variables
#'   site_formula=~ (x1 + x1.2):species + x1.SLA,
#'   n_latent=2,
#'   site_effect="random",
#'   # Starting values
#'   alpha_start=0,
#'   gamma_start=0,
#'   beta_start=0,
#'   lambda_start=0,
#'   W_start=0,
#'   V_alpha=1,
#'   # Priors
#'   shape_Valpha=0.5, rate_Valpha=0.0005,
#'   mu_gamma=0, V_gamma=10,
#'   mu_beta=0, V_beta=10,
#'   mu_lambda=0, V_lambda=10,
#'   seed=1234, verbose=1)
#' 
#' #= Parameter estimates
#' 
#' # gamma 
#' par(mfrow=c(2,2))
#' for(d in 1:nd){
#'  coda::traceplot(mod$mcmc.gamma[,d])
#'  coda::densplot(mod$mcmc.gamma[,d],
#'                 main = colnames(mod$mcmc.gamma)[d])
#'  abline(v=gamma.target[d],col='red')
#' }
#' ## beta_j
#' # summary(mod$mcmc.sp$sp_1[,1:ncol(X)])
#' mean_beta <- matrix(0,nsp,ncol(X))
#' pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_probit.pdf"))
#' par(mfrow=c(ncol(X),2))
#' for (j in 1:nsp) {
#'   mean_beta[j,] <- apply(mod$mcmc.sp[[j]]
#'                          [,1:ncol(X)], 2, mean)
#'   for (p in 1:ncol(X)){
#'     coda::traceplot(mod$mcmc.sp[[j]][,p])
#'     coda::densplot(mod$mcmc.sp[[j]][,p],
#'       main = paste0(colnames(mod$mcmc.sp[[j]])[p],"_sp",j))
#'     abline(v=beta.target[p,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#' ## lambda_j
#' # summary(mod$mcmc.sp$sp_1[,(ncol(X)+1):(ncol(X)+n_latent)])
#' # summary(mod$mcmc.sp$sp_2[,(ncol(X)+1):(ncol(X)+n_latent)])
#' mean_lambda <- matrix(0,nsp,n_latent)
#' pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_probit.pdf"))
#' par(mfrow=c(n_latent*2,2))
#' for (j in 1:nsp) {
#'   mean_lambda[j,] <- apply(mod$mcmc.sp[[j]]
#'                            [,(ncol(X)+1):(ncol(X)+n_latent)], 2, mean)
#'   for (l in 1:n_latent) {
#'     coda::traceplot(mod$mcmc.sp[[j]][,ncol(X)+l])
#'     coda::densplot(mod$mcmc.sp[[j]][,ncol(X)+l],
#'                    main=paste0(colnames(mod$mcmc.sp[[j]])[ncol(X)+l], "_sp",j))
#'     abline(v=lambda.target[l,j],col='red')
#'   }
#' }
#' dev.off()
#' 
#' # Species effects beta and factor loadings lambda
#' par(mfrow=c(1,2))
#' plot(t(beta.target), mean_beta,
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
#' # summary(mod$probit_theta_latent)
#' par(mfrow=c(1,2))
#' plot(probit_theta[-1],mod$probit_theta_latent,
#'      main="probit(theta)",xlab="obs",ylab="fitted")
#' abline(a=0,b=1,col='red')
#' 
#' ## Z
#' # summary(mod$Z_latent)
#' plot(Z_true[-1],mod$Z_latent,
#'      main="Z_latent", xlab="obs", ylab="fitted")
#' abline(a=0,b=1,col='red')
#' ## theta
#' # summary(mod$theta_latent)
#' par(mfrow=c(1,1))
#' plot(pnorm(probit_theta[-1]),mod$theta_latent,
#'      main="theta",xlab="obs",ylab="fitted")
#' abline(a=0,b=1,col='red')
#' 
#' @keywords Binomial probit regression biodiversity JSDM hierarchical Bayesian models MCMC Markov Chains Monte Carlo Gibbs Sampling
#' @export 

jSDM_binomial_probit_long_format <- function(#Iteration
                                             burnin=5000, mcmc=10000, thin=10,
                                             #Data and suitability process
                                             data, site_formula, n_latent=0,
                                             site_effect="none",
                                             # Starting values
                                             alpha_start=0, gamma_start=0,
                                             beta_start=0,
                                             lambda_start=0, W_start=0,
                                             V_alpha=1,
                                             # Priors 
                                             shape_Valpha=0.5, rate_Valpha=0.0005,
                                             mu_gamma=0, V_gamma=10,
                                             mu_beta=0, V_beta=10,
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
  }
  if(sum(Id_sp)==0){
    for (j in 1:nsp) {
      Id_sp[grepl(unique(as.character(data$species))[j], 
                  as.character(data$species))] <- j-1
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
  }
  if(sum(Id_site)==0){
    data$site <- as.character(data$site)
    for (i in 1:nsite) {
      Id_site[grepl(unique(data$site)[i],data$site)] <- i-1
    }
  }
  
  #= Suitability
  suitability <- site_formula 
  if(site_formula==~.) suitability <- ~. - site - Y 
  mf.suit <- model.frame(formula=suitability, data=data)
  # design matrix X for species effects beta
  Xterms <- stringi::stri_remove_empty(gsub(":?species:?", "", 
                                            grep("species", attr(attr(mf.suit,"terms"),"term.labels"), value=T)))
  Xformula <- paste0("~",paste0(Xterms, collapse="+"))
  mf.suit.X <- model.frame(formula=Xformula, data=data)
  attr(attr(mf.suit.X,"terms"),"intercept") <- ifelse(grepl("- *species", suitability[2]) | grepl("- *1", suitability[2]),0,1)
  X <- model.matrix(attr(mf.suit.X,"terms"), data=mf.suit.X)
  np <- ncol(X)
  
  if((attr(attr(mf.suit,"terms"),"intercept")==0) & (unique(X[,1])!=1)) {
    cat("Error: The model must include a species intercept to be interpretable.\n")
    stop("Please respecify the site_formula and call ", calling.function(), " again.",
         call.=FALSE)
  }
  # design matrix D for parameters gamma 
  Dterms <- grep("species", grep("site", attr(attr(mf.suit,"terms"),"term.labels"),
                                 value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)
  if(length(Dterms)!=0){
    Dformula <- paste0("~", paste0(Dterms, collapse="+"),"-1")
    mf.suit.D <- model.frame(formula=Dformula, data=data)
    D <- model.matrix(attr(mf.suit.D,"terms"), data=mf.suit.D)
    nd <- ncol(D)
  } else{
    nd <- 0 
  }
  # Common covariables between X and D
  if(nd!=0){
    Id_common_var <- rep(0,ncol(X))
    for(p in 1:ncol(X)){
      if(length(which(colnames(D)==colnames(X)[p]))!=0){
        Id_common_var[p] <- which(colnames(D)==colnames(X)[p])
      }
    }
    if(sum(Id_common_var)!=0){
      cat("Error: the general design matrix and the species effects design matrix should not contain the same variables
        to ensure the identifiability and interpretability of the model.\n")
      stop("Please respecify the site_formula and call ", calling.function(), " again.",
           call.=FALSE)
    }
  }
  if(site_formula==~. - site - Y) site_formula <- ~.
  
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                   Y=Y, X=as.matrix(X),
                                                   Id_sp=Id_sp, Id_site=Id_site,
                                                   beta_start=beta_start,
                                                   V_beta=Vbeta, mu_beta = mubeta,
                                                   seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj,
                                                                         start=nburn+1,
                                                                         end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         data=data, n_latent=n_latent, 
                         site_formula=site_formula,
                         beta_start=beta_start, 
                         mu_beta=mubeta, V_beta=Vbeta,
                         site_effect=site_effect,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.sp = MCMC.sp, 
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    } else{
      #========
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                          Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                          Id_sp=Id_sp, Id_site=Id_site,
                                                          gamma_start=gamma_start,
                                                          V_gamma=Vgamma, mu_gamma=mugamma,
                                                          beta_start=beta_start,
                                                          V_beta=Vbeta, mu_beta=mubeta,
                                                          seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## betaj
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X),"_sp.",j)
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj,start=nburn+1,
                                                                                end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         data=data, n_latent=n_latent, 
                         site_formula=site_formula,
                         gamma_start=gamma_start, 
                         beta_start=beta_start, 
                         mu_gamma=mugamma, V_gamma=Vgamma,
                         mu_beta=mubeta, V_beta=Vbeta,
                         site_effect=site_effect,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma = MCMC.gamma, 
                     mcmc.sp = MCMC.sp, 
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                      Y=Y, X=as.matrix(X),
                                                      Id_sp=Id_sp, Id_site=Id_site,
                                                      beta_start= beta_start,
                                                      V_beta=Vbeta, mu_beta = mubeta,
                                                      lambda_start= lambda_start,
                                                      V_lambda=Vlambda, mu_lambda = mulambda,
                                                      W_start=W_start, V_W=V_W,
                                                      seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j),
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         site_effect=site_effect,
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
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    } else {
      #========
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                             Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                             Id_sp=Id_sp, Id_site=Id_site,
                                                             gamma_start=gamma_start,
                                                             V_gamma=Vgamma, mu_gamma=mugamma,
                                                             beta_start= beta_start,
                                                             V_beta=Vbeta, mu_beta = mubeta,
                                                             lambda_start= lambda_start,
                                                             V_lambda=Vlambda, mu_lambda = mulambda,
                                                             W_start=W_start, V_W=V_W,
                                                             seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j), 
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         site_effect=site_effect,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         gamma_start=gamma_start, mu_gamma=mugamma, V_gamma=Vgamma,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.sp=MCMC.sp, mcmc.latent=MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_fixed_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                              Y=Y, X=as.matrix(X),
                                                              Id_sp=Id_sp, Id_site=Id_site,
                                                              beta_start=beta_start, V_beta=Vbeta, mu_beta=mubeta,
                                                              alpha_start=alpha_start, V_alpha=V_alpha,
                                                              seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj, 
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, alpha_start=alpha_start,
                         V_alpha=V_alpha, site_effect=site_effect,
                         mu_beta=mubeta, V_beta=Vbeta, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    } else {
      #========
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_fixed_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                     Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                     Id_sp=Id_sp, Id_site=Id_site,
                                                                     gamma_start=gamma_start,
                                                                     V_gamma=Vgamma, mu_gamma=mugamma,
                                                                     beta_start=beta_start,
                                                                     V_beta=Vbeta, mu_beta=mubeta,
                                                                     alpha_start=alpha_start, V_alpha=V_alpha,
                                                                     seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,1:np]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj, 
                                                                         start=nburn+1,
                                                                         end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         gamma_start=gamma_start,
                         beta_start=beta_start, alpha_start=alpha_start,
                         V_alpha=V_alpha, site_effect=site_effect,
                         mu_gamma=mugamma, V_gamma=Vgamma, 
                         mu_beta=mubeta, V_beta=Vbeta, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_rand_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                             Y=Y, X=as.matrix(X),
                                                             Id_sp=Id_sp, Id_site=Id_site,
                                                             beta_start=beta_start, V_beta=Vbeta, mu_beta=mubeta,
                                                             alpha_start=alpha_start, V_alpha_start=V_alpha,
                                                             shape_Valpha = shape_Valpha, rate_Valpha = rate_Valpha,
                                                             seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj,
                                                                         start=nburn+1,
                                                                         end=ngibbs, thin=nthin)
      }
      
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, alpha_start=alpha_start,
                         V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, mu_beta=mubeta, V_beta=Vbeta, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    } else {
      #========
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_rand_site_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                    Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                    Id_sp=Id_sp, Id_site=Id_site,
                                                                    gamma_start=gamma_start,
                                                                    V_gamma=Vgamma, mu_gamma=mugamma,
                                                                    beta_start=beta_start, 
                                                                    V_beta=Vbeta, mu_beta=mubeta,
                                                                    alpha_start=alpha_start, V_alpha_start=V_alpha,
                                                                    shape_Valpha = shape_Valpha,
                                                                    rate_Valpha = rate_Valpha,
                                                                    seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,1:np]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X),"_sp.",j)
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(MCMC.betaj, 
                                                                         start=nburn+1,
                                                                         end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula,
                         n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         gamma_start=gamma_start,
                         beta_start=beta_start, alpha_start=alpha_start,
                         V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha,
                         rate_Valpha=rate_Valpha,
                         site_effect=site_effect,
                         mu_gamma=mugamma, V_gamma=Vgamma, 
                         mu_beta=mubeta, V_beta=Vbeta, 
                         family="binomial", link="probit",
                         seed=seed, verbose=verbose)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma=MCMC.gamma, 
                     mcmc.alpha=MCMC.alpha,
                     mcmc.V_alpha=MCMC.V_alpha,
                     mcmc.sp = MCMC.sp,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      V_alpha_start <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_rand_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                Y=Y, X=as.matrix(X),
                                                                Id_sp=Id_sp, Id_site=Id_site,
                                                                beta_start= beta_start,
                                                                V_beta=Vbeta, mu_beta = mubeta,
                                                                lambda_start= lambda_start,
                                                                V_lambda=Vlambda, mu_lambda = mulambda,
                                                                W_start=W_start, V_W=V_W,
                                                                alpha_start=alpha_start,
                                                                V_alpha_start=V_alpha_start,
                                                                shape_Valpha = shape_Valpha,
                                                                rate_Valpha = rate_Valpha,
                                                                seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j), 
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha_start=V_alpha_start, 
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
      
    } else {
      #========
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      V_alpha_start <- check.Valpha(V_alpha)
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_rand_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                       Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                       Id_sp=Id_sp, Id_site=Id_site,
                                                                       gamma_start=gamma_start,
                                                                       V_gamma=Vgamma, mu_gamma=mugamma,
                                                                       beta_start= beta_start,
                                                                       V_beta=Vbeta, mu_beta = mubeta,
                                                                       lambda_start= lambda_start,
                                                                       V_lambda=Vlambda, mu_lambda = mulambda,
                                                                       W_start=W_start, V_W=V_W,
                                                                       alpha_start=alpha_start,
                                                                       V_alpha_start=V_alpha_start,
                                                                       shape_Valpha = shape_Valpha,
                                                                       rate_Valpha = rate_Valpha,
                                                                       seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j),
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         gamma_start=gamma_start, mu_gamma=mugamma, V_gamma=Vgamma,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha_start=V_alpha_start, 
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
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
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      V_alpha <- check.Valpha(V_alpha)
      
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_fixed_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                 Y=Y, X=as.matrix(X),
                                                                 Id_sp=Id_sp, Id_site=Id_site,
                                                                 beta_start= beta_start,
                                                                 V_beta=Vbeta, mu_beta = mubeta,
                                                                 lambda_start= lambda_start,
                                                                 V_lambda=Vlambda, mu_lambda = mulambda,
                                                                 W_start=W_start, V_W=V_W,
                                                                 alpha_start=alpha_start, V_alpha=V_alpha,
                                                                 seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j),
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha=V_alpha, 
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    } else {
      #=======
      # Initial starting values for M-H
      #========
      gamma_start <- form.beta.start(gamma_start, nd)
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      #========
      # Form and check priors
      #========
      mugamma <- check.mubeta(mu_gamma,nd)
      Vgamma <- check.Vbeta.mat(V_gamma,nd)
      mubeta <- check.mubeta(mu_beta,np)
      Vbeta <- check.Vbeta.mat(V_beta,np)
      mulambda <- check.mulambda(mu_lambda,n_latent)
      Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
      V_W <- diag(rep(1,n_latent))
      V_alpha <- check.Valpha(V_alpha)
      
      
      #========
      # call Rcpp function
      #========
      mod <- Rcpp_jSDM_binomial_probit_traits_fixed_site_lv_long_format(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                                        Y=Y, X=as.matrix(X), D=as.matrix(D),
                                                                        Id_sp=Id_sp, Id_site=Id_site,
                                                                        gamma_start=gamma_start,
                                                                        V_gamma=Vgamma, mu_gamma=mugamma,
                                                                        beta_start= beta_start,
                                                                        V_beta=Vbeta, mu_beta = mubeta,
                                                                        lambda_start= lambda_start,
                                                                        V_lambda=Vlambda, mu_lambda = mulambda,
                                                                        W_start=W_start, V_W=V_W,
                                                                        alpha_start=alpha_start, V_alpha=V_alpha,
                                                                        seed=seed, verbose=verbose)
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_", unique(data$site))
      MCMC.gamma <- coda::mcmc(mod$gamma,start=nburn+1,end=ngibbs,thin=nthin)   
      colnames(MCMC.gamma) <- paste0("gamma_",colnames(D))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.betaj <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.betaj) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(as.matrix(mod$lambda[,j,]), start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[ifelse(is.character(data$species),
                        unique(data$species)[j],
                        paste0("sp_",
                               unique(data$species)[j]))]] <- coda::mcmc(cbind(MCMC.betaj, MCMC.lambda_j),
                                                                         start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(data=data,
                         site_formula=site_formula, n_latent=n_latent,
                         burnin=burnin, mcmc=mcmc, thin=thin,
                         gamma_start=gamma_start, mu_gamma=mugamma, V_gamma=Vgamma,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=Vbeta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         alpha_start=alpha_start, V_alpha=V_alpha, 
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, W_start=W_start, V_W=V_W,
                         family="binomial", link="probit",
                         seed=seed)
      
      #= Output
      output <- list(mcmc.Deviance=MCMC.Deviance,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.alpha = MCMC.alpha,
                     mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                     Z_latent=mod$Z_latent, 
                     probit_theta_latent=mod$probit_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec) 
    }
  }
  
  class(output) <- "jSDM"
  # return S3 object output belonging to class jSDM
  # acting like list
  return(output) 
}

# End
