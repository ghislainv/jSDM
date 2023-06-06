## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name jSDM_poisson_log
#' @aliases jSDM_poisson_log 
#' @title Poisson regression with log link function 
#' @description The \code{jSDM_poisson_log} function performs a Poisson regression with log link function in a Bayesian framework. 
#' The function calls a Gibbs sampler written in 'C++' code which uses an adaptive Metropolis algorithm to estimate the conditional posterior distribution of model's parameters.
#' @param burnin The number of burnin iterations for the sampler.
#' @param mcmc The number of Gibbs iterations for the sampler. Total number of Gibbs iterations is equal to \code{burnin+mcmc}.\cr
#'  \code{burnin+mcmc} must be divisible by 10 and superior or equal to 100 so that the progress bar can be displayed.
#' @param thin The thinning interval used in the simulation. The number of mcmc iterations must be divisible by this value.
#' @param count_data  A matrix \eqn{n_{site} \times n_{species}}{n_site x n_species} indicating the abundance of each species on each site.
#' @param site_data A data frame containing the model's explanatory variables by site.
#' @param site_formula A one-sided formula of the form '~x1+...+xp' specifying the explanatory variables for the suitability process of the model,\cr
#' used to form the design matrix \eqn{X} of size \eqn{n_{site} \times np}{n_site x np}.
#' @param trait_data A data frame containing the species traits which can be included as part of the model.\cr
#'  Default to \code{NULL} to fit a model without species traits.
#' @param trait_formula A one-sided formula of the form '~ t1 + ... + tk + x1:t1 + ... + xp:tk' specifying the interactions between the environmental variables and the species traits to be considered in the model,\cr
#' used to form the trait design matrix \eqn{Tr} of size \eqn{n_{species} \times nt}{n_species x nt} \cr 
#' and to set to \eqn{0} the \eqn{\gamma} parameters corresponding to interactions not taken into account according to the formula. 
#' Default to \code{NULL} to fit a model with all possible interactions between species traits found in \code{trait_data} and environmental variables defined by \code{site_formula}.
#' @param n_latent An integer which specifies the number of latent variables to generate. Defaults to \eqn{0}.
#' @param site_effect A string indicating whether row effects are included as fixed effects (\code{"fixed"}), as random effects (\code{"random"}),\cr
#'  or not included (\code{"none"}) in the model. 
#'  If fixed effects, then for parameter identifiability the first row effect is set to zero, which analogous to acting as a reference level when dummy variables are used.
#'  If random effects, they are drawn from a normal distribution with mean zero and unknown variance, analogous to a random intercept in mixed models. Defaults to \code{"none"}.
#' @param beta_start Starting values for \eqn{\beta} parameters of the suitability process for each species must be either a scalar or a \eqn{np \times n_{species}}{np x n_species} matrix. 
#' If \code{beta_start} takes a scalar value, then that value will serve for all of the \eqn{\beta} parameters.
#' @param gamma_start Starting values for \eqn{\gamma} parameters that represent the influence of species-specific traits on species' responses \eqn{\beta},
#' \code{gamma_start} must be either a scalar, a vector of length \eqn{nt}, \eqn{np} or \eqn{nt.np} or a \eqn{nt \times np}{nt x np} matrix.
#'  If \code{gamma_start} takes a scalar value, then that value will serve for all of the \eqn{\gamma} parameters.
#'  If \code{gamma_start} is a vector of length \eqn{nt} or \eqn{nt.np} the resulting \eqn{nt \times np}{nt x np} matrix is filled by column with specified values,
#'   if a \eqn{np}-length vector is given, the matrix is filled by row. 
#' @param lambda_start Starting values for \eqn{\lambda} parameters corresponding to the latent variables for each species must be either a scalar
#'  or a \eqn{n_{latent} \times n_{species}}{n_latent x n_species} upper triangular matrix with strictly positive values on the diagonal,
#'  ignored if \code{n_latent=0}.
#'  If \code{lambda_start} takes a scalar value, then that value will serve for all of the \eqn{\lambda} parameters except those concerned by the constraints explained above.
#' @param W_start Starting values for latent variables must be either a scalar or a \eqn{nsite \times n_latent}{n_site x n_latent} matrix, ignored if \code{n_latent=0}.
#'  If \code{W_start} takes a scalar value, then that value will serve for all of the \eqn{W_{il}}{W_il} with \eqn{i=1,\ldots,n_{site}}{i=1,...,n_site} and \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}.
#' @param alpha_start Starting values for random site effect parameters must be either a scalar or a \eqn{n_{site}}{n_site}-length vector, ignored if \code{site_effect="none"}.
#'  If \code{alpha_start} takes a scalar value, then that value will serve for all of the \eqn{\alpha} parameters.
#' @param V_alpha Starting value for variance of random site effect if \code{site_effect="random"} or constant variance of the Gaussian prior distribution for the fixed site effect if 
#' \code{site_effect="fixed"}. Must be a strictly positive scalar, ignored if \code{site_effect="none"}.
#' @param shape_Valpha Shape parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}. 
#' Must be a strictly positive scalar. Default to 0.5 for weak informative prior.
#' @param rate_Valpha Rate parameter of the Inverse-Gamma prior for the random site effect variance \code{V_alpha}, ignored if \code{site_effect="none"} or \code{site_effect="fixed"}
#' Must be a strictly positive scalar. Default to 0.0005 for weak informative prior.
#' @param mu_beta Means of the Normal priors for the \eqn{\beta}{\beta} parameters of the suitability process. \code{mu_beta} must be either a scalar or a \eqn{np}-length vector.
#'  If \code{mu_beta} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\beta} parameters.
#'   The default value is set to 0 for an uninformative prior, ignored if \code{trait_data} is specified.
#' @param V_beta Variances of the Normal priors for the \eqn{\beta}{\beta} parameters of the suitability process.
#'  \code{V_beta} must be either a scalar or a \eqn{np \times np}{np x np} symmetric positive semi-definite square matrix.
#'  If \code{V_beta} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\beta} parameters,
#'   so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal. 
#' The default variance is large and set to 10 for an uninformative flat prior.
#' @param mu_gamma Means of the Normal priors for the \eqn{\gamma}{\gamma} parameters.
#'  \code{mu_gamma} must be either a scalar, a vector of length \eqn{nt}, \eqn{np} or \eqn{nt.np} or a \eqn{nt \times np}{nt x np} matrix.
#'  If \code{mu_gamma} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\gamma} parameters.
#'  If \code{mu_gamma} is a vector of length \eqn{nt} or \eqn{nt.np} the resulting \eqn{nt \times np}{nt x np} matrix is filled by column with specified values,
#'   if a \eqn{np}-length vector is given, the matrix is filled by row. 
#'  The default value is set to 0 for an uninformative prior, ignored if \code{trait_data=NULL}.
#' @param V_gamma Variances of the Normal priors for the \eqn{\gamma}{\gamma} parameters.
#'  \code{V_gamma} must be either a scalar, a vector of length \eqn{nt}, \eqn{np} or \eqn{nt.np} or a \eqn{nt \times np}{nt x np} positive matrix.
#'  If \code{V_gamma} takes a scalar value, then that value will serve as the prior variance for all of the \eqn{\gamma} parameters.
#'  If \code{V_gamma} is a vector of length \eqn{nt} or \eqn{nt.np} the resulting \eqn{nt \times np}{nt x np} matrix is filled by column with specified values,
#'   if a \eqn{np}-length vector is given, the matrix is filled by row. 
#'  The default variance is large and set to 10 for an uninformative flat prior, ignored if \code{trait_data=NULL}.
#' @param mu_lambda Means of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables. 
#' \code{mu_lambda} must be either a scalar or a \eqn{n_{latent}}{n_latent}-length vector. 
#' If \code{mu_lambda} takes a scalar value, then that value will serve as the prior mean for all of the \eqn{\lambda}{\lambda} parameters.
#'  The default value is set to 0 for an uninformative prior.
#' @param V_lambda Variances of the Normal priors for the \eqn{\lambda}{\lambda} parameters corresponding to the latent variables.
#'  \code{V_lambda} must be either a scalar or a \eqn{n_{latent} \times n_{latent}}{n_latent x n_latent} symmetric positive semi-definite square matrix.
#' If \code{V_lambda} takes a scalar value, then that value will serve as the prior variance for all of \eqn{\lambda}{\lambda} parameters,
#'  so the variance covariance matrix used in this case is diagonal with the specified value on the diagonal.
#' The default variance is large and set to 10 for an uninformative flat prior.
#' @param ropt Target acceptance rate for the adaptive Metropolis algorithm. Default to 0.44.
#' @param seed The seed for the random number generator. Default to 1234.
#' @param verbose A switch (0,1) which determines whether or not the progress of the sampler is printed to the screen. Default is 1: a progress bar is printed, indicating the step (in \%) reached by the Gibbs sampler.
#' @return An object of class \code{"jSDM"} acting like a list including : 
#'  \item{mcmc.alpha}{An mcmc object that contains the posterior samples for site effects \eqn{\alpha_i}, not returned if \code{site_effect="none"}.}
#'  \item{mcmc.V_alpha}{An mcmc object that contains the posterior samples for variance of random site effect, not returned if \code{site_effect="none"} or  \code{site_effect="fixed"}.}
#'  \item{mcmc.latent}{A list by latent variable of mcmc objects that contains the posterior samples for latent variables \eqn{W_l} with \eqn{l=1,\ldots,n_{latent}}{l=1,...,n_latent}, not returned if \code{n_latent=0}.}
#'  \item{mcmc.sp}{A list by species of mcmc objects that contains the posterior samples for species effects \eqn{\beta_j} and \eqn{\lambda_j} if \code{n_latent>0}.}
#'  \item{mcmc.gamma}{A list by covariates of mcmc objects that contains the posterior samples for \eqn{\gamma_p} parameters with \eqn{p=1,\ldots,np}{p=1,...,np} if \code{trait_data} is specified.}
#'  \item{mcmc.Deviance}{The posterior sample of the deviance (\eqn{D}) is also provided, with \eqn{D} defined as : \eqn{D=-2\log(\prod_{ij} P(y_{ij}|\beta_j,\lambda_j, \alpha_i, W_i))}{D=-2log(\prod_ij P(y_ij|\beta_j,\lambda_j, \alpha_i, W_i))}.} 
#'  \item{log_theta_latent}{Predictive posterior mean of the probability to each species to be present on each site, transformed by log link function.}
#'  \item{theta_latent}{Predictive posterior mean of the probability to each species to be present on each site.}
#'  \item{model_spec}{Various attributes of the model fitted, including the response and model matrix used, distributional assumptions as link function, family and number of latent variables, hyperparameters used in the Bayesian estimation and mcmc, burnin and thin.}
#' The \code{mcmc.} objects can be summarized by functions provided by the \code{coda} package. 
#' 
#' @details We model an ecological process where the presence or absence of species \eqn{j} on site \eqn{i} is explained by habitat suitability.
#'
#' \bold{Ecological process : }
#' \deqn{y_{ij} \sim \mathcal{P}oisson(\theta_{ij})}{y_ij ~ Poisson(\theta_ij),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j}{(\theta_ij) = X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j}{(\theta_ij) =  X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab log\eqn{(\theta_{ij}) = X_i \beta_j  + \alpha_i}{(\theta_ij) = X_i \beta_j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab log\eqn{(\theta_{ij}) = X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) = X_i  \beta_j +  W_i \lambda_j + \alpha_i}  \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab log\eqn{(\theta_{ij}) = X_i \beta_j  + \alpha_i}{(\theta_ij) = X_i \beta_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab log\eqn{(\theta_{ij}) = X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) = X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }
#' 
#' 
#' In the absence of data on species traits (\code{trait_data=NULL}), the effect of species \eqn{j}: \eqn{\beta_j};
#' follows the same \emph{a priori} Gaussian distribution such that \eqn{\beta_j \sim \mathcal{N}_{np}(\mu_{\beta},V_{\beta})}{\beta_j ~ N_np(\mu_\beta,V_\beta)},
#' for each species. 
#' 
#' If species traits data are provided, the effect of species \eqn{j}: \eqn{\beta_j};
#' follows an \emph{a priori} Gaussian distribution such that \eqn{\beta_j \sim \mathcal{N}_{np}(\mu_{\beta_j},V_{\beta})}{\beta_j ~ N_np(\mu_\betaj,V_\beta)},
#' where \eqn{\mu_{\beta_jp} = \sum_{k=1}^{nt} t_{jk}.\gamma_{kp}}{mu_\beta.jp=\Sigma_(k=1,...,nt) t_jk.\gamma_kp}, takes different values for each species.
#'
#' We assume that \eqn{\gamma_{kp} \sim \mathcal{N}(\mu_{\gamma_{kp}},V_{\gamma_{kp}})}{\gamma.kp ~ N(\mu_\gamma.kp,V_\gamma.kp)} as prior distribution.   
#' 
#' We define the matrix \eqn{\gamma=(\gamma_{kp})_{k=1,...,nt}^{p=1,...,np}}{\gamma=(\gamma_kp) for k=1,...,nt and p=1,...,np} such as : 
#' 
#' \tabular{rccccccl}{
#' \tab \strong{\eqn{x_0}} \tab \strong{\eqn{x_1}} \tab... \tab \strong{\eqn{x_p}} \tab ... \tab \strong{\eqn{x_{np}}{x_np}} \tab  \cr 
#' \tab__________\tab________\tab________\tab________\tab________\tab________\tab \cr 
#' \strong{\eqn{t_0}} | \tab  \strong{\eqn{\gamma_{0,0}}{\gamma_(0,0)}} \tab \emph{\eqn{\gamma_{0,1}}{\gamma_(0,1)}} \tab ... \tab \emph{\eqn{\gamma_{0,p}}{\gamma_(0,p)}} \tab... \tab \emph{\eqn{\gamma_{0,np}}{\gamma_(0,np)}} \tab \{ \emph{effect of } \cr
#'   | \tab \strong{intercept} \tab \tab \tab \tab \tab \tab \emph{environmental} \cr 
#'   | \tab \tab \tab \tab \tab \tab \tab \emph{   variables} \cr 
#'   \strong{\eqn{t_1}} | \tab  \eqn{\gamma_{1,0}}{\gamma_(1,0)} \tab \emph{\eqn{\gamma_{1,1}}{\gamma_(1,1)}} \tab ... \tab \emph{\eqn{\gamma_{1,p}}{\gamma_(1,p)}} \tab... \tab \emph{\eqn{\gamma_{1,np}}{\gamma_(1,np)}} \tab \cr
#'   ... | \tab ... \tab ... \tab ... \tab ... \tab ... \tab ... \tab \cr 
#'   \strong{\eqn{t_k}} | \tab \emph{\eqn{\gamma_{k,0}}{\gamma_(k,0)}} \tab \emph{\eqn{\gamma_{k,1}}{\gamma_(k,1)}} \tab... \tab \eqn{\gamma_{k,p}}{\gamma_(k,p)} \tab... \tab \eqn{\gamma_{k,np}}{\gamma_(k,np)} \tab \cr
#'   ... | \tab ... \tab  ... \tab ... \tab ... \tab ... \tab ... \tab \cr 
#'   \strong{\eqn{t_{nt}}{t_nt}} | \tab \emph{\eqn{\gamma_{nt,0}}{\gamma_(nt,0)}} \tab \eqn{\gamma_{nt,1}}{\gamma_(nt,1)}\tab ... \tab \eqn{\gamma_{nt,p}}{\gamma_(nt,p)} \tab... \tab \eqn{\gamma_{nt,np}}{\gamma_(nt,np)} \tab \cr
#'   \tab \emph{average} \tab \tab \tab \tab \tab \tab \cr 
#'   \tab \emph{trait effect} \tab \tab interaction \tab traits \tab environment \tab \tab \cr 
#' }
#'  
#' @importFrom stats na.pass
#' @examples 
#' #==============================================
#' # jSDM_poisson_log()
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
#' #= Number of species
#' nsp <- 10
#' #= Set seed for repeatability
#' seed <- 1234
#' 
#' #= Ecological process (suitability)
#' set.seed(seed)
#' x1 <- rnorm(nsite,0,1)
#' set.seed(2*seed)
#' x2 <- rnorm(nsite,0,1)
#' X <- cbind(rep(1,nsite),x1,x2)
#' np <- ncol(X)
#' set.seed(3*seed)
#' W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
#' n_latent <- ncol(W)
#' l.zero <- 0
#' l.diag <- runif(2,0,1)
#' l.other <- runif(nsp*2-3,-1,1)
#' lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],
#'                           l.diag[2],l.other[-1]),
#'                         byrow=TRUE, nrow=nsp)
#' beta.target <- matrix(runif(nsp*np,-1,1), byrow=TRUE, nrow=nsp)
#' V_alpha.target <- 0.5
#' alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
#' log.theta <- X %*% t(beta.target) + W %*% t(lambda.target) + alpha.target
#' theta <- exp(log.theta)
#' Y <- apply(theta, 2, rpois, n=nsite)
#' 
#' #= Site-occupancy model
#' # Increase number of iterations (burnin and mcmc) to get convergence
#' mod <- jSDM_poisson_log(# Chains
#'                         burnin=200,
#'                         mcmc=200,
#'                         thin=1,
#'                         # Response variable
#'                         count_data=Y,
#'                         # Explanatory variables
#'                         site_formula=~x1+x2,
#'                         site_data=X,
#'                         n_latent=n_latent,
#'                         site_effect="random",
#'                         # Starting values
#'                         beta_start=0,
#'                         lambda_start=0,
#'                         W_start=0,
#'                         alpha_start=0,
#'                         V_alpha=1,
#'                         # Priors
#'                         shape_Valpha=0.5,
#'                         rate_Valpha=0.0005,
#'                         mu_beta=0,
#'                         V_beta=10,
#'                         mu_lambda=0,
#'                         V_lambda=10,
#'                         # Various
#'                         seed=1234,
#'                         ropt=0.44,
#'                         verbose=1)
#' #==========
#' #== Outputs
#' 
#' oldpar <- par(no.readonly = TRUE)
#' 
#' #= Parameter estimates
#' 
#' ## beta_j
#' mean_beta <- matrix(0,nsp,np)
#' pdf(file=file.path(tempdir(), "Posteriors_beta_jSDM_log.pdf"))
#' par(mfrow=c(ncol(X),2))
#' for (j in 1:nsp) {
#'   mean_beta[j,] <- apply(mod$mcmc.sp[[j]][,1:ncol(X)],
#'                          2, mean)
#'   for (p in 1:ncol(X)) {
#'     coda::traceplot(mod$mcmc.sp[[j]][,p])
#'     coda::densplot(mod$mcmc.sp[[j]][,p],
#'       main = paste(colnames(
#'         mod$mcmc.sp[[j]])[p],
#'         ", species : ",j))
#'     abline(v=beta.target[j,p],col='red')
#'   }
#' }
#' dev.off()
#' 
#' ## lambda_j
#' mean_lambda <- matrix(0,nsp,n_latent)
#' pdf(file=file.path(tempdir(), "Posteriors_lambda_jSDM_log.pdf"))
#' par(mfrow=c(n_latent*2,2))
#' for (j in 1:nsp) {
#'   mean_lambda[j,] <- apply(mod$mcmc.sp[[j]]
#'                            [,(ncol(X)+1):(ncol(X)+n_latent)], 2, mean)
#'   for (l in 1:n_latent) {
#'     coda::traceplot(mod$mcmc.sp[[j]][,ncol(X)+l])
#'     coda::densplot(mod$mcmc.sp[[j]][,ncol(X)+l],
#'                    main=paste(colnames(mod$mcmc.sp[[j]])
#'                               [ncol(X)+l],", species : ",j))
#'     abline(v=lambda.target[j,l],col='red')
#'   }
#' }
#' dev.off()
#' 
#' # Species effects beta and factor loadings lambda
#' par(mfrow=c(1,2))
#' plot(beta.target, mean_beta,
#'      main="species effect beta",
#'      xlab ="obs", ylab ="fitted")
#' abline(a=0,b=1,col='red')
#' plot(lambda.target, mean_lambda,
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
#' par(mfrow=c(1,3))
#' plot(alpha.target, summary(mod$mcmc.alpha)[[1]][,"Mean"],
#'      xlab ="obs", ylab ="fitted", main="site effect alpha")
#' abline(a=0,b=1,col='red')
#' ## Valpha
#' coda::traceplot(mod$mcmc.V_alpha)
#' coda::densplot(mod$mcmc.V_alpha)
#' abline(v=V_alpha.target,col='red')
#' 
#' ## Deviance
#' summary(mod$mcmc.Deviance)
#' plot(mod$mcmc.Deviance)
#' 
#' #= Predictions
#' par(mfrow=c(1,2))
#' plot(log.theta, mod$log_theta_latent,
#'      main="log(theta)",
#'      xlab="obs", ylab="fitted")
#' abline(a=0 ,b=1, col="red")
#' plot(theta, mod$theta_latent,
#'      main="Expected abundance theta",
#'      xlab="obs", ylab="fitted")
#' abline(a=0 ,b=1, col="red")
#' par(oldpar)
#' 
#' @references
#' Gelfand, A. E.; Schmidt, A. M.; Wu, S.; Silander, J. A.; Latimer, A. and Rebelo, A. G. (2005) Modelling species diversity through species level hierarchical modelling. \emph{Applied Statistics}, 54, 1-20.
#' 
#' Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006) Building statistical models to analyze species distributions. \emph{Ecological Applications}, 16, 33-50.
#'
#' Ovaskainen, O., Tikhonov, G., Norberg, A., Blanchet, F. G., Duan, L., Dunson, D., Roslin, T. and Abrego, N. (2017) How to make more out of community data? A conceptual framework and its implementation as models and software. \emph{Ecology Letters}, 20, 561-576.
#'
#' @author 
#'  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#'  
#' Jeanne Clément <jeanne.clement16@laposte.net>
#' 
#' @seealso \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}} \code{\link{jSDM_binomial_probit}}  \code{\link{jSDM_binomial_logit}} 
#' @keywords multivariate logistic regression model poisson biodiversity MCMC Metropolis algorithm 
#' @export
#' 

jSDM_poisson_log  <- function(# Iteration
                              burnin=5000, mcmc=10000, thin=10,
                              # Data and suitability process
                              count_data, site_data, site_formula,
                              trait_data=NULL, trait_formula=NULL, 
                              n_latent=0, site_effect="none",
                              # Starting values
                              beta_start=0, 
                              gamma_start=0,
                              lambda_start=0,
                              W_start=0,
                              alpha_start=0, 
                              V_alpha=1,
                              # Priors 
                              shape_Valpha=0.5,
                              rate_Valpha=0.0005,
                              mu_beta=0, V_beta=10,
                              mu_gamma=0, V_gamma=10, 
                              mu_lambda=0, V_lambda=10,
                              # Various 
                              ropt=0.44, seed=1234, verbose=1)

{   
  #==== Basic checks ======
  check.mcmc.parameters(burnin, mcmc, thin)
  check.verbose(verbose)
  
  #====== Form response, covariate matrices and model parameters
  #=======  Response =================
  Y <- as.matrix(count_data)
  nsp <- ncol(Y)
  nsite <- nrow(Y)
  nobs <- nsite*nsp
  if(is.null(colnames(Y))){
    colnames(Y) <- paste0("sp_",1:ncol(Y))
  }
  if(is.null(rownames(Y))){
    rownames(Y) <- 1:nrow(Y)
  }
  #==== Site formula suitability ====
  mf.suit <- model.frame(formula=site_formula, data=as.data.frame(site_data),
                         na.action=na.pass) # X will contain NA's in rows corresponding to site_data.
  X <- model.matrix(attr(mf.suit,"terms"), data=mf.suit)
  np <- ncol(X)
  n_Xint <- sum(sapply(apply(X,2,unique), FUN=function(x){all(x==1)}))
  col_Xint <- which(sapply(apply(X,2,unique), FUN=function(x){all(x==1)}))
  if(n_Xint!=1){
    message("Error: The model must include one species intercept to be interpretable.\n")
    stop("Please respecify the site_formula formula and call ", calling.function(), " again.",
         call.=FALSE)
  }
  
  #===== Iterations =====
  ngibbs <- mcmc+burnin
  nthin <- thin
  nburn <- burnin
  nsamp <- mcmc/thin
  
  #==========  Check data ==========
  check.Y.poisson(Y)
  check.X(as.matrix(X), nsite)
  
  #======== function without traits ========
  if(is.null(trait_data)){
    ##======== without latent variables and site effect ========
    if(n_latent==0 && site_effect=="none"){
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      
      # Form and check priors
      mu_beta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log (ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                    Y=Y, X=as.matrix(X),
                                    beta_start=beta_start, mu_beta = mu_beta, V_beta=V_beta,
                                    ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j, start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data, n_latent=n_latent,
                         beta_start=beta_start, mu_beta=mu_beta, V_beta=V_beta,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, mcmc.Deviance=MCMC.Deviance,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with latent variables ======
    if(n_latent>0 && site_effect=="none"){
      
      if (nsp==1) {
        message("Error: Unable to adjust latent variables from data about only one species.\n n_latent must be equal to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      # Form and check priors
      mubeta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                      Y=Y, X=as.matrix(X),
                                      beta_start=beta_start, mu_beta = mubeta, V_beta=V_beta,
                                      lambda_start=lambda_start, mu_lambda = mulambda, V_lambda=Vlambda,
                                      W_start = W_start, V_W = V_W,
                                      ropt=ropt, seed=seed, verbose=verbose)
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data, n_latent=n_latent,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=V_beta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with random site effect ======
    if(n_latent==0 && site_effect=="random"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to 'none' with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      # Form and check priors
      mubeta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_rand_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                             Y=Y, X=as.matrix(X),
                                             beta_start=beta_start, mu_beta = mubeta, V_beta=V_beta,
                                             alpha_start=alpha_start, V_alpha_start=V_alpha,
                                             shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                                             ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data,  n_latent=n_latent,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=V_beta,
                         alpha_start=alpha_start, V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with fixed site effect ======
    if(n_latent==0 && site_effect=="fixed"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to 'none' with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      
      # Form and check priors
      mubeta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_fixed_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                              Y=Y, X=as.matrix(X),
                                              beta_start=beta_start, mu_beta = mubeta, V_beta=V_beta,
                                              alpha_start=alpha_start, V_alpha=V_alpha, 
                                              ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j,start=nburn+1, end=ngibbs, thin=nthin)
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data,  n_latent=n_latent,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=V_beta,
                         alpha_start=alpha_start, V_alpha=V_alpha, 
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with lv and fixed site effect======
    if(n_latent>0 && site_effect=="fixed"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect and latent variables from data about only one species.\n site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      # Form and check priors
      mubeta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_fixed_site_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                 Y=Y, X=as.matrix(X),
                                                 beta_start=beta_start, mu_beta = mubeta, V_beta=V_beta,
                                                 lambda_start=lambda_start, mu_lambda = mulambda, V_lambda=Vlambda,
                                                 W_start = W_start, V_W = V_W,
                                                 alpha_start=alpha_start, V_alpha=V_alpha,
                                                 ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data, n_latent=n_latent,
                         beta_start=beta_start, mu_beta=mubeta, V_beta=V_beta,
                         lambda_start=lambda_start, mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         alpha_start=alpha_start, V_alpha=V_alpha, site_effect=site_effect,
                         family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     mcmc.alpha = MCMC.alpha, 
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    ##======== with lv and random site effect======
    if(n_latent>0 && site_effect=="random"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect and latent variables from data about only one species.\n site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      
      # Form and check priors
      mubeta <- check.mubeta(mu_beta,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_rand_site_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                Y=Y, X=as.matrix(X),
                                                beta_start=beta_start,
                                                mu_beta = mubeta, V_beta=V_beta,
                                                lambda_start=lambda_start,
                                                mu_lambda = mulambda, V_lambda=Vlambda,
                                                W_start = W_start, V_W = V_W,
                                                alpha_start=alpha_start,
                                                V_alpha_start=V_alpha,
                                                shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                                                ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"    
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y, 
                         site_formula=site_formula,
                         site_data=site_data, n_latent=n_latent,
                         beta_start=beta_start,
                         mu_beta=mubeta, V_beta=V_beta,
                         lambda_start=lambda_start,
                         mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         alpha_start=alpha_start,
                         V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
  }
  #======== function with traits ========
  if(!is.null(trait_data)){
    if (nsp==1) {
      message("Error: Unable to estimate the influence of species-specific traits on species' responses from data about only one species.\n
          trait_data should not be specified with a single species.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    #==== Trait formula ======
    # used to fix zeros in gamma matrix (np x nt) 
    if((trait_formula==~.) || is.null(trait_formula)){
      mf.suit.tr <- model.frame(formula=~., data=as.data.frame(trait_data))
      Tr <- model.matrix(attr(mf.suit.tr,"terms"), data=mf.suit.tr)
      nt <- ncol(Tr)
      gamma_zeros <- matrix(1,nt,np)
    }
    else{
      data <- as.data.frame(trait_data)
      # add column of 1 with names of covariables in site_data 
      data[,colnames(X)] <- 1
      mf.suit.tr <- model.frame(formula=trait_formula, data=data)
      # full design matrix corresponding to formula
      mod.mat <- model.matrix(attr(mf.suit.tr,"terms"), data=mf.suit.tr)
      # Remove duplicated columns to get design matrix for traits 
      Tr <- as.matrix(mod.mat[,!duplicated(mod.mat,MARGIN=2)])
      colnames(Tr) <- colnames(mod.mat)[!duplicated(mod.mat,MARGIN=2)]
      for(p in 1:np){
        if(sum(colnames(Tr)==colnames(X)[p])==0){
          colnames(Tr) <- gsub(pattern=paste0(":",colnames(X)[p]), replacement="",
                               x=colnames(Tr), fixed=TRUE)
          colnames(Tr) <- gsub(pattern=paste0(colnames(X)[p],":"), replacement="",
                               x=colnames(Tr), fixed=TRUE)
        }
      }
      nt <- ncol(Tr)
      n_Tint <- sum(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
      col_Tint <- which(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
      if(n_Tint!=1) {
        message("Error: The model must include one trait intercept to be interpretable.\n")
        stop("Please respecify the trait_formula formula and call ", calling.function(), " again.",
             call.=FALSE)
      }
      gamma_zeros <- matrix(0,nt,np)
      rownames(gamma_zeros) <- colnames(Tr)
      colnames(gamma_zeros) <- colnames(X)
      for(t in 1:nt){
        for(p in 1:np){
          term <-  c(grep(paste0(colnames(X)[p],":"), colnames(mod.mat), value=TRUE, fixed=TRUE),grep(paste0(":",colnames(X)[p]), colnames(mod.mat), value=TRUE, fixed=TRUE))
          if(length(term)==0) next
          # fixed=TRUE pattern is a string to be matched as is 
          # not a regular expression because of special characters in formula
          gamma_zeros[t,p] <- length(c(grep(paste0(":",colnames(Tr)[t]), term, fixed=TRUE),grep(paste0(colnames(Tr)[t],":"), term, fixed=TRUE)))
        }
        gamma_zeros[t,col_Xint] <- length(which(colnames(mod.mat)==colnames(Tr)[t]))  
      }
      gamma_zeros[col_Tint,] <- 1
    }
    #===== Check traits data ======= 
    check.X(Tr, nsp)
    
    ##======== without latent variables and site effect ========
    if(n_latent==0 && site_effect=="none"){
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                          Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                          beta_start=beta_start,
                                          gamma_start=gamma_start,
                                          V_gamma=Vgamma, mu_gamma=mugamma,
                                          gamma_zeros=gamma_zeros, 
                                          V_beta=V_beta,
                                          ropt=ropt, seed=seed,
                                          verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j, start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data, 
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with latent variables ======
    if(n_latent>0 && site_effect=="none"){
      
      if (nsp==1) {
        message("Error: Unable to adjust latent variables from data about only one species.\n n_latent must be equal to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                             Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                             beta_start=beta_start,
                                             gamma_start=gamma_start,
                                             V_gamma=Vgamma, mu_gamma=mugamma,
                                             gamma_zeros=gamma_zeros, 
                                             V_beta=V_beta,
                                             lambda_start=lambda_start,
                                             mu_lambda = mulambda, V_lambda=Vlambda,
                                             W_start = W_start, V_W = V_W,
                                             ropt=ropt, seed=seed, verbose=verbose)
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data,
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         lambda_start=lambda_start,
                         mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with random site effect ======
    if(n_latent==0 && site_effect=="random"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to 'none' with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits_rand_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                    Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                                    beta_start=beta_start,
                                                    V_beta=V_beta,
                                                    gamma_start=gamma_start,
                                                    V_gamma=Vgamma, mu_gamma=mugamma,
                                                    gamma_zeros=gamma_zeros,
                                                    alpha_start=alpha_start,
                                                    V_alpha_start=V_alpha,
                                                    shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                                                    ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j,start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data, 
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         alpha_start=alpha_start,
                         V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp=MCMC.sp,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with fixed site effect ======
    if(n_latent==0 && site_effect=="fixed"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect from data about only one species.\n site_effect must be equal to 'none' with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits_fixed_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                     Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                                     beta_start=beta_start,
                                                     gamma_start=gamma_start,
                                                     V_gamma=Vgamma, mu_gamma=mugamma,
                                                     gamma_zeros=gamma_zeros,
                                                     V_beta=V_beta,
                                                     alpha_start=alpha_start, V_alpha=V_alpha, 
                                                     ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(MCMC.beta_j,start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data,
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         alpha_start=alpha_start, V_alpha=V_alpha, 
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp,
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.alpha = MCMC.alpha,
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    
    ##======== with lv and fixed site effect======
    if(n_latent>0 && site_effect=="fixed"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect and latent variables from data about only one species.\n site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits_fixed_site_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                        Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                                        beta_start=beta_start,
                                                        gamma_start=gamma_start,
                                                        V_gamma=Vgamma, mu_gamma=mugamma,
                                                        gamma_zeros=gamma_zeros,
                                                        V_beta=V_beta,
                                                        lambda_start=lambda_start,
                                                        mu_lambda = mulambda, V_lambda=Vlambda,
                                                        W_start = W_start, V_W = V_W,
                                                        alpha_start=alpha_start, V_alpha=V_alpha,
                                                        ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y,
                         site_formula=site_formula,
                         site_data=site_data,
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         lambda_start=lambda_start,
                         mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         alpha_start=alpha_start,
                         V_alpha=V_alpha, site_effect=site_effect,
                         family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     mcmc.alpha = MCMC.alpha, 
                     log_theta_latent=mod$log_theta_latent,
                     theta_latent=mod$theta_latent,
                     model_spec=model_spec)
    }
    ##======== with lv and random site effect======
    if(n_latent>0 && site_effect=="random"){
      if (nsp==1) {
        message("Error: Unable to adjust site effect and latent variables from data about only one species.\n site_effect must be equal to 'none' and n_latent to 0 with a single species.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
      
      # Initial starting values for M-H
      beta_start <- form.beta.start.sp(beta_start, np, nsp)
      lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
      alpha_start <- form.alpha.start.sp(alpha_start, nsite)
      W_start <-form.W.start.sp(W_start, nsite, n_latent)
      gamma_start <- form.gamma.start.mat(gamma_start, nt, np)
      
      # Form and check priors
      mugamma <- check.mugamma.mat(mu_gamma,nt,np)
      Vgamma <- check.Vgamma.mat(V_gamma,nt,np)
      V_beta <- check.Vbeta(V_beta,np)
      mulambda <- check.mubeta(mu_lambda,n_latent)
      Vlambda <- check.Vlambda(V_lambda,n_latent)
      V_W <- rep(1,n_latent)
      V_alpha <- check.Valpha(V_alpha)
      
      # call Rcpp function
      mod <- Rcpp_jSDM_poisson_log_traits_rand_site_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                                       Y=Y, X=as.matrix(X), Tr=as.matrix(Tr),
                                                       beta_start=beta_start,
                                                       gamma_start=gamma_start,
                                                       V_gamma=Vgamma, mu_gamma=mugamma,
                                                       gamma_zeros=gamma_zeros,
                                                       V_beta=V_beta,
                                                       lambda_start=lambda_start,
                                                       mu_lambda = mulambda, V_lambda=Vlambda,
                                                       W_start = W_start, V_W = V_W,
                                                       alpha_start=alpha_start,
                                                       V_alpha_start=V_alpha,
                                                       shape_Valpha=shape_Valpha,
                                                       rate_Valpha=rate_Valpha,
                                                       ropt=ropt, seed=seed, verbose=verbose)
      
      
      
      
      #= Transform Sample list in an MCMC object
      MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)     
      colnames(MCMC.Deviance) <- "Deviance"
      MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.alpha) <- paste0("alpha_",rownames(Y))
      MCMC.V_alpha <- coda::mcmc(mod$V_alpha,start=nburn+1,end=ngibbs,thin=nthin)
      colnames(MCMC.V_alpha) <- "V_alpha"    
      MCMC.sp <- list()
      for (j in 1:nsp) {
        ## beta_j
        MCMC.beta_j <- coda::mcmc(as.matrix(mod$beta[,j,]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
        ## lambda_j
        MCMC.lambda_j <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)	
        colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
        
        MCMC.sp[[colnames(Y)[j]]] <- coda::mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
      }
      ## gamma 
      MCMC.gamma <- list()
      for (p in 1:np) {
        MCMC.gamma_p <- coda::mcmc(as.matrix(mod$gamma[,,p]), start=nburn+1, end=ngibbs, thin=nthin)
        colnames(MCMC.gamma_p) <- paste0("gamma_",colnames(X)[p],".",colnames(Tr))
        MCMC.gamma[[p]] <- MCMC.gamma_p
      }
      ## W latent variables 
      MCMC.latent <- list()
      for (l in 1:n_latent) {
        MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
        MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
      }
      
      #= Model specification, site_formula,
      model_spec <- list(burnin=burnin, mcmc=mcmc, thin=thin,
                         count_data=Y, 
                         site_formula=site_formula,
                         site_data=site_data,
                         trait_data=trait_data,
                         trait_formula=trait_formula,
                         n_latent=n_latent, 
                         beta_start=beta_start, 
                         V_beta=V_beta,
                         gamma_start=gamma_start,
                         V_gamma=Vgamma, mu_gamma=mugamma,
                         gamma_zeros=gamma_zeros,
                         lambda_start=lambda_start,
                         mu_lambda=mulambda, V_lambda=Vlambda,
                         W_start=W_start, V_W=V_W,
                         alpha_start=alpha_start,
                         V_alpha_start=V_alpha,
                         shape_Valpha=shape_Valpha, rate_Valpha=rate_Valpha,
                         site_effect=site_effect, family="poisson", link="log",
                         ropt=ropt, seed=seed, verbose=verbose)
      colnames(mod$log_theta_latent) <- colnames(mod$theta_latent) <- colnames(Y)
      rownames(mod$log_theta_latent) <- rownames(mod$theta_latent) <- rownames(Y)
      #= Output
      output <- list(mcmc.sp= MCMC.sp, 
                     mcmc.gamma=MCMC.gamma,
                     mcmc.Deviance=MCMC.Deviance,
                     mcmc.latent = MCMC.latent,
                     mcmc.alpha = MCMC.alpha, mcmc.V_alpha = MCMC.V_alpha,
                     log_theta_latent=mod$log_theta_latent,
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
