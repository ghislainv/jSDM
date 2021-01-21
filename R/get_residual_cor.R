## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name get_residual_cor
#' @aliases get_residual_cor
#' @title Calculate the residual correlation matrix from a latent variable model (LVM)
#' @description This function use coefficients \eqn{(\lambda_{jl} with j=1,\dots,n_species and l=1,\dots,n_latent)}{(\lambda_jl with j=1,...,n_species and l=1,...,n_latent)}, corresponding to latent variables fitted using \code{jSDM} package, to calculate the variance-covariance matrix which controls correlation between species.
#' @param mod An object of class \code{"jSDM"}
#' @return results A list including : \tabular{ll}{
#' cov.mean \tab Average over the MCMC samples of the variance-covariance matrix \cr 
#' cov.median \tab Median over the MCMC samples of the variance-covariance matrix \cr
#' cor.mean \tab Average over the MCMC samples of the residual correlation matrix \cr
#' cor.median \tab Median over the MCMC samples of the residual correlation matrix \cr
#' }
#' @details  After fitting the jSDM with latent variables, the \bold{fullspecies residual correlation matrix} : \eqn{R=(R_ij) avec i=1,\ldots, nspecies et j=1,\ldots, nspecies}{R=(R_ij) avec i=1,..., nspecies et j=1,..., nspecies} can be derived from the covariance in the latent variables such as : 
#' \tabular{lll}{
#' \eqn{\Sigma_{ij}}{Sigma_ij} \tab \eqn{= \lambda_i .\lambda_j' + 1}{= \lambda_i . \lambda_j' + 1} \tab if i=j \cr
#'          \tab \eqn{= \lambda_i .\lambda_j'}{= \lambda_i . \lambda_j'} \tab else, \cr}
#' then we compute correlations from covariances :
#'\deqn{R_{ij} = \frac{\Sigma_{ij}}{\sqrt{\Sigma_ii\Sigma _jj}}}{R_ij = Sigma_ij / sqrt(Sigma_ii.Sigma _jj)}.
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr }
#' @references  \tabular{l}{
#' Hui FKC (2016). “boral: Bayesian Ordination and Regression Analysis of Multivariate Abundance Data in R.” Methods in Ecology and Evolution, 7, 744–750. \cr
#' Ovaskainen et al. (2016). Using latent variable models to identify large networks of species-to-species associations at different spatial scales. Methods in Ecology and Evolution, 7, 549-555.\cr
#' Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). Methods in Ecology and Evolution, 5, 397-406.\cr }
#' @seealso \code{\link{jSDM-package}} \code{\link{jSDM_binomial_probit_block}} \code{\link{jSDM_binomial_logit}} 
#' @examples 
#' library(jSDM)
#' # frogs data
#'  data(frogs, package="jSDM")
#'  # Arranging data
#'  PA_frogs <- frogs[,4:12]
#'  # Normalized continuous variables
#'  Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
#'  colnames(Env_frogs) <- colnames(frogs[,1:3])
#'  Env_frogs <- as.data.frame(Env_frogs)
#'  # Parameter inference
#'# Increase the number of iterations to reach MCMC convergence
#' mod <- jSDM_binomial_probit_block(# Response variable
#'                                   presence_site_sp=PA_frogs,
#'                                   # Explanatory variables
#'                                   site_suitability = ~.,
#'                                   site_data = Env_frogs,
#'                                   n_latent=2,
#'                                   site_effect="random",
#'                                   # Chains
#'                                   burnin=100,
#'                                   mcmc=100,
#'                                   thin=1,
#'                                   # Starting values
#'                                   alpha_start=0,
#'                                   beta_start=0,
#'                                   lambda_start=0,
#'                                   W_start=0,
#'                                   V_alpha=1,
#'                                   # Priors
#'                                   shape=0.5, rate=0.0005,
#'                                   mu_beta=0, V_beta=10,
#'                                   mu_lambda=0, V_lambda=10,
#'                                   # Various
#'                                   seed=1234, verbose=1)
#' # Calcul of residual correlation between species 
#'  result <- get_residual_cor(mod)
#'  result$cov.mean
#'  result$cor.mean
#' @keywords stats::cov2cor
#' @seealso \code{\link{get_enviro_cor}} \code{\link[stats]{cov2cor}}  \code{\link{jSDM_binomial_logit}}  \code{\link{jSDM_poisson_log}}  \code{\link{jSDM_binomial_probit_block}} 
#' @importFrom stats cov2cor
#' @export


# Calculate the residual correlation matrix from a LVM
get_residual_cor <- function(mod) {
  #== Check
  if (!class(mod)=="jSDM"){
    stop("Please provide an object of class jSDM in", calling.function(), call.=FALSE)
  }
  if (mod$model_spec$n_latent==0) {
    cat("Error: The jSDM class object provided is not a latent variable model (LVM).\n
        The lambdas parameters needed to compute the residual correlation matrix have not been estimated. \n")
    stop("Please fit a LVM on your data and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.null(mod$model_spec$presences)){
    n.species <- ncol(mod$model_spec$presences)
  }
  if(!is.null(mod$model_spec$data)){
    n.species <- length(unique(mod$model_spec$data$species))
  }
  n.mcmc <- nrow(mod$mcmc.latent[[1]])
  Tau.arr <- matrix(NA,n.mcmc,n.species^2)
  Tau.cor.arr <- matrix(NA,n.mcmc,n.species^2)
  
  for(t in 1:n.mcmc) { 
    lv.coefs <- t(sapply(mod$mcmc.sp, "[", t, grep("lambda",colnames(mod$mcmc.sp$sp_1))))
    Tau.mat <- lv.coefs%*%t(lv.coefs) + diag(n.species)
    Tau.arr[t,] <- as.vector(Tau.mat) 
    Tau.cor.mat <- cov2cor(Tau.mat)
    Tau.cor.arr[t,] <- as.vector(Tau.cor.mat) 
  }
  
  ## Average/Median over the MCMC samples
  Tau.mat.mean <- sig.Tau.mat.mean <- matrix(apply(Tau.arr,2,mean),n.species,byrow=F)
  Tau.mat.median <- sig.Tau.mat.median <- matrix(apply(Tau.arr,2,median),n.species,byrow=F)
  Tau.cor.mean <- sig.Tau.cor.mean <- matrix(apply(Tau.cor.arr,2,mean),n.species,byrow=F)
  Tau.cor.median <- sig.Tau.cor.median <- matrix(apply(Tau.cor.arr,2,median),n.species,byrow=F)
  results <- list(cov.mean = Tau.mat.mean, cov.median = Tau.mat.median, cor.mean = Tau.cor.mean, cor.median = Tau.cor.median)
  return(results)
}	

