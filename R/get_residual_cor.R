## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name get_residual_cor
#' @aliases get_residual_cor
#' @title Calculate the residual correlation matrix from a latent variable model (LVM)
#' @description This function use coefficients \eqn{(\lambda_{jl} with j=1,\dots,n_{species} and l=1,\dots,n_{latent})}{(\lambda_jl with j=1,...,n_species and l=1,...,n_latent)}, corresponding to latent variables fitted using \code{jSDM} package, to calculate the variance-covariance matrix which controls correlation between species.
#' @param mod An object of class \code{"jSDM"}
#' @param prob A numeric scalar in the interval \eqn{(0,1)} giving the target probability coverage of the highest posterior density (HPD) intervals, by which to determine whether the correlations are "significant". Defaults to 0.95.
#' @return results A list including : 
#' \item{cov.mean}{Average over the MCMC samples of the variance-covariance matrix.} 
#' \item{cov.median}{Median over the MCMC samples of the variance-covariance matrix.}
#' \item{cov.lower}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the lower limits of the  (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval of variance-covariance matrices over the MCMC samples.}
#' \item{cov.upper}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the upper limits of the  (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval of variance-covariance matrices over the MCMC samples.}
#' \item{cov.sig}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the value 1 corresponding to the “significant" co-variances and the value 0 corresponding to "non-significant" co-variances, whose (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval over the MCMC samples contain zero.}
#' \item{cor.mean}{Average over the MCMC samples of the residual correlation matrix.}
#' \item{cor.median}{Median over the MCMC samples of the residual correlation matrix.}
#' \item{cor.lower}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the lower limits of the  (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval of correlation matrices over the MCMC samples.}
#' \item{cor.upper}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the upper limits of the  (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval of correlation matrices over the MCMC samples.}
#' \item{cor.sig}{A \eqn{n_{species} \times n_{species}}{n_species x n_species} matrix containing the value \eqn{1} corresponding to the “significant" correlations and the value \eqn{0} corresponding to "non-significant" correlations,
#'  whose (\eqn{100 \times prob \%}{100 x prob \%}) HPD interval over the MCMC samples contain zero.}
#'
#' @details  After fitting the jSDM with latent variables, the \bold{fullspecies residual correlation matrix} : \eqn{R=(R_{ij})}{R=(R_ij)} with \eqn{i=1,\ldots, n_{species}}{i=1,..., n_species} and \eqn{j=1,\ldots, n_{species}}{j=1,..., n_species} can be derived from the covariance in the latent variables such as : 
#' \tabular{lll}{
#' \eqn{\Sigma_{ij}}{Sigma_ij} \tab \eqn{= \lambda_i .\lambda_j' + 1}{= \lambda_i . \lambda_j' + 1} \tab if i=j \cr
#'          \tab \eqn{= \lambda_i .\lambda_j'}{= \lambda_i . \lambda_j'} \tab else, \cr}
#' then we compute correlations from covariances :
#'\deqn{R_{ij} = \frac{\Sigma_{ij}}{\sqrt{\Sigma_ii\Sigma _jj}}}{R_ij = Sigma_ij / sqrt(Sigma_ii.Sigma _jj)}.
#' @author
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#' 
#' Jeanne Clément <jeanne.clement16@laposte.net>
#' 
#' @references  
#' Hui FKC (2016). “boral: Bayesian Ordination and Regression Analysis of Multivariate Abundance Data in R.” \emph{Methods in Ecology and Evolution}, 7, 744–750.
#' 
#' Ovaskainen et al. (2016). Using latent variable models to identify large networks of species-to-species associations at different spatial scales. \emph{Methods in Ecology and Evolution}, 7, 549-555.
#' 
#' Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). \emph{Methods in Ecology and Evolution}, 5, 397-406.
#' 
#' @seealso \code{\link{get_enviro_cor}} \code{\link[stats]{cov2cor}} \code{\link{jSDM-package}} \code{\link{jSDM_binomial_probit}} \cr 
#'  \code{\link{jSDM_binomial_logit}}  \code{\link{jSDM_poisson_log}} 
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
#' # Increase the number of iterations to reach MCMC convergence
#'  mod <- jSDM_binomial_probit(# Response variable
#'                              presence_data=PA_frogs,
#'                              # Explanatory variables
#'                              site_formula = ~.,
#'                              site_data = Env_frogs,
#'                              n_latent=2,
#'                              site_effect="random",
#'                              # Chains
#'                              burnin=100,
#'                              mcmc=100,
#'                              thin=1,
#'                              # Starting values
#'                              alpha_start=0,
#'                              beta_start=0,
#'                              lambda_start=0,
#'                              W_start=0,
#'                              V_alpha=1,
#'                              # Priors
#'                              shape=0.5, rate=0.0005,
#'                              mu_beta=0, V_beta=10,
#'                              mu_lambda=0, V_lambda=10,
#'                              # Various
#'                              seed=1234, verbose=1)
#' # Calcul of residual correlation between species 
#' result <- get_residual_cor(mod, prob=0.95)
#' # Residual variance-covariance matrix
#' result$cov.mean
#' ## All non-significant co-variances are set to zero.
#' result$cov.mean * result$cov.sig
#' # Residual correlation matrix
#' result$cor.mean
#' ## All non-significant correlations are set to zero.
#' result$cor.mean * result$cor.sig
#' @keywords stats::cov2cor
#' @importFrom stats cov2cor
#' @importFrom coda HPDinterval
#' @export


# Calculate the residual correlation matrix from a LVM
get_residual_cor <- function(mod, prob=0.95) {
  #== Check
  if(!class(mod)=="jSDM"){
    stop("Please provide an object of class jSDM in", calling.function(), call.=FALSE)
  }
  if(mod$model_spec$n_latent==0) {
    cat("Error: The jSDM class object provided is not a latent variable model (LVM).\n
        The factor loadings needed to compute the residual correlation matrix have not been estimated. \n")
    stop("Please fit a LVM on your data and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(mod$model_spec$n_latent==1) {
    cat("Error: Residual correlation matrix is reliably modeled only with two or more latent variables. \n")
    stop("Please fit a LVM with n_latent > 1 on your data and call ", calling.function(), " again.",
         call.=FALSE)
  }
  
  if(prob>1 || prob<0) {
    cat("Error: The target probability coverage of the intervals must be in the interval ]0,1[. \n")
    stop("Please specify a probability for prob and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.null(mod$model_spec$presence_data)){
    n.species <- ncol(mod$model_spec$presence_data)
  }
  if(!is.null(mod$model_spec$count_data)){
    n.species <- ncol(mod$model_spec$count_data)
  }
  if(!is.null(mod$model_spec$data)){
    n.species <- length(unique(mod$model_spec$data$species))
  }
  n.mcmc <- nrow(mod$mcmc.latent[[1]])
  Tau.arr <- matrix(NA,n.mcmc,n.species^2)
  Tau.cor.arr <- matrix(NA,n.mcmc,n.species^2)
  
  for(t in 1:n.mcmc) { 
    lv.coefs <- t(sapply(mod$mcmc.sp, "[", t, grep("lambda",colnames(mod$mcmc.sp[[1]]))))
    Tau.mat <- lv.coefs%*%t(lv.coefs) + diag(n.species)
    Tau.arr[t,] <- as.vector(Tau.mat) 
    Tau.cor.mat <- cov2cor(Tau.mat)
    Tau.cor.arr[t,] <- as.vector(Tau.cor.mat) 
  }
  # Lower/upper limits of HPD interval over the MCMC samples
  cov.lower <- cov.upper <- cor.lower <- cor.upper <-  matrix(NA, n.species, n.species)
  sig.Tau.cor <- sig.Tau.mat <- matrix(NA, n.species, n.species)
  for(j in 1:n.species){
    for(jprim in 1:n.species){
      ## Residual variance-covariance matrices
      get.hpd.covs <- coda::HPDinterval(coda::as.mcmc(Tau.arr[,(j-1)*n.species + jprim]),
                                        prob = prob)
      cov.lower[j, jprim] <- get.hpd.covs[1]
      cov.upper[j, jprim] <- get.hpd.covs[2]
      ## Significant values whose HPD interval does not contain zero
      sig.Tau.mat[j, jprim] <- ifelse((0 > get.hpd.covs[1]) & (0 < get.hpd.covs[2]), 0, 1)
      ## Residual correlations matrices
      get.hpd.cors <- coda::HPDinterval(coda::as.mcmc(Tau.cor.arr[,(j-1)*n.species + jprim]),
                                        prob = prob)
      cor.lower[j, jprim] <- get.hpd.cors[1]
      cor.upper[j, jprim] <- get.hpd.cors[2]
      ## Significant values whose HPD interval does not contain zero
      sig.Tau.cor[j, jprim] <-  ifelse((0 > get.hpd.cors[1]) & (0 < get.hpd.cors[2]), 0, 1)
    }
  }
  ## Average/Median over the MCMC samples
  Tau.mat.mean <-  matrix(apply(Tau.arr,2,mean),n.species,byrow=F)
  Tau.mat.median <-  matrix(apply(Tau.arr,2,median),n.species,byrow=F)
  Tau.cor.mean <- matrix(apply(Tau.cor.arr,2,mean),n.species,byrow=F)
  Tau.cor.median <- matrix(apply(Tau.cor.arr,2,median),n.species,byrow=F)
  # Species names in results 
  if(!is.null(mod$model_spec$presence_data)){
    names.sp <- colnames(mod$model_spec$presence_data)
  }
  if(!is.null(mod$model_spec$count_data)){
    names.sp <-  colnames(mod$model_spec$count_data)
  }
  if(!is.null(mod$model_spec$data)){
    names.sp <- unique(mod$model_spec$data$species)
  }
  dimnames(Tau.mat.mean) <- dimnames(Tau.mat.median) <- dimnames(cov.lower) <- dimnames(cov.upper) <- dimnames(sig.Tau.mat) <- list(names.sp, names.sp)
  dimnames(Tau.cor.mean) <- dimnames(Tau.cor.median) <- dimnames(cor.lower) <- dimnames(cor.upper) <- dimnames(sig.Tau.cor) <- list(names.sp, names.sp)
  # Results 
  results <- list(cov.mean = Tau.mat.mean, cov.median = Tau.mat.median,
                  cov.lower= cov.lower, cov.upper = cov.upper, cov.sig = sig.Tau.mat,
                  cor.mean = Tau.cor.mean, cor.median = Tau.cor.median,
                  cor.lower=cor.lower, cor.upper=cor.upper, cor.sig=sig.Tau.cor)
  return(results)
}	

