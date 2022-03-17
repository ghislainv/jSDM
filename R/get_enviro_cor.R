## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name get_enviro_cor
#' @aliases get_enviro_cor
#' @title Extract covariances and correlations due to shared environmental responses
#' @description Calculates the correlation between columns of the response matrix,
#'  due to similarities in the response to explanatory variables i.e., shared environmental response.
#' @param mod An object of class \code{"jSDM"}
#' @param type A choice of either the posterior median (\code{type = "median"}) or posterior mean (\code{type = "mean"}), which are then treated as estimates and the fitted values are calculated from.
#' Default is posterior mean.
#' @param prob A numeric scalar in the interval \eqn{(0,1)} giving the target probability coverage of the intervals, by which to determine whether the correlations are "significant".
#' Defaults to 0.95.
#' @return results A list including : 
#' \item{cor, cor.lower, cor.upper}{A set of \eqn{np \times np}{np x np} correlation matrices, containing either the posterior median or mean estimate  over the MCMC samples plus lower and upper limits of the corresponding 95 \% highest posterior interval.} 
#' \item{sig.cor}{A \eqn{np \times np}{np x np} correlation matrix containing only the “significant" correlations whose 95 \% highest posterior density (HPD) interval does not contain zero. All non-significant correlations are set to zero.}
#' \item{cov}{Average over the MCMC samples of the \eqn{np \times np}{np x np} covariance matrix.}
#' 
#' @details In both independent response and correlated response models, where each of the columns of the response matrix \eqn{Y} are fitted to a set of explanatory variables given by \eqn{X},
#' the covariance between two columns \eqn{j} and \eqn{j'}, due to similarities in their response to the model matrix, is thus calculated based on the linear predictors \eqn{X \beta_j} and \eqn{X \beta_j'}, where \eqn{\beta_j} are species effects relating to the explanatory variables.
#  For multivariate abundance data, the correlation calculated by this function can be interpreted as the correlation attributable to similarities in the environmental response between species. 
#' Such correlation matrices are discussed and found in Ovaskainen et al., (2010), Pollock et al., (2014). 
#' @author 
#' 
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#' 
#' Jeanne Clément <jeanne.clement16@laposte.net> 
#' 
#'@references
#' Hui FKC (2016). “boral: Bayesian Ordination and Regression Analysis of Multivariate Abundance Data in R.” \emph{Methods in Ecology and Evolution}, 7, 744–750.  
#' 
#' Ovaskainen et al. (2010). Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. \emph{Ecology}, 91, 2514-2521. 
#' 
#' Pollock et al. (2014). Understanding co-occurrence by modelling species simultaneously with a Joint Species Distribution Model (JSDM). \emph{Methods in Ecology and Evolution}, 5, 397-406. 
#' 
#'@seealso \code{\link[stats]{cov2cor}} \code{\link{get_residual_cor}} \code{\link{jSDM-package}} \code{\link{jSDM_binomial_probit}} \cr
#' \code{\link{jSDM_binomial_logit}} \code{\link{jSDM_poisson_log}} 
#' @examples 
#' library(jSDM)
#' # frogs data
#'  data(frogs, package="jSDM")
#'  # Arranging data
#'  PA_frogs <- frogs[,4:12]
#'  # Normalized continuous variables
#'  Env_frogs <- cbind(scale(frogs[,1]),frogs[,2], 
#'                     scale(frogs[,3]))
#'  colnames(Env_frogs) <- colnames(frogs[,1:3])
#'  Env_frogs <- as.data.frame(Env_frogs)
#'  # Parameter inference
#'# Increase the number of iterations to reach MCMC convergence
#' mod <- jSDM_binomial_probit(# Response variable
#'                              presence_data=PA_frogs,
#'                              # Explanatory variables
#'                              site_formula = ~.,
#'                              site_data = Env_frogs,
#'                              n_latent=0,
#'                              site_effect="random",
#'                              # Chains
#'                              burnin=100,
#'                              mcmc=100,
#'                              thin=1,
#'                              # Starting values
#'                              alpha_start=0,
#'                              beta_start=0,
#'                              V_alpha=1,
#'                              # Priors
#'                              shape=0.5, rate=0.0005,
#'                              mu_beta=0, V_beta=10,
#'                              # Various
#'                              seed=1234, verbose=1)
#' # Calcul of residual correlation between species 
#'  enviro.cors <- get_enviro_cor(mod)
# corrplot::corrplot(enviro.cors$sig.cor, title = "Shared response correlations",
#                    type = "lower", diag = FALSE, mar = c(3,0.5,2,1), tl.srt = 45)
#' @keywords stats::cov2cor
#' @importFrom stats cov cor
#' @importFrom coda as.mcmc HPDinterval
#' @importFrom stringi stri_remove_empty
#' @export
get_enviro_cor <- function(mod, type = "mean", prob = 0.95) {
  ##= Check
  if (!class(mod)=="jSDM"){
    stop("Please provide an object of class jSDM in", calling.function(),call.=FALSE)
  }
  if (!(type %in% c("mean","median"))) {stop("type must be \"mean\" or \"median\"")}
  if (prob<0 | prob>1) {stop("prob must be a probability between (0,1)")}
  model.spec <- mod$model_spec
  n.X.coeff <- nrow(model.spec$beta_start)
  if (n.X.coeff==0) {
    cat("Error: The object provided is not a latent variable model (LVM).\n
        The lambdas parameters needed to compute the residual correlation matrix have not been estimated. \n")
    stop("Please fit a LVM on your data and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.null(model.spec$presence_data)){
    y <- model.spec$presence_data
    species <- colnames(y)
    n.species <- ncol(model.spec$presence_data)
    n.sites <- nrow(model.spec$presence_data)
    # Suitability process
    suitability <- model.spec$site_formula
    mf.suit <- model.frame(formula=suitability,data=as.data.frame(model.spec$site_data))
    X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  }
  if(!is.null(model.spec$count_data)){
    y <- model.spec$count_data
    species <- colnames(y)
    n.species <- ncol(model.spec$count_data)
    n.sites <- nrow(model.spec$count_data)
    # Suitability process
    suitability <- model.spec$site_formula
    mf.suit <- model.frame(formula=suitability,data=as.data.frame(model.spec$site_data))
    X <- model.matrix(attr(mf.suit,"terms"),data=mf.suit)
  }
  if(!is.null(model.spec$data)){
    y <- model.spec$data$Y
    species <- unique(model.spec$data$species)
    n.species <- length(unique(model.spec$data$species))
    n.sites <- length(unique(model.spec$data$site))
    # Suitability process
    suitability <- model.spec$site_formula 
    if(model.spec$site_formula==~.) suitability <- ~. - site - Y
    mf.suit <- model.frame(formula=suitability, data=as.data.frame(model.spec$data))
    # design matrix X for species effects beta
    Xterms <- stringi::stri_remove_empty(gsub(":?species:?", "", 
                                              grep("species", attr(attr(mf.suit,"terms"),"term.labels"), value=T)))  
    Xformula <- paste0("~",paste0(Xterms, collapse="+"))
    mf.suit.X <- model.frame(formula=Xformula, data=as.data.frame(model.spec$data))
    attr(attr(mf.suit.X,"terms"),"intercept") <- ifelse(grepl("- *species", suitability[2]),0,1)
    X <- model.matrix(attr(mf.suit.X,"terms"), data=mf.suit.X)
  }
  n.mcmc <- nrow(mod$mcmc.sp[[1]])
  
  if (length(grep("beta",colnames(mod$mcmc.sp[[1]]))) == 0) {
    stop("Cannot find MCMC sample corresponding to coefficients for X")
  }
  enviro_cor_mat <- enviro_cor_mat_cilower <- enviro_cor_mat_ciupper <- enviro_cov_mat <- matrix(0, n.species, n.species)
  sig_enviro_cor_mat <- matrix(0, n.species, n.species)
  rownames(enviro_cor_mat) <- rownames(enviro_cor_mat_cilower) <- rownames(enviro_cor_mat_ciupper) <- rownames(enviro_cov_mat) <- rownames(sig_enviro_cor_mat) <- species
  colnames(enviro_cor_mat) <- colnames(enviro_cor_mat_cilower) <- colnames(enviro_cor_mat_ciupper) <- colnames(enviro_cov_mat) <- colnames(sig_enviro_cor_mat) <- species
  all_enviro_cov_mat <- all_enviro_cor_mat <- array(0, dim = c(n.mcmc, n.species, n.species))
  
  
  for (t in 1:n.mcmc)
  {
    beta <- sapply(mod$mcmc.sp, "[", t, grep("beta",colnames(mod$mcmc.sp[[1]])))
    enviro.linpreds <- as.matrix(X) %*% as.matrix(beta)
    all_enviro_cov_mat[t, , ] <- cov(enviro.linpreds)
    all_enviro_cor_mat[t, , ] <- cor(enviro.linpreds)
  }
  
  for (j in 1:n.species) {
    for (j2 in 1:n.species)
    { ## Average/Median over the MCMC samples
      if (type == "median") {
        enviro_cov_mat[j, j2] <- median(all_enviro_cov_mat[, j, j2])
        enviro_cor_mat[j, j2] <- median(all_enviro_cor_mat[, j, j2])
      }
      if (type == "mean") {
        enviro_cov_mat[j, j2] <- mean(all_enviro_cov_mat[, j, j2])
        enviro_cor_mat[j, j2] <- mean(all_enviro_cor_mat[, j, j2])
      }
      
      sig_enviro_cor_mat[j, j2] <- enviro_cor_mat[j, j2]
      get.hpd.cors <- coda::HPDinterval(coda::as.mcmc(all_enviro_cor_mat[, j, j2]), prob = prob)
      enviro_cor_mat_cilower[j, j2] <- get.hpd.cors[1]
      enviro_cor_mat_ciupper[j, j2] <- get.hpd.cors[2]
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) {
        sig_enviro_cor_mat[j, j2] <- 0
      }
    }
  }
  
  return(list(cor = enviro_cor_mat, cor.lower = enviro_cor_mat_cilower,
              cor.upper = enviro_cor_mat_ciupper,
              sig.cor = sig_enviro_cor_mat,
              cov = enviro_cov_mat))
}