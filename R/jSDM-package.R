#' @docType package
#' @name jSDM-package
#' @aliases jSDM-package
#' @title joint species distribution models
#' @description \code{jSDM} is an R package for fitting joint species distribution models (jSDM) in a hierarchical Bayesian framework. The Gibbs sampler is written in C++. It uses Rcpp, Armadillo and GSL to maximize computation efficiency.
#'
#'   The package includes the following functions to fit various species distribution models : \tabular{llc}{
#'   function \tab data type \tab fitted model \cr
#'   jSDM_binomial_logit_one_species \tab presence-absence \tab  \eqn{y_n \sim \mathcal{B}inomial(\theta_n, t_n) / where, / logit(\theta_n) = X_n \beta}{y_n ~ Binomial(\theta_n, t_n) where, logit(\theta_n) = X_n \beta} \cr
#'   jSDM_binomial_logit \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / logit(\theta_{ij}) = X_i \beta_j}{y_ij ~ Binomial(\theta_ij, t_i) where, logit(\theta_ij) = X_i \beta_j} \cr
#'   jSDM_binomial_logit_rand_site \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / logit(\theta_{ij}) = X_i \beta_j   / and / \alpha_i \sim \mathcal{N}(0,V_{\alpha})}{y_ij ~ Binomial(\theta_ij, t_i) where, logit(\theta_ij) = X_i \beta_j + \alpha_i and \alpha_i ~ N(0,V_\alpha)} \cr
#'   jSDM_binomial_logit_lv \tab presence-absence \tab   \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / logit(\theta_{ij}) = X_i \beta_j + W_i \lambda_j}{y_ij ~ Binomial(\theta_ij, t_i) where, logit(\theta_ij) = X_i \beta_j + W_i \lambda_j} \cr
#'   jSDM_binomial_logit_rand_site_lv \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_ij, t_i) / where, / logit(\theta_{ij}) = X_i \beta_j +  W_i \lambda_j + \alpha_i / and /  \alpha_i \sim \mathcal{N}(0,V_{\alpha}) }{y_ij ~ Binomial(\theta_ij, t_i) where, logit(\theta_ij) = X_i \beta_j +  W_i \lambda_j + \alpha_i and \alpha_i ~ N(0,V_\alpha)} \cr
#'   jSDM_binomial_probit_block_one_species \tab presence-absence \tab  \eqn{y_n \sim \mathcal{B}inomial(\theta_n, t_n) \text{where, } probit(\theta_n) = X_n \beta}{y_n ~ Binomial(\theta_n, t_n) where, probit(\theta_n) = X_n \beta} \cr
#'   jSDM_binomial_probit_block \tab presence-absence \tab   \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / probit(\theta_{ij}) = X_i \beta_j}{y_ij ~ Binomial(\theta_ij, t_i) where, probit(\theta_ij) = X_i \beta_j} \cr
#'   jSDM_binomial_probit_block_rand_site \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / probit(\theta_{ij}) = X_i \beta_j + \alpha_i / and /  \alpha_i \sim \mathcal{N}(0,V_{\alpha})}{y_ij ~ Binomial(\theta_ij, t_i) where, probit(\theta_ij) = X_i \beta_j + \alpha_i and \alpha_i ~ N(0,V_\alpha)} \cr
#'   jSDM_binomial_probit_block_lv \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / probit(\theta_{ij}) = X_i \beta_j + W_i \lambda_j}{y_ij ~ Binomial(\theta_ij, t_i) where, probit(\theta_ij) = X_i \beta_j + W_i \lambda_j} \cr
#'   jSDM_binomial_probit_block_rand_site_lv \tab presence-absence \tab  \eqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) / where, / probit(\theta_{ij}) = X_i \beta_j +  W_i \lambda_j + \alpha_i / and /  \alpha_i \sim \mathcal{N}(0,V_{\alpha})}{y_ij ~ Binomial(\theta_ij, t_i) where, probit(\theta_ij) = X_i \beta_j +  W_i \lambda_j + \alpha_i and \alpha_i ~ N(0,V_\alpha)} \cr }
#' @details \tabular{ll}{
#'    Package: \tab jSDM\cr
#'    Type: \tab Package\cr
#'    Version: \tab 0.1.0\cr
#'    Date: \tab 2019-01-11\cr
#'    License: \tab GPL-3 \cr
#'    LazyLoad: \tab yes\cr }
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr
#'}
#' @keywords MCMC, Metropolis, algorithm, binomial, biodiversity, logistic model multivariate regression
#' @importFrom Rcpp evalCpp
#' @importFrom graphics par
#' @importFrom stats median model.frame model.matrix pnorm quantile nobs
#' @useDynLib jSDM , .registration=TRUE
NULL