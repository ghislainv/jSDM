#' @docType package
#' @name jSDM-package
#' @aliases jSDM-package
#' @title joint species distribution models
#' @description \code{jSDM} is an R package for fitting joint species distribution models (jSDM) in a hierarchical Bayesian framework.
#' The Gibbs sampler is written in C++. It uses Rcpp, Armadillo and GSL to maximize computation efficiency.
#' #'
#' The package includes the following functions to fit various species distribution models :
#' \itemize{
#'   
#'   \item{ \code{\link{jSDM_binomial_probit_block}} : 
#'
#' \bold{Ecological process:}
#' \deqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij},t_i)}{y_ij ~ Binomial(\theta_ij,t_i),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j}{(\theta_ij) = \beta_0j + X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j}{(\theta_ij) = \beta_0j + X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j  + \alpha_i}{(\theta_ij) = \beta_0j + X_i \beta_j + \alpha_i}  and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) = \beta_0j + X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' 
#'   \item{ \code{\link{jSDM_binomial_logit}} : 
#'
#' \bold{Ecological process : }
#' \deqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij},t_i)}{y_ij ~ Binomial(\theta_ij,t_i),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab logit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j}{(\theta_ij) = \beta_0j + X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab logit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j}{(\theta_ij) = \beta_0j + X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab logit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j  + \alpha_i}{(\theta_ij) = \beta_0j + X_i \beta_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab logit\eqn{(\theta_{ij}) = \beta_{0j} + X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) = \beta_0j + X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' }
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