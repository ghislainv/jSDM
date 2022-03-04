#' @docType package
#' @name jSDM-package
#' @aliases jSDM-package
#' @title joint species distribution models
#' @description \code{jSDM} is an R package for fitting joint species distribution models (JSDM) in a hierarchical Bayesian framework.
#'
#' The Gibbs sampler is written in C++. It uses Rcpp, Armadillo and GSL to maximize computation efficiency.
#'\tabular{ll}{
#'    Package: \tab jSDM\cr
#'    Type: \tab Package\cr
#'    Version: \tab 0.1.0\cr
#'    Date: \tab 2019-01-11\cr
#'    License: \tab GPL-3 \cr
#'    LazyLoad: \tab yes\cr }
#'    
#' @details The package includes the following functions to fit various species distribution models :
#' 
#' \tabular{ll}{
#'    function \tab data-type \cr
#'    \code{\link{jSDM_binomial_logit}} \tab presence-absence \cr
#'    \code{\link{jSDM_binomial_probit}} \tab presence-absence \cr
#'    \code{\link{jSDM_binomial_probit_sp_constrained}}  \tab presence-absence \cr
#'    \code{\link{jSDM_binomial_probit_long_format}} \tab presence-absence \cr
#'    \code{\link{jSDM_poisson_log}} \tab abundance \cr }
#'    
#' \itemize{
#'   
#' \item{ \code{\link{jSDM_binomial_probit}} : \cr
#' \cr
#' \bold{Ecological process:}
#' \deqn{y_{ij} \sim \mathcal{B}ernoulli(\theta_{ij})}{y_ij ~ Bernoulli(\theta_ij),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{ij}) = X_i \beta_j}{(\theta_ij) =  X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab probit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j}{(\theta_ij) =  X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i}  and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab probit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' 
#' \item{ \code{\link{jSDM_binomial_probit_sp_constrained}} : \cr
#' \cr
#' This function allows to fit the same models than the function \code{\link{jSDM_binomial_probit}} except for models not including latent variables, indeed \code{n_latent} must be greater than zero in this function. 
#' It aims to improve the convergence of latent variable models fitting by selecting the species constrained to positive values on the diagonal of the \eqn{\Lambda} matrix of factor loadings as the ones that structure themselves most clearly on each latent axis.
#' This function returns the fitted JSDM considering as constrained species those that maximize on each latent axis the Gelman–Rubin convergence diagnostic (\eqn{\hat{R}}{Rhat}) computed from the factor loadings \eqn{\lambda }
#' of a first JSDM adjusted with the constrained species arbitrarily chosen as the first ones in the presence-absence data-set.
#'}
#'
#'\item{ \code{\link{jSDM_binomial_logit}} : \cr
#' \cr
#' \bold{Ecological process : }
#' \deqn{y_{ij} \sim \mathcal{B}inomial(\theta_{ij},t_i)}{y_ij ~ Binomial(\theta_ij,t_i),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j}{(\theta_ij) =  X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j}{(\theta_ij) =  X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i}  \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab logit\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' 
#' \item{ \code{\link{jSDM_poisson_log}} : \cr
#' \cr
#' \bold{Ecological process : }
#' \deqn{y_{ij} \sim \mathcal{P}oisson(\theta_{ij})}{y_ij ~ Poisson(\theta_ij),}
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j}{(\theta_ij) =  X_i \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j}{(\theta_ij) =  X_i \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i}  \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j  + \alpha_i}{(\theta_ij) =  X_i \beta_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab log\eqn{(\theta_{ij}) =  X_i \beta_j + W_i \lambda_j + \alpha_i}{(\theta_ij) =  X_i  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' 
#' \item{ \code{\link{jSDM_binomial_probit_long_format}} : \cr
#' \cr
#' \bold{Ecological process:}
#' \deqn{y_n \sim \mathcal{B}ernoulli(\theta_n)}{y_n ~ Bernoulli(\theta_n)} such as \eqn{species_n=j} and \eqn{site_n=i},
#' where \tabular{ll}{
#'  if \code{n_latent=0} and \code{site_effect="none"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j} \cr
#'  if \code{n_latent>0} and \code{site_effect="none"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j +  W_i \lambda_j} \cr
#'  if \code{n_latent=0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j + \alpha_i}  and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#'  if \code{n_latent>0} and \code{site_effect="fixed"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n  \beta_j +  W_i \lambda_j + \alpha_i} \cr
#'  if \code{n_latent=0} and \code{site_effect="random"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n \beta_j + \alpha_i} \cr
#'  if \code{n_latent>0} and \code{site_effect="random"} \tab probit\eqn{(\theta_n) = D_n \gamma + X_n  \beta_j +  W_i \lambda_j + \alpha_i} and \eqn{\alpha_i \sim \mathcal{N}(0,V_\alpha)}{\alpha_i ~ N(0,V_\alpha)} \cr
#' }}
#' 
#' }
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr
#' Frédéric Gosselin <frederic.gosselin@inrae.fr>\cr
#'}
#' @keywords MCMC Metropolis algorithm binomial biodiversity logistic model multivariate regression
#' @importFrom Rcpp evalCpp
#' @importFrom graphics par
#' @importFrom stats median model.frame model.matrix pnorm quantile nobs
#' @useDynLib jSDM , .registration=TRUE
NULL