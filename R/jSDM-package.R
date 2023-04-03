#' @docType package
#' @name jSDM-package
#' @aliases jSDM-package
#' @title joint species distribution models
#' @description \code{jSDM} is an R package for fitting joint species distribution models (JSDM) in a hierarchical Bayesian framework.
#'
#' The Gibbs sampler is written in 'C++'. It uses 'Rcpp', 'Armadillo' and 'GSL' to maximize computation efficiency.
#'\tabular{ll}{
#'    Package: \tab jSDM\cr
#'    Type: \tab Package\cr
#'    Version: \tab 0.2.1\cr
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
#' At first, the function fit a JSDM with the constrained species arbitrarily chosen as the first ones in the presence-absence data-set.
#' Then, the function evaluates the convergence of MCMC \eqn{\lambda} chains using the Gelman-Rubin convergence diagnostic (\eqn{\hat{R}}).
#' It identifies the species (\eqn{\hat{j}_l}) having the higher \eqn{\hat{R}} for \eqn{\lambda_{\hat{j}_l}}.
#' These species drive the structure of the latent axis \eqn{l}.
#' The \eqn{\lambda} corresponding to this species are constrained to be positive and placed on the diagonal of the \eqn{\Lambda} matrix for fitting a second model.
#' This usually improves the convergence of the latent variables and factor loadings. The function returns the parameter posterior distributions for this second model.
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
#' @author 
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#' 
#' Jeanne Clément <jeanne.clement16@laposte.net>
#' 
#' Frédéric Gosselin <frederic.gosselin@inrae.fr>
#'
#'@references
#' Chib, S. and Greenberg, E. (1998) Analysis of multivariate probit models. \emph{Biometrika}, 85, 347-361. 
#' 
#' Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.
#'
#' Ovaskainen, O., Tikhonov, G., Norberg, A., Blanchet, F. G., Duan, L., Dunson, D., Roslin, T. and Abrego, N. (2017) How to make more out of community data? A conceptual framework and its implementation as models and software. \emph{Ecology Letters}, 20, 561-576.
#' 
#' @keywords MCMC Metropolis algorithm binomial biodiversity logistic model multivariate regression
#' @importFrom Rcpp evalCpp
#' @importFrom graphics par
#' @importFrom stats median model.frame model.matrix pnorm quantile nobs
#' @useDynLib jSDM , .registration=TRUE
NULL