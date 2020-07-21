#' @docType package
#' @name jSDM-package
#' @aliases jSDM-package
#' @title joint species distribution models
#' @description \code{jSDM} is an R package for fitting joint species distribution models (jSDM) in a hierarchical Bayesian framework.
#' The Gibbs sampler is written in C++. It uses Rcpp, Armadillo and GSL to maximize computation efficiency.
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