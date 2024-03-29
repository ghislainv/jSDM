## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#' @name plot_residual_cor
#' @aliases plot_residual_cor
#' @title Plot the residual correlation matrix from a latent variable model (LVM).
#' @description Plot the posterior mean estimator of residual correlation matrix reordered by first principal component using \code{\link[corrplot]{corrplot}} function from the package of the same name.
#' @param mod An object of class \code{"jSDM"}.
#' @param prob A numeric scalar in the interval \eqn{(0,1)} giving the target probability coverage of the intervals, by which to determine whether the correlations are "significant".
#'   If \code{prob=0.95} is specified only significant correlations, whose \eqn{95\%} HPD interval does not contain zero, are represented. 
#'   Defaults to \code{prob=NULL} to represent all correlations significant or not.
#' @param main Character, title of the graph.
#' @param cex.main Numeric, title's size. 
#' @param diag Logical, whether display the correlation coefficients on the principal diagonal.
#' @param type Character, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.
#' @param method Character, the visualization method of correlation matrix to be used. Currently, it supports seven methods, named "circle" (default), "square", "ellipse", "number", "pie", "shade" and "color". 
#' @param mar See \code{\link[graphics]{par}}
#' @param tl.cex Numeric, for the size of text label (variable names).
#' @param tl.srt Numeric, for text label string rotation in degrees, see \code{\link[graphics]{text}}.
#' @param ... Further arguments passed to \code{\link[corrplot]{corrplot}} function
#' @return No return value. Displays a reordered correlation matrix.
#' @author 
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
#' 
#' Jeanne Clément <jeanne.clement16@laposte.net>
#' 
#' @seealso \code{\link[corrplot]{corrplot}} \code{\link{jSDM-package}} \code{\link{jSDM_binomial_probit}} \cr
#'          \code{\link{jSDM_binomial_logit}} \code{\link{jSDM_poisson_log}} 
#' @examples 
#' library(jSDM)
#' # frogs data
#' data(frogs, package="jSDM")
#' # Arranging data
#' PA_frogs <- frogs[,4:12]
#' # Normalized continuous variables
#'  Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
#'  colnames(Env_frogs) <- colnames(frogs[,1:3])
#' # Parameter inference
#' # Increase the number of iterations to reach MCMC convergence
#' mod<-jSDM_binomial_probit(# Response variable
#'                            presence_data = PA_frogs,
#'                            # Explanatory variables
#'                            site_formula = ~.,
#'                            site_data = Env_frogs,
#'                            n_latent=2,
#'                            site_effect="random",
#'                            # Chains
#'                            burnin=100,
#'                            mcmc=100,
#'                            thin=1,
#'                            # Starting values
#'                            alpha_start=0,
#'                            beta_start=0,
#'                            lambda_start=0,
#'                            W_start=0,
#'                            V_alpha=1,
#'                            # Priors
#'                            shape=0.1, rate=0.1,
#'                            mu_beta=0, V_beta=1,
#'                            mu_lambda=0, V_lambda=1,
#'                            # Various
#'                            seed=1234, verbose=1)
#' # Representation of residual correlation between species
#' plot_residual_cor(mod)
#' plot_residual_cor(mod, prob=0.95)
#' @references
#' Taiyun Wei and Viliam Simko (2017). R package "corrplot": Visualization of a Correlation Matrix (Version 0.84)   
#' 
#' Warton, D. I.; Blanchet, F. G.; O'Hara, R. B.; O'Hara, R. B.; Ovaskainen, O.; Taskinen, S.; Walker, S. C. and Hui, F. K. C. (2015) So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution}, 30, 766-779.
#'
#' @keywords corrplot
#' @importFrom corrplot corrplot
#' @export 

## Plot of residual correlation matrix (posterior mean estimator, and reordered by first principal component)
plot_residual_cor <- function(mod, prob=NULL,
                              main = "Residual Correlation Matrix from LVM",
                              cex.main= 1.5,
                              diag = FALSE, type = "lower",
                              method = "color", 
                              mar = c(1,1,3,1),
                              tl.srt = 45, tl.cex = 0.5, ...) {
  lv2.cor <- get_residual_cor(mod, prob=ifelse(is.null(prob), 0.95, prob))
  # All non-significant correlations are set to zero, 
  # according to the (100 x prob)% HPD interval.
  if(!is.null(prob)){
  lv2.cor$cor.mean <- lv2.cor$cor.mean * lv2.cor$cor.sig
  }
  lv2.cor$reorder.cor.mean <- corrplot::corrMatOrder(lv2.cor$cor.mean, order = "FPC", hclust.method = "average")
  if(!is.null(mod$model_spec$presence_data)){
  rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- colnames(mod$model_spec$presence_data)
  }
  if(!is.null(mod$model_spec$count_data)){
    rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- colnames(mod$model_spec$count_data)
  }
  if(!is.null(mod$model_spec$data)){
    rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- rownames(lv2.cor$cor.mean) <- colnames(lv2.cor$cor.mean) <- unique(mod$model_spec$data$species)
  }
  oldpar <- par(no.readonly = TRUE) 
  on.exit(par(oldpar))
  par(cex=1, cex.main=cex.main)
  corrplot::corrplot(lv2.cor$cor.mean[lv2.cor$reorder.cor.mean,lv2.cor$reorder.cor.mean], diag = diag,
                     type = type, title = main, mar = mar, method = method , tl.srt = tl.srt, tl.cex = tl.cex, ...)
  }

# End