## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================
#'@name logit
#'@aliases  logit 
#'@title Generalized logit function
#'@description Compute generalized logit function.
#' @param x value(s) to be transformed
#' @param min Lower end of logit interval
#' @param max Upper end of logit interval 
#' @details The generalized logit function takes values on \eqn{[min, max]} and transforms them to span \eqn{[ -\infty, +\infty ]}{[-Inf,+Inf]} it is defined as:
#' \deqn{y = log(\frac{p}{(1-p)})}{y = log(p/(1-p))} 
#' \deqn{where}
#'\deqn{p=\frac{(x-min)}{(max-min)}}{p=(x-min)/(max-min)}
#' @return y Transformed value(s).
#'@author{Gregory R. Warnes <greg@warnes.net>}
#' @examples x <- seq(0,10, by=0.25)
#' xt <- jSDM::logit(x, min=0, max=10)
#' cbind(x,xt)
#' y <- jSDM::inv_logit(xt, min=0, max=10)
#' cbind(x,xt,y)  
#' @keywords math logistic logit
#' @export


logit <- function(x, min=0, max=1)
  {
    p <- (x-min)/(max-min)
    log(p/(1-p))
}
#'@name inv_logit
#'@aliases inv_logit
#'@title Generalized inverse logit function
#'@description Compute generalized inverse logit function.
#' @param x value(s) to be transformed
#' @param min Lower end of logit interval
#' @param max Upper end of logit interval 
#' @details The generalized inverse logit function takes values on [-Inf,Inf] and transforms them to span [min, max] :
#' \deqn{y = p' (max-min) + min}{y = p * (max-min) + min} 
#' \deqn{where}
#'\deqn{p =\frac{exp(x)}{(1+exp(x))}}{p =exp(x)/(1+exp(x))}
#' @return y Transformed value(s).
#'@author{Gregory R. Warnes <greg@warnes.net>}
#' @examples x <- seq(0,10, by=0.25)
#' xt <- jSDM::logit(x, min=0, max=10)
#' cbind(x,xt)
#' y <- jSDM::inv_logit(xt, min=0, max=10)
#' cbind(x,xt,y)  
#' @keywords math logistic logit
#' @export
inv_logit <- function(x, min=0, max=1)
  {
    p <- exp(x)/(1+exp(x))
    p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
    p * (max-min) + min
  }
                 
# End