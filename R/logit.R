## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================


logit <- function(x, min=0, max=1)
  {
    p <- (x-min)/(max-min)
    log(p/(1-p))
  }

inv_logit <- function(x, min=0, max=1)
  {
    p <- exp(x)/(1+exp(x))
    p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
    p * (max-min) + min
  }
                 
# End