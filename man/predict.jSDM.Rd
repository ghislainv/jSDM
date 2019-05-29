\name{predict.jSDM}
\alias{predict.jSDM}

\title{Predict method for models fitted with jSDM}

\description{Predicted values for models fitted with jSDM}

\usage{\method{predict}{jSDM}(object,newdata=NULL,type="mean",probs=c(0.025,0.975),...)}

\arguments{
  
  \item{object}{An object of class \code{"jSDM"}.}

  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}

  \item{type}{Type of prediction. Can be \code{"mean"} for predictive posterior mean, \code{"quantile"} for producing sample quantiles from the predictive posterior corresponding to the given probabilities (see \code{probs} argument) or \code{"posterior"} for the full predictive posterior for each prediction. Using \code{"quantile"} or \code{"posterior"} might lead to memory problem depending on the number of predictions and the number of samples for the jSDM model's parameters.}
  
  \item{probs}{Numeric vector of probabilities with values in [0,1] and used when \code{type="quantile"}.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{Return a vector for the predictive posterior mean when \code{type="mean"}, a data-frame with the mean and quantiles when \code{type="quantile"} or an \code{mcmc} object (see \code{coda} package) with posterior distribution for each prediction when \code{type="posterior"}.} 

\author{
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}
}

\seealso{
  \code{\link{jSDM}}
}

\keyword{prediction}
\keyword{predictive posterior}
\keyword{credible interval}