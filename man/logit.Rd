\name{logit}
\alias{logit}
\alias{inv_logit}
\title{Generalized logit and inverse logit function}
\description{
  Compute generalized logit and generalized inverse logit functions.
}
\usage{
logit(x, min = 0, max = 1)
inv_logit(x, min = 0, max = 1)
}
\arguments{
  \item{x}{value(s) to be transformed}
  \item{min}{Lower end of logit interval}
  \item{max}{Upper end of logit interval}
}
\details{

  The generalized logit function takes values on [min, max] and
  transforms them to span [-Inf,Inf] it is defined as:

    \deqn{y = log(\frac{p}{(1-p)})}{y = log(p/(1-p))}

  where
    
    \deqn{p=\frac{(x-min)}{(max-min)}}{p=(x-min)/(max-min)}

  The generalized inverse logit function provides the inverse
  transformation:
    
    \deqn{x = p' (max-min) + min}{x = p * (max-min) + min}

  where

    \deqn{p'=\frac{exp(y)}{(1+exp(y))}}{exp(y)/(1+exp(y))}


}
\value{
  Transformed value(s).
}
\author{Gregory R. Warnes <greg@warnes.net>}
\examples{
  x <- seq(0,10, by=0.25)
  xt <- jSDM::logit(x, min=0, max=10)
  cbind(x,xt)
  y <- jSDM::inv_logit(xt, min=0, max=10)
  cbind(x,xt,y)  
}
\keyword{math}
