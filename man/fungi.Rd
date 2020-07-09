\name{fungi}

\alias{fungi}

\docType{data}

\title{fungi dataset}

\description{
Presence or absence of 11 species of fungi on dead-wood objects at 800 sites and 12 covariates collected at these sites.}

\usage{data("fungi")}
\format{
  A data frame with 800 observations on the following 23 variables :
  \describe{
  
  11 fungi species which presence on sites is indicated by a 1 and absence by a 0 :
    \item{\code{antser}}{a binary vector}
    \item{\code{antsin}}{a binary vector}
    \item{\code{astfer}}{a binary vector}
    \item{\code{fompin}}{a binary vector}
    \item{\code{hetpar}}{a binary vector}
    \item{\code{junlut}}{a binary vector}
    \item{\code{phefer}}{a binary vector}
    \item{\code{phenig}}{a binary vector}
    \item{\code{phevit}}{a binary vector}
    \item{\code{poscae}}{a binary vector}
    \item{\code{triabi}}{a binary vector}
    
   12 covariates collected on the 800 sites : 

    \item{\code{diam}}{a numeric vector indicating the diameter of dead-wood object}
    \item{\code{dc1}}{a binary vector indicating if the decay class is 1 measured in the scale 1, 2, 3, 4, 5 (from freshly  decayed to almost completely decayed)}
    \item{\code{dc2}}{a binary vector indicating if the decay class is 2}
    \item{\code{dc3}}{a binary vector indicating if the decay class is 3}
    \item{\code{dc4}}{a binary vector indicating if the decay class is 4}
    \item{\code{dc5}}{a binary vector indicating if the decay class is 5}
    \item{\code{quality3}}{a binary vector indicating if the quality is level 3}
    \item{\code{quality4}}{a binary vector indicating if the quality is level 4}
    \item{\code{ground3}}{a binary vector indicating if the ground contact is level 3 as 2 = no ground contact, 3 = less than half of the log in ground contact and 4 = more than half of the log in ground contact}
    \item{\code{ground4}}{a binary vector a binary vector indicating if the ground contact is level 4}
    \item{\code{epi}}{a numeric vector indicating the epiphyte cover}
    \item{\code{bark}}{a numeric vector indicating the bark cover}
  }
}

\source{
Wilkinson, D. P.; Golding, N.; Guillera-Arroita, G.; Tingley, R. and McCarthy, M. A. (2018) A comparison of joint species distribution models for presence-absence data. \emph{Methods in Ecology and Evolution}.}

\examples{
data(fungi, package="jSDM")
head(fungi)
}

\keyword{datasets}
