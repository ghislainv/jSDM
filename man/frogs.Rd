\name{frogs}

\alias{frogs}

\docType{data}

\title{frogs dataset}

\description{Presence or absence of 9 species of frogs on 104 sites, 3 covariates   collected on those site and their coordinates.
}

\usage{frogs}

\format{
  \code{frogs} is a data frame with 104 observations on the following 14 variables.
  \describe{
    \item{\code{Species_}}{1 to 9 indicate by a 0 the absence of the species on one site and by a 1 its presence} 
    \item{\code{Covariates_}}{1 and 3 continuous variables}
    \item{\code{Covariates_}}{2 discrete variables}
    \item{\code{x}}{a numeric vector of first coordinates corresponding to each site}
    \item{\code{y}}{a numeric vector of second coordinates corresponding to each site}
  }
}

\source{
Wilkinson, D. P.; Golding, N.; Guillera-Arroita, G.; Tingley, R. and McCarthy, M. A. (2018) \emph{A comparison of joint species distribution models for presence-absence data.} Methods in Ecology and Evolution.}

\examples{
data(frogs, package="jSDM")
head(frogs)
}

\keyword{datasets}
