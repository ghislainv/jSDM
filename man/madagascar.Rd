\name{madagascar}

\alias{madagascar}

\docType{data}

\title{Madagascar's forest inventory}

\description{Dataset compiled from the national forest inventories carried out on 753 sites on the island of Madagascar, listing the presence or absence of 555 plant species on each of these sites between 1994 and 1996.
We use these forest inventories to calculate a matrix indicating the presence by a 1 and the absence by a 0 of the species at each site by removing observations for which the species is not identified. This presence-absence matrix therefore records the occurrences of 483 species at 751 sites.
}

\format{
  \code{madagascar} is a data frame with 751 rows corresponding to the inventory sites and 483 columns corresponding to the species whose presence or absence has been recorded on the sites.
  \describe{
    \item{\code{sp_}}{1 to 483 indicate by a 0 the absence of the species on one site and by a 1 its presence} 
    \item{\code{site}}{"1" to "753" inventory sites identifiers.}
  }
}

\examples{
data(madagascar, package="jSDM")
head(madagascar)
}

\keyword{datasets}
