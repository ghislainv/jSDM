\name{mites}
\alias{mites}
\docType{data}

\title{mites dataset}

\description{
This example data set is composed of 70 cores of mostly Sphagnum mosses collected on the territory of the Station de biologie des Laurentides of University of  Montreal, Quebec, Canada in June 1989.   

The whole sampling area was 2.5 m x 10 m in size and thirty-five taxa were recognized as species, though many were not given a species name, owing to the incomplete stage of systematic knowledge of the North American Oribatid fauna.

The data set comprises the abundances of 35 morphospecies, 5 substrate and micritopographic variables, and the x-y Cartesian coordinates of the 70 sampling sites. 

See Borcard et al. (1992, 1994) for details.
}

\details{ 
Oribatid mites (Acari: Oribatida) are a very diversified group of small (0.2-1.2 mm) soil-dwelling, mostly microphytophagous and detritivorous arthropods. A well aerated soil or a complex substrate like Sphagnum mosses present in bogs and wet forests can harbour up to several hundred thousand individuals per square metre.   

Local assemblages are sometimes composed of over a hundred species, including many rare ones. This diversity makes oribatid mites an interesting target group to study community-environment relationships at very local scales.
}

\usage{data("mites")}
\format{ A data frame with 70 observations on the following 42 variables.
    \describe{
    Abundance of \bold{35 Oribatid mites morphospecies} named : 
     \item{\code{Brachy}}{a vector of integers}
     \item{\code{PHTH}}{a vector of integers}
     \item{\code{HPAV}}{a vector of integers}
     \item{\code{RARD}}{a vector of integers}
     \item{\code{SSTR}}{a vector of integers}
     \item{\code{Protopl}}{a vector of integers}
     \item{\code{MEGR}}{a vector of integers}
     \item{\code{MPRO}}{a vector of integers}
     \item{\code{TVIE}}{a vector of integers}
     \item{\code{HMIN}}{a vector of integers}
     \item{\code{HMIN2}}{a vector of integers}
     \item{\code{NPRA}}{a vector of integers}
     \item{\code{TVEL}}{a vector of integers}
     \item{\code{ONOV}}{a vector of integers}
     \item{\code{SUCT}}{a vector of integers}
     \item{\code{LCIL}}{a vector of integers}
     \item{\code{Oribatul1}}{a vector of integers}
     \item{\code{Ceratoz1}}{a vector of integers}
     \item{\code{PWIL}}{a vector of integers}
     \item{\code{Galumna1}}{a vector of integers}
     \item{\code{Steganacarus2}}{a vector of integers}
     \item{\code{HRUF}}{a vector of integers}
     \item{\code{Trhypochth1}}{a vector of integers}
     \item{\code{PPEL}}{a vector of integers}
     \item{\code{NCOR}}{a vector of integers}
     \item{\code{SLAT}}{a vector of integers}
     \item{\code{FSET}}{a vector of integers}
     \item{\code{Lepidozetes}}{a vector of integers}
     \item{\code{Eupelops}}{a vector of integers}
     \item{\code{Minigalumna}}{a vector of integers}
     \item{\code{LRUG}}{a vector of integers}
     \item{\code{PLAG2}}{a vector of integers}
     \item{\code{Ceratoz3}}{a vector of integers}
     \item{\code{Oppia.minus}}{a vector of integers}
     \item{\code{Trimalaco2}}{a vector of integers}  
     
     
     
     \bold{5 covariates} collected on the 70 sites and their coordinates :  
     \item{\code{substrate}}{a categorical vector indicating substrate type using a 7-level unordered factor :  \code{sph1}, \code{sph2}, \code{sph3}, \code{sph4}, \code{litter}, \code{peat} and \code{inter} for interface.}
     \item{\code{shrubs}}{a categorical vector indicating shrub density using a 3-level ordered factor : \code{None}, \code{Few} and \code{Many}}
     \item{\code{topo}}{a categorical vector indicating microtopography using a 2-level factor: \code{blanket} or \code{hummock}}
     \item{\code{density}}{a numeric vector indicating the substrate density (g/L)}
     \item{\code{water}}{a numeric vector indicating the water content of the substrate (g/L)}
     \item{\code{x}}{a numeric vector indicating first coordinates of sampling sites}
     \item{\code{y}}{a numeric vector indicating second coordinates of sampling sites}
    }
}

\source{ Pierre Legendre  }

\references{
Borcard, D.; Legendre, P. and Drapeau, P. (1992) Partialling out the spatial component of ecological variation. \emph{Ecology} 73: 1045-1055.

Borcard, D. and Legendre, P. (1994) Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). \emph{Environmental and Ecological Statistics} 1: 37-61. 

Borcard, D. and Legendre, P. (2002) All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. \emph{Ecological Modelling} 153: 51-68.
}


\examples{
data(mites, package="jSDM")
head(mites)
}

\keyword{datasets}