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

\details{ Oribatid mites (Acari: Oribatida) are a very diversified group of small (0.2-1.2 mm) soil-dwelling, mostly microphytophagous and detritivorous arthropods. A well aerated soil or a complex substrate like Sphagnum mosses present in bogs and wet forests can harbour up to several hundred thousand individuals per square metre. Local assemblages are sometimes composed of over a hundred species, including many rare ones. This diversity makes oribatid mites an interesting target group to study community-environment relationships at very local scales.
}

\usage{data("mites")}
\format{ A data frame with 70 observations on the following 42 variables.

Abundance of 35 Oribatid mites morphospecies named : 

\code{Brachy}, \code{PHTH}, \code{HPAV}, \code{RARD}, \code{SSTR}, \code{Protopl}, \code{MEGR}, \code{MPRO}, \code{TVIE}, \code{HMIN},\code{HMIN2}, \code{NPRA}, \code{TVEL}, \code{ONOV}, \code{SUCT}, \code{LCIL}, \code{Oribatul1}, \code{Ceratoz1}, \code{PWIL}, \code{Galumna1}, \code{Steganacarus2}, \code{HRUF}, \code{Trhypochth1}, \code{PPEL}, \code{NCOR}, \code{SLAT}, \code{FSET}, \code{Lepidozetes}, \code{Eupelops}, \code{Minigalumna}, \code{LRUG}, \code{PLAG2}, \code{Ceratoz3}, \code{Oppia.minus} and \code{Trimalaco2}.
    
    5 covariates collected on the 70 sites and their coordinates : 
    
    \tabular{ll}{
    \bold{\code{substrate}} \tab  a categorical vector indicating substrate type using a 7-level unordered factor : \code{sph1}, \code{sph2}, \code{sph3}, \code{sph4}, \code{litter}, \code{peat} and \code{inter} for interface. \cr
    \tab  \cr
    \bold{\code{shrubs}} \tab a categorical vector indicating shrub density using a 3-level ordered factor : \code{None}, \code{Few} and \code{Many} \cr
    \tab  \cr
    \bold{\code{topo}} \tab  a categorical vector indicating microtopography using a factor with levels \code{blanket} and \code{hummock} \cr
    \tab  \cr
    \bold{\code{density}} \tab   a numeric vector indicating the substrate density (g/L)\cr
    \tab \cr
    \bold{\code{water}} \tab  a numeric vector indicating the water content of the substrate (g/L) \cr
    \tab \cr
    \bold{\code{x}} \tab  a numeric vector indicating first coordinates of sampling sites \cr 
            \tab \cr
    \bold{\code{y}} \tab  a numeric vector indicating second coordinates of sampling sites  \cr 
    \tab \cr 
}
}

\source{ Pierre Legendre  }

\references{
\tabular{l}{
Borcard, D.; Legendre, P. and Drapeau, P. (1992) Partialling out the spatial component of ecological variation. \emph{Ecology} 73: 1045-1055.\cr
\cr
Borcard, D. and Legendre, P. (1994) Environmental control and spatial structure in ecological communities: an example using Oribatid mites (Acari, Oribatei). \emph{Environmental and Ecological Statistics} 1: 37-61. \cr
\cr
Borcard, D. and Legendre, P. (2002) All-scale spatial analysis of ecological data by means of principal coordinates of neighbour matrices. \emph{Ecological Modelling} 153: 51-68. \cr
}}


\examples{
data(mites, package="jSDM")
head(mites)
}

\keyword{datasets}

