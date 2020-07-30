\name{mite}
\alias{mite}
\docType{data}

\title{mite dataset}

\description{
This example data set is composed of 70 cores of mostly Sphagnum mosses collected on the territory of the Station de biologie des Laurentides of Université de Montréal, Québec, Canada in June 1989. 

The whole sampling area was 2.5 m x 10 m in size and thirty-five taxa were recognized as species, though many were not given a species name, owing to the incomplete stage of systematic knowledge of the North American Oribatid fauna.

The data set comprises the abundances of 35 morphospecies, 5 substrate and micritopographic variables, and the x-y Cartesian coordinates of the of the 70 sampling sites. 

See Borcard et al. (1992, 1994) for details.
}

\details{ Oribatid mites (Acari: Oribatida) are a very diversified group of small (0.2–1.2 mm) soil-dwelling, mostly microphytophagous and detritivorous arthropods. A well aerated soil or a complex substrate like Sphagnum mosses present in bogs and wet forests can harbour up to several hundred thousand individuals per square metre. Local assemblages are sometimes composed of over a hundred species, including many rare ones. This diversity makes oribatid mites an interesting target group to study community–environment relationships at very local scales.
}

\usage{data("mite")}
\format{
  A data frame with 70 observations on the following 42 variables.
  
  \describe{
\item{\code{Abundance of 35 Oribatid mites morphospecies named:}}{Brachy, PHTH, HPAV, RARD, SSTR, Protopl, MEGR, MPRO, TVIE, HMIN, HMIN2, NPRA, TVEL, ONOV, SUCT, LCIL, Oribatl1, Ceratoz1, PWIL, Galumna1, Stgncrs2, HRUF, Trhypch1, PPEL, NCOR, SLAT, FSET, Lepidzts, Eupelops, Miniglmn, LRUG, PLAG2, Ceratoz3, Oppiminu, Trimalc2.}
    
    5 covariates collected on the 70 sites and their coordinates : 

    \item{\code{SubsDens}}{a numeric vector indicating the substrate density (g/L)}
    \item{\code{WatrCont}}{a numeric vector indicating the water content of the substrate (g/L)}
    \item{\code{Substrate}}{a categorical vector indicating substrate type using a 7-level unordered factor:Sphagn1, Sphagn2, Sphagn3, Sphagn, Litter, Barepeat and Interface.}
    \item{\code{Shrub}}{a categorical vector indicating shrub density using a 3-level ordered factor: None, Few and Many}
    \item{\code{Topo}}{a categorical vector indicating microtopography using a factor with levels Blanket and Hummock}
    \item{\code{x}}{a numeric vector indicating first coordinates of sampling sites}
    \item{\code{y}}{a numeric vector indicating second coordinates of sampling sites}
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
data(mite, package="jSDM")
head(mite)
}

\keyword{datasets}

