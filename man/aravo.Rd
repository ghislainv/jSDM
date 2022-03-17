\name{aravo}
\alias{aravo}
\docType{data}
\title{Distribution of Alpine plants in Aravo (Valloire, France)}
\description{This dataset describe the distribution of 82 species of Alpine plants in 75 sites. Species traits and environmental variables are also measured.}
\usage{data("aravo")}
\format{
aravo is a list containing the following objects :
  \describe{
    \item{\code{spe}}{is a data.frame with the abundance values of 82 species (columns) in 75 sites (rows).}
    \item{\code{env}}{is a data.frame with the measurements of 6 environmental variables for the sites.}
    \item{\code{traits}}{is data.frame with the measurements of 8 traits for the species.}
    \item{\code{spe.names}}{is a vector with full species names.}
    }
}

\details{
  \describe{
The environmental variables are : 
 \tabular{ll}{
 \code{Aspect} \tab Relative south aspect (opposite of the sine of aspect with flat coded 0) \cr
 \code{Slope} \tab Slope inclination (degrees) \cr
    \code{Form} \tab Microtopographic landform index: 1 (convexity); 2 (convex slope); 3 (right slope); \cr
    \tab 4 (concave slope); 5 (concavity) \cr
    \code{Snow} \tab Mean snowmelt date (Julian day) averaged over 1997-1999 \cr
    \code{PhysD} \tab Physical disturbance, i.e., percentage of unvegetated soil due to physical processes \cr
    \code{ZoogD} \tab Zoogenic disturbance, i.e., quantity of unvegetated soil due to marmot activity: no; some; high \cr }
    
The species traits for the plants are:
 \tabular{ll}{
    \code{Height} \tab Vegetative height (cm) \cr
    \code{Spread} \tab Maximum lateral spread of clonal plants (cm) \cr
    \code{Angle} \tab Leaf elevation angle estimated at the middle of the lamina \cr
    \code{Area} \tab Area of a single leaf \cr
    \code{Thick} \tab Maximum thickness of a leaf cross section (avoiding the midrib) \cr
    \code{SLA} \tab Specific leaf area \cr
    \code{Nmass} \tab Mass-based leaf nitrogen content \cr
    \code{Seed} \tab Seed mass \cr } }
}
\source{
Choler, P. (2005) Consistent shifts in Alpine plant traits along a mesotopographical gradient. Arctic, Antarctic, and Alpine Research 37,444-453.

Dray S, Dufour A (2007). The ade4 Package: Implementing the Duality Diagram for Ecologists. \emph{Journal of Statistical Software}, 22(4), 1-20. \doi{10.18637/jss.v022.i04}. 
}

\examples{
data(aravo, package="jSDM")
summary(aravo)
}

\keyword{datasets}
