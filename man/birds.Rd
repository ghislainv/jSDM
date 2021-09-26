\name{birds}
\alias{birds}
\docType{data}
\title{birds dataset}

\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of 158 common species since 1999. The MHB sample from \code{data(MHB2014, package="AHMbook")} consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol, where each square is surveyed up to three times during the breeding season (only twice above the tree line). Surveys are conducted along a transect that does not change over the years. 

The \code{birds} dataset has the data for 2014, except one quadrat not surveyed in 2014. 
}


\usage{data("birds")}
\format{ A data frame with 266 observations on the following 166 variables.

158 bird species named in latin and whose occurrences are indicated as the number of visits to each site during which the species was observed, including 13 species not recorded in the year 2014 and 5 covariates collected on the 266 1x1 km quadrat as well as their identifiers and coordinates : 
    
    \tabular{ll}{
    \bold{\code{siteID}} \tab an alphanumeric site identifier \cr
    \tab \cr
    \bold{\code{coordx}} \tab a numeric vector indicating the x coordinate of the centre of the quadrat; the coordinate reference system intentionally not specified. \cr
    \tab \cr
    \bold{\code{coordy}} \tab a numeric vector indicating the y coordinate of the centre of the quadrat.\cr
        \tab  \cr
    \bold{\code{elev}} \tab a numeric vector indicating the mean elevation of the quadrat (m).\cr
        \tab  \cr
    \bold{\code{rlength}} \tab  the length of the route walked in the quadrat (km).\cr
    \tab  \cr
    \bold{\code{nsurvey}} \tab  a numeric vector indicating the number of replicate surveys planned in the quadrat; above the tree-line 2, otherwise 3. \cr 
            \tab \cr
    \bold{\code{forest}} \tab a numeric vector indicating the percentage of forest cover in the quadrat.\cr 
    \tab \cr 
    \bold{\code{obs14}} \tab a categorical vector indicating the identifying number of the observer.\cr 
}
}

\source{Swiss Ornithological Institute}

\details{
Only the Latin names of bird species are given in this dataset but you can find the corresponding English names in the original dataset : code{data(MHB2014, package="AHMbook")}. 
}
\references{
\tabular{l}{
Kéry and Royle (2016) Applied Hierarachical Modeling in Ecology Section 11.3 \cr
}}


\examples{
data(birds, package="jSDM")
head(birds)
# find species not recorded in 2014
which(colSums(birds[,1:158])==0)
}

\keyword{datasets}