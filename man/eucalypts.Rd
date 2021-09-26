\name{eucalypts}

\alias{eucalypts}

\docType{data}

\title{eucalypts dataset}

\description{
The Eucalyptus data set includes 12 taxa recorded in 458 plots spanning elevation gradients in the Grampians National Park, Victoria, which is known for high species diversity and endemism. The park has three mountain ranges interspersed with alluvial valleys and sand sheet and has a semi-Mediterranean climate with warm, dry summers and cool, wet winters.

This dataset records presence or absence at 458 sites of 12 eucalypts species, 7 covariates collected at these sites as well as their longitude and latitude.}

\usage{data("eucalypts")}
\format{
  A data frame with 458 observations on the following 21 variables.
  
  \describe{
  
    12 eucalypts species which presence on sites is indicated by a 1 and absence by a 0 : 
    
    \item{\code{ALA}}{a binary vector indicating the occurrence of the species E. alaticaulis}
    \item{\code{ARE}}{a binary vector indicating the occurrence of the species E. arenacea}
    \item{\code{BAX}}{a binary vector indicating the occurrence of the species E. baxteri}
    \item{\code{CAM}}{a binary vector indicating the occurrence of the species E. camaldulensis}
    \item{\code{GON}}{a binary vector indicating the occurrence of the species E. goniocalyx}
    \item{\code{MEL}}{a binary vector indicating the occurrence of the species E. melliodora}
    \item{\code{OBL}}{a binary vector indicating the occurrence of the species E. oblique}
    \item{\code{OVA}}{a binary vector indicating the occurrence of the species E. ovata}
    \item{\code{WIL}}{a binary vector indicating the occurrence of the species E. willisii subsp. Falciformis}
    \item{\code{ALP}}{a binary vector indicating the occurrence of the species E. serraensis, E. verrucata and E. victoriana}
    \item{\code{VIM}}{a binary vector indicating the occurrence of the species E. viminalis subsp. Viminalis and Cygnetensis}
    \item{\code{ARO.SAB}}{a binary vector indicating the occurrence of the species E. aromaphloia and E. sabulosa}
    
    7 covariates collected on the 458 sites and their coordinates : 

    \item{\code{Rockiness}}{a numeric vector taking values from 0 to 95 corresponding to the rock cover of the site in percent estimated in 5 \% increments in field plots}
    \item{\code{Sandiness}}{a binary vector indicating if soil texture categorie is sandiness based on soil texture classes from field plots and according to relative amounts of sand, silt, and clay particles}
    \item{\code{VallyBotFlat}}{a numeric vector taking values from 0 to 6 corresponding to the valley bottom flatness GIS-derived variable defining flat areas relative to surroundings likely to accumulate sediment (units correspond to the percentage of slope e.g. 0.5 = 16 \%slope, 4.5 = 1 \%slope, 5.5 = 0.5 \%slope)}
    \item{\code{PPTann}}{a numeric vector taking values from 555 to 1348 corresponding to annual precipitation in millimeters measured as the sum of monthly precipitation estimated using BIOCLIM based on 20m grid cell Digital Elevation Model}
    \item{\code{Loaminess}}{a binary vector indicating if soil texture categorie is loaminess based on soil texture classes from field plots and according to relative amounts of sand, silt, and clay particles}
    \item{\code{cvTemp}}{a numeric vector taking values from 136 to 158 corresponding to coefficient of variation of temperature seasonality in percent measured as the standard deviation of weekly mean temperatures as a percentage of the annual mean temperature from BIOCLIM}
    \item{\code{T0}}{a numeric vector corresponding to solar radiation in \eqn{WH/m^2} measured as the amount of incident solar energy based on the visible sky and the sun's position. Derived from Digital Elevation Model in ArcGIS 9.2 Spatial Analyst for the summer solstice (December 22)}
    \item{\code{latitude}}{a numeric vector indicating the latitude of the studied site}
    \item{\code{longitude}}{a numeric vector indicating the longitude of the studied site}
  }
}

\source{
Wilkinson, D. P.; Golding, N.; Guillera-Arroita, G.; Tingley, R. and McCarthy, M. A. (2018) A comparison of joint species distribution models for presence-absence data. \emph{Methods in Ecology and Evolution}.}

\examples{
data(eucalypts, package="jSDM")
head(eucalypts)
}

\keyword{datasets}

