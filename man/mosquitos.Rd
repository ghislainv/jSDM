\name{mosquitos}

\alias{mosquitos}

\docType{data}

\title{mosquitos dataset}

\description{
Presence or absence at 167 sites of 16 species that constitute the aquatic faunal community studied, 13 covariates collected at these sites and their coordinates.}

\usage{data("mosquitos")}

\format{
  A data frame with 167 observations on the following 31 variables : 
  \describe{
  16 aquatic species including larvae of four mosquito species (all potential vectors of human disease), which presence on sites is indicated by a 1 and absence by a 0 : 
    \item{\code{Culex_pipiens_sl}}{a binary vector (mosquito species)}
    \item{\code{Culex_modestus}}{a binary vector (mosquito species)}
    \item{\code{Culiseta_annulata}}{a binary vector (mosquito species)}
    \item{\code{Anopheles_maculipennis_sl}}{a binary vector (mosquito species)}
    \item{\code{waterboatmen__Corixidae}}{a binary vector}
    \item{\code{diving_beetles__Dysticidae}}{a binary vector}
    \item{\code{damselflies__Zygoptera}}{a binary vector}
    \item{\code{swimming_beetles__Haliplidae}}{a binary vector}
    \item{\code{opossum_shrimps__Mysidae}}{a binary vector}
    \item{\code{ditch_shrimp__Gammarus}}{a binary vector}
    \item{\code{beetle_larvae__Coleoptera}}{a binary vector}
    \item{\code{dragonflies__Anisoptera}}{a binary vector}
    \item{\code{mayflies__Ephemeroptera}}{a binary vector}
    \item{\code{newts__Pleurodelinae}}{a binary vector}
    \item{\code{fish}}{a binary vector}
    \item{\code{saucer_bugs__Ilyocoris}}{a binary vector}
    
    
    13 covariates collected on the 167 sites and their coordinates : 
    
    \item{\code{depth__cm}}{a numeric vector corresponding to the water depth in cm recorded as the mean of the depth at the edge and the centre of each dip site}
    \item{\code{temperature__C}}{a numeric vector corresponding to the temperature in °C}
    \item{\code{oxidation_reduction_potential__Mv}}{a numeric vector corresponding to the redox potential of the water in millivolts (mV)}
    \item{\code{salinity__ppt}}{a numeric vector corresponding to the salinity of the water in parts per thousand (ppt)}
        
    High-resolution digital photographs were taken of vegetation at the edge and centre dip points and the presence or absence of different vegetation types at each dipsite was determined from these photographs using field guides :
    \item{\code{water_crowfoot__Ranunculus}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of the plant species Ranunculus aquatilis which common name is water-crowfoot}
    \item{\code{rushes__Juncus_or_Scirpus}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of rushes from the Juncus or Scirpus genus}
    \item{\code{filamentous_algae}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of filamentous algae}
    \item{\code{emergent_grass}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of emergent grass}
    \item{\code{ivy_leafed_duckweed__Lemna_trisulca}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of ivy leafed duckweed of species Lemna trisulca}
    \item{\code{bulrushes__Typha}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of bulrushes from the Typha genus}
    \item{\code{reeds_Phragmites}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of reeds from the Phragmites genus}
    \item{\code{marestail__Hippuris}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of plants from the Hippuris genus known as mare's-tail}
    \item{\code{common_duckweed__Lemna_minor}}{a binary vector indicating presence on sites by a 1 and absence by a 0 of common duckweed of species Lemna minor}
    \item{\code{x}}{a numeric vector of first coordinates corresponding to each site}
    \item{\code{y}}{a numeric vector of second coordinates corresponding to each site}
  }
}

\source{
Wilkinson, D. P.; Golding, N.; Guillera-Arroita, G.; Tingley, R. and McCarthy, M. A. (2018) A comparison of joint species distribution models for presence-absence data. \emph{Methods in Ecology and Evolution}.}

\examples{
data(mosquitos, package="jSDM")
head(mosquitos)
}

\keyword{datasets}
