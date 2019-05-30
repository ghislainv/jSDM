\name{get_residual_cor}

\alias{get_residual_cor}

\title{Calculate the residual correlation matrix from a LVM.}

\description{This function use coefficients \eqn{(\lambda_j^l) with j=1,\dots,n_species and l=1,\dots,n_latent}$ corresponding to latent variables to calculate the variance-covariance matrix which controls correlation between species.}

\usage{get_residual_cor(mod)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{ 
    An object of class \code{"jSDM"}
  }
}

\value{
\item{cov.mean}{Average over the MCMC samples of the variance-covariance matrix.}
\item{cov.median}{Median over the MCMC samples of the variance-covariance matrix.}
\item{cor.mean}{Average over the MCMC samples of the residual correlation matrix.}
\item{cor.median}{Median over the MCMC samples of the residual correlation matrix.}
}

\author{
  Jeanne Clement \email{jeanne.clement16@laposte.net}
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{boral}% use one of  RShowDoc("KEYWORDS")
\keyword{stats::cov2cor}
