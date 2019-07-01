\name{predict.jSDM}
\alias{predict.jSDM}

\title{Predict method for models fitted with jSDM}

\description{Predicted values for models fitted with jSDM}

\usage{\method{predict}{jSDM}(object, newdata=NULL, Id_species, Id_sites, type="mean",
  probs=c(0.025,0.975), ...)}

\arguments{
  
  \item{object}{An object of class \code{"jSDM"}.}

  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  
  \item{Id_species}{An vector of character or integer indicating for which species the probabilities of presence on chosen sites will be predicted.}
  
  \item{Id_sites}{An vector of integer indicating for which sites the probabilities of presence of specified species will be predicted.}

  \item{type}{Type of prediction. Can be \code{"mean"} for predictive posterior mean, \code{"quantile"} for producing sample quantiles from the predictive posterior corresponding to the given probabilities (see \code{probs} argument) or \code{"posterior"} for the full predictive posterior for each prediction. Using \code{"quantile"} or \code{"posterior"} might lead to memory problem depending on the number of predictions and the number of samples for the jSDM model's parameters.}
  
  \item{probs}{Numeric vector of probabilities with values in [0,1] and used when \code{type="quantile"}.}
  
  \item{...}{Further arguments passed to or from other methods.}

}

\value{Return a vector for the predictive posterior mean when \code{type="mean"}, a data-frame with the mean and quantiles when \code{type="quantile"} or an \code{mcmc} object (see \code{coda} package) with posterior distribution for each prediction when \code{type="posterior"}.} 

\author{
  Jeanne Clement \email{jeanne.clement16@laposte.net}
  Ghislain Vieilledent \email{ghislain.vieilledent@cirad.fr}}

\seealso{
  \code{\link{jSDM}}
}

\examples{

# frogs data
data(frogs, package="jSDM")

# Arranging data
PA_frogs <- frogs[,4:12]

# Normalized continuous variables
Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
colnames(Env_frogs) <- colnames(frogs[,1:3])

# Parameter inference
# Increase the number of iterations to reach MCMC convergence
mod_jSDM_block_frogs <- jSDM::jSDM_probit_block(
  # Response variable 
  presence_site_sp = as.matrix(PA_frogs), 
  # Explanatory variables 
  site_suitability = ~.,   
  site_data = as.data.frame(Env_frogs), n_latent=2,
  # Chains
  burnin=1000, mcmc=1000, thin=1,
  # Starting values
  alpha_start=0, beta_start=0,
  lambda_start=0, W_start=0,
  V_alpha_start=1, 
  # Priors
  shape=0.5, rate=0.0005,
  mu_beta=0, V_beta=1.0E6,
  mu_lambda=0, V_lambda=10,
  # Various 
  seed=1234, verbose=1)
  
# Select site and species for predictions
## 30 sites
Id_sites <- sample.int(nrow(PA_frogs), 30)
## 5 species
Id_species <- sample(colnames(PA_frogs), 5)

# Predictions 
theta_pred <- jSDM::predict.jSDM(mod_jSDM_block_frogs,
  Id_species=Id_species, Id_sites=Id_sites, type="mean")
hist(theta_pred, main="Predicted theta with simulated covariates")

}
 
\keyword{prediction}
\keyword{predictive posterior}
\keyword{credible interval}
