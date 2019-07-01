\name{plot_residual_cor}

\alias{plot_residual_cor}

\title{Plot the residual correlation matrix from a LVM.}

\description{
Plot the posterior mean estimator of residual correlation matrix reordered by first principal component.}

\usage{
plot_residual_cor(mod)
}

\arguments{
  \item{mod}{
    An object of class \code{"jSDM"}.
  }
}

\author{
  Jeanne Clement \email{jeanne.clement16@laposte.net}
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
mod_jSDM_block_frogs <- jSDM::jSDM_probit_block (
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
  
# Representation of residual correlation between species 
jSDM::plot_residual_cor(mod_jSDM_block_frogs)

}

\keyword{corrplot}
