# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# license         :GPLv3
# ==============================================================================

jSDM_binomial_fit <- function (X, Y, T,
						   burnin=5000, mcmc=10000, thin=10,
						   beta_start,
						   mubeta=0, Vbeta=1.0E6,
						   seed=1234, ropt=0.44, verbose=1)
	
{   
	#========
	# Basic checks
	#========
	check.mcmc.parameters(burnin, mcmc, thin)
	check.verbose(verbose)
	
	#======== 
	# Form response, covariate matrices and model parameters
	#========
	
	#= Response
	Y <- as.numeric(Y)
	T <- as.numeric(T)
	#= Suitability
	X <- as.matrix(X)
	np <- ncol(X)
	#= Iterations
	ngibbs <- mcmc+burnin
	nthin <- thin
	nburn <- burnin
	nsamp <- mcmc/thin
	
	#========== 
	# Check data
	#==========
	check.T.binomial(T,nobs)
	check.Y.binomial(Y,T)
	check.X(X,nobs)
	
	#========
	# Initial starting values for M-H
	#========
	beta_start <- form.beta.start(beta_start,np)
	
	#========
	# Form and check priors
	#========
	mubeta <- check.mubeta(mubeta,np)
	Vbeta <- check.Vbeta(Vbeta,np)
	
	#========
	# call Rcpp function
	#========
	mod <- Rcpp_jSDM_binomial(ngibbs, nthin, nburn,
							  Y, T, X, beta_start, mubeta, Vbeta,
							  seed, ropt, verbose)
	
	#= Matrix of MCMC samples
	Matrix <- cbind(mod$beta, mod$Deviance)
	names.fixed <- paste("beta_",colnames(X),sep="")
	colnames(Matrix) <- c(names.fixed,"Deviance")
	
	#= Transform Sample list in an MCMC object
	MCMC <- coda::mcmc(Matrix,start=nburn+1,end=ngibbs,thin=nthin)
	
	#= Output
	output <- list(mcmc=MCMC, theta_latent=mod$theta_latent)
	
	return(output)
	
}

# End