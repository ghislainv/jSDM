## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

jSDM_probit_block <- function (presence_site_sp, site_suitability,
							   site_data, n_latent=2,
							   burnin=5000, mcmc=10000, thin=10,
							   alpha_start=0, beta_start=0, lambda_start=0, W_start=0,
							   V_alpha_start=1, shape=0.5, rate=0.0005,
							   mu_beta=0, V_beta=1.0E6,
							   mu_lambda=0, V_lambda=10,
							   seed=1234, verbose=1)

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
	Y <- presence_site_sp
	nsp <- ncol(Y)
	nsite <- nrow(Y)
	nobs <- nsite*nsp
	T <- matrix(1, nsite, nsp)
	#= Suitability
	mf.suit <- model.frame(formula=site_suitability, data=site_data)
	X <- model.matrix(attr(mf.suit,"terms"), data=mf.suit)
	np <- ncol(X)
	#= Iterations
	ngibbs <- mcmc+burnin
	nthin <- thin
	nburn <- burnin
	nsamp <- mcmc/thin
	
	#========== 
	# Check data
	#==========
	check.T.binomial(c(T), nobs)
	check.Y.binomial(c(Y), c(T))
	check.X(X, nsite)
	
	#========
	# Initial starting values for M-H
	#========
	beta_start <- form.beta.start.sp(beta_start, np, nsp)
	lambda_start <- form.lambda.start.sp(lambda_start, n_latent, nsp)
	
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
	
	#= Model specification
	model_spec <- list(presences=presences, trials=trials,
					   suitability=suitability,
					   data=data,
					   burnin=burnin, mcmc=mcmc, thin=thin,
					   beta_start=beta_start, mubeta=mubeta, Vbeta=Vbeta,
					   seed=seed, ropt=ropt, verbose=verbose)
	
	#= Output
	output <- list(mcmc=MCMC, theta_latent=mod$theta_latent,
				   family="binomial",
				   model_spec=model_spec)
	
	class(output) <- "jSDM"
	return(output)
	
}

# End
