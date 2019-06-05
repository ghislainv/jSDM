## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

jSDM_probit_block <- function (presence_site_sp, site_suitability,
                               site_data, n_latent=2,
                               burnin=5000, mcmc=15000, thin=10,
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
  alpha_start <- form.alpha.start.sp(alpha_start, nsite)
  W_start <-form.W.start.sp(W_start, nsite, n_latent)
  param_start = rbind(beta_start,lambda_start)
  
  #========
  # Form and check priors
  #========
  mubeta <- check.mubeta(mu_beta,np)
  Vbeta <- check.Vbeta.mat(V_beta,np)
  mulambda <- check.mubeta(mu_lambda,n_latent)
  Vlambda <- check.Vlambda.mat(V_lambda,n_latent)
  Vparam <- diag(c(diag(Vbeta),diag(Vlambda)))
  muparam <- c(mubeta,mulambda)
  VW <- diag(rep(1,n_latent))
  V_alpha_start <- check.Valpha(V_alpha_start)
  
  #========
  # call Rcpp function
  #========
  mod <- Rcpp_jSDM_probit_block(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
                                Y=Y,T=T,X=X,
                                param_start= param_start, Vparam=Vparam, muparam = muparam,
                                W_start=W_start, VW=VW,
                                alpha_start=alpha_start, Valpha_start=V_alpha_start,
                                shape = shape, rate = rate,
                                seed=seed, verbose=verbose)

  
  
  
  #= Transform Sample list in an MCMC object
  MCMC.Deviance <- coda::mcmc(mod$Deviance,start=nburn+1,end=ngibbs,thin=nthin)
  MCMC.alpha <- coda::mcmc(mod$alpha,start=nburn+1,end=ngibbs,thin=nthin)
  MCMC.Valpha <- coda::mcmc(mod$Valpha,start=nburn+1,end=ngibbs,thin=nthin)
  MCMC.sp <- list()
  for (j in 1:nsp) {
    ## beta_j
    MCMC.beta_j <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
    colnames(MCMC.beta_j) <- paste0("beta_",colnames(X))
    ## lambda_j
    MCMC.lambda_j <- coda::mcmc(mod$param[,j,(np+1):(n_latent+np)], start=nburn+1, end=ngibbs, thin=nthin)	
    colnames(MCMC.lambda_j) <- paste0("lambda_",1:n_latent)
    
    MCMC.sp[[paste0("sp_",j)]] <- coda::as.mcmc(cbind(MCMC.beta_j, MCMC.lambda_j),start=nburn+1, end=ngibbs, thin=nthin)
  }
  ## W latent variables 
  MCMC.latent <- list()
  for (l in 1:n_latent) {
   MCMC.lv_l <- coda::mcmc(mod$W[,,l], start=nburn+1, end=ngibbs, thin=nthin)
   MCMC.latent[[paste0("lv_",l)]] <- MCMC.lv_l
    }
  
  #= Model specification, site_suitability,
  model_spec <- list(presences=presence_site_sp,
                     site_suitability=site_suitability,
                     site_data=site_data, n_latent=n_latent,
                     burnin=burnin, mcmc=mcmc, thin=thin,
                     beta_start=beta_start, mubeta=mubeta, Vbeta=Vbeta,
                     lambda_start=lambda_start, mulambda=mulambda, Vlambda=Vlambda,
                     alpha_start=alpha_start, V_alpha_start=V_alpha_start, 
                     W_start=W_start, VW=VW,
                     seed=seed, verbose=verbose)
  
  #= Output
  output <- list(mcmc.Deviance=MCMC.Deviance,
                 mcmc.alpha = MCMC.alpha, mcmc.Valpha = MCMC.Valpha,
                 mcmc.sp = MCMC.sp, mcmc.latent = MCMC.latent,
                 Z_latent=mod$Z_latent, 
                 probit_theta_pred=mod$probit_theta_pred,
                 model_spec=model_spec)
  
  class(output) <- "jSDM"
  return(output)
  
}

# End
