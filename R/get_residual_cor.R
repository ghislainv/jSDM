## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# Calculate the residual correlation matrix from a LVM
get_residual_cor <- function(mod) {
  
  n.species <- ncol(mod$model_spec$presences)
  n.sites <- nrow(mod$model_spec$presences)
  n.mcmc <- nrow(mod$mcmc.latent[[1]])
  n.lv <- length(mod$mcmc.latent)
  n.X.coeff <- nrow(mod$model_spec$beta_start)
  Tau.arr <- matrix(NA,n.mcmc,n.species^2)
  Tau.cor.arr <- matrix(NA,n.mcmc,n.species^2)

  for(t in 1:n.mcmc) { 
    lv.coefs <- mod$mcmc.sp[["sp_1"]][t,(n.X.coeff+1):(n.X.coeff+n.lv)]
    for(j in 2:n.species) { 
    lv.coefs <- rbind(lv.coefs,mod$mcmc.sp[[paste0("sp_",j)]][t,(n.X.coeff+1):(n.X.coeff+n.lv)])
    }
    Tau.mat <- lv.coefs%*%t(lv.coefs) + diag(n.species)
    Tau.arr[t,] <- as.vector(Tau.mat) 
    Tau.cor.mat <- cov2cor(Tau.mat)
    Tau.cor.arr[t,] <- as.vector(Tau.cor.mat) 
  }
  
  ## Average/Median over the MCMC samples
  Tau.mat.mean <- sig.Tau.mat.mean <- matrix(apply(Tau.arr,2,mean),n.species,byrow=F)
  Tau.mat.median <- sig.Tau.mat.median <- matrix(apply(Tau.arr,2,median),n.species,byrow=F)
  Tau.cor.mean <- sig.Tau.cor.mean <- matrix(apply(Tau.cor.arr,2,mean),n.species,byrow=F)
  Tau.cor.median <- sig.Tau.cor.median <- matrix(apply(Tau.cor.arr,2,median),n.species,byrow=F)
  
  return(list(cov.mean = Tau.mat.mean, cov.median = Tau.mat.median, cor.mean = Tau.cor.mean, cor.median = Tau.cor.median))
}	

