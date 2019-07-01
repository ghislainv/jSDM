## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

predict.jSDM <- function(object, newdata=NULL, Id_species, Id_sites, type="mean", probs=c(0.025,0.975), ...) {
  
  ##= Check
  if (!(type %in% c("mean","quantile","posterior"))) {stop("type must be \"mean\", \"quantile\" or \"posterior\"")}
  if (sum(probs<0)>0 | sum(probs>1)>0) {stop("probs must be a vector of probabilities between (0,1)")}
  
  ##= Link function probit for family=="binomial"
  inv.link <- function (x) {pnorm(x)}
  
  ##= Model specifications
  model.spec <- object$model_spec
  species <- colnames(model.spec$presences)
  
  ##= Matrix for predictions
  if (is.null(newdata)) {
    newdata <- model.spec$site_data[Id_sites,]
  }
  suitability <- model.spec$site_suitability
  mf.pred <- model.frame(formula=suitability,data=newdata)
  X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
  npred <- nrow(X.pred)
  
  ##= Model parameters
  np <- ncol(X.pred)
  nsp <- length(Id_species)
  species <- colnames(model.spec$presences)
  if (is.character(Id_species)) {
    num_species <- rep(0,nsp)
    for(j in 1:nsp) {
      num_species[j] <- which(species == Id_species[j])
    }
  }
  if (is.numeric(Id_species)) {
    num_species <- Id_species
  }
  nl <- model.spec$n_latent
  ngibbs <- model.spec$mcmc + model.spec$burnin
  nthin <- model.spec$thin
  nburn <- model.spec$burnin
  nsamp <- model.spec$mcmc/nthin
  ##= Posterior mean
  if (type=="mean") {
    term.pred <- matrix(0, npred, nsp)
    colnames(term.pred) <- species[num_species]
    rownames(term.pred) <- Id_sites
    ##= Loop on species
    for(j in 1:nsp) {
      ##= Matrix of MCMC parameters
      beta_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,c(1:np)])
      lambda_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,(np+1):(np+nl)])
      term <- rep(0, npred)
      
      ##= Loop on samples
      for (t in 1:nsamp) {
        
        W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,Id_sites])
        for(l in 2:nl) {
          W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,Id_sites]))
        }
        beta_j <- beta_j.mat[t,]
        lambda_j <- lambda_j.mat[t,]
        
        link.term <- X.pred %*% beta_j + W.mat %*% lambda_j + as.matrix(object$mcmc.alpha[t,Id_sites])
        term <- term + inv.link(link.term)
      }
      term.pred[,j]<- as.numeric(term/nsamp)
    }
  }
  
  ##= Full posterior
  if (type %in% c("quantile","posterior")) {
    term.pred <- list()
    ##= Loop on species
    for(j in 1:nsp) {
      term <- matrix(0,nsamp,npred)
      ##= Matrix of MCMC parameters
      beta_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,c(1:np)])
      lambda_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,(np+1):(np+nl)])
      
      ##= Loop on samples
      for (t in 1:nsamp) {
        W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,Id_sites])
        for(l in 2:nl) {
          W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,Id_sites]))
        }
        beta_j <- beta_j.mat[t,]
        lambda_j <- lambda_j.mat[t,]
        
        link.term <- X.pred %*% beta_j + W.mat %*% lambda_j + as.matrix(object$mcmc.alpha[t,Id_sites])
        term[t,] <-  inv.link(link.term)
      }
      if (type=="quantile") {
        term.mean <- apply(term,2,mean)
        term.quant <- apply(term,2,quantile,probs)
        term.pred[[Id_species[j]]] <- as.data.frame(t(rbind(term.mean,term.quant)))
        names(term.pred[[Id_species[j]]])[1] <- c("mean")
        rownames(term.pred[[Id_species[j]]]) <- Id_sites
      }
      if (type=="posterior") {
        term.pred[[Id_species[j]]] <- coda::mcmc(term,start=nburn+1,end=ngibbs,thin=nthin)
        colnames(term.pred[[Id_species[j]]]) <- Id_sites
      }
    }
  }
  
  ## Output
  return(term.pred)
}

# End