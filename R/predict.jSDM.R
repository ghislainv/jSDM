## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

#'@name predict.jSDM
#' @aliases predict.jSDM
#' @title Predict method for models fitted with jSDM
#' @description Prediction of species probabilities of occurrence from models fitted using the jSDM package
#'@param object An object of class \code{"jSDM"}.
#'@param newdata An optional data frame in which explanatory variables can be searched for prediction. If omitted, the adjusted values are used.
#'@param Id_species An vector of character or integer indicating for which species the probabilities of presence on chosen sites will be predicted.
#'@param Id_sites An vector of integer indicating for which sites the probabilities of presence of specified species will be predicted.
#'@param type Type of prediction. Can be \code{"mean"} for predictive posterior mean,
#' \code{"quantile"} for producing sample quantiles from the predictive posterior corresponding to the given probabilities (see \code{probs} argument)
#' or \code{"posterior"} for the full predictive posterior for each prediction.
#' Using \code{"quantile"} or \code{"posterior"} might lead to memory problem depending on the number of predictions and the number of samples for the jSDM model's parameters.
#'@param probs Numeric vector of probabilities with values in \eqn{[0,1]}, used when \code{type="quantile"}.
#'@param ... Further arguments passed to or from other methods.
#' @return Return a vector for the predictive posterior mean when \code{type="mean"}, a data-frame with the mean and quantiles when \code{type="quantile"} or an \code{mcmc} object (see \code{coda} package) with posterior distribution for each prediction when \code{type="posterior"}.
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr }
#' @seealso \code{\link{jSDM-package}} \code{\link{jSDM_binomial_logit}} \code{\link{jSDM_binomial_logit_one_species}}
#' @examples 
#' library(jSDM)
#' # frogs data
#' data(frogs, package="jSDM")
#'# Arranging data
#' PA_frogs <- frogs[,4:12]
#'# Normalized continuous variables
#'Env_frogs <- cbind(scale(frogs[,1]),frogs[,2],scale(frogs[,3]))
#' colnames(Env_frogs) <- colnames(frogs[,1:3])
#'# Parameter inference
#'# Increase the number of iterations to reach MCMC convergence
#'mod<-jSDM_binomial_probit_block_rand_site_lv(# Response variable 
#'                                            presence_site_sp = 
#'                                            PA_frogs, 
#'                                            # Explanatory variables 
#'                                            site_suitability = ~.,   
#'                                            site_data = Env_frogs,
#'                                            n_latent=2,
#'                                            # Chains
#'                                            burnin=100,
#'                                            mcmc=100,
#'                                            thin=1,
#'                                            # Starting values
#'                                            alpha_start=0,
#'                                            beta_start=0,
#'                                            lambda_start=0,
#'                                            W_start=0,
#'                                            V_alpha_start=1, 
#'                                            # Priors
#'                                            shape=0.5, rate=0.0005,
#'                                            mu_beta=0, V_beta=1.0E6,
#'                                            mu_lambda=0, V_lambda=10,
#'                                            # Various 
#'                                            seed=1234, verbose=1)
#'                                                
#'# Select site and species for predictions
#'## 30 sites
#'Id_sites <- sample.int(nrow(PA_frogs), 30)
#' ## 5 species
#' Id_species <- sample(colnames(PA_frogs), 5)
#' 
#'# Predictions 
#'theta_pred <- predict(mod,
#'                      Id_species=Id_species,
#'                      Id_sites=Id_sites,
#'                      type="mean")
#'hist(theta_pred, main="Predicted theta with simulated covariates")
#' @keywords prediction predictive posterior credible interval
#' @export predict.jSDM
#' @export 
#' 
#' 
predict.jSDM <- function(object, newdata=NULL, Id_species, Id_sites, type="mean", probs=c(0.025,0.975), ...) {
  
  ##= Check
  if (!(type %in% c("mean","quantile","posterior"))) {stop("type must be \"mean\", \"quantile\" or \"posterior\"")}
  if (sum(probs<0)>0 | sum(probs>1)>0) {stop("probs must be a vector of probabilities between (0,1)")}
  
  ##= Model specifications
  model.spec <- object$model_spec
  species <- colnames(model.spec$presences)
  
  ##= Link function probit for family=="binomial"
  if( model.spec$link=="probit") inv.link <- function (x) {pnorm(x)}
  if( model.spec$link=="logit") inv.link <- inv_logit
  
  ##= Matrix for predictions
  if (is.null(newdata)) {
    newdata <- model.spec$site_data[Id_sites,]
  }
  suitability <- model.spec$site_suitability
  mf.pred <- model.frame(formula=suitability,data=as.data.frame(newdata))
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
  if(!is.null(model.spec$n_latent)){
    nl <- model.spec$n_latent
  }
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
      term <- rep(0, npred)
      ##= Matrix of MCMC parameters
      if(is.null(model.spec$n_latent)){
        beta_j.mat <- as.matrix(object$mcmc.betas[[paste0("sp_",num_species[j])]])
      }
      if(!is.null(model.spec$n_latent)){
        beta_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,c(1:np)])
        lambda_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,(np+1):(np+nl)])
      }
      
      ##= Loop on samples
      for (t in 1:nsamp) {
        if(!is.null(model.spec$n_latent)){
          W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,Id_sites])
          for(l in 2:nl) {
            W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,Id_sites]))
          }
          lambda_j <- lambda_j.mat[t,]
        }
        
        beta_j <- beta_j.mat[t,]
        link.term <- X.pred %*% beta_j
        
        if(!is.null(model.spec$n_latent)){
          link.term <- link.term + W.mat %*% lambda_j 
        }
        
        if(!is.null(model.spec$alpha_start)){
          link.term <- link.term + as.matrix(object$mcmc.alpha[t,Id_sites])
        }
        
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
      
      if(is.null(model.spec$n_latent)){
        beta_j.mat <- as.matrix(object$mcmc.betas[[paste0("sp_",num_species[j])]])
      }
      
      if(!is.null(model.spec$n_latent)){
        ##= Matrix of MCMC parameters
        beta_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,c(1:np)])
        lambda_j.mat <- as.matrix(object$mcmc.sp[[paste0("sp_",num_species[j])]][,(np+1):(np+nl)])
      }
      
      ##= Loop on samples
      for (t in 1:nsamp) {
        if(!is.null(model.spec$n_latent)){
          W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,Id_sites])
          for(l in 2:nl) {
            W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,Id_sites]))
          }
          lambda_j <- lambda_j.mat[t,]
        }
        
        beta_j <- beta_j.mat[t,]
        link.term <- X.pred %*% beta_j
        if(!is.null(model.spec$n_latent)){
          link.term <- link.term + W.mat %*% lambda_j
        }
        if(!is.null(model.spec$alpha_start)){
          link.term <- link.term + as.matrix(object$mcmc.alpha[t,Id_sites])
        }
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