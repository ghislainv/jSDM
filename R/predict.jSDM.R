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
#'@param type Type of prediction. Can be : \tabular{l}{
#' \code{"mean"} for predictive posterior mean, \cr
#' \code{"quantile"} for producing sample quantiles from the predictive posterior corresponding to the given probabilities (see \code{probs} argument), \cr
#' \code{"posterior"} for the full predictive posterior for each prediction. \cr }
#' Using \code{"quantile"} or \code{"posterior"} might lead to memory problem depending on the number of predictions and the number of samples for the jSDM model's parameters.
#'@param probs Numeric vector of probabilities with values in [0,1], \cr
#' used when \code{type="quantile"}.
#'@param ... Further arguments passed to or from other methods.
#' @return Return a vector for the predictive posterior mean when \code{type="mean"}, a data-frame with the mean and quantiles when \code{type="quantile"} or an \code{mcmc} object (see \code{coda} package) with posterior distribution for each prediction when \code{type="posterior"}.
#' @author \tabular{l}{
#' Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>\cr
#' Jeanne Clément <jeanne.clement16@laposte.net>\cr }
#' @seealso \code{\link{jSDM-package}} \code{\link{jSDM_binomial_logit}}  \code{\link{jSDM_binomial_probit}} \code{\link{jSDM_poisson_log}}
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
#' # Increase the number of iterations to reach MCMC convergence
#' mod<-jSDM_binomial_probit(# Response variable
#'                           presence_data=PA_frogs,
#'                           # Explanatory variables
#'                           site_formula = ~.,
#'                           site_data = Env_frogs,
#'                           n_latent=2,
#'                           site_effect="random",
#'                           # Chains
#'                           burnin=100,
#'                           mcmc=100,
#'                           thin=1,
#'                           # Starting values
#'                           alpha_start=0,
#'                           beta_start=0,
#'                           lambda_start=0,
#'                           W_start=0,
#'                           V_alpha=1,
#'                           # Priors
#'                           shape=0.5, rate=0.0005,
#'                           mu_beta=0, V_beta=10,
#'                           mu_lambda=0, V_lambda=10,
#'                           # Various
#'                           seed=1234, verbose=1)
#' 
#'# Select site and species for predictions
#'## 30 sites
#'Id_sites <- sample.int(nrow(PA_frogs), 30)
#'## 5 species
#'Id_species <- sample(colnames(PA_frogs), 5)
#' 
#'# Predictions 
#'theta_pred <- predict(mod,
#'                      Id_species=Id_species,
#'                      Id_sites=Id_sites,
#'                      type="mean")
#'hist(theta_pred, main="Predicted theta with simulated covariates")
#' @keywords prediction predictive posterior credible interval
#' @importFrom stringi stri_remove_empty
#' @export 
#' 
#' 
predict.jSDM <- function(object, newdata=NULL, Id_species, Id_sites, type="mean", probs=c(0.025,0.975), ...) {
  ##= Check
  if (!class(object)=="jSDM"){
    stop("Please provide an object of class jSDM in", calling.function(),
         call.=FALSE)
  }
  if (!(type %in% c("mean","quantile","posterior"))) {stop("type must be \"mean\", \"quantile\" or \"posterior\"")}
  if (sum(probs<0)>0 | sum(probs>1)>0) {stop("probs must be a vector of probabilities between (0,1)")}
  
  ##= Model specifications
  model.spec <- object$model_spec
  if(!is.null(model.spec$presence_data)){
    species <- colnames(model.spec$presence_data)
    sites <- rownames(model.spec$presence_data)
  }
  if(!is.null(model.spec$count_data)){
    species <- colnames(model.spec$count_data)
    sites <- rownames(model.spec$count_data)
  }
  if(!is.null(model.spec$data)){
    species <- unique(model.spec$data$species)
    sites <- unique(model.spec$data$site)
  }
  
  ##= Link function probit for family=="binomial"
  if( model.spec$link=="probit") inv.link <- function (x) {pnorm(x)}
  if( model.spec$link=="logit") inv.link <- function (x) {inv_logit(x)}
  if( model.spec$link=="log") inv.link <- function (x) {exp(x)}
  
  ##= Matrix for predictions
  if(!is.null(model.spec$site_data)){
    if(is.null(newdata)) newdata <- model.spec$site_data[Id_sites,]
    if (!all(colnames(newdata) %in% colnames(model.spec$site_data)) && !(ncol(newdata)==ncol(model.spec$site_data))) {stop("newdata must have the same number of columns as object$model_spec$site_data and identical columns names\n")}
  }
  if(!is.null(model.spec$data)){
    nobs <- length(Id_sites)
    if (is.null(newdata)){
      for (i in 1:nobs){
        rowId_site <- which(model.spec$data$site==Id_sites[i] & model.spec$data$species==Id_species[i])
        newdata <- rbind(newdata, model.spec$data[rowId_site,!c(grepl("Y",colnames(model.spec$data)) | grepl("species",colnames(model.spec$data)) | grepl("site",colnames(model.spec$data)))])
      }
    }
    if (!all(colnames(newdata) %in% colnames(model.spec$data)) &&
        !(ncol(newdata)==(nrow(model.spec$beta_start) + ifelse(sum(grepl("(Intercept)", colnames(object$mcmc.sp[[1]])))==1,-1,0)))) 
      {stop("newdata must have columns names corresponding to names of all covariables in object$model_spec$data \n")}
  }
  newdata <- data.frame(newdata)
  
  if(!is.null(model.spec$data)){
    #= Suitability
    suitability <- model.spec$site_formula 
    if(model.spec$site_formula==~.) suitability <- ~. - site - Y
    mf.suit <- model.frame(formula=suitability, data=as.data.frame(newdata))
    # design matrix X for species effects beta
    Xterms <- stringi::stri_remove_empty(gsub(":?species:?", "", 
                                              grep("species", attr(attr(mf.suit,"terms"),"term.labels"), value=T)))  
    Xformula <- paste0("~",paste0(Xterms, collapse="+"))
    mf.suit.X <- model.frame(formula=Xformula, data=as.data.frame(newdata))
    attr(attr(mf.suit.X,"terms"),"intercept") <- ifelse(grepl("- *species", suitability[2]),0,1)
    X.pred <- model.matrix(attr(mf.suit.X,"terms"), data=mf.suit.X)
    # design matrix D for parameters gamma
    Dterms <- grep("species", grep("site", attr(attr(mf.suit,"terms"),"term.labels"), value=T, invert=T), value=T, invert=T)
    if(length(Dterms)!=0){
      Dformula <- paste0("~", paste0(Dterms, collapse="+"),"-1")
      mf.suit.D <- model.frame(formula=Dformula, data=as.data.frame(newdata))
      D.pred <- model.matrix(attr(mf.suit.D,"terms"), data=mf.suit.D)
      nd <- ncol(D.pred)
    }
  }
  
  if(!is.null(model.spec$presence_data) | !is.null(model.spec$count_data)){
    # Suitability process
    suitability <- model.spec$site_formula
    mf.pred <- model.frame(formula=suitability,data=as.data.frame(newdata))
    X.pred <- model.matrix(attr(mf.pred,"terms"),data=mf.pred)
  }
  
  ##= Model parameters
  np <- ncol(X.pred)
  npred <- nrow(X.pred)
  nsp <- length(unique(Id_species))
  nsite <- length(unique(Id_sites))
  
  if(model.spec$n_latent > 0){
    nl <- model.spec$n_latent
  }
  if (is.character(Id_species) || is.factor(Id_species)) {
    if(!is.null(model.spec$presence_data) | !is.null(model.spec$count_data)){
      num_species <- rep(0,nsp)
      for(j in 1:nsp) {
        num_species[j] <- which(species == Id_species[j])
      }
    }
    if(!is.null(model.spec$data)){
      num_species <- rep(0,npred)
      for(i in 1:npred) {
        num_species[i] <- which(species == Id_species[i])
      }
    }
  }
  
  if (is.numeric(Id_species)) {
    num_species <- Id_species
  }
  
  if (is.character(Id_sites) || is.factor(Id_sites)) {
    if(!is.null(model.spec$presence_data)| !is.null(model.spec$count_data)){
      num_sites <- rep(0,nsite)
      for(i in 1:nsite) {
        num_sites[i] <- which(sites == Id_sites[i])
      }
    }
    if(!is.null(model.spec$data)){
      num_sites <- rep(0,npred)
      for(i in 1:npred) {
        num_sites[i] <- which(sites == Id_sites[i])
      }
    }
  }
  
  if (is.numeric(Id_sites)) {
    num_sites <- Id_sites
  }
  
  # Chains parameters 
  ngibbs <- model.spec$mcmc + model.spec$burnin
  nthin <- model.spec$thin
  nburn <- model.spec$burnin
  nsamp <- model.spec$mcmc/nthin
  
  if(!is.null(model.spec$presence_data) | !is.null(model.spec$count_data)){
    ##= Posterior mean
    if (type=="mean") {
      term.pred <- matrix(0, npred, nsp)
      colnames(term.pred) <- species[num_species]
      rownames(term.pred) <- Id_sites
      
      ##= Loop on species
      for(j in 1:nsp) {
        term <- rep(0, npred)
        ##= Matrix of MCMC parameters
        if(model.spec$n_latent==0){
          if(length(model.spec$beta_start)==np){
            betaj.mat <- as.matrix(object$mcmc[,grepl("beta", colnames(object$mcmc))])
          }
          else{
            betaj.mat <- as.matrix(object$mcmc.sp[[num_species[j]]])
          }
        }
        if(model.spec$n_latent > 0){
          betaj.mat <- as.matrix(object$mcmc.sp[[num_species[j]]][,c(1:np)])
          lambda_j.mat <- as.matrix(object$mcmc.sp[[num_species[j]]][,(np+1):(np+nl)])
        }
        
        ##= Loop on samples
        for (t in 1:nsamp) {
          
          betaj <- betaj.mat[t,]
          link.term <- X.pred %*% as.vector(betaj)
          
          if(model.spec$n_latent > 0){
            W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,num_sites])
            for(l in 2:nl) {
              W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,num_sites]))
            }
            lambda_j <- lambda_j.mat[t,]
            link.term <- link.term + W.mat %*% lambda_j 
          }
          
          if(!is.null(model.spec$alpha_start)){
            link.term <- link.term + as.matrix(object$mcmc.alpha[t,num_sites])
          }
          
          term <- term + inv.link(link.term)
        }
        term.pred[,j] <- as.numeric(term/nsamp)
      }
    }
    
    ##= Full posterior
    if (type %in% c("quantile","posterior")) {
      term.pred <- list()
      ##= Loop on species
      for(j in 1:nsp) {
        term <- matrix(0,nsamp,npred)
        
        if(model.spec$n_latent==0){
          if(length(model.spec$beta_start)==np){
            betaj.mat <- as.matrix(object$mcmc[,grepl("beta", colnames(object$mcmc))])
          }
          else{
            betaj.mat <- as.matrix(object$mcmc.sp[[num_species[j]]])
          }
        }
        
        if(model.spec$n_latent > 0){
          ##= Matrix of MCMC parameters
          betaj.mat <- as.matrix(object$mcmc.sp[[num_species[j]]][,c(1:np)])
          lambda_j.mat <- as.matrix(object$mcmc.sp[[num_species[j]]][,(np+1):(np+nl)])
        }
        
        ##= Loop on samples
        for (t in 1:nsamp) {
          
          betaj <- betaj.mat[t,]
          link.term <- X.pred %*% as.vector(betaj)
          
          if(model.spec$n_latent > 0){
            W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,num_sites])
            for(l in 2:nl) {
              W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,num_sites]))
            }
            lambda_j <- lambda_j.mat[t,]
            link.term <- link.term + W.mat %*% lambda_j
          }
          
          if(!is.null(model.spec$alpha_start)){
            link.term <- link.term + as.matrix(object$mcmc.alpha[t,num_sites])
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
  }
  
  if(!is.null(model.spec$data)){
    ##= Posterior mean
    if (type=="mean") {
      term.pred <- rep(0,npred)
      ##= Loop on species
      for(i in 1:npred) {
        term <- 0
        ##= Matrix of MCMC parameters
        if(length(Dterms)!=0){
          gamma.mat <- as.matrix(object$mcmc.gamma)
        }
        if(model.spec$n_latent==0){
          if(length(model.spec$beta_start)==np){
            betaj.mat <- as.matrix(object$mcmc[,grepl("beta", colnames(object$mcmc))])
          }
          else{
            betaj.mat <- as.matrix(object$mcmc.sp[[num_species[i]]])
          }
        }
        if(model.spec$n_latent > 0){
          betaj.mat <- as.matrix(object$mcmc.sp[[num_species[i]]][,c(1:np)])
          lambda_j.mat <- as.matrix(object$mcmc.sp[[num_species[i]]][,(np+1):(np+nl)])
        }
        
        ##= Loop on samples
        for (t in 1:nsamp) {
  
          betaj <- betaj.mat[t,]
          link.term <- X.pred[i,] %*% as.vector(betaj)
          
          if(length(Dterms)!=0){
            gamma <- gamma.mat[t,]
            link.term <- link.term + D.pred[i,] %*% as.vector(gamma) 
          }
          
          if(model.spec$n_latent > 0){
            W_i <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,num_sites[i]])
            for(l in 2:nl) {
              W_i <- cbind(W_i, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,num_sites[i]]))
            }
            lambda_j <- lambda_j.mat[t,]
            link.term <- link.term + W_i %*% lambda_j 
          }
          
          if(!is.null(model.spec$alpha_start)){
            link.term <- link.term + object$mcmc.alpha[t,num_sites[i]]
          }
          
          term <- term + inv.link(link.term)
        }
        term.pred[i]<- as.numeric(term/nsamp)
      }
    }
    
    ##= Full posterior
    if (type %in% c("quantile","posterior")) {
      ##= Loop on species
      for(i in 1:npred) {
        term <- rep(0,nsamp)
        
        if(length(Dterms)!=0){
          gamma.mat <- as.matrix(object$mcmc.gamma)
        }
        
        if(model.spec$n_latent==0){
          if(length(model.spec$beta_start)==np){
            betaj.mat <- as.matrix(object$mcmc[,grepl("beta", colnames(object$mcmc))])
          }
          else{
            betaj.mat <- as.matrix(object$mcmc.sp[[num_species[i]]])
          }
        }
        
        if(model.spec$n_latent > 0){
          ##= Matrix of MCMC parameters
          betaj.mat <- as.matrix(object$mcmc.sp[[num_species[i]]][,c(1:np)])
          lambda_j.mat <- as.matrix(object$mcmc.sp[[num_species[i]]][,(np+1):(np+nl)])
        }
        
        ##= Loop on samples
        for (t in 1:nsamp) {
          
          betaj <- betaj.mat[t,]
          link.term <- X.pred[i,] %*% as.vector(betaj)
          
          if(length(Dterms)!=0){
            gamma <- gamma.mat[t,]
            link.term <- link.term + D.pred[i,] %*% as.vector(gamma) 
          }
          
          if(model.spec$n_latent > 0){
            W.mat <- as.matrix(object$mcmc.latent[[paste0("lv_",1)]][t,num_sites[i]])
            for(l in 2:nl) {
              W.mat <- cbind(W.mat, as.matrix(object$mcmc.latent[[paste0("lv_",l)]][t,num_sites[i]]))
            }
            lambda_j <- lambda_j.mat[t,]
            link.term <- link.term + W.mat %*% lambda_j
          }
          
          if(!is.null(model.spec$alpha_start)){
            link.term <- link.term + as.matrix(object$mcmc.alpha[t,num_sites[i]])
          }
          
          term[t] <-  inv.link(link.term)
        }
        
        if (type=="quantile") {
          term.mean <- mean(term)
          term.quant <- quantile(term,probs)
          term.pred[[i]] <- data.frame(theta_pred=c(term.mean,term.quant), row.names=c("mean", names(term.quant)))
        }
        
        if (type=="posterior") {
          term.pred[[i]] <- coda::mcmc(term,start=nburn+1,end=ngibbs,thin=nthin)
        }
      }
    }
  }
  ## Output
  return(term.pred)
}

# End