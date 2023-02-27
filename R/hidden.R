## ==============================================================================
## author          :Ghislain Vieilledent, Jeanne Clément
## email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# =======================================================================
#
# calling.function
# return name of the calling function
#
# =======================================================================

calling.function <- function(parentheses=TRUE) {
  calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
  if (parentheses){
    calling.function <- paste(calling.function, "()", sep="")
  }
  return(calling.function)
}

# =======================================================================
#
# Check mcmc parameters
#
# =======================================================================

check.mcmc.parameters <- function(burnin, mcmc, thin) {
  
  if(mcmc %% thin != 0) {
    cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(mcmc <= 0) {
    cat("Error: MCMC iterations must be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
  }
  if(burnin < 0) {
    cat("Error: Burnin iterations must be positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(((burnin+mcmc) %% 10 != 0) || (burnin+mcmc)<100) {
    cat("Error: Value 'burnin+mcmc' should be divisible by 10 and >= 100.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(thin < 1) {
    cat("Error: Thinning interval must be superior or equal to 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

# =======================================================================
#
# Check verbose
#
# =======================================================================

check.verbose <- function (verbose) {
  if (!(verbose %in% c(0,1))) {
    cat("Error: verbose must take value 0 or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

# =======================================================================
#
# Check data
#
# =======================================================================

check.Y.poisson <- function (Y) {
  if(!is.numeric(Y)) {
    cat("Error: 'counts' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(Y))>0) {
    cat("Error: 'counts' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(Y<0 | Y%%1!=0)>0) {
    cat("Error: 'counts' must be a vector of positive integers.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.T.binomial <- function (T, nobs) { # = This version of the function assumes T>0
  if(length(T)!=nobs) {
    cat("Error: 'trials' must have the same length as the response variable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if(!is.numeric(T)) {
    cat("Error: 'trials' must be a vector of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(T))>0) {
    cat("Error: 'trials' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(T<=0 | T%%1!=0)>0) {
    cat("Error: 'trials' must be a vector of integers superior to zero.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.Y.binomial <- function (Y,T) {
  if(!is.numeric(Y)) {
    cat("Error: 'presences' must be a vector or a matrix of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(Y))>0) {
    cat("Error: 'presences' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(Y<0 | Y%%1!=0)>0) {
    cat("Error: 'presences' must be a vector of positive integers.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  } 
  if (sum(Y>T)>0) {
    cat("Error: 'presences' must be less than or equal to 'trials'.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

# check.U <- function (U,nobs) {
#   if(length(U)!=nobs) {
#     cat("Error: 'alteration' must have the same length as the response variable.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if(!is.numeric(U)) {
#     cat("Error: 'alteration' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(U))>0) {
#     cat("Error: 'alteration' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (!all(U>=0 && U<=1)) {
#     cat("Error: 'alteration' must be in the interval [0,1].\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(0)
# }

check.X <- function (X,n) {
  if(!is.numeric(c(X))) {
    cat("Error: 'suitability' only accept vectors of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(c(X)))>0) {
    cat("Error: 'suitability' do not accept vectors with missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (dim(X)[1]!=n) {
    cat("Error: Incorrect vector length for the 'suitability' argument.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}

check.Y.gaussian <- function(Y) {
  if(!is.numeric(Y)) {
    cat("Error: 'response_data' must be a vector or a matrix of numeric values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(is.na(Y))>0) {
    cat("Error: 'response_data' must not contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(0)
}
# check.W <- function (W,n) {
#   if(!is.numeric(c(W))) {
#     cat("Error: 'observability' only accept vectors of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(c(W)))>0) {
#     cat("Error: 'observability' do not accept vectors with missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (dim(W)[1]!=n) {
#     cat("Error: Incorrect vector length for the 'observability' argument.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(0)
# }

# check.neighbors <- function (n.neighbors,ncell,neighbors) {
#   # Length of n.neighbors=ncell
#   if(length(n.neighbors)!=ncell) {
#     cat("Error: 'n.neighbors' must have a length equal to the number of cells.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Numeric values
#   if(!is.numeric(neighbors)) {
#     cat("Error: 'neighbors' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if(!is.numeric(n.neighbors)) {
#     cat("Error: 'n.neighbors' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # NA
#   if (sum(is.na(neighbors))>0) {
#     cat("Error: 'neighbors' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(n.neighbors))>0) {
#     cat("Error: 'n.neighbors' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Positive integer
#   if (sum(neighbors<=0 | neighbors%%1!=0)>0) {
#     cat("Error: 'neighbors' must be a vector of integers superior to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(n.neighbors<=0 | n.neighbors%%1!=0)>0) {
#     cat("Error: 'n.neighbors' must be a vector of integers superior to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Number of neighbors inferior to ncell-1
#   if (!all(n.neighbors < (ncell-1))) {
#     cat("Error: 'n.neighbors' must not contain values superior to ncell-1.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Number of neighbors and length of neighbors
#   if (sum(n.neighbors)!=length(neighbors)) {
#     cat("Error: 'neighbors' must be a vector of length equal to sum(n.neighbors).\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Check values in neighbors
#   if (sum(!(neighbors %in% c(1:ncell)))>0) {
#     cat("Error: 'neighbors' must be a vector of integers between 1 and ncell.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   # Check that the target cell is not in the list of neighbors: --> ToDoList
#   return(0)
# }
# 
# check.sites <- function (sites,nobs) {
#   if(length(sites)!=nobs) {
#     cat("Error: 'sites' must have the same length as the response variable.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if(!is.numeric(sites)) {
#     cat("Error: 'sites' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(sites))>0) {
#     cat("Error: 'sites' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(sites<=0 | sites%%1!=0)>0) {
#     cat("Error: 'sites' must be a vector of integers superior to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(0)
# }
# 
# check.cells <- function (cells,nsite) {
#   if(length(cells)!=nsite) {
#     cat("Error: 'spatial.entity' must be of length equals to the number of sites.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if(!is.numeric(cells)) {
#     cat("Error: 'spatial.entity' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(cells))>0) {
#     cat("Error: 'spatial.entity' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(cells<=0 | cells%%1!=0)>0) {
#     cat("Error: 'spatial.entity' must be a vector of integers superior to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(0)
# }
# 
# check.cells.pred <- function (cells.pred,npred) {
#   if(length(cells.pred)!=npred) {
#     cat("Error: 'spatial.entity.pred' must be of length equals to the number of predictions.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if(!is.numeric(cells.pred)) {
#     cat("Error: 'spatial.entity.pred' must be a vector of numeric values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(is.na(cells.pred))>0) {
#     cat("Error: 'spatial.entity.pred' must not contain missing values.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (sum(cells.pred<=0 | cells.pred%%1!=0)>0) {
#     cat("Error: 'spatial.entity.pred' must be a vector of integers superior to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(0)
# }
# 
# # =======================================================================
# #
# # Check and form starting parameters
# #
# # =======================================================================
# 
#
# check.Vrho.start <- function (Vrho.start) {
#   if (length(Vrho.start)!=1) {
#     cat("Error: Vrho.start should be a scalar.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (!is.numeric(Vrho.start)) {
#     cat("Error: Vrho.start should be a numeric.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (Vrho.start<=0) {
#     cat("Error: Vrho.start should be strictly positive.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(Vrho.start)
# }
#
# =======================================================================
#
# Check and form priors
#
# =======================================================================

# check.Vbeta <- function(Vbeta, np) {
#   if (!all(Vbeta>0)) {
#     cat("Error: Vbeta should be strictly positive.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (is.null(dim(Vbeta))) {
#     Vbeta <- rep(Vbeta,np) 
#   }
#   else if (length(Vbeta)!=np) {
#     cat("Error: Vbeta not conformable.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(Vbeta)
# }

# check.mugamma <- function(mugamma, nq) {
#   if (is.null(dim(mugamma))) {
#     mugamma <- rep(mugamma,nq) 
#   }
#   else if (length(mugamma)!=nq) {
#     cat("Error: mugamma not conformable.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(mugamma)
# }

# check.Vgamma <- function(Vgamma, nq) {
#   if (!all(Vgamma>0)) {
#     cat("Error: Vgamma should be strictly positive.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (is.null(dim(Vgamma))) {
#     Vgamma <- rep(Vgamma,nq) 
#   }
#   else if (length(Vgamma)!=nq) {
#     cat("Error: Vgamma not conformable.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(Vgamma)
#}
#
# check.ig.prior <- function(nu, delta) {
#   
#   if(nu <= 0) {
#     cat("Error: in IG(nu,delta) prior, nu less than or equal to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.\n",
#          call.=FALSE)
#   }
#   if(delta <= 0) {
#     cat("Error: in IG(nu,delta) prior, delta less than or equal to zero.\n")
#     stop("Please respecify and call ", calling.function(), " again.\n",
#          call.=FALSE)    
#   }
#   return(0)
# }
# 
# check.Vrho.max <- function (Vrho.max) {
#   if (length(Vrho.max)!=1) {
#     cat("Error: Vrho.max should be a scalar.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (!is.numeric(Vrho.max)) {
#     cat("Error: Vrho.max should be a numeric.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   if (Vrho.max<=0) {
#     cat("Error: Vrho.max should be strictly positive.\n")
#     stop("Please respecify and call ", calling.function(), " again.",
#          call.=FALSE)
#   }
#   return(Vrho.max)
# }
# 
# form.priorVrho <- function (priorVrho) {
#   if (is.numeric(priorVrho[1]) && priorVrho[1] > 0.0) {
#     priorVrho <- priorVrho[1]
#   }
#   else if (priorVrho=="Uniform") {
#     priorVrho <- -2.0
#   }
#   else if (priorVrho=="1/Gamma") {
#     priorVrho <- -1.0
#   }
#   else {
#     priorVrho <- -1.0
#     cat("priorVrho has been set to \"1/Gamma\" \n")
#   }
#   return(priorVrho)
# }

# =======================================================================
#
# Check and form starting parameters 
#
# =======================================================================

is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1L
}

form.beta.start <- function (beta.start, np) {
  if (sum(is.na(beta.start))>0) { 
    beta.start <- rep(0,np)
  }
  else if(!is.na(beta.start)[1] && length(beta.start)!=np) {
    beta.start <- rep(beta.start[1],np) 
  }
  else if(length(beta.start)!=np) {
    cat("Error: beta.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(beta.start)
}

form.gamma.start <- function (gamma.start,nq) {
  if (sum(is.na(gamma.start))>0) {
    gamma.start <- rep(0,nq)
  }
  else if(!is.na(gamma.start)[1] && length(gamma.start)!=nq) {
    gamma.start <- rep(gamma.start[1],nq)
  }
  else if(length(gamma.start)!=nq) {
    cat("Error: gamma.start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(gamma.start)
}


form.gamma.start.mat <- function(gamma.start,nt,np){
  if(!all(!is.na(gamma.start))){
    cat("Error: gamma_start must no contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  else if(is.null(dim(gamma.start))){
    if(is.vector(gamma.start)){
      if(length(gamma.start)==(np*nt) || length(gamma.start)==nt || length(gamma.start)==1){
        gamma.start.mat <- matrix(gamma.start,nt,np) 
      }
      else if(length(gamma.start)==np){
        gamma.start.mat <- matrix(gamma.start,nt,np, byrow=TRUE) 
      }
      else if(sum(c(np,nt,np*nt,1)==length(gamma.start))==0){
        cat("Error: gamma_start not conformable, you can specify a vector of length np=", np,
            "(number of covariates plus intercept), nt=", nt,
            "(number of traits plus intercept) or np.nt=", np*nt, "to fill matrix gamma_start.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
    }
  }
  else if(sum(dim(gamma.start)==c(nt,np))==2){
    gamma.start.mat <- gamma.start
  }
  else if(sum(dim(gamma.start)==c(nt,np))!=2){
    cat("Error: gamma_start not conformable, should form a matrix of size (number of traits plus intercept) nt x np (number of covariates plus intercept) :"
        , nt," x", np,". \n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(gamma.start.mat)
}

form.beta.start.sp <- function (beta.start, np, nsp) {
  if (sum(dim(beta.start) != c(np,nsp))==0 && sum(is.na(beta.start))==0) { 
    beta.start.mat <- beta.start 
  }
  if (sum(is.na(beta.start))>0) { 
    beta.start.mat <- matrix(0, np, nsp)
  }
  else if(is.scalar(beta.start) && sum(is.na(beta.start))==0) {
    beta.start.mat <- matrix(beta.start, np, nsp) 
  }
  else if(sum(dim(beta.start) != c(np, nsp)) > 0) {
    stop("Error: beta.start not conformable. \n Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
  }
  return(beta.start.mat)
}

form.b.start <- function (b.start, nd) {
  if (sum(is.na(b.start))>0) { 
    b.start.vec <- rep(0,nd)
  }
  else if(is.scalar(b.start) && !is.na(b.start)) {
    b.start.vec <- rep(b.start, nd) 
  }
  else if(sum(length(b.start) != nd) > 0) {
    stop("Error: b.start not conformable.\n")
  }
  return(b.start.vec)
}

form.lambda.start.sp <- function (lambda.start, n_latent, nsp) {
  if (sum(is.na(lambda.start))>0) { 
    lambda.start.mat <- matrix(0, n_latent, nsp)
    for (i in 1:n_latent) {
      lambda.start.mat[i, i] <- 1
    }
  }
  else if(is.scalar(lambda.start) && sum(is.na(lambda.start))==0) {
    lambda.start.mat <- matrix(lambda.start, n_latent, nsp)
    for (i in 1:n_latent) {
      if (lambda.start > 0) {
        lambda.start.mat[i, i] <- lambda.start
      } else {
        lambda.start.mat[i, i] <- 1
      }
      for (j in 1:n_latent) {
        if (i > j) {
          lambda.start.mat[i, j] <- 0
        }
      }
    }
  }
  else if(sum(dim(lambda.start) != c(n_latent, nsp)) > 0) {
    stop("Error: lambda.start not conformable.\n")
  }
  else if(sum(dim(lambda.start) != c(n_latent, nsp)) == 0 && sum(is.na(lambda.start))==0) {
    lambda.start.mat <- lambda.start
    for (i in 1:n_latent) {
      if (lambda.start.mat[i, i]<=0) {
        stop("Error: lambda must be positive on the diagonal.\n")
      }
      for (j in 1:n_latent) {
        if (i > j && lambda.start.mat[i, j] != 0) {
          stop("Error: lambda_start must be an upper triangular matrix, values should be constrained to zero on lower diagonal.\n")
        }
      }
    }
  }
  return(lambda.start.mat)
}

form.alpha.start.sp <- function (alpha.start, nsite) {
  if (sum(is.na(alpha.start))>0) { 
    alpha.start <- rep(0, nsite)
  }
  else if(is.scalar(alpha.start) && !is.na(alpha.start)) {
    alpha.start <- rep(alpha.start, nsite) 
  }
  else if(length(alpha.start) != nsite ) {
    stop("Error: alpha.start not conformable.\n")
  }
  return(alpha.start)
}

form.W.start.sp <- function (W.start, nsite, n_latent) {
  if (sum(dim(W.start) == c(nsite, n_latent))==2 && sum(is.na(W.start))==0) { 
    W.start.mat <- W.start
  }
  if (sum(is.na(W.start))>0) { 
    W.start.mat <- matrix(0, nsite, n_latent)
  }
  else if(is.scalar(W.start) && !is.na(W.start)) {
    W.start.mat <- matrix(W.start, nsite, n_latent) 
  }
  else if(sum(dim(W.start) != c(nsite, n_latent)) > 0) {
    stop("Error: W.start not conformable.\n")
  }
  return(W.start.mat)
}

#======================
# Check and form prior
#======================

check.mugamma.mat <- function(mugamma,nt,np){
  if(!all(!is.na(mugamma))){
    cat("Error: mu_gamma must no contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  else if(is.null(dim(mugamma))){
    if(is.vector(mugamma)){
      if(length(mugamma)==(np*nt) || length(mugamma)==nt || length(mugamma)==1){
        mugamma.mat <- matrix(mugamma,nt,np) 
      }
      else if(length(mugamma)==np){
        mugamma.mat <- matrix(mugamma,nt,np, byrow=TRUE) 
      }
      else if(sum(c(np,nt,np*nt,1)==length(mugamma))==0){
        cat("Error: mu_gamma not conformable, you can specify a vector of length np=", np,
            "(number of covariates plus intercept), nt=", nt,
            "(number of traits plus intercept) or np.nt=", np*nt, "to fill matrix mu_gamma.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
    }
  }
  else if(sum(dim(mugamma)==c(nt,np))==2){
    mugamma.mat <- mugamma
  }
  else if(sum(dim(mugamma)==c(nt,np))!=2){
    cat("Error: mu_gamma not conformable, should form a matrix of size (number of traits plus intercept) nt x np (number of covariates plus intercept) :"
        , nt," x", np,". \n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mugamma.mat)
}

check.Vgamma.mat <- function(Vgamma,nt,np){
  if(!all(!is.na(Vgamma))){
    cat("Error: V_gamma must no contain missing values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (!all(Vgamma>=0)) {
    cat("Error: V_gamma should be positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  else if (is.null(dim(Vgamma))) {
    if(is.vector(Vgamma)){
      if(length(Vgamma)==(np*nt) || length(Vgamma)==nt || length(Vgamma)==1){
        Vgamma.mat <- matrix(Vgamma,nt,np)
      }
      else if(length(Vgamma)==np){
        Vgamma.mat <- matrix(Vgamma,nt,np, byrow=TRUE)
      }
      else if(sum(c(np,nt,np*nt,1)==length(Vgamma))==0){
        cat("Error: V_gamma not conformable, you can specify a vector of length np=", np,
            "(number of covariates plus intercept), nt=", nt,
            "(number of traits plus intercept) or nt.np=", nt*np, "to fill matrix V_gamma.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)
      }
    }
  }
  else if(sum(dim(Vgamma)==c(nt,np))==2){
    Vgamma.mat <- Vgamma
  }
  else if(sum(dim(Vgamma)==c(nt,np))!=2){
    cat("Error: V_gamma not conformable, should form a matrix of size (number of traits plus intercept) nt x np (number of covariates plus intercept) :",
        nt,"x", np,". \n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vgamma.mat)
}

check.mubeta <- function(mubeta, np){
  if (is.null(dim(mubeta))){
    if(is.scalar(mubeta)){
      mubeta <- rep(mubeta,np) 
    }
    else if(is.vector(mubeta) && length(mubeta)==np){
      mubeta <- mubeta
    }
  }
  else if (length(mubeta)!=np) {
    cat("Error: mubeta not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mubeta)
}

check.Vbeta <- function(Vbeta, np) {
  if (!all(Vbeta>0) && sum(is.na(Vbeta))!=0) {
    cat("Error: Vbeta should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.", call.=FALSE)
  }
  if (is.scalar(Vbeta)) {
    Vbeta <- rep(Vbeta,np)
  }
  else if (length(Vbeta)!=np){
    cat("Error: Vbeta not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",call.=FALSE)
  }
  if(!is.null(dim(Vbeta)) && length(Vbeta)==np){
    Vbeta <- as.vector(Vbeta)
  }
  return(Vbeta)
}

check.Vbeta.mat <- function(Vbeta, np) {
  if (sum(dim(Vbeta)==c(np,np))==2){
    if (!all(diag(Vbeta)>0) && sum(is.na(Vbeta))!=0) {
      cat("Error: V_beta should be strictly positive on diagonal.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    Vbeta.mat <- Vbeta
  }
  else if (!all(Vbeta>0) && sum(is.na(Vbeta))!=0) {
    cat("Error: V_beta should be strictly positive on diagonal.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (is.null(dim(Vbeta))){
    if(is.scalar(Vbeta)){
      Vbeta.mat <- matrix(0,np,np)
      diag(Vbeta.mat) <- Vbeta 
    }
    else if(is.vector(Vbeta)){
      if(length(Vbeta)==np){
        Vbeta.mat <- matrix(0,np,np)
        diag(Vbeta.mat) <- Vbeta 
      } else{
        cat("Error: V_beta not conformable, you must specify a ", np,"-length vector to fill the diagonal of the square matrix", np,"x",np,".\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE) 
      }
    }
  }
  else if (sum(dim(Vbeta) != c(np, np)) > 0) {
    cat("Error: V_beta not conformable, should form a square matrix", np,"x",np,".\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vbeta.mat)
}

check.mub <- function(mub, nd) {
  if (is.null(dim(mub))) {
    mub <- rep(mub,nd) 
  }
  else if (length(mub)!=nd) {
    cat("Error: mu_b not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mub)
}

check.Vb.mat <- function(Vb, nd) {
  if (!all(Vb>0)) {
    cat("Error: V_b should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (is.null(dim(Vb))) {
    Vb <- diag(rep(Vb,nd))
  }
  else if (sum(dim(Vb) != c(nd, nd)) > 0) {
    cat("Error: V_b not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vb)
}

check.mulambda <- function(mulambda, n_latent) {
  if (is.null(dim(mulambda))) {
    if(is.scalar(mulambda)){
      mulambda <- rep(mulambda,n_latent)
    }
    else if(is.vector(mulambda) && length(mulambda)==n_latent){
      mulambda <- mulambda
    }
  }
  else if (length(mulambda)!=n_latent) {
    cat("Error: mulambda not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(mulambda)
}

check.Vlambda.mat <- function(Vlambda, n_latent) {
  
  if (sum(dim(Vlambda)==c(n_latent,n_latent))==2) {
    if (!all(diag(Vlambda)>0) && sum(is.na(Vlambda))!=0) {
      cat("Error: V_lambda should be strictly positive on diagonal.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)
    }
    Vlambda.mat <- Vlambda
  }
  else if (!all(Vlambda>0) && sum(is.na(Vlambda))!=0) {
    cat("Error: Vlambda should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (is.null(dim(Vlambda))) {
    if(is.scalar(Vlambda)){
      if(n_latent==1){
        Vlambda.mat <- as.matrix(Vlambda)
      }else {
        Vlambda.mat <- diag(rep(Vlambda,n_latent))
      }
    }
    else if(is.vector(Vlambda) && length(Vlambda)==n_latent){
      Vlambda.mat <- diag(Vlambda)
    }
  }
  else if (sum(dim(Vlambda) != c(n_latent, n_latent)) > 0) {
    cat("Error: Vlambda not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(Vlambda.mat)
}

check.Vlambda <- function(Vlambda, n_latent) {
  if (!all(Vlambda>0) && sum(is.na(Vlambda))!=0) {
    cat("Error: Vlambda should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.", call.=FALSE)
  }
  if (is.scalar(Vlambda)) {
    Vlambda <- rep(Vlambda,n_latent)
  }
  else if (length(Vlambda)!=n_latent){
    cat("Error: Vlambda not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",call.=FALSE)
  }
  if(!is.null(dim(Vlambda)) && length(Vlambda)==n_latent){
    Vlambda <- as.vector(Vlambda)
  }
  return(Vlambda)
}

check.Valpha <- function(V_alpha_start) {
  if (!(V_alpha_start>0)) {
    cat("Error: V_alpha_start should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  else if (!is.null(dim(V_alpha_start))) {
    cat("Error: V_alpha_start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(V_alpha_start)
}

check.V <- function(V_start) {
  if (!(V_start>0)) {
    cat("Error: The variance of residuals V should be strictly positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  else if (!is.null(dim(V_start))) {
    cat("Error: V_start not conformable.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  return(V_start)
}
# End
