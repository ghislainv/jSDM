// ==============================================================================
// author          :Ghislain Vieilledent, Jeanne Clement
// email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <cmath>
#include "Rcpp_jSDM_useful.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

using namespace arma;
using namespace std;

/* ************************************************************ */
/* Gibbs sampler function */

// [[Rcpp::export]]
Rcpp::List Rcpp_jSDM_binomial_probit_block_lv(const int ngibbs,const int nthin,const int nburn, 
                                              const arma::umat &Y, 
                                              const arma::mat &X,
                                              const arma::mat &param_start,
                                              const arma::mat &V_param,
                                              const arma::vec &mu_param,
                                              arma::mat V_W,
                                              arma::mat W_start,
                                              const int seed,
                                              const int verbose) {
  
  ////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Defining and initializing objects
  
  ////////////////////////////////////////
  // Initialize random number generator //
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(s, seed);
  
  ///////////////////////////
  // Redefining constants //
  const int NGIBBS = ngibbs;
  const int NTHIN = nthin;
  const int NBURN = nburn;
  const int NSAMP = (NGIBBS-NBURN)/NTHIN;
  const int NSITE = Y.n_rows;
  const int NP = X.n_cols;
  const int NSP = Y.n_cols;
  const int NL = W_start.n_cols; 
  
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> param; param.zeros(NSAMP, NSP, NP+NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  /* Latent variable */
  arma::mat probit_theta_pred; probit_theta_pred.zeros(NSITE, NSP);
  arma::mat Z_latent; Z_latent.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat param_run = param_start;
  // w latent variables (nsite*nl)
  arma::mat W_run = W_start;
  // Z latent (nsite*nsp)
  arma::mat Z_run; Z_run.zeros(NSITE,NSP);
  // probit_theta_ij = X_i*beta_j + W_i*lambda_j
  arma::mat probit_theta_run; probit_theta_run.zeros(NSITE,NSP);
  // data 
  arma::mat data = arma::join_rows(X,W_run);
  
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Gibbs sampler 
  
  for (int g=0; g < NGIBBS; g++) {
    
    ////////////////////////////////////////////////
    // latent variable Z // 
    
    for (int j=0; j<NSP; j++) {
      for (int i=0; i<NSITE; i++) {
        // Actualization
        if (Y(i,j) == 0) {
          Z_run(i,j) = rtnorm(s, R_NegInf, 0, probit_theta_run(i,j), 1);
        } else {
          Z_run(i,j) = rtnorm(s, 0, R_PosInf, probit_theta_run(i,j), 1);
        }
      }
    }
    
    /////////////////////////////////////////////
    // mat latent variable W: Gibbs algorithm //
    
    // Loop on sites
    for (int i=0; i<NSITE; i++) {
      arma::mat beta_run = param_run.submat(0,0,NP-1,NSP-1);
      arma::mat lambda_run = param_run.submat(NP,0,NP+NL-1,NSP-1);
      // big_V
      arma::mat big_V = inv(inv(V_W)+lambda_run*lambda_run.t());
      
      // small_v
      arma::vec small_v =lambda_run*(Z_run.row(i)-X.row(i)*beta_run).t();
      
      // Draw in the posterior distribution
      arma::vec W_i = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      W_run.row(i) = W_i.t();
    }
    
    data = arma::join_rows(X, W_run);
    
    //////////////////////////////////
    // mat param: Gibbs algorithm //
    
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // small_v
      arma::vec small_v = inv(V_param)*mu_param + data.t()*(Z_run.col(j));
      // big_V
      arma::mat big_V = inv(inv(V_param)+data.t()*data);
      
      // Draw in the posterior distribution
      arma::vec param_prop = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      
      // constraints on lambda
      for (int l=0; l<NL; l++) {
        if (l > j) {
          param_prop(NP+l) = 0;
        }
        if ((l==j) & (param_prop(NP+l) < 0)) {
          param_prop(NP+l) = param_run(NP+l,j);
        }
      }
      param_run.col(j) = param_prop;
    }
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int i = 0; i < NSITE; i++ ) {
      for ( int j = 0; j < NSP; j++ ) {
        // probit(theta_ij) = X_i*beta_j + W_i*lambda_j 
        probit_theta_run(i,j) = arma::as_scalar(data.row(i)*param_run.col(j));
        // link function probit is the inverse of N(0,1) repartition function 
        double theta = gsl_cdf_ugaussian_P(probit_theta_run(i,j));
        
        /* log Likelihood */
        logL += R::dbinom(Y(i,j), 1, theta, 1);
      } // loop on species
    } // loop on sites
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        param.tube(isamp-1,j) = param_run.col(j);
        for ( int i=0; i<NSITE; i++ ) {
          W.tube(isamp-1,i) = W_run.row(i);
          Z_latent(i,j) += Z_run(i,j) / NSAMP; // We compute the mean of NSAMP values
          probit_theta_pred(i,j) += probit_theta_run(i,j)/NSAMP;        
        }
      }
      Deviance(isamp-1) = Deviance_run;
    }
    
    //////////////////////////////////////////////////
    // Progress bar
    double Perc=100*(g+1)/(NGIBBS);
    if (((g+1)%(NGIBBS/100))==0 && (verbose==1)) {  
      Rprintf("*");
      R_FlushConsole();
      //R_ProcessEvents(); for windows
      if (((g+1)%(NGIBBS/10))==0) {
        Rprintf(":%.1f%% \n",Perc);
        R_FlushConsole();
        //R_ProcessEvents(); for windows
      }
    } 
    
    
    //////////////////////////////////////////////////
    // User interrupt
    R_CheckUserInterrupt(); // allow user interrupts 	    
    
  } // Gibbs sampler
  
  // Free memory
  gsl_rng_free(s);
  
  // Return results as a Rcpp::List
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("param") = param,
                                          Rcpp::Named("W") = W,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("probit_theta_pred") = probit_theta_pred );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_block_lv

// Test
/*** R
# #===================================================
# #Data
# #===================================================
# 
# nsp<- 50
# nsite <- 200
# np <- 3
# nl <- 2
# seed <- 123
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# colnames(X) <- c("Int","x1","x2")
# W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
# data = cbind (X,W)
# beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
# l.zero <- 0
# l.diag <- runif(2,0,2)
# l.other <- runif(nsp*2-3,-2,2)
# lambda.target <- t(matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=T, nrow=nsp))
# param.target <- rbind(beta.target,lambda.target)
# probit_theta <- X %*% beta.target + W %*% lambda.target
# e <- matrix(rnorm(nsp*nsite,0,1),nsite,nsp)
# Z_true <- probit_theta + e
# 
# Y <- matrix (NA, nsite,nsp)
# for (i in 1:nsite){
#   for (j in 1:nsp){
#     if ( Z_true[i,j] > 0) {Y[i,j] <- 1}
#     else {Y[i,j] <- 0}
#   }
# }
# 
# param_start=matrix(0,np+nl,nsp)
# for (i in 1:nl){
#   param_start[np+i,i] = 1
# }
# 
# # Call to C++ function
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# mod_block <- Rcpp_jSDM_binomial_probit_block_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                           Y=Y, X=X,
#                                           param_start=param_start, W_start=matrix(0,nsite,nl),
#                                           V_param=diag(c(rep(10,np),rep(1,nl))),
#                                           mu_param = rep(0,np+nl), V_W=diag(rep(1,nl)),
#                                           seed=123, verbose=1)
# 
# # ===================================================
# # Result analysis
# # ===================================================
# 
# # Parameter estimates
# ## beta_j
# par(mfrow=c(np,2))
# mean_beta <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta[j,] <-apply(mod_block$param[,j,1:np],2,mean)
#   if(j<5){
#     for (p in 1:np) {
#       MCMC.betaj <- coda::mcmc(mod_block$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.betaj)
#       coda::traceplot(MCMC.betaj[,p])
#       coda::densplot(MCMC.betaj[,p], main = paste0("beta",p,j))
#       abline(v=beta.target[p,j],col='red')
#     }
#   }
# }
# ## lambda_j
# par(mfrow=c(nl*2,2))
# mean_lambda <- matrix(0,nsp,nl)
# for (j in 1:nsp) {
#   mean_lambda[j,] <- apply(mod_block$param[,j,(np+1):(nl+np)],2,mean)
#   if(j<5){
#     for (l in 1:nl) {
#       MCMC.lambdaj <- coda::mcmc(mod_block$param[,j,(np+1):(nl+np)], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.lambdaj)
#       coda::traceplot(MCMC.lambdaj[,l])
#       coda::densplot(MCMC.lambdaj[,l],main = paste0("lambda",l,j))
#       abline(v=lambda.target[l,j],col='red')
#     }
#   }
# }
# 
# ## Species effect beta and loading factors lambda
# par(mfrow=c(1,2),oma=c(1, 0, 1, 0))
# plot(t(beta.target),mean_beta, xlab="obs", ylab="fitted",main="beta")
# title("Fixed species effects", outer = T)
# abline(a=0,b=1,col='red')
# plot(t(lambda.target),mean_lambda, xlab="obs", ylab="fitted",main="lambda")
# abline(a=0,b=1,col='red')
# 
# # W latent variables
# par(mfrow=c(1,2),oma=c(1, 0, 1, 0))
# mean_W <- apply(mod_block$W, c(2,3), mean)
# plot(W[,1],mean_W[,1], main="W1",xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# title("Variables latentes", outer = T)
# plot(W[,2],mean_W[,2], main="W2", xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# 
# # lambda * W
# par(mfrow=c(1,1))
# plot(W %*% lambda.target,mean_W %*%t(mean_lambda),
#      xlab="obs", ylab= "fitted", main="W_i.lambda_j")
# abline(a=0,b=1,col='red')
# 
# ## Deviance
# mean(mod_block$Deviance)
# ## Prediction
# # probit_theta
# plot(probit_theta,mod_block$probit_theta_pred,xlab="obs", ylab="fitted",main="probit(theta)")
# abline(a=0,b=1,col='red')
# # Z
# plot(Z_true,mod_block$Z_latent, xlab="obs", ylab="fitted",main="Z_latent" )
# abline(a=0,b=1,col='red')
*/