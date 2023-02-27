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
Rcpp::List Rcpp_jSDM_gaussian_lv(const int ngibbs,const int nthin, const int nburn, 
                                 const arma::mat &Y, 
                                 const arma::mat &X,
                                 const arma::mat &beta_start,
                                 const arma::mat &V_beta,
                                 const arma::vec &mu_beta,
                                 const arma::mat &lambda_start,
                                 const arma::mat &W_start,
                                 const arma::vec &mu_lambda,
                                 const arma::mat &V_lambda,
                                 const arma::mat &V_W,
                                 const double &V_start,
                                 const double &shape_V,
                                 const double & rate_V,
                                 const int &seed,
                                 const int &verbose) {
  
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
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::Cube<double> lambda; lambda.zeros(NSAMP, NSP, NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  arma::vec V; V.zeros(NSAMP);
  /* Latent variable */
  arma::mat Y_pred; Y_pred.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters (np*nsp)
  arma::mat beta_run = beta_start;
  // factor loadings coefficients for latent variables (nl*nsp)
  arma::mat lambda_run = lambda_start;
  // w latent variables (nsite*nl)
  arma::mat W_run = W_start;
  // Residuals standard deviation 
  double V_run = V_start; 
  // Residuals
  arma::mat e; e.zeros(NSITE,NSP);
  // Y_hat = X_i*beta_j 
  arma::mat Y_hat ; Y_hat.zeros(NSITE,NSP);
  // inverse of matrix fo conjugate priors formula
  arma::mat inv_Vbeta = inv(V_beta);
  arma::mat inv_VW = inv(V_W);
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Gibbs sampler 
  
  for (int g=0; g < NGIBBS; g++) {
    
    // Loop on sites
    for (int i=0; i<NSITE; i++) {
      /////////////////////////////////////////////
      // mat latent variable W: Gibbs algorithm //
      
      // big_V
      arma::mat big_V = inv(inv_VW + lambda_run*lambda_run.t()/V_run);
      
      // small_v
      arma::vec small_v = lambda_run*(Y.row(i)-X.row(i)*beta_run).t()/V_run;
      
      // Draw in the posterior distribution
      arma::vec W_i = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      W_run.row(i) = W_i.t();
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    // Centering and reducing W_i
    for ( int q = 0; q < NL; q++ ) {
      W_run.col(q) = W_run.col(q) - arma::mean(W_run.col(q));
      W_run.col(q) = W_run.col(q)*1.0/arma::stddev(W_run.col(q));
    }
    
    ////////////////////////////////////////////
    // Variance of residuals V : Gibbs algorithm 
    e = Y - X*beta_run - W_run*lambda_run; // Y-Y_hat
    double sum = arma::accu(e % e); 
    // Parameters
    double shape_posterior = shape_V + 0.5*(NSITE*NSP); //shape
    double rate_posterior = rate_V + 0.5*sum; //rate
    V_run = rate_posterior/gsl_ran_gamma_mt(s, shape_posterior, 1.0);
    
    
    // Loop on species
    for (int j=0; j<NSP; j++){
      
      //////////////////////////////////
      // mat beta: Gibbs algorithm //
      // small_v
      arma::vec small_v = inv_Vbeta*mu_beta + X.t()*(Y.col(j)-W_run*lambda_run.col(j))/V_run;
      // big_V
      arma::mat big_V = inv(inv_Vbeta + X.t()*X/V_run);
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      //////////////////////////////////////////////////
      // User interrupt
      R_CheckUserInterrupt(); // allow user interrupts 	 
      
      //////////////////////////////////
      // mat lambda : Gibbs algorithm //
      for (int l=0; l<NL; l++){
        if (l > j) {
          lambda_run(l,j) = 0;
        } else {
          arma::vec lambda_prop = lambda_run.col(j);
          lambda_prop(l) = 0.0; 
          // small_v
          double small_v = arma::as_scalar(1.0/(V_lambda(l,l))*mu_lambda(l) + W_run.col(l).t()*(Y.col(j)-X*beta_run.col(j)-W_run*lambda_prop)/V_run);
          // big_V
          double big_V = arma::as_scalar(1.0/(1.0/V_lambda(l,l)+ W_run.col(l).t()*W_run.col(l)/V_run));
          if (l!=j){
            // Draw in the posterior distribution
            lambda_run(l,j) = big_V*small_v + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V));
          } else {
            lambda_run(l,j) = rtnorm(s,0,R_PosInf, small_v*big_V, std::sqrt(big_V));
          }
        }
        R_CheckUserInterrupt(); // allow user interrupt
      }
    }
    
  
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for (int i = 0; i < NSITE; i++ ) {
      for (int j = 0; j < NSP; j++ ) {
        // Y_hat_ij = X_i*beta_j 
        Y_hat(i,j) = arma::as_scalar(X.row(i)*beta_run.col(j) + W_run.row(i)*lambda_run.col(j));
        /* log Likelihood */
        logL += R::dnorm(Y(i,j), Y_hat(i,j), sqrt(V_run), 1);
        R_CheckUserInterrupt(); // allow user interrupts 	 
      } // loop on species
      R_CheckUserInterrupt(); // allow user interrupts 	 
    } // loop on sites
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1-NBURN)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        beta.tube(isamp-1,j) = beta_run.col(j);
        lambda.tube(isamp-1,j) = lambda_run.col(j);
      }
      for (int i=0; i<NSITE; i++) {
        W.tube(isamp-1,i) = W_run.row(i);
      }
      Deviance(isamp-1) = Deviance_run;
      V(isamp-1) = V_run;    
      // We compute the mean of NSAMP values
      Y_pred += Y_hat/NSAMP;    
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
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                          Rcpp::Named("lambda") = lambda,
                                          Rcpp::Named("W") = W,
                                          Rcpp::Named("V") = V,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Y_pred") = Y_pred);  
  return results;
  
} // end Rcpp_jSDM_gaussian_lv

// Test
/*** R
# #===================================================
# #Data
# #===================================================
# 
# nsp <- 70
# nsite <- 210
# np <- 3
# nl <- 2
# seed <- 1234
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# colnames(X) <- c("Int","x1","x2")
# beta.target <- t(matrix(runif(nsp*np,-1, 1), byrow=TRUE, nrow=nsp))
# W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
# mat <- t(matrix(runif(nsp*nl,-1,1), byrow=TRUE, nrow=nsp))
# diag(mat) <- runif(nl,0,1)
# lambda.target <- matrix(0,nl,nsp)
# lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# V.target <- 0.2
# mu.target <- X %*% beta.target + W %*% lambda.target
# Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)
# hist(Y)
# 
# # Call to C++ function
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# lambda_start=matrix(0,nl,nsp)
# for (i in 1:nl){
#   lambda_start[i,i] = 1
# }
# mod <- Rcpp_jSDM_gaussian_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                              Y=Y, X=X,
#                              beta_start=matrix(0,np,nsp),
#                              V_beta=diag(rep(1,np)), mu_beta=rep(0,np),
#                              lambda_start=lambda_start,
#                              V_lambda=diag(rep(1,nl)), mu_lambda = rep(0,nl),
#                              W_start=matrix(0,nsite,nl), V_W=diag(rep(1,nl)),
#                              V_start=1 , shape_V=0.001, rate_V=0.001,
#                              seed=123, verbose=1)
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
#   mean_beta[j,] <-apply(mod$beta[,j,], 2, mean)
#   if(j<5){
#     for (p in 1:np) {
#       MCMC.betaj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.betaj)
#       coda::traceplot(MCMC.betaj[,p])
#       coda::densplot(MCMC.betaj[,p], main = paste0("beta",p,j))
#       abline(v=beta.target[p,j],col='red')
#     }
#   }
# }
# 
# ## lambda_j
# par(mfrow=c(nl*2,2))
# mean_lambda <- matrix(0,nsp,nl)
# for (j in 1:nsp) {
#   mean_lambda[j,] <- apply(mod$lambda[,j,],2,mean)
#   if(j<5){
#     for (l in 1:nl) {
#       MCMC.lambdaj <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.lambdaj)
#       coda::traceplot(MCMC.lambdaj[,l])
#       coda::densplot(MCMC.lambdaj[,l],main = paste0("lambda",l,j))
#       abline(v=lambda.target[l,j],col='red')
#     }
#   }
# }
# 
# ## Variance of residuals
# par(mfrow=c(1,2))
# MCMC.V <- coda::mcmc(mod$V, start=nburn+1, end=ngibbs, thin=nthin)
# coda::traceplot(MCMC.V)
# coda::densplot(MCMC.V, main="Variance of residuals")
# abline(v=V.target, col='red')
# 
# ## Species effect beta and factor loadings
# par(mfrow=c(1,1),oma=c(1, 0, 1, 0))
# plot(t(beta.target),mean_beta, xlab="obs", ylab="fitted",main="beta")
# abline(a=0,b=1,col='red')
# plot(t(lambda.target),mean_lambda, xlab="obs", ylab="fitted",main="lambda")
# abline(a=0,b=1,col='red')
# title(main="Fixed species effects", outer=TRUE)
# 
# # W latent variables
# par(mfrow=c(1,2),oma=c(1, 0, 1, 0))
# mean_W <- apply(mod$W, c(2,3), mean)
# plot(W[,1],mean_W[,1], main="W1",xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# title("Variables latentes", outer=TRUE)
# plot(W[,2],mean_W[,2], main="W2", xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# 
# ## Deviance
# mean(mod$Deviance)
# ## Prediction
# # theta
# plot(Y, mod$Y_pred, xlab="obs", ylab="fitted", main="Y")
# abline(a=0,b=1,col='red')
*/