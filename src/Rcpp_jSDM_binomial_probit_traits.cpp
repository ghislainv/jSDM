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
Rcpp::List Rcpp_jSDM_binomial_probit_traits(const int ngibbs,const int nthin,const int nburn, 
                                            const arma::umat &Y, 
                                            const arma::mat &X,
                                            const arma::mat &Tr,
                                            const arma::mat &gamma_start,
                                            const arma::mat &gamma_zeros,
                                            const arma::mat &beta_start,
                                            const arma::mat &V_gamma,
                                            const arma::mat &mu_gamma,
                                            const arma::mat &V_beta,
                                            const int &seed,
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
  const int NT = Tr.n_cols;
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::Cube<double> gamma; gamma.zeros(NSAMP,NT,NP);
  /* Latent variable */
  arma::mat probit_theta_latent; probit_theta_latent.zeros(NSITE, NSP);
  arma::mat theta_latent; theta_latent.zeros(NSITE, NSP);
  arma::mat Z_latent; Z_latent.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat beta_run = beta_start;
  arma::mat gamma_run=gamma_start;
  arma::mat mu_beta_run = Tr*gamma_run;
  // Z latent (nsite*nsp)
  arma::mat Z_run; Z_run.zeros(NSITE,NSP);
  // probit_theta_ij = X_i*beta_j 
  arma::mat probit_theta_run; probit_theta_run.zeros(NSITE,NSP);
  arma::mat theta_run; theta_run.zeros(NSITE,NSP);
  // inverse of matrix fo conjugate priors formula
  arma::mat inv_Vbeta = inv(V_beta);
  arma::mat inv_VbetaXtX = inv(inv_Vbeta+X.t()*X);
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
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    /////////////////////////////////////////////
    //mat gamma: Gibbs algorithm //
    for (int t= 0; t < NT; t++ ) {
      for (int p= 0; p < NP; p++ ) {
        if(gamma_zeros(t,p)!=0.0){
          arma::vec gamma_prop=gamma_run.col(p);
          gamma_prop(t)=0.0;
          // small_v
          double small_v = arma::as_scalar(1.0/(V_gamma(t,p))*mu_gamma(t,p) +
                                           Tr.col(t).t()*(beta_run.row(p).t()-Tr*gamma_prop));
          // big_V
          double big_V = 1.0/arma::as_scalar((1.0/V_gamma(t,p) + Tr.col(t).t()*Tr.col(t)));
          // Draw in the posterior distribution
          gamma_run(t,p)= big_V*small_v + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V));
        }
      }
    }
    
    //////////////////////////////////
    // mat beta: Gibbs algorithm //
    mu_beta_run = Tr*gamma_run;
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // small_v
      arma::vec small_v = inv_Vbeta*mu_beta_run.row(j).as_col() + X.t()*(Z_run.col(j));
      // big_V
      arma::mat big_V = inv_VbetaXtX;
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for (int i = 0; i < NSITE; i++ ) {
      for (int j = 0; j < NSP; j++ ) {
        // probit(theta_ij) = X_i*beta_j
        probit_theta_run(i,j) = arma::as_scalar(X.row(i)*beta_run.col(j));
        // link function probit is the inverse of N(0,1) repartition function 
        theta_run(i,j) = gsl_cdf_ugaussian_P(probit_theta_run(i,j));
        
        /* log Likelihood */
        logL += R::dbinom(Y(i,j), 1, theta_run(i,j), 1);
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on species
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on sites
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1-NBURN)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        beta.tube(isamp-1,j) = beta_run.col(j);
      }
      for(int t=0; t<NT; t++){
        gamma.tube(isamp-1,t) = gamma_run.row(t);
      }
      Deviance(isamp-1) = Deviance_run;
      // We compute the mean of NSAMP values
      Z_latent += Z_run / NSAMP;
      probit_theta_latent += probit_theta_run/NSAMP;    
      theta_latent += theta_run/NSAMP;    
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
                                          Rcpp::Named("gamma") = gamma,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("theta_latent") = theta_latent,
                                          Rcpp::Named("probit_theta_latent") = probit_theta_latent);  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_traits

// Test
/*** R
# #===================================================
# #Data
# #===================================================
# 
# nsp <- 70
# nsite <- 210
# np <- 3
# seed <- 123
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# colnames(X) <- c("Int","x1","x2")
# Tr <- data.frame(Int=1, WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
# nt <- ncol(Tr)
# gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
# mu_beta <- as.matrix(Tr) %*% gamma.target
# V_beta <- diag(1,np)
# beta.target <- matrix(NA,nrow=np,ncol=nsp)
# library(MASS)
# for(j in 1:nsp){
#   beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
# }
# probit_theta <- X %*% beta.target 
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
# # Call to C++ function
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# mod <- Rcpp_jSDM_binomial_probit_traits(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                         Y=Y, X=X, Tr=as.matrix(Tr),
#                                         gamma_start=matrix(0,nt,np),
#                                         gamma_zeros=matrix(1,nt,np),
#                                         V_gamma=matrix(10,nt,np), mu_gamma = matrix(0,nt,np),
#                                         beta_start=matrix(0,np,nsp), V_beta=diag(rep(1,np)),
#                                         seed=1234, verbose=1)
# 
# # ===================================================
# # Result analysis
# # ===================================================
# 
# # Parameter estimates
# ## gamma
# par(mfrow=c(2,2))
# for(p in 1:np){
#   MCMC.gamma_p <- coda::mcmc(mod$gamma[,,p], start=nburn+1, end=ngibbs, thin=nthin)
#   for(t in 1:nt){
#     coda::traceplot(MCMC.gamma_p[,t])
#     coda::densplot(MCMC.gamma_p[,t], main = paste0("gamma_",colnames(X)[p],".",colnames(Tr)[t]))
#     abline(v=gamma.target[t,p],col='red')
#   }
# }
# 
# ## beta_j
# par(mfrow=c(np,2))
# mean_beta <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta[j,] <-apply(mod$beta[,j,],2,mean)
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
# ## Species effect beta 
# par(mfrow=c(1,1),oma=c(1, 0, 1, 0))
# plot(t(beta.target),mean_beta, xlab="obs", ylab="fitted",main="beta")
# abline(a=0,b=1,col='red')
# title(main="Fixed species effects", outer=TRUE)
# 
# ## Deviance
# mean(mod$Deviance)
# ## Prediction
# # theta
# plot(pnorm(probit_theta),mod$theta_latent,xlab="obs", ylab="fitted",main="theta")
# abline(a=0,b=1,col='red')
# # probit_theta
# plot(probit_theta,mod$probit_theta_latent,xlab="obs", ylab="fitted",main="probit(theta)")
# abline(a=0,b=1,col='red')
# # Z
# plot(Z_true,mod$Z_latent, xlab="obs", ylab="fitted",main="Z_latent" )
# abline(a=0,b=1,col='red')
*/