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
Rcpp::List Rcpp_jSDM_gaussian_traits_fixed_site(const int ngibbs,const int nthin, const int nburn, 
                                                const arma::mat &Y, 
                                                const arma::mat &X,
                                                const arma::mat &Tr,
                                                const arma::mat &gamma_start,
                                                const arma::mat &gamma_zeros,
                                                const arma::mat &beta_start,
                                                const arma::mat &V_beta,
                                                const arma::mat &mu_gamma,
                                                const arma::mat &V_gamma,
                                                const arma::vec &alpha_start,
                                                const double &V_alpha,
                                                const double &V_start,
                                                const double &shape_V,
                                                const double &rate_V,
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
  const int NT = Tr.n_cols;
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  arma::Cube<double> gamma; gamma.zeros(NSAMP,NT,NP);
  arma::vec V; V.zeros(NSAMP);
  /* Latent variable */
  arma::mat Y_pred; Y_pred.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  // mat of species effects parameters 
  arma::mat beta_run = beta_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  // constraints of identifiability on alpha
  alpha_run(0) = 0.0;
  // coefficients for species traits and environment interaction
  arma::mat gamma_run=gamma_start;
  arma::mat mu_beta_run = Tr*gamma_run;
  // Variance of residuals 
  double V_run = V_start; 
  // Residuals
  arma::mat e; e.zeros(NSITE,NSP);
  // Y_hat = X_i*beta_j 
  arma::mat Y_hat ; Y_hat.zeros(NSITE,NSP);
  // inverse of matrix fo conjugate priors formula
  arma::mat inv_Vbeta = inv(V_beta);
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Gibbs sampler 
  
  for (int g=0; g < NGIBBS; g++) {
    
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
      arma::vec small_v = inv_Vbeta*mu_beta_run.row(j).as_col() + X.t()*(Y.col(j) - alpha_run)/V_run;
      // big_V
      arma::mat big_V = inv(inv_Vbeta + X.t()*X/V_run);
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      //////////////////////////////////////////////////
      // User interrupt
      R_CheckUserInterrupt(); // allow user interrupts 	 
    }
    
    // Loop on sites
    for (int i=1; i<NSITE; i++) {
      ///////////////////////////////
      // vec alpha : Gibbs algorithm //
      
      // small_v
      double small_v2 = arma::sum(Y.row(i)-X.row(i)*beta_run)/V_run;
      
      // big_V
      double big_V2 = 1.0/(1.0/V_alpha + NSP*1.0/V_run);
      
      // Draw in the posterior distribution
      alpha_run(i) = big_V2*small_v2 + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V2));
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    ////////////////////////////////////////////
    // Variance of residuals V : Gibbs algorithm 
    e = Y - X*beta_run; // Y-Y_hat
    e.each_col() -= alpha_run;
    double sum_V = arma::accu(e % e); 
    // Parameters
    double shape_posterior_V = shape_V + 0.5*(NSITE*NSP); //shape
    double rate_posterior_V = rate_V + 0.5*sum_V; //rate
    V_run = rate_posterior_V/gsl_ran_gamma_mt(s, shape_posterior_V, 1.0);
    
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for (int i = 0; i < NSITE; i++ ) {
      for (int j = 0; j < NSP; j++ ) {
        // Y_hat_ij = X_i*beta_j 
        Y_hat(i,j) = arma::as_scalar(X.row(i)*beta_run.col(j) + alpha_run(i));
        /* log Likelihood */
        logL += R::dnorm(Y(i,j), Y_hat(i,j), sqrt(V_run), 1.0);
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
      }
      for(int t=0; t<NT; t++){
        gamma.tube(isamp-1,t) = gamma_run.row(t);
      }
      alpha.row(isamp-1) = alpha_run.t();
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
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("gamma") = gamma,
                                          Rcpp::Named("V") = V,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Y_pred") = Y_pred);  
  return results;
  
} // end Rcpp_jSDM_gaussian_traits_fixed_site

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
# gamma.target <- matrix(runif(nt*np,-1,1), byrow=TRUE, nrow=nt)
# mu_beta <- as.matrix(Tr) %*% gamma.target
# V_beta <- diag(1,np)
# beta.target <- matrix(NA,nrow=np,ncol=nsp)
# library(MASS)
# for(j in 1:nsp){
#   beta.target[,j] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
# }
# alpha.target <- runif(nsite,-1,1)
# alpha.target[1] <- 0
# mu.target <- X %*% beta.target + alpha.target
# # Variance of residuals
# V.target <- 0.2
# Y <- matrix(rnorm(nsite*nsp, mu.target, sqrt(V.target)), nrow=nsite)
# hist(Y)
# # Call to C++ function
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# mod <- Rcpp_jSDM_gaussian_traits_fixed_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                             Y=Y, X=X, Tr=as.matrix(Tr),
#                                             gamma_start=matrix(0,nt,np),
#                                             gamma_zeros=matrix(1,nt,np),
#                                             V_gamma=matrix(1,nt,np), mu_gamma = matrix(0,nt,np),
#                                             beta_start=matrix(0,np,nsp),
#                                             V_beta=diag(rep(1,np)),
#                                             alpha_start=rep(0,nsite), V_alpha=1,
#                                             V_start=1 , shape_V=0.001, rate_V=0.001,
#                                             seed=123, verbose=1)
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
# ## alpha
# par(mfrow=c(1,1))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"],
#      ylab ="alpha.estimated", main="Fixed site effect")
# abline(a=0,b=1,col='red')
# 
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
# ## Variance of residuals
# par(mfrow=c(1,2))
# MCMC.V <- coda::mcmc(mod$V, start=nburn+1, end=ngibbs, thin=nthin)
# coda::traceplot(MCMC.V)
# coda::densplot(MCMC.V, main="Variance of residuals")
# abline(v=V.target, col='red')
# 
# ## Species effect beta
# par(mfrow=c(1,1),oma=c(1, 0, 1, 0))
# plot(t(beta.target), mean_beta, xlab="obs", ylab="fitted", main="beta")
# abline(a=0,b=1,col='red')
# title(main="Fixed species effects", outer=TRUE)
# 
# ## Deviance
# mean(mod$Deviance)
# ## Prediction
# # theta
# plot(Y, mod$Y_pred, xlab="obs", ylab="fitted", main="Y")
# abline(a=0,b=1,col='red')
*/