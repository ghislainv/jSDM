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
Rcpp::List Rcpp_jSDM_binomial_probit_fixed_site_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::mat& X,
    const arma::mat& beta_start,
    const arma::mat& V_beta,
    const arma::vec& mu_beta,
    const arma::vec& alpha_start,
    double V_alpha,
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
  const int NSITE = alpha_start.n_elem;
  const int NP = X.n_cols;
  const int NSP = beta_start.n_cols;
  const int NOBS = Y.n_elem;
  
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  /* Latent variable */
  arma::vec theta_latent; theta_latent.zeros(NOBS);
  arma::vec probit_theta_latent; probit_theta_latent.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat beta_run = beta_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit(theta_ij) = X_i*beta_j + alpha_i
  arma::mat probit_theta_run; probit_theta_run.zeros(NOBS);
  arma::mat theta_run; theta_run.zeros(NOBS);
  // Small fixed vectors indexed on i (site) or j (species) for data access
  arma::field<arma::uvec> rowId_site(NSITE); 
  arma::field<arma::uvec> rowId_sp(NSP); 
  for (int i=0; i<NSITE; i++) {
    rowId_site[i] = arma::find(Id_site==i);
  }
  for (int j=0; j<NSP; j++) {
    rowId_sp[j] = arma::find(Id_sp==j);
  }
  // inv(Vbeta)
  arma::mat inv_Vbeta=inv(V_beta);
  
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
    
    for (int n=0; n<NOBS; n++) {
      // Actualization
      if (Y(n) == 0) {
        Z_run(n) = rtnorm(s, R_NegInf, 0, probit_theta_run(n), 1);
      } else {
        Z_run(n) = rtnorm(s, 0, R_PosInf, probit_theta_run(n), 1);
      }
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    
    ///////////////////////////////
    // vec alpha : Gibbs algorithm //
    
    // Loop on sites 
    for (int i=0; i<NSITE; i++) {
      if(i==0){
        // constraints of identifiability on alpha
        alpha_run(i) = 0.0;
      } else {
        int nobs_i = rowId_site[i].n_elem;
        double small_v=0.0;
        // small_v
        for (int n=0; n<nobs_i; n++) {
          small_v += arma::as_scalar(Z_run.row(rowId_site[i].at(n))-X.row(rowId_site[i].at(n))*beta_run.col(Id_sp(rowId_site[i].at(n))));
        }
        // big_V
        double big_V = 1.0/(1.0/V_alpha + nobs_i);
        
        // Draw in the posterior distribution
        alpha_run(i) = big_V*small_v + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V));
      }
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    //////////////////////////////////
    // mat beta: Gibbs algorithm //
    
    // Loop on species
    for (int j=0; j<NSP; j++) {
      
      // small_v
      arma::vec small_v = inv_Vbeta*mu_beta + X.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) - alpha_run(Id_site(rowId_sp[j])));
      
      // big_V
      arma::mat big_V = inv(inv_Vbeta + X.rows(rowId_sp[j]).t()*X.rows(rowId_sp[j]));
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_j + alpha_i 
      probit_theta_run(n) = arma::as_scalar(X.row(n)*beta_run.col(Id_sp(n)) + alpha_run(Id_site(n)));
      // link function probit is the inverse of N(0,1) repartition function 
      theta_run(n) = gsl_cdf_ugaussian_P(probit_theta_run(n));
      
      /* log Likelihood */
      logL += R::dbinom(Y(n), 1, theta_run(n), 1);
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on observations
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    //////////////////////////////////////////////////
    // Output
    if (((g+1)>NBURN) && (((g+1-NBURN)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        beta.tube(isamp-1,j) = beta_run.col(j);
      }
      // We compute the mean of NSAMP values
      Z_latent += Z_run/NSAMP; 
      probit_theta_latent += probit_theta_run/NSAMP;        
      theta_latent += theta_run/NSAMP;        
      alpha.row(isamp-1) = alpha_run.t();
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
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("probit_theta_latent") = probit_theta_latent,
                                          Rcpp::Named("theta_latent") = theta_latent);  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_fixed_site_long_format

// Test
/*** R
# # ===================================================
# # Data
# # ===================================================
# 
# nsp <- 70
# nsite <- 210
# seed <- 123
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# colnames(X) <- c("Int","x1","x2")
# np <- ncol(X)
# beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
# alpha.target <- runif(nsite,-2,2)
# alpha.target[1] <- 0
# probit_theta <- X %*% beta.target + alpha.target
# X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
# probit_theta_supObs <- X_supObs%*%beta.target + alpha.target
# probit_theta <- c(probit_theta, probit_theta_supObs)
# nobs <- length(probit_theta)
# e <- rnorm(nobs,0,1)
# Z_true <- probit_theta + e
# Y<-rep(0,nobs)
# for (n in 1:nobs){
#   if ( Z_true[n] > 0) {Y[n] <- 1}
# }
# 
# Id_site <- rep(1:nsite,nsp)-1
# Id_sp <- rep(1:nsp,each=nsite)-1
# data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
#                    intercept=rep(1,nobs), x1=c(rep(x1,nsp),rep(X_supObs[,2],nsp)),
#                    x2=c(rep(x2,nsp),rep(X_supObs[,3],nsp)))
# # missing observation
# data <- data[-1,]
# X=as.matrix(data[,c("intercept","x1","x2")])
# # Call to C++ function
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# T1 <- Sys.time()
# mod <- Rcpp_jSDM_binomial_probit_fixed_site_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, Id_site=data$site, Id_sp=data$species,
#   beta_start=matrix(0,np,nsp), V_beta=diag(rep(100,np)),
#   mu_beta = rep(0,np),
#   alpha_start=rep(0,nsite), V_alpha=10,
#   seed=123, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# 
# # Parameter estimates
# ## theta
# par(mfrow=c(1,1))
# plot(pnorm(probit_theta[-1]), mod$theta_latent, xlab="obs", ylab ="fitted", main="theta")
# abline(a=0,b=1,col='red')
# ## probit_theta
# par(mfrow=c(1,2))
# plot(probit_theta[-1],mod$probit_theta_latent, xlab="obs", ylab ="fitted", main="probit(theta)")
# abline(a=0,b=1,col='red')
# ## Z
# plot(Z_true[-1],mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# abline(a=0,b=1,col='red')
# 
# ## beta_j
# par(mfrow=c(np,2))
# mean_beta <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta[j,] <-apply(mod$beta[,j,1:np],2,mean)
#   if(j<5){
#     for (p in 1:np) {
#       MCMC.betaj <- coda::mcmc(mod$beta[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.betaj)
#       coda::traceplot(MCMC.betaj[,p])
#       coda::densplot(MCMC.betaj[,p], main = paste0("beta",p,j))
#       abline(v=beta.target[p,j],col='red')
#     }
#   }
# }
# 
# ## Fixed pecies effect beta and site effect alpha
# par(mfrow=c(1,2))
# ## beta
# plot(t(beta.target),mean_beta, xlab="obs", ylab="fitted", main="Fixed species effect beta")
# abline(a=0,b=1,col='red')
# ## alpha
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], xlab="obs", ylab ="fitted", main="Fixed site effect alpha")
# abline(a=0,b=1,col='red')
# 
# mean(mod$Deviance)
*/
