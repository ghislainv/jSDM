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
Rcpp::List Rcpp_jSDM_binomial_probit_fixed_site_lv_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::mat& X,
    const arma::mat& beta_start,
    const arma::mat& lambda_start,
    arma::vec alpha_start,
    const arma::mat& V_beta,
    const arma::vec& mu_beta,
    const arma::mat& V_lambda,
    const arma::vec& mu_lambda,
    const arma::mat& V_W,
    const arma::mat& W_start,
    const double& V_alpha,
    const int &seed,
    const int &verbose){
  
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
  const int NP = X.n_cols;
  const int NSP = beta_start.n_cols;
  const int NOBS = Y.n_elem;
  const int NSITE = W_start.n_rows;
  const int NL = W_start.n_cols; 
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::Cube<double> lambda; lambda.zeros(NSAMP, NSP, NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  /* Latent variable */
  arma::vec probit_theta_latent; probit_theta_latent.zeros(NOBS);
  arma::vec theta_latent; theta_latent.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat beta_run = beta_start;
  arma::mat lambda_run = lambda_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  // W latent variables (nsite*nl)
  arma::mat W_run = W_start;
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit(theta_ij) = X_i*beta_j + alpha_i + W_i*lambda_j 
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
  // compute matrix to use at each iteration 
  arma::mat inv_Vbeta=inv(V_beta);
  arma::mat inv_VW=inv(V_W);
  
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
    
    
    /////////////////////////////////////////////
    // mat latent variable W: Gibbs algorithm //
    // Loop on sites
    for (int i=0; i<NSITE; i++) {
      int nobs_i = rowId_site[i].n_elem;
      arma::rowvec Zprim; Zprim.zeros(nobs_i);
      for (int n=0; n<nobs_i; n++) {
        Zprim(n)=arma::as_scalar(Z_run.row(rowId_site[i].at(n)) -
          X.row(rowId_site[i].at(n))*beta_run.col(Id_sp(rowId_site[i].at(n))) -
          alpha_run(Id_site(rowId_site[i].at(n))));
      }
      // small_v
      arma::vec small_v=lambda_run.cols(Id_sp(rowId_site[i]))*Zprim.t();
      // big_V
      arma::mat big_V = inv(inv_VW + lambda_run.cols(Id_sp(rowId_site[i]))*lambda_run.cols(Id_sp(rowId_site[i])).t());
      
      // Draw in the posterior distribution
      arma::vec W_i = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      W_run.row(i) = W_i.t();
      
      ///////////////////////////////
      // vec alpha : Gibbs algorithm //
      if(i==0){
        // constraints of identifiability on alpha
        alpha_run(i) = 0.0;
      } else {
        double small_v2=0.0;
        // small_v
        for (int n=0; n<nobs_i; n++) {
          small_v2 += arma::as_scalar(Z_run.row(rowId_site[i].at(n))
                                        -X.row(rowId_site[i].at(n))*beta_run.col(Id_sp(rowId_site[i].at(n)))
                                        -W_run.row(i)*lambda_run.col(Id_sp(rowId_site[i].at(n))));
        }
        // big_V
        double big_V2 = 1.0/(1.0/V_alpha + nobs_i);
        
        // Draw in the posterior distribution
        alpha_run(i) = big_V2*small_v2 + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V2));
      }
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on sites
    
    // Centering and reducing W_i
    for ( int q = 0; q < NL; q++ ) {
      W_run.col(q) = W_run.col(q) - arma::mean(W_run.col(q));
      W_run.col(q) = W_run.col(q)*1.0/arma::stddev(W_run.col(q));
    }
    
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // latent variables W in long format 
      arma::mat W_run_long=W_run.rows(Id_site(rowId_sp[j]));
      //////////////////////////////////
      // mat beta: Gibbs algorithm //
      // small_v
      arma::vec small_v = inv_Vbeta*mu_beta + X.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) -
        W_run_long*lambda_run.col(j) -
        alpha_run(Id_site(rowId_sp[j])));
      
      // big_V
      arma::mat big_V = inv(inv_Vbeta + X.rows(rowId_sp[j]).t()*X.rows(rowId_sp[j]));
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));  
      
      //////////////////////////////////
      // mat lambda : Gibbs algorithm //
      for (int l=0; l<NL; l++) {
        if (l > j) {
          lambda_run(l,j) = 0.0;
        } else {
          arma::vec lambda_prop = lambda_run.col(j);
          lambda_prop(l) = 0.0; 
          // small_v
          double small_v = arma::as_scalar(1.0/(V_lambda(l,l))*mu_lambda(l) +
                                           W_run_long.col(l).t()*(Z_run(rowId_sp[j]) -
                                           X.rows(rowId_sp[j])*beta_run.col(j) -
                                           W_run_long*lambda_prop -
                                           alpha_run(Id_site(rowId_sp[j]))));
          // big_V
          double big_V = arma::as_scalar(1.0/(1.0/V_lambda(l,l)+
                                         W_run_long.col(l).t()*W_run_long.col(l)));
          if (l!=j){
            // Draw in the posterior distribution
            lambda_run(l,j) = big_V*small_v + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V));
          } else {
            lambda_run(l,j) = rtnorm(s,0,R_PosInf, small_v*big_V, std::sqrt(big_V));
          }
        }
      }
      R_CheckUserInterrupt(); // allow user interrupt
    } 
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_j + W_i*lambda_j 
      probit_theta_run(n) = arma::as_scalar(X.row(n)*beta_run.col(Id_sp(n))+
        W_run.row(Id_site(n))*lambda_run.col(Id_sp(n))+
        alpha_run(Id_site(n)));
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
        lambda.tube(isamp-1,j) = lambda_run.col(j);
      }
      for ( int i=0; i<NSITE; i++ ) {
        W.tube(isamp-1,i) = W_run.row(i);
      }
      alpha.row(isamp-1) = alpha_run.t();
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
                                          Rcpp::Named("lambda") = lambda,
                                          Rcpp::Named("W") = W,
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("theta_latent") = theta_latent,
                                          Rcpp::Named("probit_theta_latent") = probit_theta_latent );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_fixed_site_lv_long_format

// Test
/*** R
# # ===================================================
# # Data
# # ===================================================
# nsp <- 100
# nl <- 2
# nsite <- 300
# seed <- 1234
# set.seed(seed)
# 
# # Ecological process (suitability)
# ## X
# x1 <- rnorm(nsite,0,1)
# x1.2 <- scale(x1^2)
# X <- cbind(rep(1,nsite),x1,x1.2)
# colnames(X) <- c("Int","x1","x1.2")
# np <- ncol(X)
# ## W
# W <- matrix(rnorm(nsite*nl,0,1),nrow=nsite,byrow=TRUE)
# ## parameters
# beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
# mat <- t(matrix(runif(nsp*nl,-2,2), byrow=TRUE, nrow=nsp))
# diag(mat) <- runif(nl,0,2)
# lambda.target <- matrix(0,nl,nsp)
# lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# alpha.target <- runif(nsite,-2,2)
# alpha.target[1] <- 0
# ## probit_theta
# probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target)+ rep(alpha.target, nsp)
# nobs <- length(probit_theta)
# hist(probit_theta)
# e <- rnorm(nobs,0,1)
# 
# Z_true <- probit_theta + e
# Y<-rep(0,nobs)
# for (n in 1:nobs){
#   if ( Z_true[n] > 0) {Y[n] <- 1}
# }
# Id_site <- rep(1:nsite,nsp)-1
# Id_sp <- rep(1:nsp,each=nsite)-1
# data <- data.frame(site=Id_site, species=Id_sp, Y=Y,
#                    intercept=rep(1,nobs), x1=x1,
#                    x1.2=x1.2)
# lambda_start=matrix(0,nl,nsp)
# for (i in 1:nl){
#   lambda_start[i,i] = 1
# }
# X=as.matrix(data[,c("intercept","x1","x1.2")])
# # Call to C++ function
# # Iterations
# nsamp <- 10000
# nburn <- 10000
# nthin <- 10
# ngibbs <- nsamp+nburn
# T1 <- Sys.time()
# mod <- Rcpp_jSDM_binomial_probit_fixed_site_lv_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, Id_site=data$site, Id_sp=data$species,
#   beta_start=matrix(0,np,nsp),
#   V_beta=diag(c(10,rep(1,np-1))),
#   mu_beta = rep(0,np),
#   lambda_start=lambda_start,
#   V_lambda=diag(rep(1,nl)),
#   mu_lambda = rep(0,nl),
#   W_start=matrix(0,nsite,nl), V_W=diag(rep(1,nl)),
#   alpha_start=rep(0,nsite), V_alpha=10,
#   seed=1234, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# # Parameter estimates
# par(mfrow=c(1,2))
# hist(probit_theta)
# ## theta
# #plot(pnorm(probit_theta[-1]),mod$theta_latent)
# plot(pnorm(probit_theta),mod$theta_latent)
# abline(a=0,b=1,col='red')
# ## probit_theta
# par(mfrow=c(1,2))
# #plot(probit_theta[-1],mod$probit_theta_latent, xlab="obs", ylab ="fitted", main="probit(theta)")
# plot(probit_theta,mod$probit_theta_latent, xlab="obs", ylab ="fitted", main="probit(theta)")
# abline(a=0,b=1,col='red')
# ## Z
# #plot(Z_true[-1],mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# plot(Z_true,mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# abline(a=0,b=1,col='red')
# ## alpha
# par(mfrow=c(1,1))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="alpha.estimated")
# abline(a=0,b=1,col='red')
# ## beta_j
# par(mfrow=c(np,2))
# mean_beta <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta[j,] <-apply(mod$beta[,j,],2,mean)
#   if(j<5){
#     MCMC.betaj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#     for (p in 1:np) {
#       summary(MCMC.betaj)
#       coda::traceplot(MCMC.betaj[,p])
#       coda::densplot(MCMC.betaj[,p], main = paste0("beta_",j,p))
#       abline(v=beta.target[p,j],col='red')
#     }
#   }
# }
# 
# ## lambda_j
# par(mfrow=c(nl,2))
# mean_lambda <- matrix(0,nsp,nl)
# for (j in 1:nsp) {
#   mean_lambda[j,] <-apply(mod$lambda[,j,],2,mean)
#   MCMC.lambdaj <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#   if(j<5){
#     for (l in 1:nl) {
#       summary(MCMC.lambdaj)
#       coda::traceplot(MCMC.lambdaj[,l])
#       coda::densplot(MCMC.lambdaj[,l], main = paste0("lambda",j,l))
#       abline(v=lambda.target[l,j],col='red')
#     }
#   }
# }
# 
# ## Fixed species effect beta and loading factors lambda
# par(mfrow=c(1,2))
# ## beta
# plot(t(beta.target), mean_beta, xlab="obs", ylab="fitted", main="Fixed species effect beta")
# abline(a=0,b=1,col='red')
# ## lambda
# plot(t(lambda.target), mean_lambda, xlab="obs", ylab="fitted", main="Loading factors lambda")
# abline(a=0,b=1,col='red')
# 
# # W latent variables
# par(mfrow=c(1,2))
# mean_W <- apply(mod$W, c(2,3), mean)
# plot(W[,1],mean_W[,1])
# abline(a=0,b=1,col='red')
# plot(W[,2],mean_W[,2])
# abline(a=0,b=1,col='red')
# 
# # lambda * W
# par(mfrow=c(1,1))
# plot(W %*% lambda.target, mean_W %*%t(mean_lambda))
# abline(a=0,b=1,col='red')
# 
# # Deviance
# colnames(mod$Deviance) <- "Deviance"
# mcmc.Deviance <- coda::mcmc(mod$Deviance, start=nburn+1, end=ngibbs, thin=nthin)
# plot(mcmc.Deviance)
*/
