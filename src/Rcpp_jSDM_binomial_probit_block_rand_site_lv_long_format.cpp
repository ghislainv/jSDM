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
Rcpp::List Rcpp_jSDM_binomial_probit_block_rand_site_lv_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::mat& X,
    const arma::mat& param_start,
    const arma::mat& V_param,
    const arma::vec& mu_param,
    const arma::mat& V_W,
    const arma::mat& W_start,
    const arma::vec& alpha_start,
    double V_alpha_start,
    double shape,
    double rate,
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
  const int NSITE = W_start.n_rows;
  const int NP = X.n_cols;
  const int NSP = param_start.n_cols;
  const int NL = W_start.n_cols; 
  const int NOBS = Y.n_elem;
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> param; param.zeros(NSAMP, NSP, NP+NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  arma::vec V_alpha; V_alpha.zeros(NSAMP);
  /* Latent variable */
  arma::vec probit_theta_pred; probit_theta_pred.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat param_run = param_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  double V_alpha_run = V_alpha_start;
  // w latent variables (nsite*nl)
  arma::mat W_run = W_start;
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit_theta_ij = X_i*beta_j + W_i*lambda_j + alpha_i
  arma::mat probit_theta_run; probit_theta_run.zeros(NOBS);
  // Data 
  // data-set with covariables and latent variables 
  arma::mat data; data.zeros(NOBS,NP+NL);
  data = arma::join_rows(X,W_run.rows(Id_site));
  // Small fixed vectors indexed on i (site) or j (species) for data access
  arma::field<arma::uvec> rowId_site(NSITE); 
  arma::field<arma::uvec> rowId_sp(NSP); 
  
  for (int i=0; i<NSITE; i++) {
    rowId_site[i] = arma::find(Id_site==i);
  }
  for (int j=0; j<NSP; j++) {
    rowId_sp[j] = arma::find(Id_sp==j);
  }
  
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
    }
    
    arma::mat beta_run; beta_run = param_run.submat(0,0,NP-1,NSP-1);
    arma::mat lambda_run; lambda_run = param_run.submat(NP,0,NP+NL-1,NSP-1);
    // Loop on sites
    for (int i=0; i<NSITE; i++) {
      /////////////////////////////////////////////
      // mat latent variable W: Gibbs algorithm //
      int nobs_i = rowId_site[i].n_elem;
      arma::rowvec Zprim; Zprim.zeros(nobs_i);
      for (int n=0; n<nobs_i; n++) {
        Zprim(n)=arma::as_scalar(Z_run.row(rowId_site[i].at(n))
                                   -X.row(rowId_site[i].at(n))*beta_run.col(Id_sp(rowId_site[i].at(n)))
                                   -alpha_run(Id_site(rowId_site[i].at(n))));
      }
      // small_v
      arma::vec small_v=lambda_run.cols(Id_sp(rowId_site[i]))*Zprim.t();
      // big_V
      arma::mat big_V = inv(inv(V_W)+lambda_run.cols(Id_sp(rowId_site[i]))*lambda_run.cols(Id_sp(rowId_site[i])).t());
      
      // Draw in the posterior distribution
      arma::vec W_i = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      W_run.row(i) = W_i.t();
      data = arma::join_rows(X,W_run.rows(Id_site));
      
      ///////////////////////////////
      // vec alpha : Gibbs algorithm //
      double small_v2=0.0;
      // small_v
      for (int n=0; n<nobs_i; n++) {
        small_v2 += arma::as_scalar(Z_run.row(rowId_site[i].at(n))-data.row(rowId_site[i].at(n))*param_run.col(Id_sp(rowId_site[i].at(n))));
      }
      // big_V
      double big_V2 = 1/(1/V_alpha_run + nobs_i);
      
      // Draw in the posterior distribution
      alpha_run(i) = big_V2*small_v2 + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V2));
      }
    ////////////////////////////////////////////////
    // V_alpha
    double sum = arma::as_scalar(alpha_run.t()*alpha_run);
    double shape_posterior = shape + 0.5*NSITE;
    double rate_posterior = rate + 0.5*sum;
    
    V_alpha_run = rate_posterior/gsl_ran_gamma_mt(s, shape_posterior, 1.0);
    
    //////////////////////////////////
    // mat param: Gibbs algorithm //
    
    // Loop on species
    for (int j=0; j<NSP; j++) {
      
      // small_v
      arma::vec small_v = inv(V_param)*mu_param + data.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) - alpha_run(Id_site(rowId_sp[j])));
      
      // big_V
      arma::mat big_V = inv(inv(V_param) + data.rows(rowId_sp[j]).t()*data.rows(rowId_sp[j]));
      
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
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_j + W_i*lambda_j + alpha_i 
      probit_theta_run(n) = arma::as_scalar(data.row(n)*param_run.col(Id_sp(n)) + alpha_run(Id_site(n)));
      // link function probit is the inverse of N(0,1) repartition function 
      double theta = gsl_cdf_ugaussian_P(probit_theta_run(n));
      
      /* log Likelihood */
      logL += R::dbinom(Y(n), 1, theta, 1);
    } // loop on observations
    
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
        }
      }
      for ( int n=0; n<NOBS; n++ ) {
        Z_latent(n) += Z_run(n) / NSAMP; // We compute the mean of NSAMP values
        probit_theta_pred(n) += probit_theta_run(n)/NSAMP;        
      }
      alpha.row(isamp-1) = alpha_run.t();
      V_alpha(isamp-1) = V_alpha_run;
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
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("V_alpha") = V_alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("probit_theta_pred") = probit_theta_pred );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_block_rand_site_lv_long_format

// Test
/*** R
# # ===================================================
# # Data
# # ===================================================
# 
# nsp <- 70
# nsite <- 210
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
# beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
# l.zero <- 0
# l.diag <- runif(2,0,2)
# l.other <- runif(nsp*2-3,-2,2)
# lambda.target <- t(matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=T, nrow=nsp))
# param.target <- rbind(beta.target,lambda.target)
# V_alpha.target <- 0.5
# alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
# probit_theta <- X %*% beta.target + W %*% lambda.target + alpha.target
# X_supObs <- cbind(rep(1,nsite),rnorm(nsite),rnorm(nsite))
# probit_theta_supObs <- X_supObs%*%beta.target + W%*%lambda.target + alpha.target
# probit_theta <- c(probit_theta, probit_theta_supObs)
# nobs <- length(probit_theta)
# e <- rnorm(nobs,0,1)
# Z_true <- probit_theta + e
# Y<-rep(0,nobs)
# for (n in 1:nobs){
#   if ( Z_true[n] > 0) {Y[n] <- 1}
# }
# 
# param_start=matrix(0,np+nl,nsp)
# for (i in 1:nl){
#   param_start[np+i,i] = 1
# }
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
# mod <- Rcpp_jSDM_binomial_probit_block_rand_site_lv_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, Id_site=data$site, Id_sp=data$species,
#   param_start=param_start, V_param=diag(c(rep(1.0E6,np),rep(10,nl))),
#   mu_param = rep(0,np+nl), W_start=matrix(0,nsite,nl), V_W=diag(rep(1,nl)),
#   alpha_start=rep(0,nsite), V_alpha_start=1, shape=0.5, rate=0.0005,
#   seed=123, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# 
# # Parameter estimates
# ## probit_theta
# par(mfrow=c(1,1))
# plot(probit_theta[-1],mod$probit_theta_pred)
# abline(a=0,b=1,col='red')
# ## Z
# plot(Z_true[-1],mod$Z_latent)
# abline(a=0,b=1,col='red')
# ## alpha
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="alpha.estimated")
# abline(a=0,b=1,col='red')
# ## V_alpha
# MCMC_V_alpha <- coda::mcmc(mod$V_alpha, start=nburn+1, end=ngibbs, thin=nthin)
# summary(MCMC_V_alpha)
# par(mfrow=c(1,2))
# coda::traceplot(MCMC_V_alpha)
# coda::densplot(MCMC_V_alpha, main ="V_alpha" )
# abline(v=V_alpha.target,col='red')
# 
# ## beta_j
# par(mfrow=c(np,2))
# mean_beta <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta[j,] <-apply(mod$param[,j,1:np],2,mean)
#   if(j<5){
#     for (p in 1:np) {
#       MCMC.betaj <- coda::mcmc(mod$param[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
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
#   mean_lambda[j,] <- apply(mod$param[,j,(np+1):(nl+np)],2,mean)
#   if(j<5){
#     for (l in 1:nl) {
#       MCMC.lambdaj <- coda::mcmc(mod$param[,j,(np+1):(nl+np)], start=nburn+1, end=ngibbs, thin=nthin)
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
# ## Deviance
# mean(mod$Deviance)
*/
