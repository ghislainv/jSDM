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
Rcpp::List Rcpp_jSDM_binomial_probit_block_traits_rand_site_lv_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::vec& Id_common_var,
    const arma::mat& X,
    const arma::mat& D,
    const arma::mat& param_sp_start,
    arma::vec alpha_start,
    double V_alpha_start,
    const arma::mat& V_param_sp,
    const arma::vec& mu_param_sp,
    const arma::vec& beta_start,
    const arma::mat& V_beta,
    const arma::vec& mu_beta,
    const arma::mat& V_W,
    const arma::mat& W_start,
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
  const int NP = X.n_cols;
  const int ND = D.n_cols;
  const int NSP = param_sp_start.n_cols;
  const int NOBS = Y.n_elem;
  const int NSITE = W_start.n_rows;
  const int NL = W_start.n_cols; 
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> param_sp; param_sp.zeros(NSAMP, NSP, NP+NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  arma::vec V_alpha; V_alpha.zeros(NSAMP);
  arma::mat beta; beta.zeros(NSAMP, ND);
  /* Latent variable */
  arma::vec probit_theta_pred; probit_theta_pred.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters and coefficients for latent variables (nl+np,nsp)
  arma::mat param_sp_run = param_sp_start;
  arma::mat beta_run = beta_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  double V_alpha_run = V_alpha_start;
  // W latent variables (nsite*nl)
  arma::mat W_run = W_start;
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit_theta_ij = X_i*beta_j + alpha_i + W_i*lambda_j + D_ij*beta
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
    
    
    /////////////////////////////////////////////
    // mat latent variable W: Gibbs algorithm //
    // Loop on sites
    arma::mat beta_sp_run; beta_sp_run = param_sp_run.submat(0,0,NP-1,NSP-1);
    arma::mat lambda_sp_run; lambda_sp_run = param_sp_run.submat(NP,0,NP+NL-1,NSP-1);
    
    for (int i=0; i<NSITE; i++) {
      int nobs_i = rowId_site[i].n_elem;
      arma::rowvec Zprim; Zprim.zeros(nobs_i);
      for (int n=0; n<nobs_i; n++) {
        Zprim(n)=arma::as_scalar(Z_run.row(rowId_site[i].at(n))
                                   -X.row(rowId_site[i].at(n))*beta_sp_run.col(Id_sp(rowId_site[i].at(n)))
                                   -D.row(rowId_site[i].at(n))*beta_run
                                   -alpha_run(Id_site(rowId_site[i].at(n))));
      }
      // small_v
      arma::vec small_v=lambda_sp_run.cols(Id_sp(rowId_site[i]))*Zprim.t();
      // big_V
      arma::mat big_V = inv(inv(V_W)+lambda_sp_run.cols(Id_sp(rowId_site[i]))*lambda_sp_run.cols(Id_sp(rowId_site[i])).t());
      
      // Draw in the posterior distribution
      arma::vec W_i = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      W_run.row(i) = W_i.t();
      data = arma::join_rows(X,W_run.rows(Id_site));
    
    ///////////////////////////////
    // vec alpha : Gibbs algorithm //
    double small_v2=0.0;
    // small_v
    for (int n=0; n<nobs_i; n++) {
      small_v2 += arma::as_scalar(Z_run.row(rowId_site[i].at(n))
                                    -data.row(rowId_site[i].at(n))*param_sp_run.col(Id_sp(rowId_site[i].at(n)))
                                    -D.row(rowId_site[i].at(n))*beta_run);
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
    // vec beta: Gibbs algorithm //
    
    // Loop on species
    arma::vec small_v; small_v.zeros(ND);
    for (int j=0; j<NSP; j++) {
      // small_v
      small_v += D.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) 
                                            - data.rows(rowId_sp[j])*param_sp_run.col(j)
                                            -alpha_run(Id_site(rowId_sp[j])));
    }
    // small_v
    small_v += inv(V_beta)*mu_beta;
    
    // big_V
    arma::mat big_V = inv(inv(V_beta) + D.t()*D);
    
    // Draw in the posterior distribution
    beta_run = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
    
    //////////////////////////////////
    // mat param_sp: Gibbs algorithm //
    
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // small_v
      arma::vec small_v = inv(V_param_sp)*mu_param_sp + data.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j])
                                                                                      -D.rows(rowId_sp[j])*beta_run
                                                                                      -alpha_run(Id_site(rowId_sp[j])));
      
      // big_V
      arma::mat big_V = inv(inv(V_param_sp) + data.rows(rowId_sp[j]).t()*data.rows(rowId_sp[j]));
      
      // Draw in the posterior distribution
      arma::vec param_sp_prop = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      if(j==0 && !Id_common_var.has_nan()){
        // constraint of identifiability on beta_sp
        // Id of covariables in common between D and X 
        param_sp_prop(conv_to<uvec>::from(Id_common_var)).fill(0.0);
      }
      // constraints on lambda
      for (int l=0; l<NL; l++) {
        if (l > j) {
          param_sp_prop(NP+l) = 0;
        }
        if ((l==j) & (param_sp_prop(NP+l) < 0)) {
          param_sp_prop(NP+l) = param_sp_run(NP+l,j);
        }
      }
      param_sp_run.col(j) = param_sp_prop;
    }
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_sp_j + W_i*lambda_sp_j + D_ij*beta
      probit_theta_run(n) = arma::as_scalar(data.row(n)*param_sp_run.col(Id_sp(n))
                                              + D.row(n)*beta_run 
                                              + alpha_run(Id_site(n)));
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
        param_sp.tube(isamp-1,j) = param_sp_run.col(j);
      }
      for ( int i=0; i<NSITE; i++ ) {
        W.tube(isamp-1,i) = W_run.row(i);
      }
      beta.row(isamp-1) = beta_run.t();
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
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("param_sp") = param_sp,
                                          Rcpp::Named("beta") = beta,
                                          Rcpp::Named("W") = W,
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("V_alpha") = V_alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("probit_theta_pred") = probit_theta_pred );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_block_traits_rand_site_lv_long_format

// Test
/*** R
# # ===================================================
# # Data
# # ===================================================
# nsp <- 100
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
# ## D
# SLA <- runif(nsp,-1,1)
# D <- data.frame(Int=1, x1=x1, x1.2=x1.2, x1.SLA= scale(c(x1 %*% t(SLA))))
# nd <- ncol(D)
# ## parameters
# beta_sp.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
# # constraint of identifiability
# beta_sp.target[,1] <- 0.0
# mat <- t(matrix(runif(nsp*nl,-2,2), byrow=TRUE, nrow=nsp))
# diag(mat) <- runif(nl,0,2)
# lambda_sp.target <- matrix(0,nl,nsp)
# lambda_sp.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
# beta.target <-runif(nd,-1,1)
# V_alpha.target <- 0.5
# alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
# ## probit_theta
# probit_theta <- c(X %*% beta_sp.target) + c(W %*% lambda_sp.target) + as.matrix(D) %*% beta.target + rep(alpha.target, nsp)
# # x1_supObs <- rnorm(nsite)
# # X_supObs <- cbind(rep(1,nsite),x1_supObs,scale(x1_supObs^2))
# # D_supObs <- data.frame(Int=1, x1=x1_supObs, x1.2=scale(x1_supObs^2), x1.SLA= scale(c(x1_supObs %*% t(SLA))))
# # probit_theta_supObs <- c(X_supObs%*%beta_sp.target) + as.matrix(D_supObs) %*% beta.target
# # probit_theta <- c(probit_theta, probit_theta_supObs)
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
# # data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
# #                    intercept=rep(1,nobs), x1=c(D$x1,D_supObs$x1),
# #                    x1.2=c(D$x1.2,D_supObs$x1.2),x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# data <- data.frame(site=Id_site, species=Id_sp, Y=Y,
#                    intercept=rep(1,nobs), x1=D$x1,
#                    x1.2=D$x1.2,x1.SLA=D$x1.SLA, SLA=rep(SLA,each=nsite),
#                    W1=W[,1], W2=W[,2])
# # missing observation
# #data <- data[-1,]
# # # Remove species with less than 5 presences
# # rowId_rare_sp <- NULL
# # for(j in 0:(nsp-1)){
# #   if(sum(data[data$species==j,"Y"]) < 5){
# #     print(j)
# #     rowId_rare_sp <- c(rowId_rare_sp, which(data$species==j))
# #   }
# # }
# # if(!is.null(rowId_rare_sp)){
# #   data <- data[-rowId_rare_sp, ]
# #   nobs <- nrow(data)
# #   nsp <- length(unique(data$species))
# #   nsite <- length(unique(data$site))
# #   probit_theta <- probit_theta[-rowId_rare_sp]
# #   Z_true <- Z_true[-rowId_rare_sp]
# # }
# param_sp_start=matrix(0,np+nl,nsp)
# for (i in 1:nl){
#   param_sp_start[np+i,i] = 1
# }
# X=as.matrix(data[,c("intercept","x1","x1.2")])
# D=as.matrix(data[,c("intercept","x1","x1.2","x1.SLA")])
# # Call to C++ function
# # Iterations
# nsamp <- 10000
# nburn <- 10000
# nthin <- 10
# ngibbs <- nsamp+nburn
# T1 <- Sys.time()
# mod <- Rcpp_jSDM_binomial_probit_block_traits_rand_site_lv_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, D=D, Id_site=data$site, Id_sp=data$species,
#   Id_common_var=c(0,1,2),
#   beta_start=rep(0,nd), V_beta=diag(rep(2,nd)),
#   mu_beta = rep(0,nd),
#   param_sp_start=param_sp_start,
#   V_param_sp=diag(rep(2,np+nl)),
#   mu_param_sp = rep(0,np+nl),
#   W_start=matrix(0,nsite,nl), V_W=diag(rep(1,nl)),
#   alpha_start=rep(0,nsite), V_alpha_start=1, shape=0.5, rate=0.0005,
#   seed=1234, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# pdf(file="~/Documents/jSDM.test/results_binomial_probit_block_traits_rand_site_lv_long_format.pdf")
# par(mfrow=c(1,1))
# hist(probit_theta)
# # Parameter estimates
# ## probit_theta
# par(mfrow=c(1,2))
# #plot(probit_theta[-1],mod$probit_theta_pred, xlab="obs", ylab ="fitted", main="probit(theta)")
# plot(probit_theta,mod$probit_theta_pred, xlab="obs", ylab ="fitted", main="probit(theta)")
# abline(a=0,b=1,col='red')
# ## Z
# #plot(Z_true[-1],mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# plot(Z_true,mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# abline(a=0,b=1,col='red')
# ## alpha
# par(mfrow=c(1,3))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="alpha.estimated")
# abline(a=0,b=1,col='red')
# ## V_alpha
# MCMC_V_alpha <- coda::mcmc(mod$V_alpha, start=nburn+1, end=ngibbs, thin=nthin)
# summary(MCMC_V_alpha)
# coda::traceplot(MCMC_V_alpha)
# coda::densplot(MCMC_V_alpha, main ="V_alpha" )
# abline(v=V_alpha.target,col='red')
# ## beta
# par(mfrow=c(2,2))
# MCMC.beta <- coda::mcmc(mod$beta, start=nburn+1, end=ngibbs, thin=nthin)
# for(d in 1:nd){
# coda::traceplot(MCMC.beta[,d])
# coda::densplot(MCMC.beta[,d], main = paste0("beta_",colnames(D)[d]))
# abline(v=beta.target[d],col='red')
# }
# 
# ## beta_sp_j
# par(mfrow=c(np,2))
# mean_beta_sp <- matrix(0,nsp,np)
# for (j in 1:nsp) {
#   mean_beta_sp[j,] <-apply(mod$param_sp[,j,1:np],2,mean)
#   if(j<5){
#     MCMC.beta_spj <- coda::mcmc(mod$param_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
#     for (p in 1:np) {
#       summary(MCMC.beta_spj)
#       coda::traceplot(MCMC.beta_spj[,p])
#       coda::densplot(MCMC.beta_spj[,p], main = paste0("beta_sp",j,p))
#       abline(v=beta_sp.target[p,j],col='red')
#     }
#   }
# }
# 
# ## lambda_sp_j
# par(mfrow=c(nl,2))
# mean_lambda_sp <- matrix(0,nsp,nl)
# for (j in 1:nsp) {
#   mean_lambda_sp[j,] <-apply(mod$param_sp[,j,(np+1):(np+nl)],2,mean)
#   MCMC.lambda_spj <- coda::mcmc(mod$param_sp[,j,(np+1):(np+nl)], start=nburn+1, end=ngibbs, thin=nthin)
#   if(j<5){
#     for (l in 1:nl) {
#       summary(MCMC.lambda_spj)
#       coda::traceplot(MCMC.lambda_spj[,l])
#       coda::densplot(MCMC.lambda_spj[,l], main = paste0("lambda_sp",j,l))
#       abline(v=lambda_sp.target[l,j],col='red')
#     }
#   }
# }
# 
# ## Fixed species effect beta_sp and loading factors lambda_sp
# par(mfrow=c(1,2))
# ## beta_sp
# plot(t(beta_sp.target), mean_beta_sp, xlab="obs", ylab="fitted", main="Fixed species effect beta")
# abline(a=0,b=1,col='red')
# ## lambda_sp
# plot(t(lambda_sp.target), mean_lambda_sp, xlab="obs", ylab="fitted", main="Loading factors lambda")
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
# plot(W %*% lambda_sp.target, mean_W %*%t(mean_lambda_sp))
# abline(a=0,b=1,col='red')
# 
# # Deviance
# colnames(mod$Deviance) <- "Deviance"
# mcmc.Deviance <- coda::mcmc(mod$Deviance, start=nburn+1, end=ngibbs, thin=nthin)
# plot(mcmc.Deviance)
# 
# par(mfrow=c(1,1))
# colnames(mod$beta) <- paste0("beta_",colnames(D))
# boxplot(mod$beta, main="beta fitted by jSDM")
# # create data for segments
# n <- ncol(mod$beta)
# # width of each boxplot is 0.8
# x0s <- 1:n - 0.4
# x1s <- 1:n + 0.4
# # these are the y-coordinates for the horizontal lines
# # that you need to set to the desired values.
# y0s <- beta.target
# # add segments
# segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red", lwd=2)
# legend("bottomright", legend=c("beta target"), col=c('red'), lty=1, lwd=2)
# dev.off()
# #test en SDM logit binomial
# #voir obs dans IC
*/
