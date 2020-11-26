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
Rcpp::List Rcpp_jSDM_binomial_probit_block_traits_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::vec& Id_common_var,
    const arma::mat& X,
    const arma::mat& D,
    const arma::mat& beta_sp_start,
    const arma::mat& V_beta_sp,
    const arma::vec& mu_beta_sp,
    const arma::vec& beta_start,
    const arma::mat& V_beta,
    const arma::vec& mu_beta,
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
  const int NSP = beta_sp_start.n_cols;
  const int NOBS = Y.n_elem;
  
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta_sp; beta_sp.zeros(NSAMP, NSP, NP);
  arma::mat beta; beta.zeros(NSAMP, ND);
  /* Latent variable */
  arma::vec probit_theta_pred; probit_theta_pred.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters (np,nsp)
  arma::mat beta_sp_run = beta_sp_start;
  arma::mat beta_run = beta_start;
  
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit_theta_ij = X_i*beta_sp_j + D_ij*beta
  arma::mat probit_theta_run; probit_theta_run.zeros(NOBS);
  // Small fixed vectors indexed on j (species) for data access
  arma::field<arma::uvec> rowId_sp(NSP); 
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
    
    //////////////////////////////////
    // vec beta: Gibbs algorithm //
    
    // Loop on species
    arma::vec small_v; small_v.zeros(ND);
    for (int j=0; j<NSP; j++) {
      // small_v
      small_v += D.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) - X.rows(rowId_sp[j])*beta_sp_run.col(j));
    }
    // small_v
    small_v += inv(V_beta)*mu_beta;
    
    // big_V
    arma::mat big_V = inv(inv(V_beta) + D.t()*D);
    
    // Draw in the posterior distribution
    beta_run = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
    
    //////////////////////////////////
    // mat beta_sp: Gibbs algorithm //
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // small_v
      arma::vec small_v = inv(V_beta_sp)*mu_beta_sp + X.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j])-D.rows(rowId_sp[j])*beta_run);
      
      // big_V
      arma::mat big_V = inv(inv(V_beta_sp) + X.rows(rowId_sp[j]).t()*X.rows(rowId_sp[j]));
      
      // Draw in the posterior distribution
      beta_sp_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      if(j==0 && !Id_common_var.has_nan()){
        // constraint of identifiability on beta_sp
        // Id of columns in common between D and X 
        arma::uvec Id_sp0; Id_sp0.zeros(1);
        beta_sp_run.submat(conv_to<uvec>::from(Id_common_var),Id_sp0).fill(0.0);
      }
    }
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_sp_j + D_ij*beta
      probit_theta_run(n) = arma::as_scalar(X.row(n)*beta_sp_run.col(Id_sp(n))+D.row(n)*beta_run);
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
        beta_sp.tube(isamp-1,j) = beta_sp_run.col(j);
      }
      beta.row(isamp-1) = beta_run.t();
      for ( int n=0; n<NOBS; n++ ) {
        Z_latent(n) += Z_run(n) / NSAMP; // We compute the mean of NSAMP values
        probit_theta_pred(n) += probit_theta_run(n)/NSAMP;        
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
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("beta_sp") = beta_sp,
                                          Rcpp::Named("beta") = beta,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("probit_theta_pred") = probit_theta_pred );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_block_traits_long_format

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
# x1 <- rnorm(nsite,0,1)
# x1.2 <- scale(x1^2)
# X <- cbind(rep(1,nsite),x1,x1.2)
# colnames(X) <- c("Int","x1","x1.2")
# np <- ncol(X)
# beta_sp.target <- t(matrix(runif(nsp*np,-1,1), byrow=TRUE, nrow=nsp))
# # constraint of identifiability
# beta_sp.target[,1] <- 0.0
# SLA <- runif(nsp,-1,1)
# D <- data.frame(Int=1, x1=x1, x1.2=x1.2, x1.SLA= scale(c(x1 %*% t(SLA))))
# nd <- ncol(D)
# beta.target <-runif(nd,-1,1)
# probit_theta <- c(X %*% beta_sp.target) + as.matrix(D) %*% beta.target
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
#                    x1.2=D$x1.2,x1.SLA=D$x1.SLA, SLA=rep(SLA,each=nsite))
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
# X=as.matrix(data[,c("intercept","x1","x1.2")])
# D=as.matrix(data[,c("intercept","x1","x1.2","x1.SLA")])
# # Call to C++ function
# # Iterations
# nsamp <- 10000
# nburn <- 10000
# nthin <- 10
# ngibbs <- nsamp+nburn
# T1 <- Sys.time()
# mod <- Rcpp_jSDM_binomial_probit_block_traits_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, D=D, Id_site=data$site, Id_sp=data$species,
#   Id_common_var=c(0,1,2),
#   beta_start=rep(0,nd), V_beta=diag(rep(10,nd)),
#   mu_beta = rep(0,nd),
#   beta_sp_start=matrix(0,np,nsp), V_beta_sp=diag(rep(10,np)),
#   mu_beta_sp = rep(0,np),
#   seed=1234, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# pdf(file="~/Documents/jSDM.test/results_binomial_probit_block_traits_long_format.pdf")
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
# 
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
#   mean_beta_sp[j,] <-apply(mod$beta_sp[,j,1:np],2,mean)
#   if(j<5){
#     for (p in 1:np) {
#       MCMC.beta_spj <- coda::mcmc(mod$beta_sp[,j,1:np], start=nburn+1, end=ngibbs, thin=nthin)
#       summary(MCMC.beta_spj)
#       coda::traceplot(MCMC.beta_spj[,p])
#       coda::densplot(MCMC.beta_spj[,p], main = paste0("beta_sp",j,p))
#       abline(v=beta_sp.target[p,j],col='red')
#     }
#   }
# }
# 
# ## Fixed pecies effect beta_sp
# par(mfrow=c(1,1))
# ## beta_sp
# plot(t(beta_sp.target),mean_beta_sp, xlab="obs", ylab="fitted", main="Fixed species effect beta")
# abline(a=0,b=1,col='red')
# 
# # Deviance
# colnames(mod$Deviance) <- "Deviance"
# mcmc.Deviance <- coda::mcmc(mod$Deviance, start=nburn+1, end=ngibbs, thin=nthin)
# plot(mcmc.Deviance)
# 
# ## Fitting Generalized Linear Model
# data$species <- as.factor(data$species)
# glm <- glm( Y ~  1 + x1 + x1.2 + x1.SLA + species + x1:species + x1.2:species, family=binomial(link="probit"), data=data)
# par(mfrow=c(1,2))
# plot(pnorm(mod$probit_theta_pred),glm$fitted.values, main="theta",
#      xlab="fitted by jSDM", ylab="fitted by glm")
# abline(a=0,b=1,col='red')
# plot(pnorm(probit_theta),glm$fitted.values, main="theta",
#      xlab="obs", ylab="fitted by glm")
# abline(a=0,b=1,col='red')
# par(mfrow=c(1,3))
# for(p in 1:np){
#   plot(mean_beta_sp[-1,p], glm$coefficients[(nd+1+(nsp-1)*(p-1)):(nd+(nsp-1)*p)])
#   abline(a=0,b=1,col='red')
# }
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
# segments(x0 = x0s, x1 = x1s, y0=glm$coefficients[c("(Intercept)", "x1","x1.2", "x1.SLA")], col = "blue", lwd=2)
# legend("topright", legend=c("beta target", "beta fitted by glm"), col=c('red','blue'), lty=1, lwd=2)
# dev.off()
# test en SDM logit binomial
# voir obs dans IC
*/
