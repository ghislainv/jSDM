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
Rcpp::List Rcpp_jSDM_binomial_probit_traits_fixed_site_long_format(
    const int ngibbs, int nthin, int nburn, 
    const arma::uvec& Y, 
    const arma::uvec& Id_sp,
    const arma::uvec& Id_site,
    const arma::mat& X,
    const arma::mat& D,
    const arma::mat& beta_start,
    const arma::mat& V_beta,
    const arma::vec& mu_beta,
    const arma::vec& gamma_start,
    const arma::mat& V_gamma,
    const arma::vec& mu_gamma,
    arma::vec alpha_start,
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
  const int NP = X.n_cols;
  const int ND = D.n_cols;
  const int NSITE=alpha_start.n_elem;
  const int NSP = beta_start.n_cols;
  const int NOBS = Y.n_elem;
  
  
  ///////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::mat gamma; gamma.zeros(NSAMP, ND);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  /* Latent variable */
  arma::vec probit_theta_latent; probit_theta_latent.zeros(NOBS);
  arma::vec theta_latent; theta_latent.zeros(NOBS);
  arma::vec Z_latent; Z_latent.zeros(NOBS);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
  /////////////////////////////////////
  // Initializing running parameters //
  
  //  mat of species effects parameters (np,nsp)
  arma::mat beta_run = beta_start;
  arma::mat gamma_run = gamma_start;
  // alpha vec of sites effects (nsite)
  arma::vec alpha_run = alpha_start;
  // Z latent (nobs)
  arma::vec Z_run; Z_run.zeros(NOBS);
  // probit(theta_ij) = X_i*beta_j + D_ij*gamma + alpha_i
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
  arma::mat inv_Vgamma=inv(V_gamma);
  arma::mat bigV_gamma=inv(inv_Vgamma + D.t()*D);
  arma::mat chol_decomp_bigV_gamma=chol_decomp(bigV_gamma);
  
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
    // vec alpha : Gibbs algorithm 
    
    for (int i=0; i<NSITE; i++) {
      if(i==0){
        // constraints of identifiability on alpha
        alpha_run(i) = 0.0;
      } else {
        int nobs_i = rowId_site[i].n_elem;
        double small_v2=0.0;
        // small_v
        for (int n=0; n<nobs_i; n++) {
          small_v2 += arma::as_scalar(Z_run.row(rowId_site[i].at(n))
                                        -X.row(rowId_site[i].at(n))*beta_run.col(Id_sp(rowId_site[i].at(n)))
                                        -D.row(rowId_site[i].at(n))*gamma_run);
        }
        // big_V
        double big_V2 = 1/(1/V_alpha + nobs_i);
        
        // Draw in the posterior distribution
        alpha_run(i) = big_V2*small_v2 + gsl_ran_gaussian_ziggurat(s, std::sqrt(big_V2));
      }
      R_CheckUserInterrupt(); // allow user interrupt
    }
    //////////////////////////////////
    // vec gamma: Gibbs algorithm //
    
    // Loop on species
    arma::vec small_v; small_v.zeros(ND);
    for (int j=0; j<NSP; j++) {
      R_CheckUserInterrupt(); // allow user interrupt
      // small_v
      small_v += D.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j]) 
                                            - X.rows(rowId_sp[j])*beta_run.col(j)
                                            - alpha_run(Id_site(rowId_sp[j])));
    }
    // small_v
    small_v += inv_Vgamma*mu_gamma;
    
    // Draw in the posterior distribution
    gamma_run = arma_mvgauss(s, bigV_gamma*small_v, chol_decomp_bigV_gamma);
    
    //////////////////////////////////
    // mat beta: Gibbs algorithm //
    // Loop on species
    for (int j=0; j<NSP; j++) {
      // small_v
      arma::vec small_v = inv_Vbeta*mu_beta + X.rows(rowId_sp[j]).t()*(Z_run(rowId_sp[j])-D.rows(rowId_sp[j])*gamma_run-alpha_run(Id_site(rowId_sp[j])));
      
      // big_V
      arma::mat big_V = inv(inv_Vbeta + X.rows(rowId_sp[j]).t()*X.rows(rowId_sp[j]));
      
      // Draw in the posterior distribution
      beta_run.col(j) = arma_mvgauss(s, big_V*small_v, chol_decomp(big_V));
      // if(j==0 && !Id_common_var.has_nan()){
      //   // constraint of identifiability on beta
      //   // Id of columns in common between D and X 
      //   arma::uvec Id_sp0; Id_sp0.zeros(1);
      //   beta_run.submat(conv_to<uvec>::from(Id_common_var),Id_sp0).fill(0.0);
      // }
      R_CheckUserInterrupt(); // allow user interrupt
    }
    
    //////////////////////////////////////////////////
    //// Deviance
    
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < NOBS; n++ ) {
      // probit(theta_ij) = X_i*beta_j + D_ij*gamma
      probit_theta_run(n) = arma::as_scalar(X.row(n)*beta_run.col(Id_sp(n)) + D.row(n)*gamma_run + alpha_run(Id_site(n)));
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
      gamma.row(isamp-1) = gamma_run.t();
      alpha.row(isamp-1) = alpha_run.t();
      Deviance(isamp-1) = Deviance_run;
      // We compute the mean of NSAMP values
      Z_latent += Z_run/NSAMP;
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
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("Z_latent") = Z_latent,
                                          Rcpp::Named("theta_latent") = theta_latent,
                                          Rcpp::Named("probit_theta_latent") = probit_theta_latent );  
  return results;
  
} // end Rcpp_jSDM_binomial_probit_traits_fixed_site_long_format

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
# X <- cbind(rep(1,nsite),x1)
# colnames(X) <- c("Int","x1")
# np <- ncol(X)
# beta.target <- t(matrix(runif(nsp*np,-1,1), byrow=TRUE, nrow=nsp))
# # constraint of identifiability
# beta.target[,1] <- 0.0
# SLA <- runif(nsp,-2,2)
# D <- data.frame(x1.2=x1.2, x1.SLA= scale(c(x1 %*% t(SLA))))
# nd <- ncol(D)
# gamma.target <-runif(nd,-1,1)
# alpha.target <- runif(nsite,-1,1)
# alpha.target[1] <- 0
# probit_theta <- c(X %*% beta.target) + as.matrix(D) %*% gamma.target + alpha.target
# # x1_supObs <- rnorm(nsite)
# # X_supObs <- cbind(rep(1,nsite),x1_supObs)
# # D_supObs <- data.frame( 1.2=scale(x1_supObs^2), x1.SLA= scale(c(x1_supObs %*% t(SLA))))
# # probit_theta_supObs <- c(X_supObs%*%beta.target) + as.matrix(D_supObs) %*% gamma.target
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
# #                    intercept=rep(1,nobs), x1=c(x1,x1_supObs),
# #                    x1.2=c(D$x1.2,D_supObs$x1.2),x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# data <- data.frame(site=Id_site, species=Id_sp, Y=Y,
#                    intercept=rep(1,nobs), x1=x1,
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
# X=as.matrix(data[,c("intercept","x1")])
# D=as.matrix(data[,c("x1.2","x1.SLA")])
# # Call to C++ function
# # Iterations
# nsamp <- 10000
# nburn <- 10000
# nthin <- 10
# ngibbs <- nsamp+nburn
# T1 <- Sys.time()
# mod <- Rcpp_jSDM_binomial_probit_traits_fixed_site_long_format(
#   ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#   Y=data$Y, X=X, D=D, Id_site=data$site, Id_sp=data$species,
# #   Id_common_var=c(0,1,2),
#   gamma_start=rep(0,nd), V_gamma=diag(rep(10,nd)),
#   mu_gamma = rep(0,nd),
#   beta_start=matrix(0,np,nsp), V_beta=diag(rep(10,np)),
#   mu_beta = rep(0,np),
#   alpha_start=rep(0,nsite), V_alpha=10,
#   seed=1234, verbose=1)
# T2 <- Sys.time()
# T <- difftime(T2,T1)
# # ===================================================
# # Result analysis
# # ===================================================
# pdf(file="~/Documents/jSDM.test/results_binomial_probit_traits_fixed_site_long_format.pdf")
# par(mfrow=c(1,2))
# hist(probit_theta)
# # Parameter estimates
# ## theta
# #plot(pnorm(probit_theta[-1]),mod$theta_latent)
# plot(pnorm(probit_theta),mod$theta_latent)
# abline(a=0,b=1,col='red')
# ## probit_theta
# par(mfrow=c(1,2))
# #plot(probit_theta[-1],mod$probit_theta_latent, xlab="obs", ylab ="fitted", main="probit(theta)")
# plot(probit_theta,mod$probit_theta_latent, xlab="obs", 
#      ylab ="fitted", main="probit(theta)")
# abline(a=0,b=1,col='red')
# ## Z
# #plot(Z_true[-1],mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# plot(Z_true,mod$Z_latent, xlab="obs", ylab ="fitted", main="Latent variable Z")
# abline(a=0,b=1,col='red')
# 
# ## alpha
# par(mfrow=c(1,1))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="alpha.estimated")
# abline(a=0,b=1,col='red')
# 
# ## gamma
# par(mfrow=c(2,2))
# MCMC.gamma <- coda::mcmc(mod$gamma, start=nburn+1, end=ngibbs, thin=nthin)
# for(d in 1:nd){
#   coda::traceplot(MCMC.gamma[,d])
#   coda::densplot(MCMC.gamma[,d], main = paste0("gamma_",colnames(D)[d]))
#   abline(v=gamma.target[d],col='red')
# }
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
#       coda::densplot(MCMC.betaj[,p], main = paste0("beta_",j,p))
#       abline(v=beta.target[p,j],col='red')
#     }
#   }
# }
# 
# ## Fixed pecies effect beta
# par(mfrow=c(1,1))
# ## beta
# plot(t(beta.target),mean_beta, xlab="obs", ylab="fitted",
#      main="Fixed species effect beta")
# abline(a=0,b=1,col='red')
# 
# # Deviance
# colnames(mod$Deviance) <- "Deviance"
# mcmc.Deviance <- coda::mcmc(mod$Deviance, start=nburn+1, end=ngibbs, thin=nthin)
# plot(mcmc.Deviance)
# 
# # ## Fitting Generalized Linear Model
# # data$species <- as.factor(data$species)
# # data$site <- as.factor(data$site)
# # Model failed to converge in 10000 evaluations with 50 species and 150 sites
# # glm <- lme4::glmer( Y ~  1 + x1 + x1.2 + x1.SLA + species + (1|site) + x1:species + x1.2:species,
# #                     family=binomial(link="probit"), data=data)
# # par(mfrow=c(1,2))
# # plot(pnorm(mod$probit_theta_latent),glm$fitted.values, main="theta",
# #      xlab="fitted by jSDM", ylab="fitted by glm")
# # abline(a=0,b=1,col='red')
# # plot(pnorm(probit_theta),glm$fitted.values, main="theta",
# #      xlab="obs", ylab="fitted by glm")
# # abline(a=0,b=1,col='red')
# # par(mfrow=c(1,3))
# # for(p in 1:np){
# #   plot(mean_beta[-1,p], glm$coefficients[(nd+1+(nsp-1)*(p-1)):(nd+(nsp-1)*p)])
# #   abline(a=0,b=1,col='red')
# # }
# # par(mfrow=c(1,1))
# colnames(mod$gamma) <- paste0("gamma_",colnames(D))
# boxplot(mod$gamma, main="gamma fitted by jSDM")
# # create data for segments
# n <- ncol(mod$gamma)
# # width of each boxplot is 0.8
# x0s <- 1:n - 0.4
# x1s <- 1:n + 0.4
# # these are the y-coordinates for the horizontal lines
# # that you need to set to the desired values.
# y0s <- gamma.target
# # add segments
# segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red", lwd=2)
# legend("topright", legend=c("gamma target"), col=c('red'), lty=1, lwd=2)
# dev.off()
*/
