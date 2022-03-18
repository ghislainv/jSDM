#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Rcpp_jSDM_useful.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

/* ************************ */
/* Gibbs sampler function */

// [[Rcpp::export]]
Rcpp::List  Rcpp_jSDM_binomial_logit_traits_fixed_site(
    const int ngibbs, int nthin, int nburn, // Number of iterations, burning and samples
    const arma::umat &Y, // Number of successes (presences)
    const arma::uvec &T, // Number of trials
    const arma::mat &X, // Suitability covariates
    const arma::mat &Tr, // Species traits
    const arma::mat &gamma_zeros, // Interactions between covariates and traits to consider
    const arma::mat &beta_start,
    const arma::mat &gamma_start,
    const arma::vec &alpha_start,//alpha
    const double &V_alpha,
    const arma::mat &V_gamma, // Priors 
    const arma::mat &mu_gamma, 
    const arma::vec &V_beta,
    const int &seed, // Various 
    const double &ropt,
    const int &verbose) {
  
  ////////////////////////////////////////
  // Defining and initializing objects //
  
  ////////////////////////////////////////
  // Initialize random number generator //
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, seed);
  
  ///////////////////////////
  // Redefining constants //
  const int NGIBBS = ngibbs;
  const int NTHIN = nthin;
  const int NBURN = nburn;
  const int NSAMP = (NGIBBS-NBURN)/NTHIN;
  const int NSITE = X.n_rows;
  const int NP = X.n_cols;
  const int NSP = Y.n_cols;
  const int NL=0;
  const int NT = Tr.n_cols;
  
  ////////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::Cube<double> gamma; gamma.zeros(NSAMP,NT,NP);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  /* Latent variable */
  arma::mat logit_theta_run; logit_theta_run.zeros(NSITE, NSP);
  arma::mat logit_theta_latent; logit_theta_latent.zeros(NSITE, NSP);
  arma::mat theta_run; theta_run.zeros(NSITE, NSP);
  arma::mat theta_latent; theta_latent.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  // gamma 
  arma::mat gamma_run = gamma_start;
  
  //////////////////////////////////////////////////////////
  // Set up and initialize structure for density function //
  dens_par dens_data;
  // Data 
  dens_data.NSITE = NSITE;
  dens_data.NSP = NSP;
  dens_data.NL = NL;
  // Y
  dens_data.Y = Y;
  // T
  dens_data.T = T;
  // Suitability process 
  dens_data.NP = NP;
  dens_data.X = X;
  // beta
  dens_data.pos_beta = 0;
  dens_data.sp_beta = 0;
  dens_data.mu_beta = arma::zeros(NP);
  dens_data.V_beta = V_beta;
  dens_data.beta_run = beta_start;
  // alpha
  dens_data.site_alpha = 0;
  dens_data.V_alpha_run = V_alpha;
  dens_data.alpha_run = alpha_start.t();
  dens_data.shape = 0.0;
  dens_data.rate = 0.0;
  // lambda 
  dens_data.pos_lambda = 0;
  dens_data.sp_lambda = 0;
  dens_data.mu_lambda = arma::zeros(NL);
  dens_data.V_lambda = arma::zeros(NL);
  dens_data.lambda_run = arma::zeros(NL,NSP);
  // W
  dens_data.site_W = 0;
  dens_data.pos_W = 0;
  dens_data.V_W = arma::zeros(NL);
  dens_data.W_run = arma::zeros(NSITE,NL);
  
  ////////////////////////////////////////////////////////////
  // Proposal variance and acceptance for adaptive sampling //
  
  // beta
  arma::mat sigmap_beta; sigmap_beta.ones(NSP,NP);
  arma::mat nA_beta; nA_beta.zeros(NSP,NP);
  arma::mat Ar_beta; Ar_beta.zeros(NSP,NP); // Acceptance rate
  
  // alpha
  arma::vec sigma_alpha; sigma_alpha.ones(NSITE);
  arma::vec nA_alpha; nA_alpha.zeros(NSITE);
  arma::vec Ar_alpha; Ar_alpha.zeros(NSITE); // Acceptance rate
  
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  
  ////////////////////
  // Gibbs sampler //
  
  for ( int g = 0; g < NGIBBS; g++ ) {
    
    for ( int i = 0; i < NSITE; i++ ) {
      // alpha
      if(i==0){
        // constraints of identifiability on alpha
        dens_data.alpha_run(i) = 0.0;
      } else {
        dens_data.site_alpha = i; // Specifying the site 
        double x_now = dens_data.alpha_run(i);
        double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigma_alpha(i));
        double p_now = alphadens_logit(x_now, &dens_data);
        double p_prop = alphadens_logit(x_prop, &dens_data);
        double ratio = std::exp(p_prop - p_now); // ratio
        double z = gsl_rng_uniform(r);
        // Actualization
        if ( z < ratio ) {
          dens_data.alpha_run(i) = x_prop;
          nA_alpha(i)++;
        } 
      }
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on sites 
    
    /////////////////////////////////////////////
    //mat gamma: Gibbs algorithm //
    for (int t= 0; t < NT; t++ ) {
      for (int p= 0; p < NP; p++ ) {
        if(gamma_zeros(t,p)!=0.0){
          arma::vec gamma_prop=gamma_run.col(p);
          gamma_prop(t)=0.0;
          // small_v
          double small_v = arma::as_scalar(1.0/(V_gamma(t,p))*mu_gamma(t,p) +
                                           Tr.col(t).t()*(dens_data.beta_run.row(p).t()-Tr*gamma_prop));
          // big_V
          double big_V = 1.0/arma::as_scalar((1.0/V_gamma(t,p) + Tr.col(t).t()*Tr.col(t)));
          // Draw in the posterior distribution
          gamma_run(t,p)= big_V*small_v + gsl_ran_gaussian_ziggurat(r, std::sqrt(big_V));
        }
      }
    }
    
    for ( int j=0; j < NSP; j++ ) {
      dens_data.mu_beta = (Tr.row(j)*gamma_run).t(); // compute mu_j 
      // beta
      dens_data.sp_beta = j; // Specifying the species
      for ( int p = 0; p < NP; p++ ) {
        dens_data.pos_beta = p; // Specifying the rank of the parameter of interest
        double x_now = dens_data.beta_run(p,j);
        double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigmap_beta(j,p));
        double p_now = betadens_logit(x_now, &dens_data);
        double p_prop = betadens_logit(x_prop, &dens_data);
        double ratio = std::exp(p_prop - p_now); // ratio
        double z = gsl_rng_uniform(r);
        // Actualization
        if ( z < ratio ) {
          dens_data.beta_run(p,j) = x_prop;
          nA_beta(j,p)++;
        } 
      } // loop on rank of parameters
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on species
    
    
    ///////////////
    // Deviance //
    
    // logLikelihood
    double logL = 0.0;
    for ( int i = 0; i < NSITE; i++ ) {
      for ( int j = 0; j < NSP; j++ ) {
        /* theta */
        double Xpart_theta = 0.0;
        for ( int p = 0; p < NP; p++ ) {
          Xpart_theta += dens_data.X(i,p) * dens_data.beta_run(p,j);
        }
        Xpart_theta += dens_data.alpha_run(i);
        logit_theta_run(i,j) = Xpart_theta;
        theta_run(i,j) = invlogit(Xpart_theta);
        /* log Likelihood */
        logL += R::dbinom(dens_data.Y(i,j), dens_data.T(i), theta_run(i,j), 1);
      } // loop on species
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on sites
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    
    /////////////
    // Output //
    if (((g+1)>NBURN) && (((g+1-NBURN)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        beta.tube(isamp-1,j) = dens_data.beta_run.col(j);
      }// loop on species
      for(int t=0; t<NT; t++){
        gamma.tube(isamp-1,t) = gamma_run.row(t);
      }
      alpha.row(isamp-1) = dens_data.alpha_run;
      Deviance(isamp-1) = Deviance_run;
      // We compute the mean of NSAMP values
      logit_theta_latent+= logit_theta_run/NSAMP; 
      theta_latent+=theta_run/NSAMP; 
    }
    
    ///////////////////////////////////////////////
    // Adaptive sampling (on the burnin period) //
    const double ROPT = ropt;
    int DIV = 0;
    if ( NGIBBS >= 1000 ) DIV=100;
    else DIV = NGIBBS / 10;
    /* During the burnin period */
    if ( (g+1)%DIV== 0 && (g+1)<=NBURN ) {
      for (int j=0; j<NSP; j++) {
        for ( int p=0; p<NP; p++ ) {
          Ar_beta(j,p) = ((double) nA_beta(j,p)) / DIV;
          if ( Ar_beta(j,p) >= ROPT )
            sigmap_beta(j,p) = sigmap_beta(j,p)*(2-(1-Ar_beta(j,p)) / (1-ROPT));
          else sigmap_beta(j,p) = sigmap_beta(j,p) / (2-Ar_beta(j,p) / ROPT);
          nA_beta(j,p) = 0.0; // We reinitialize the number of acceptance to zero for beta
        } // loop on rank of parameters
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on species 
      for (int i=0; i<NSITE; i++) {
        Ar_alpha(i) = ((double) nA_alpha(i)) / DIV;
        if ( Ar_alpha(i) >= ROPT ) sigma_alpha(i) = sigma_alpha(i) * (2-(1-Ar_alpha(i)) / (1-ROPT));
        else sigma_alpha(i) = sigma_alpha(i) / (2-Ar_alpha(i) / ROPT);
        nA_alpha(i) = 0.0; // We reinitialize the number of acceptance for alpha to zero
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on sites
    }
    
    /* After the burnin period */
    if ( (g+1) % DIV == 0 && (g+1) > NBURN ) {
      for (int j=0; j<NSP; j++) {
        for (int p=0; p<NP; p++) {
          Ar_beta(j,p) = ((double) nA_beta(j,p)) / DIV;
          nA_beta(j,p) = 0.0; // We reinitialize the number of acceptance to zero for beta
        } // loop on rank of parameters
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on species
      for (int i=0; i<NSITE; i++) {
        Ar_alpha(i) = ((double) nA_alpha(i)) / DIV;
        nA_alpha(i) = 0.0; // We reinitialize the number of acceptance for alpha to zero
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on sites
    }
    
    //////////////////////////////////////////////////
    // Progress bar
    double Perc = 100 * (g+1) / (NGIBBS);
    if ( (g+1) % (NGIBBS/100) == 0 && verbose == 1) {
      Rprintf("*");
      R_FlushConsole();
      if( (g+1) % (NGIBBS/10) == 0 ) {
        double mAr_beta=0; // Mean acceptance rate of beta
        for ( int j = 0; j < NSP; j++ ) {
          for ( int p = 0; p < NP; p++ ) {
            mAr_beta += Ar_beta(j,p) / (NSP*NP);
          } // loop on rank of parameters
        } // loop on species
        
        double mAr_alpha=0; // Mean acceptance rate of alphas
        for ( int i = 0; i < NSITE; i++ ) {
          mAr_alpha += Ar_alpha(i) / NSITE;
        }// loop on sites
        Rprintf(":%.1f%%, mean accept. rates= beta:%.3f alpha:%4.3f\n",
                Perc, mAr_beta, mAr_alpha);
        R_FlushConsole();
      }
    }
    
    //////////////////////////////////////////////////
    // User interrupt
    R_CheckUserInterrupt(); // allow user interrupt
    
  } // Gibbs sampler
  
  // Free memory
  gsl_rng_free(r);
  
  // Return results as a Rcpp::List
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                          Rcpp::Named("gamma") = gamma,
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("logit_theta_latent") = logit_theta_latent,
                                          Rcpp::Named("theta_latent") = theta_latent);
  
  return results;
  
}// end Rcpp_jSDM_binomial_logit_traits_fixed_site function

// Test
/*** R
# library(coda)
# library(jSDM)
# 
# nsp <- 70
# nsite <- 210
# seed <- 1234
# set.seed(seed)
# visits<- rpois(nsite,3)
# visits[visits==0] <- 1
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# np <- ncol(X)
# colnames(X) <- c("Int","x1","x2")
# Tr <- data.frame(Int=1, WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
# nt <- ncol(Tr)
# gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
# mu_beta <- as.matrix(Tr) %*% gamma.target
# V_beta <- diag(1,np)
# beta.target <- matrix(NA,ncol=np,nrow=nsp)
# library(MASS)
# for(j in 1:nsp){
#   beta.target[j,] <- mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
# }
# alpha.target <- runif(nsite,-2,2)
# alpha.target[1] <- 0
# logit.theta <- X %*% t(beta.target) + alpha.target
# theta <- inv_logit(logit.theta)
# Y <- apply(theta, 2, rbinom, n=nsite, size=visits)
# 
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# 
# # Call to C++ function
# mod <- Rcpp_jSDM_binomial_logit_traits_fixed_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                            Y=Y, T=visits, X=X, Tr=as.matrix(Tr),
#                                            beta_start=matrix(0,np,nsp),
#                                            alpha_start=rep(0,nsite),
#                                            gamma_start=matrix(0,nt,np),
#                                            gamma_zeros=matrix(1,nt,np),
#                                            mu_gamma=matrix(0,nt,np),
#                                            V_gamma=matrix(10,nt,np),
#                                            V_alpha=10,
#                                            V_beta=rep(10,np),
#                                            seed=1234, ropt=0.44, verbose=1)
# 
# # Parameter estimates
# # ## gamma
# par(mfrow=c(2,2))
# for(p in 1:np){
#   MCMC.gamma_p <- coda::mcmc(mod$gamma[,,p], start=nburn+1, end=ngibbs, thin=nthin)
#   for(t in 1:nt){
#     coda::traceplot(MCMC.gamma_p[,t])
#     coda::densplot(MCMC.gamma_p[,t], main = paste0("gamma_",colnames(X)[p],".",colnames(Tr)[t]))
#     abline(v=gamma.target[t,p],col='red')
#   }
# }
# ##alpha
# par(mfrow=c(1,1),oma=c(1, 0, 1.4, 0))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="fitted",
#      xlab="obs",main="Fixed site effect alpha",cex.main=1.4)
# abline(a=0,b=1,col='red')
# 
# ## species effect beta
# par(mfrow=c(np,2))
# for (j in 1:4) {
#   for (p in 1:np) {
#     MCMC.betaj <- coda::mcmc(mod$beta[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#     summary(MCMC.betaj)
#     coda::traceplot(MCMC.betaj[,p])
#     coda::densplot(MCMC.betaj[,p], main = paste0("beta",p, " sp", j))
#     abline(v=beta.target[j,p],col='red')
#   }
# }
# 
# mean_beta <- matrix(0,nsp,np)
# mean_beta <-apply(mod$beta,c(2,3),mean)
# par(mfrow=c(1,1),oma=c(1, 0, 1, 0))
# plot(beta.target,mean_beta, xlab="obs", ylab="fitted",main="beta")
# title("Fixed species effects", outer = T)
# abline(a=0,b=1,col='red')
# 
# ## Deviance
# mean(mod$Deviance)
# 
# # Predictions
# ##logit_theta
# par(mfrow=c(1,2))
# plot(logit.theta,mod$logit_theta_latent, ylab ="fitted",
#      xlab="obs", main="logit(theta)")
# abline(a=0,b=1,col='red')
# ##theta
# plot(theta,mod$theta_latent,ylab ="fitted",
#      xlab="obs", main="Probabilities of occurrence theta")
# abline(a=0,b=1,col='red')
*/

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////