#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Rcpp_jSDM_useful.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

/* ************************ */
/* Gibbs sampler function */

// [[Rcpp::export]]
Rcpp::List  Rcpp_jSDM_poisson_log_rand_site(
    const int ngibbs, int nthin, int nburn, // Number of iterations, burning and samples
    const arma::umat &Y, // Number of successes (presences)
    const arma::mat &X, // Suitability covariates
    arma::mat beta_start,
    arma::vec alpha_start,//alpha
    double V_alpha_start,
    arma::vec mu_beta, // Priors 
    arma::vec V_beta,
    double shape,
    double rate,
    const int seed, // Various 
    const double ropt,
    const int verbose) {
  
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
  
  ////////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::mat alpha; alpha.zeros(NSAMP, NSITE);
  arma::vec V_alpha; V_alpha.zeros(NSAMP);
  /* Latent variable */
  arma::mat theta_run; theta_run.zeros(NSITE, NSP);
  arma::mat theta_latent; theta_latent.zeros(NSITE, NSP);
  /* Deviance */
  arma::vec Deviance; Deviance.zeros(NSAMP);
  
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
  dens_data.T = NA_REAL;
  // Suitability process 
  dens_data.NP = NP;
  dens_data.X = X;
  // beta
  dens_data.pos_beta = 0;
  dens_data.sp_beta = 0;
  dens_data.mu_beta = mu_beta;
  dens_data.V_beta = V_beta;
  dens_data.beta_run = beta_start;
  // lambda 
  dens_data.pos_lambda =  NA_REAL;
  dens_data.sp_lambda =  NA_REAL;
  dens_data.mu_lambda = NA_REAL;
  dens_data.V_lambda = NA_REAL;
  dens_data.lambda_run = NA_REAL;
  // W
  dens_data.site_W = NA_REAL;
  dens_data.pos_W = NA_REAL;
  dens_data.V_W = NA_REAL;
  dens_data.W_run = NA_REAL;
  // alpha
  dens_data.site_alpha = 0;
  dens_data.shape = shape;
  dens_data.rate = rate;
  dens_data.V_alpha_run = V_alpha_start;
  dens_data.alpha_run = alpha_start.t();
  
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
      dens_data.site_alpha = i; // Specifying the site 
      double x_now = dens_data.alpha_run(i);
      double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigma_alpha(i));
      double p_now = alphadens_pois(x_now, &dens_data);
      double p_prop = alphadens_pois(x_prop, &dens_data);
      double ratio = std::exp(p_prop - p_now); // ratio
      double z = gsl_rng_uniform(r);
      // Actualization
      if ( z < ratio ) {
        dens_data.alpha_run(i) = x_prop;
        nA_alpha(i)++;
      } 
    } // loop on sites
    
    // Centering alpha 
    dens_data.alpha_run = dens_data.alpha_run - arma::mean(dens_data.alpha_run);
    
    // V_alpha
    double sum = arma::as_scalar(dens_data.alpha_run*dens_data.alpha_run.t());
    double shape_posterior = dens_data.shape + 0.5*NSITE;
    double rate_posterior = dens_data.rate + 0.5*sum;
    
    dens_data.V_alpha_run = rate_posterior/gsl_ran_gamma_mt(r, shape_posterior, 1.0);
    
    for ( int j = 0; j < NSP; j++ ) {
      // beta
      dens_data.sp_beta = j; // Specifying the species
      for ( int p = 0; p < NP; p++ ) {
        dens_data.pos_beta = p; // Specifying the rank of the parameter of interest
        double x_now = dens_data.beta_run(p,j);
        double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigmap_beta(j,p));
        double p_now = betadens_pois(x_now, &dens_data);
        double p_prop = betadens_pois(x_prop, &dens_data);
        double ratio = std::exp(p_prop - p_now); // ratio
        double z = gsl_rng_uniform(r);
        // Actualization
        if ( z < ratio ) {
          dens_data.beta_run(p,j) = x_prop;
          nA_beta(j,p)++;
        } 
      } // loop on rank of parameters
    } // loop on species
    
    
    ///////////////
    // Deviance //
    
    // logLikelihood
    double logL = 0.0;
    for ( int i = 0; i < NSITE; i++ ) {
      for ( int j = 0; j < NSP; j++ ) {
        /* theta */
        double log_theta = 0.0;
        for ( int p = 0; p < NP; p++ ) {
          log_theta += dens_data.X(i,p) * dens_data.beta_run(p,j);
        }
        log_theta += dens_data.alpha_run(i);
        theta_run(i,j) = exp(log_theta);
        /* log Likelihood */
        logL += R::dpois(dens_data.Y(i,j), theta_run(i,j), 1);
      } // loop on species
    } // loop on sites
    
    // Deviance
    double Deviance_run = -2 * logL;
    
    
    /////////////
    // Output //
    if (((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
      int isamp=((g+1)-NBURN)/(NTHIN);
      for ( int j=0; j<NSP; j++ ) {
        beta.tube(isamp-1,j) = dens_data.beta_run.col(j);
        for ( int i=0; i<NSITE; i++ ) {
          theta_latent(i,j) += theta_run(i,j) / NSAMP; // We compute the mean of NSAMP values
        }//loop on sites
      }// loop on species
      alpha.row(isamp-1) = dens_data.alpha_run;
      V_alpha(isamp-1) = dens_data.V_alpha_run;
      Deviance(isamp-1) = Deviance_run;
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
      } // loop on species 
      for (int i=0; i<NSITE; i++) {
        Ar_alpha(i) = ((double) nA_alpha(i)) / DIV;
        if ( Ar_alpha(i) >= ROPT ) sigma_alpha(i) = sigma_alpha(i) * (2-(1-Ar_alpha(i)) / (1-ROPT));
        else sigma_alpha(i) = sigma_alpha(i) / (2-Ar_alpha(i) / ROPT);
        nA_alpha(i) = 0.0; // We reinitialize the number of acceptance for alpha to zero
      } // loop on sites
    }
    
    /* After the burnin period */
    if ( (g+1) % DIV == 0 && (g+1) > NBURN ) {
      for (int j=0; j<NSP; j++) {
        for (int p=0; p<NP; p++) {
          Ar_beta(j,p) = ((double) nA_beta(j,p)) / DIV;
          nA_beta(j,p) = 0.0; // We reinitialize the number of acceptance to zero for beta
        } // loop on rank of parameters
      } // loop on species
      for (int i=0; i<NSITE; i++) {
        Ar_alpha(i) = ((double) nA_alpha(i)) / DIV;
        nA_alpha(i) = 0.0; // We reinitialize the number of acceptance for alpha to zero
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
                                          Rcpp::Named("alpha") = alpha,
                                          Rcpp::Named("V_alpha") = V_alpha,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("theta_latent") = theta_latent);
  
  return results;
  
}// end Rcpp_jSDM_poisson_log_rand_site function

// Test
/*** R
# library(coda)
# 
# nsp <- 50
# nsite <- 150
# seed <- 1234
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# np <- ncol(X)
# beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
# V_alpha.target <- 0.5
# alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
# log.theta <- X %*% t(beta.target) + alpha.target
# theta <- exp(log.theta)
# Y <- apply(theta, 2, rpois, n=nsite)
# 
# # Iterations
# nsamp <- 2000
# nburn <- 2000
# nthin <- 2
# ngibbs <- nsamp+nburn
# 
# # Call to C++ function
# mod <- Rcpp_jSDM_poisson_log_rand_site(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                        Y=Y, X=X,
#                                        beta_start=matrix(0,np,nsp),
#                                        alpha_start=rep(0,nsite),
#                                        V_alpha_start=1, shape = 0.5, rate = 0.0005,
#                                        mu_beta=rep(0,np), V_beta=rep(1.0E6,np),
#                                        seed=1234, ropt=0.44, verbose=1)
# 
# # Parameter estimates
# ##alpha
# par(mfrow=c(1,3),oma=c(1, 0, 1.4, 0))
# MCMC_alpha <- coda::mcmc(mod$alpha, start=nburn+1, end=ngibbs, thin=nthin)
# plot(alpha.target,summary(MCMC_alpha)[[1]][,"Mean"], ylab ="fitted",
#      xlab="obs",main="alpha",cex.main=1.4)
# abline(a=0,b=1,col='red')
# ##V_alpha
# title(main="Random site effect and its variance",outer=T,cex.main=1.8)
# MCMC_V_alpha <- coda::mcmc(mod$V_alpha, start=nburn+1, end=ngibbs, thin=nthin)
# coda::traceplot(MCMC_V_alpha,main="Trace V_alpha",cex.main=1.4)
# coda::densplot(MCMC_V_alpha,main="Density V_alpha",cex.main=1.4)
# abline(v=V_alpha.target,col='red')
# legend("topright", lty=1, col='red',legend="V_alpha obs",cex=0.8)
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
# ##log_theta
# par(mfrow=c(1,2),oma=c(1, 0, 1, 0))
# log_theta_pred <- apply(mod$theta_latent,c(1,2),log)
# plot(log.theta,log_theta_pred, ylab ="fitted",
#      xlab="obs", main="logit(theta)")
# title(main="Probabilities of occurrence",outer=T)
# abline(a=0,b=1,col='red')
# ##theta
# plot(theta,mod$theta_latent,ylab ="fitted",
#      xlab="obs", main="theta")
# abline(a=0,b=1,col='red')

*/

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////