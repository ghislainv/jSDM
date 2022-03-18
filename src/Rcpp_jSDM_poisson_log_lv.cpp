#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Rcpp_jSDM_useful.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

/* ************************ */
/* Gibbs sampler function */

// [[Rcpp::export]]
Rcpp::List  Rcpp_jSDM_poisson_log_lv(
    const int ngibbs, int nthin, int nburn, // Number of iterations, burning and samples
    const arma::umat &Y, // Number of successes (presences)
    const arma::mat &X, // Suitability covariates
    const arma::mat &W_start, // Starting values
    const arma::mat &lambda_start,
    const arma::mat &beta_start,
    const arma::vec &mu_beta, // Priors 
    const arma::vec &V_beta,
    const arma::vec &mu_lambda,
    const arma::vec &V_lambda,
    const arma::vec &V_W,
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
  const int NL = lambda_start.n_rows;
  
  ////////////////////////////////////////////
  // Declaring new objects to store results //
  /* Parameters */
  arma::Cube<double> beta; beta.zeros(NSAMP, NSP, NP);
  arma::Cube<double> lambda; lambda.zeros(NSAMP, NSP, NL);
  arma::Cube<double> W; W.zeros(NSAMP, NSITE, NL);
  // Latent variable 
  arma::mat log_theta_run; log_theta_run.zeros(NSITE, NSP);
  arma::mat log_theta_latent; log_theta_latent.zeros(NSITE, NSP);
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
  dens_data.T = arma::ones<arma::uvec>(NSITE);
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
  dens_data.pos_lambda = 0;
  dens_data.sp_lambda = 0;
  dens_data.mu_lambda = mu_lambda;
  dens_data.V_lambda = V_lambda;
  dens_data.lambda_run = lambda_start;
  // W
  dens_data.site_W = 0;
  dens_data.pos_W = 0;
  dens_data.V_W = V_W;
  dens_data.W_run = W_start;
  // alpha
  dens_data.site_alpha = 0;
  dens_data.shape = 0.0;
  dens_data.rate = 0.0;
  dens_data.V_alpha_run = 0.0;
  dens_data.alpha_run = arma::zeros<arma::rowvec>(NSITE);
  
  ////////////////////////////////////////////////////////////
  // Proposal variance and acceptance for adaptive sampling //
  
  // beta
  arma::mat sigmap_beta; sigmap_beta.ones(NSP,NP);
  arma::mat nA_beta; nA_beta.zeros(NSP,NP);
  arma::mat Ar_beta; Ar_beta.zeros(NSP,NP); // Acceptance rate
  
  // lambda
  arma::mat sigmaq_lambda; sigmaq_lambda.ones(NSP,NL);
  arma::mat nA_lambda; nA_lambda.zeros(NSP,NL);
  arma::mat Ar_lambda; Ar_lambda.zeros(NSP,NL); // Acceptance rate
  
  // W
  arma::mat sigmaq_W; sigmaq_W.ones(NSITE,NL);
  arma::mat nA_W; nA_W.zeros(NSITE,NL);
  arma::mat Ar_W; Ar_W.zeros(NSITE,NL); // Acceptance rate
  
  ////////////
  // Message//
  Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
  R_FlushConsole();
  
  ////////////////////
  // Gibbs sampler //
  
  for ( int g = 0; g < NGIBBS; g++ ) {
    
    for ( int i = 0; i < NSITE; i++ ) {
      // W
      dens_data.site_W = i; // Specifying the site
      for ( int q = 0; q < NL; q++ ) {
        dens_data.pos_W = q; // Specifying the rank of the latent variable of interest
        double x_now = dens_data.W_run(i,q);
        double x_prop = x_now + gsl_ran_gaussian_ziggurat(r,sigmaq_W(i,q));
        double p_now = Wdens_pois(x_now, &dens_data);
        double p_prop = Wdens_pois(x_prop, &dens_data);
        double ratio = std::exp(p_prop - p_now); // ratio
        double z = gsl_rng_uniform(r);
        // Actualization
        if ( z < ratio ) {
          dens_data.W_run(i,q) = x_prop;
          nA_W(i,q)++;
        } 
      } // loop on rank of latent variable 
      R_CheckUserInterrupt(); // allow user interrupt
    } // loop on sites 
    
    // Centering and reducing W_i
    for ( int i = 0; i < NSITE; i++ ) {
      for ( int q = 0; q < NL; q++ ) {
        dens_data.W_run(i,q) = dens_data.W_run(i,q) - arma::mean(dens_data.W_run.col(q));
        dens_data.W_run(i,q) = dens_data.W_run(i,q)/arma::stddev(dens_data.W_run.col(q));
      }
    }
    
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
      
      // lambda 
      dens_data.sp_lambda = j; // Specifying the species
      for ( int q = 0; q < NL; q++ ) {
        dens_data.pos_lambda = q ; // Specifying the rank of the parameter of interest
        if (q < j ) {
          double x_now = dens_data.lambda_run(q,j);
          double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigmaq_lambda(j,q));
          double p_now = lambdadens_pois(x_now, &dens_data);
          double p_prop = lambdadens_pois(x_prop, &dens_data);
          double ratio = std::exp(p_prop - p_now); // ratio
          double z = gsl_rng_uniform(r);
          // Actualization
          if ( z < ratio ) {
            dens_data.lambda_run(q,j) = x_prop;
            nA_lambda(j,q)++;
          }
        }
        if (q == j) { 
          double x_now = dens_data.lambda_run(q,j);
          double x_prop = x_now + gsl_ran_gaussian_ziggurat(r,sigmaq_lambda(j,q));
          double p_now = lambdaUdens_pois(x_now, &dens_data);
          double p_prop = lambdaUdens_pois(x_prop, &dens_data);
          double ratio = std::exp(p_prop - p_now); // ratio
          double z = gsl_rng_uniform(r);
          // Actualization
          if ( z < ratio ) {
            dens_data.lambda_run(q,j) = x_prop;
            nA_lambda(j,q)++;
          }  
        }
        if (q > j) { 
          dens_data.lambda_run(q,j) = 0;
        } 
      } // loop on rank of latent variable
      R_CheckUserInterrupt(); // allow user interrupt
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
        for ( int q = 0; q < NL; q++ ) {
          log_theta += dens_data.W_run(i,q) * dens_data.lambda_run(q,j);
        }
        log_theta_run(i,j) = log_theta;
        theta_run(i,j) = exp(log_theta);
        /* log Likelihood */
        logL += R::dpois(dens_data.Y(i,j), theta_run(i,j), 1);
        R_CheckUserInterrupt(); // allow user interrupt
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
        lambda.tube(isamp-1,j) = dens_data.lambda_run.col(j);
      }// loop on species
      for ( int i=0; i<NSITE; i++ ) {
        W.tube(isamp-1,i) = dens_data.W_run.row(i);
      }//loop on sites
      Deviance(isamp-1) = Deviance_run;
      // We compute the mean of NSAMP values
      log_theta_latent += log_theta_run/NSAMP; 
      theta_latent += theta_run/NSAMP; 
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
        for ( int q=0; q<NL; q++ ) {
          Ar_lambda(j,q) = ((double) nA_lambda(j,q)) / DIV;
          if ( Ar_lambda(j,q) >= ROPT ) 
            sigmaq_lambda(j,q) = sigmaq_lambda(j,q)*(2-(1-Ar_lambda(j,q)) / (1-ROPT));
          else sigmaq_lambda(j,q) = sigmaq_lambda(j,q) / (2-Ar_lambda(j,q) / ROPT);
          nA_lambda(j,q) = 0.0; // We reinitialize the number of acceptance to zero for lambda 
        } // loop on rank of latent variable
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on species 
      for (int i=0; i<NSITE; i++) {
        for ( int q=0; q<NL; q++ ) {
          Ar_W(i,q) = ((double) nA_W(i,q)) / DIV;
          if ( Ar_W(i,q) >= ROPT ) sigmaq_W(i,q) = sigmaq_W(i,q) * (2-(1-Ar_W(i,q)) / (1-ROPT));
          else sigmaq_W(i,q) = sigmaq_W(i,q) / (2-Ar_W(i,q) / ROPT);
          nA_W(i,q) = 0.0; // We reinitialize the number of acceptance to zero for z
        } // loop on rank of latent variable
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
        for (int q=0; q<NL; q++) {
          Ar_lambda(j,q) = ((double) nA_lambda(j,q)) / DIV;
          nA_lambda(j,q) = 0.0; // We reinitialize the number of acceptance to zero for lambda
        } // loop on rank of latent variable 
        R_CheckUserInterrupt(); // allow user interrupt
      } // loop on species
      for (int i=0; i<NSITE; i++) {
        for (int q=0; q<NL; q++) {
          Ar_W(i,q) = ((double) nA_W(i,q)) / DIV;
          nA_W(i,q) = 0.0; // We reinitialize the number of acceptance to zero for z
        } // loop on rank of latent variable
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
        double mAr_lambda=0; // Mean acceptance rate of lambda
        for ( int j = 0; j < NSP; j++ ) {
          for ( int p = 0; p < NP; p++ ) {
            mAr_beta += Ar_beta(j,p) / (NSP*NP);
          } // loop on rank of parameters
          for ( int q = 0; q < NL; q++ ) {
            mAr_lambda += Ar_lambda(j,q) / (NSP*NL-NL*(NL-1)*0.5);
          } // loop on rank of latent variable 
        } // loop on species
        
        double mAr_W=0; // Mean acceptance rate of W
        for ( int i = 0; i < NSITE; i++ ) {
          for ( int q = 0; q < NL; q++ ) {
            mAr_W += Ar_W(i,q) / (NSITE*NL);
          }// loop on rank of latent variable 
        }// loop on sites
        Rprintf(":%.1f%%, mean accept. rates= beta:%.3f lambda:%.3f W:%.3f\n",
                Perc, mAr_beta, mAr_lambda, mAr_W);
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
                                          Rcpp::Named("lambda") = lambda,
                                          Rcpp::Named("W") = W,
                                          Rcpp::Named("Deviance") = Deviance,
                                          Rcpp::Named("log_theta_latent") = log_theta_latent,
                                          Rcpp::Named("theta_latent") = theta_latent);
  
  return results;
  
}// end Rcpp_jSDM_poisson_log_lv function

// Test
/*** R
# library(coda)
# 
# nsp <- 50
# nsite <- 150
# seed <- 123
# set.seed(seed)
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# np <- ncol(X)
# set.seed(2*seed)
# W <- cbind(rnorm(nsite,0,1),rnorm(nsite,0,1))
# nl <- ncol(W)
# l.zero <- 0
# l.diag <- runif(2,0,2)
# l.other <- runif(nsp*2-3,-2,2)
# lambda.target <- matrix(c(l.diag[1],l.zero,l.other[1],l.diag[2],l.other[-1]), byrow=T, nrow=nsp)
# beta.target <- matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp)
# log.theta <- X %*% t(beta.target) + W %*% t(lambda.target)
# theta <- exp(log.theta)
# Y <- apply(theta, 2, rpois, n=nsite)
# 
# # Iterations
# nsamp <- 5000
# nburn <- 5000
# nthin <- 5
# ngibbs <- nsamp+nburn
# 
# # Call to C++ function
# mod <- Rcpp_jSDM_poisson_log_lv(ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#                                    Y=Y, X=X,
#                                    beta_start=matrix(0,np,nsp),
#                                    lambda_start=matrix(0,nl,nsp),
#                                    W_start=matrix(0,nsite,nl),
#                                    mu_beta=rep(0,np), V_beta=rep(10,np),
#                                    mu_lambda=rep(0,nl), V_lambda=rep(10,nl),
#                                    V_W=rep(1,nl),
#                                    seed=1234, ropt=0.44, verbose=1)
# 
# # Parameter estimates
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
# ## lambda_j
# par(mfrow=c(nl*2,2))
# for (j in 1:4) {
#   for (l in 1:nl) {
#     MCMC.lambdaj <- coda::mcmc(mod$lambda[,j,], start=nburn+1, end=ngibbs, thin=nthin)
#     summary(MCMC.lambdaj)
#     coda::traceplot(MCMC.lambdaj[,l])
#     coda::densplot(MCMC.lambdaj[,l],main = paste0("lambda",l,j))
#     abline(v=lambda.target[j,l],col='red')
#   }
# }
# par(mfrow=c(1,2))
# mean_beta <- matrix(0,nsp,np)
# mean_lambda <- matrix(0,nsp,nl)
# mean_beta <-apply(mod$beta,c(2,3),mean)
# mean_lambda <- apply(mod$lambda,c(2,3),mean)
# plot(beta.target,mean_beta, xlab="obs", ylab="fitted",main="beta")
# abline(a=0,b=1,col="red")
# plot(lambda.target,mean_lambda, xlab="obs", ylab="fitted",main="lambda")
# abline(a=0,b=1,col="red")
# 
# # W latent variables
# par(mfrow=c(1,2),oma=c(1, 0, 1, 0))
# mean_W <- apply(mod$W, c(2,3), mean)
# plot(W[,1],mean_W[,1], main="W1",xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# title("Variables latentes", outer = T)
# plot(W[,2],mean_W[,2], main="W2", xlab="obs", ylab= "fitted")
# abline(a=0,b=1,col='red')
# 
# # lambda * W
# par(mfrow=c(1,1))
# plot(W %*% t(lambda.target),mean_W %*%t(mean_lambda),
#      xlab="obs", ylab= "fitted", main="W_i.lambda_j")
# abline(a=0,b=1,col='red')
# 
# ## Deviance
# mean(mod$Deviance)
# 
# # Predictions
# ##log_theta
# par(mfrow=c(1,2))
# plot(log.theta, mod$log_theta_latent, ylab ="fitted",
#      xlab="obs", main="log(theta)")
# abline(a=0,b=1,col='red')
# ##theta
# plot(theta,mod$theta_latent,ylab ="fitted",
#      xlab="obs", main="Probabilities of occurrence theta")
# abline(a=0,b=1,col='red')
*/

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////