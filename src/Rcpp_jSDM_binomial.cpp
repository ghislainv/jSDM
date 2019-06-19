// ==============================================================================
// author          :Ghislain Vieilledent, Jeanne Clément
// email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Rcpp_jSDM_useful.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

/* ********************************************************************* */
/* dens_par */

struct dens_par {
    /* Data */
    int NOBS;
    arma::uvec Y;
    arma::uvec T;
    /* Suitability */
    int NP;
    int pos_beta;
    arma::mat X;
    arma::vec mubeta;
    arma::vec Vbeta;
    arma::rowvec beta_run;
};


/* ************************************************************ */
/* betadens */

double betadens (double beta_k, void *dens_data) {
    // Pointer to the structure: d
    dens_par *d;
    d = static_cast<dens_par *> (dens_data);
    // Indicating the rank of the parameter of interest
    int k = d->pos_beta;
    // logLikelihood
    double logL = 0.0;
    for ( int n = 0; n < d->NOBS; n++ ) {
        /* theta */
        double Xpart_theta = 0.0;
        for ( int p = 0; p < d->NP; p++ ) {
            if ( p != k ) {
                Xpart_theta += d->X(n, p) * d->beta_run(p);
            }
        }
        Xpart_theta += d->X(n, k) * beta_k;
        double theta = invlogit(Xpart_theta);
        /* log Likelihood */
        logL += R::dbinom(d->Y(n), d->T(n), theta, 1);
    }
    // logPosterior = logL + logPrior
    double logP = logL + R::dnorm(beta_k, d->mubeta(k), std::sqrt(d->Vbeta(k)), 1);
    return logP;
}

/* ************************************************************ */
/* Gibbs sampler function */

// [[Rcpp::export]]
Rcpp::List Rcpp_jSDM_binomial(const int ngibbs, int nthin, int nburn, // Number of iterations, burning and samples
                              arma::uvec Y, // Number of successes (presences)
                              arma::uvec T, // Number of trials
                              arma::mat X, // Suitability covariates
                              arma::vec beta_start,
                              arma::vec mubeta,
                              arma::vec Vbeta,
                              const int seed,
                              const double ropt,
                              const int verbose) {

    ////////////////////////////////////////////////////////////////////////////////
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Defining and initializing objects

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
    const int NOBS = X.n_rows;
    const int NP = X.n_cols;

    ////////////////////////////////////////////
    // Declaring new objects to store results //
    /* Parameters */
    arma::mat beta; beta.zeros(NSAMP, NP);
    /* Latent variable */
    arma::vec theta_run; theta_run.zeros(NOBS);
    arma::vec theta_latent; theta_latent.zeros(NOBS);
    /* Deviance */
    arma::vec Deviance; Deviance.zeros(NSAMP);

    //////////////////////////////////////////////////////////
    // Set up and initialize structure for density function //
    dens_par dens_data;
    /* Data */
    dens_data.NOBS = NOBS;
    // Y
    dens_data.Y = Y;
    // T
    dens_data.T = T;
    /* Suitability process */
    dens_data.NP = NP;
    dens_data.pos_beta = 0;
    dens_data.X = X;
    dens_data.mubeta = mubeta;
    dens_data.Vbeta = Vbeta;
    dens_data.beta_run = beta_start.t();

    ////////////////////////////////////////////////////////////
    // Proposal variance and acceptance for adaptive sampling //

    // beta
    arma::vec sigmap_beta; sigmap_beta.ones(NP);
    arma::vec nA_beta; nA_beta.zeros(NP);
    arma::vec Ar_beta; Ar_beta.zeros(NP); // Acceptance rate

    ////////////
    // Message//
    Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
    R_FlushConsole();

    ///////////////////////////////////////////////////////////////////////////////////////
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Gibbs sampler
    
    for ( int g = 0; g < NGIBBS; g++ ) {

        ////////////////////////////////////////////////
        // beta

        for ( int p = 0; p < NP; p++ ) {
            dens_data.pos_beta = p; // Specifying the rank of the parameter of interest
            double x_now = dens_data.beta_run(p);
            double x_prop = x_now + gsl_ran_gaussian_ziggurat(r, sigmap_beta(p));
            double p_now = betadens(x_now, &dens_data);
            double p_prop = betadens(x_prop, &dens_data);
            double ratio = std::exp(p_prop - p_now); // ratio
            double z = gsl_rng_uniform(r);
            // Actualization
            if ( z < ratio ) {
                dens_data.beta_run(p) = x_prop;
                nA_beta(p)++;
            }
        }


        //////////////////////////////////////////////////
        // Deviance

        // logLikelihood
        double logL = 0.0;
        for ( int n = 0; n < NOBS; n++ ) {
            /* theta */
            double Xpart_theta = 0.0;
            for ( int p = 0; p < NP; p++ ) {
                Xpart_theta += dens_data.X(n, p) * dens_data.beta_run(p);
            }
            theta_run(n) = invlogit(Xpart_theta);
            /* log Likelihood */
            logL += R::dbinom(dens_data.Y(n), dens_data.T(n), theta_run(n), 1);
        }

        // Deviance
        double Deviance_run = -2 * logL;


        //////////////////////////////////////////////////
        // Output
        if ( ((g+1) > NBURN) && (((g+1) % NTHIN) == 0) ) {
            int isamp = ((g+1)-NBURN) / NTHIN;
            beta.row(isamp-1) = dens_data.beta_run;
            Deviance(isamp-1) = Deviance_run;
            for ( int n=0; n<NOBS; n++ ) {
                theta_latent(n) += theta_run(n) / NSAMP; // We compute the mean of NSAMP values
            }
        }


        ///////////////////////////////////////////////////////
        // Adaptive sampling (on the burnin period)
        const double ROPT = ropt;
        int DIV = 0;
        if ( NGIBBS >= 1000 ) DIV=100;
        else DIV = NGIBBS / 10;
        /* During the burnin period */
        if ( (g+1)%DIV== 0 && (g+1)<=NBURN ) {
            for ( int p=0; p<NP; p++ ) {
                Ar_beta(p) = ((double) nA_beta(p)) / DIV;
                if ( Ar_beta(p) >= ROPT ) sigmap_beta(p) = sigmap_beta(p) * (2-(1-Ar_beta(p)) / (1-ROPT));
                else sigmap_beta(p) = sigmap_beta(p) / (2-Ar_beta(p) / ROPT);
                nA_beta(p) = 0.0; // We reinitialize the number of acceptance to zero
            }
        }
        /* After the burnin period */
        if ( (g+1) % DIV == 0 && (g+1) > NBURN ) {
            for (int p=0; p<NP; p++) {
                Ar_beta(p) = ((double) nA_beta(p)) / DIV;
                nA_beta(p) = 0.0; // We reinitialize the number of acceptance to zero
            }
        }


        //////////////////////////////////////////////////
        // Progress bar
        double Perc = 100 * (g+1) / (NGIBBS);
        if ( (g+1) % (NGIBBS/100) == 0 && verbose == 1) {
            Rprintf("*");
            R_FlushConsole();
            if( (g+1) % (NGIBBS/10) == 0 ) {
                double mAr_beta=0; // Mean acceptance rate
                for ( int p = 0; p < NP; p++ ) {
                    mAr_beta += Ar_beta(p) / NP;
                }
                Rprintf(":%.1f%%, mean accept. rates= beta:%.3f\n", Perc, mAr_beta);
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
    Rcpp::List z = Rcpp::List::create(Rcpp::Named("beta") = beta,
                                      Rcpp::Named("Deviance") = Deviance,
                                      Rcpp::Named("theta_latent") = theta_latent);
    
    return z;
    
} // end hSDM_binomial function

// Test
/*** R
# library(coda)
# 
# inv.logit <- function(x, min=0, max=1) {
#     p <- exp(x)/(1+exp(x))
#     p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
#     p * (max-min) + min
# }
# 
# nsite <- 200
# seed <- 1234
# set.seed(seed)
# visits<- rpois(nsite,3)
# visits[visits==0] <- 1
# 
# # Ecological process (suitability)
# x1 <- rnorm(nsite,0,1)
# x2 <- rnorm(nsite,0,1)
# X <- cbind(rep(1,nsite),x1,x2)
# beta.target <- c(-1,1,-1)
# logit.theta <- X %*% beta.target
# theta <- inv.logit(logit.theta)
# Y <- rbinom(nsite,visits,theta)
#     
# # Data-sets
# data.obs <- data.frame(Y,visits,x1,x2)
# 
# # Iterations
# nsamp <- 1000
# nburn <- 1000
# nthin <- 1
# ngibbs <- nsamp+nburn
# 
# mf.suit <- model.frame(formula=~x1+x2, data=data.obs)
# X <- model.matrix(attr(mf.suit,"terms"), data=mf.suit)
#     
# # Call to C++ function
# mod <- Rcpp_jSDM_binomial(
#     ngibbs=ngibbs, nthin=nthin, nburn=nburn,
#     Y=data.obs$Y,
#     T=data.obs$visits,
#     X=X,
#     beta_start=rep(0, 3),
#     mubeta=rep(0, 3), Vbeta=rep(1.0E6, 3),
#     seed=1234, ropt=0.44, verbose=1)
# 
# # Parameter estimates
# MCMC <- mcmc(mod$beta, start=nburn+1, end=ngibbs, thin=nthin)
# summary(MCMC)
# mean(mod$Deviance)
# plot(MCMC)
#            
# # GLM resolution to compare
# mod.glm <- glm(cbind(Y,visits-Y)~x1+x2,family="binomial",data=data.obs)
# summary(mod.glm)
#           
# # Predictions
# plot(theta,mod$theta_latent, main="theta",
#      xlab="obs", ylab="fitted")
# abline(a=0,b=1,col="red")
*/

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
