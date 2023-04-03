// ==============================================================================
// author          :Ghislain Vieilledent, Jeanne Clément
// email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <cmath>

/* dens_par.h */
struct dens_par {
  // Data 
  int NSITE;
  int NSP;
  arma::umat Y;
  arma::uvec N;
  // Suitability 
  // beta
  int NP;
  arma::mat X;
  int pos_beta;
  int sp_beta;
  arma::vec mu_beta;
  arma::vec V_beta;
  arma::mat beta_run;
  // lambda
  int NL; 
  int pos_lambda;
  int sp_lambda;
  arma::vec mu_lambda;
  arma::vec V_lambda;
  arma::mat lambda_run;
  // W
  int site_W;
  int pos_W;
  arma::vec V_W;
  arma::mat W_run;
  //alpha
  int site_alpha; 
  double V_alpha_run;
  double shape_Valpha;
  double rate_Valpha;
  arma::rowvec alpha_run;
};

// Prototype of useful functions
/* dens_logit */
double betadens_logit (double beta_jk, void *dens_data);
double lambdadens_logit (double lambda_jq, void *dens_data);
double lambdaUdens_logit (double lambda_jq, void *dens_data);
double alphadens_logit(double alpha_i, void *dens_data);
double Wdens_logit (double W_iq, void *dens_data);
double betadens_pois (double beta_jk, void *dens_data);
double lambdadens_pois (double lambda_jq, void *dens_data);
double lambdaUdens_pois (double lambda_jq, void *dens_data);
double alphadens_pois(double alpha_i, void *dens_data);
double Wdens_pois (double W_iq, void *dens_data);
arma::mat chol_decomp (arma::mat V);
int my_gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, 
                                      const gsl_matrix * L, gsl_vector * result);
arma::vec arma_mvgauss (const gsl_rng* r, const arma::vec mu, const arma::mat L);
double logit (double x);
double invlogit (double x);

//------------------------------------------------------------
// Pseudorandom numbers from a truncated Gaussian distribution
// The Gaussian has parameters mu (default 0) and sigma (default 1)
// and is truncated on the interval [a,b].
// Returns the random variable x and its probability p(x).
//  The code is based on the implementation by Guillaume Dollé, Vincent Mazet
//  available from http://miv.u-strasbg.fr/mazet/rtnorm/.
double rtnorm (gsl_rng *gen, double a,double b,const double mu,const double sigma);
double pdf_tnorm(double x, double mu, double sigma, double a, double b);
// End