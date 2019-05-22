// ==============================================================================
// author          :Ghislain Vieilledent, Jeanne Clément
// email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

#include <RcppArmadillo.h>
#include <gsl/gsl_rng.h>

// Prototype of useful functions
arma::mat chol_decomp (arma::mat V);
arma::vec arma_mvgauss (const gsl_rng *r, const arma::vec mu, const arma::mat L);
double left_truncated_normal_sample (double mu, double sigma, double a, gsl_rng* &seed);
double right_truncated_normal_sample (double mu, double sigma, double b, gsl_rng* &seed);
double logit (double x);
double invlogit (double x);

// End