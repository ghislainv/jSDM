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


// Prototype of useful functions
arma::mat chol_decomp (arma::mat V);
int my_gsl_ran_multivariate_gaussian (const gsl_rng * r, const gsl_vector * mu, 
                                      const gsl_matrix * L, gsl_vector * result);
arma::vec arma_mvgauss (const gsl_rng* r, const arma::vec mu, const arma::mat L);
double logit (double x);
double invlogit (double x);
//------------------------------------------------------------
// Compute y_l from y_k
double yl(int k);

//------------------------------------------------------------
// Rejection algorithm with a truncated exponential proposal
double c_rtexp(double a, double b, const gsl_rng* seed);

//------------------------------------------------------------
// One-sided rejection algorithm with exponential proposal
double c_rtexp_onesided(double a, const gsl_rng* seed);

//------------------------------------------------------------
// Chopin's two-sided ziggurat algorithm
double rtchopin_twosided(double a, double b, const gsl_rng* seed);

//------------------------------------------------------------
// Chopin's one-sided ziggurat algorithm
double rtchopin_onesided(double a, const gsl_rng* seed);

//------------------------------------------------------------
// Naive accept reject (two-sided)
double c_rtnaive(double a, double b, const gsl_rng* seed );

//------------------------------------------------------------
// Naive accept reject (one-sided)
double c_rtnaive_onesided(double a, const gsl_rng* seed );

//------------------------------------------------------------
// Pseudorandom numbers from a truncated Gaussian distribution
// The Gaussian has parameters mu (default 0) and sigma (default 1)
// and is truncated on the interval [a,b].
// Returns the random variable x and its probability p(x).
//' Pseudorandom numbers from a Gaussian distribution that is truncated to an interval.
double rtnorm(double lower, double upper, double mu, double sigma, const gsl_rng* seed );

// End