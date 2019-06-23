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
// Pseudorandom numbers from a truncated Gaussian distribution
// The Gaussian has parameters mu (default 0) and sigma (default 1)
// and is truncated on the interval [a,b].
// Returns the random variable x and its probability p(x).
//  The code is based on the implementation by Guillaume Dollé, Vincent Mazet
//  available from http://miv.u-strasbg.fr/mazet/rtnorm/.
double rtnorm (gsl_rng *gen, double a,double b,const double mu,const double sigma);

// End