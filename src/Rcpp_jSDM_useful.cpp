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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

extern "C" {
#include "mvgauss.c"
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

using namespace arma;
using namespace std;

//* ************************************************************ */
/* chol_decomp */
arma::mat chol_decomp(arma::mat V) {
	// Cholesky decomposition
	arma::mat L = arma::chol(V, "lower"); 
	return L;
} 

/* ************************************************************ */
/* arma_mvgauss */
arma::vec arma_mvgauss(const gsl_rng *r, const arma::vec mu,
                       const arma::mat L) {
	
	// gsl vector mu
	gsl_vector *gsl_mu = gsl_vector_alloc(mu.n_elem);
	int size_mu = mu.n_elem;
	for (int i=0; i < size_mu; i++) {
		gsl_vector_set(gsl_mu, i, mu(i));
	}
	
	// gsl matrix L
	gsl_matrix *gsl_L = gsl_matrix_alloc(L.n_rows, L.n_cols);
	int nrows_L =  L.n_rows;
	int ncols_L =  L.n_cols;
	for (int i=0; i < nrows_L; i++) {
		for (int j=0; j < ncols_L; j++) {
			gsl_matrix_set(gsl_L, i, j, L(i,j));
		}
	}
	
	// gsl vector R
	gsl_vector *gsl_R = gsl_vector_alloc(mu.n_elem);
	gsl_vector_set_zero(gsl_R);
	
	// Call to gsl_ran_multivariate_gaussian
	gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_L, gsl_R);
	
	// arma vec R
	arma::vec R; R.zeros(gsl_R->size);
	int size_R = gsl_R->size;
	for (int i=0; i < size_R; i++) {
		R(i) = gsl_vector_get(gsl_R, i);
	} 
	
	// free the memory
	gsl_vector_free(gsl_mu);
	gsl_matrix_free(gsl_L);
	gsl_vector_free(gsl_R);
	
	// Return result
	return R;
}

/* ************************************************************ */
/* left_truncated_normal*/
double left_truncated_normal_sample ( double mu, double sigma, double a, 
                                      gsl_rng* &seed )
                                      
{
	double alpha;
	double alpha_cdf;
	double u;
	double x;
	double xi;
	double xi_cdf;
	
	alpha = ( a - mu ) / sigma;
	
	alpha_cdf = gsl_cdf_ugaussian_P(alpha);
	u = gsl_rng_uniform ( seed );
	xi_cdf = alpha_cdf + u * ( 1.0 - alpha_cdf );
	xi = gsl_cdf_ugaussian_Pinv ( xi_cdf);
	
	x = mu + sigma * xi;
	
	return x;
}

/* ************************************************************ */
/* right_truncated_normal*/

double right_truncated_normal_sample ( double mu, double sigma, double b, gsl_rng* &seed )
{
	double beta;
	double beta_cdf;
	double u;
	double x;
	double xi;
	double xi_cdf;
	
	beta = ( b - mu ) / sigma;
	
	beta_cdf = gsl_cdf_ugaussian_P ( beta );
	
	u = gsl_rng_uniform( seed );
	xi_cdf = u * beta_cdf;
	xi = gsl_cdf_ugaussian_Pinv ( xi_cdf );
	
	x = mu + sigma * xi;
	
	return x;
}

/******************/
/* Function logit */
double logit (double x) {
	return std::log(x) - std::log(1-x);
}

/*********************/
/* Function invlogit */
double invlogit (double x) {
	if (x > 0) {
		return 1 / (1 + std::exp(-x));
	}
	else {
		return std::exp(x) / (1 + std::exp(x));
	}
}

// End