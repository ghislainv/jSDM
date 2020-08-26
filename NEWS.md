# jSDM 0.1.1

* Use of `roxygen2` for documentation and NAMESPACE 
* Rename  `jSDM_binomial_probit_block()` the function `jSDM_probit_block()`.
* New function `jSDM_poisson_log()` for fitting joint species distribution models from abundance data inspired by Hui and Francis K. C. 2016 _Methods in Ecology and Evolution_ (doi:10.1111/2041-210X.12514).
* New `jSDM_binomial_logit()` function for fitting joint species distribution models from presence-absence data at multiple-visited sites using a bayesian inference method inspired by Albert, James H. and Chib Siddhartha 1993 _Journal of the American Statistical Association_ (doi:10.1080/01621459.1993.10476321).
* Functions to fit models in which site effects are included as fixed effects, as random effects or not included and with or without latent variables.
* Add datasets (`mosquitos`, `fungi`, `eucalypts`, `birds`, `mites`, `aravo`). 
* Three new vignettes (“Bayesian inference methods”, "Poisson log-linear regression" "Bernoulli probit regression" and “Binomial logistic regression”) are available.
* Complete vignette "Comparison jSDM-boral". 
* New package website available on GitHub: <https://ecology.ghislainv.fr/jSDM>.

# jSDM 0.1.0

* First version of the jSDM R package
* Use of Rcpp and C++ code for Gibbs sampling
* Use of GSL (RcppGSL) for random draws
* Use of Armadillo (RcppArmadillo) for vector and matrix operations
* Functions to fit models from Warton et al. 2014 _Trends in Ecology and Evolution_ (doi:10.1016/j.tree.2015.09.007).
* We use `pkgdown` to build package website.
* Package website available on GitHub: <https://ecology.ghislainv.fr/jSDM>.
