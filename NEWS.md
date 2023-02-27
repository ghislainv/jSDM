# jSDM 0.2.2
* New function [`jSDM_gaussian`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_gaussian.html) for fitting joint species distribution models from continuous Gaussian data, including an overdispersion parameter.
* Use of the R package `terra` instead of `raster` and `sp`, which are depreciated, in the vignette [Estimation of Madagascar's plant biodiversity](https://ecology.ghislainv.fr/jSDM/articles/Madagascar.html). 

# jSDM 0.2.1
* Add the possibility to consider only significant correlations in the [`get_residual_cor`](https://ecology.ghislainv.fr/jSDM/reference/get_residual_cor.html) and [`plot_residual_cor`](https://ecology.ghislainv.fr/jSDM/reference/plot_residual_cor.html) functions. 
* Documentation correction

# jSDM 0.2.0
* New function [`jSDM_binomial_probit_sp_constrained`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_probit_sp_constrained.html) which aims to improve the convergence of latent variable models fitting by selecting the species constrained to have positive values of factor loadings $\lambda$ and new vignette [Bernoulli probit regression with selected constrained species](https://ecology.ghislainv.fr/jSDM/articles/jSDM_binomial_probit_sp_constrained.html) to illustrate its use. 
* New vignette [Estimation of Madagascar's plant biodiversity](https://ecology.ghislainv.fr/jSDM/articles/Madagascar.html) available. 

# jSDM 0.1.2

* Add the possibility of considering an additional hierarchical level in the Bayesian models of the [`jSDM_binomial_probit`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_probit.html), [`jSDM_binomial_logit`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_logit.html) and [`jSDM_poisson_log`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_poisson_log.html) functions to take into account interactions between species-specific traits and the environment in estimating species effects.
* New vignette [Bernoulli probit regression including species traits](https://ecology.ghislainv.fr/jSDM/articles/jSDM_with_traits.html) available.
* Separate the drawing of species effects beta and factor loading lambda in the functions `jSDM_binomial_probit_block` and `jSDM_binomial_probit_block_long_format` renamed [`jSDM_binomial_probit`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_probit.html) and [`jSDM_binomial_probit_long_format`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_probit_long_format.html).
* New function [`plot_associations`](https://ecology.ghislainv.fr/jSDM/reference/plot_associations.html) to plot species-species associations.

# jSDM 0.1.1
* Use of `roxygen2` for documentation and NAMESPACE 
* Rename  `jSDM_binomial_probit_block` the function `jSDM_probit_block`.
* New function `jSDM_binomial_probit_block_long_format` for fitting joint species distribution models from presence-absence data in long format able to handle missing observations, multiple visits at sites and to integer species traits as explanatory variables.  
* New function [`jSDM_poisson_log`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_poisson_log.html) for fitting joint species distribution models from abundance data inspired by Hui and Francis K. C. 2016 _Methods in Ecology and Evolution_ ([doi:10.1111/2041-210X.12514](https://doi.org/10.1111/2041-210X.12514)).
* New function [`jSDM_binomial_logit`](https://ecology.ghislainv.fr/jSDM/reference/jSDM_binomial_logit.html) for fitting joint species distribution models from presence-absence data at multiple-visited sites using a bayesian inference method inspired by Albert, James H. and Chib Siddhartha 1993 _Journal of the American Statistical Association_ ([doi:10.1080/01621459.1993.10476321](https://doi.org/10.1080/01621459.1993.10476321)).
* Functions to fit models in which site effects are included as fixed effects, as random effects or not included and with or without latent variables.
* New function [`get_enviro_cor`](https://ecology.ghislainv.fr/jSDM/reference/get_enviro_cor.html) to extract covariances and correlations due to shared environmental responses. 
* Complete [`jSDM-package`](https://ecology.ghislainv.fr/jSDM/reference/jSDM-package.html) documentation 
* Add datasets ([`mosquitos`](https://ecology.ghislainv.fr/jSDM/reference/mosquitos.html), [`fungi`](https://ecology.ghislainv.fr/jSDM/reference/fungi.html), [`eucalypts`](https://ecology.ghislainv.fr/jSDM/reference/eucalypts.html), [`birds`](https://ecology.ghislainv.fr/jSDM/reference/birds.html), [`mites`](https://ecology.ghislainv.fr/jSDM/reference/mites.html), [`aravo`](https://ecology.ghislainv.fr/jSDM/reference/aravo.html)). 
* Seven new vignettes ([Bayesian inference methods](https://ecology.ghislainv.fr/jSDM/articles/proof.html), [Poisson log-linear regression](https://ecology.ghislainv.fr/jSDM/articles/jSDM_poisson_log.html), [Bernoulli probit regression](https://ecology.ghislainv.fr/jSDM/articles/jSDM_binomial_probit.html), [Bernoulli probit regression with missing data and species traits](https://ecology.ghislainv.fr/jSDM/articles/jSDM_binomial_probit_long_format.html), [Binomial logistic regression](https://ecology.ghislainv.fr/jSDM/articles/jSDM_binomial_logit.html), [Running jSDM in parallel](https://ecology.ghislainv.fr/jSDM/articles/jSDM_in_parallel.html), [Comparing SDMs and JSDMs](https://ecology.ghislainv.fr/jSDM/articles/SDM_JSDM.html) and [Comparison jSDM-Hmsc](https://ecology.ghislainv.fr/jSDM/articles/jSDM_Hmsc.html)) are available.
* Complete and correct the vignette [Comparison jSDM-boral](https://ecology.ghislainv.fr/jSDM/articles/jSDM_boral.html). 
* Add Code of conduct and Contributing section  
* New package website available on GitHub: <https://ecology.ghislainv.fr/jSDM/>.

# jSDM 0.1.0

* First version of the jSDM R package
* Use of Rcpp and C++ code for Gibbs sampling
* Use of GSL (RcppGSL) for random draws
* Use of Armadillo (RcppArmadillo) for vector and matrix operations
* Functions to fit models from Warton et al. 2014 _Trends in Ecology and Evolution_ ([doi:10.1016/j.tree.2015.09.007](https://doi.org/10.1016/j.tree.2015.09.007)).
* We use `pkgdown` to build package website.
* Package website available on GitHub: <https://ecology.ghislainv.fr/jSDM/>.
