
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jSDM R Package <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/ghislainv/jSDM/workflows/R-CMD-check/badge.svg)](https://github.com/ghislainv/jSDM/actions)
[![CRAN
Status](https://www.r-pkg.org/badges/version/jSDM)](https://cran.r-project.org/package=jSDM)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3253460.svg)](https://doi.org/10.5281/zenodo.3253460)
[![Downloads](https://cranlogs.r-pkg.org/badges/jSDM)](https://cran.r-project.org/package=jSDM)

Package for fitting joint species distribution models (jSDM) in a
hierarchical Bayesian framework (Warton *et al.*
[2015](#ref-Warton2015)). The Gibbs sampler is written in C++. It uses
Rcpp, Armadillo and GSL to maximize computation efficiency.

## System requirements

Make sure the GNU Scientific Library
([GSL](https://www.gnu.org/software/gsl/)) is installed on your system.

## Installation

### Stable version from CRAN

Install the latest stable version of **jSDM** from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("jSDM")
```

### Development version

Or install the development version of **jSDM** from
[GitHub](https://github.com/ghislainv/jSDM) with:

``` r
devtools::install_github("ghislainv/jSDM")
```

Or the binary release of **jSDM**’s development version compiled with R
version 4.0.3 can be found here
:

[jSDM\_windows](https://nextcloud.fraisedesbois.net/index.php/s/bEQNBdwe2RCSK9F).

<!-- ## Available functions  -->

<!-- The package includes the following functions to fit various species distribution models :  -->

<!-- |                function                |    data type     |                      fitted model                                 | -->

<!-- |:---------------------------------------|-----------------:|-----------------------------------------------:| -->

<!-- | jSDM_binomial_logit_one_species        | presence-absence |  $y_n \sim \mathcal{B}inomial(\theta_n, t_n)  where, / logit(\theta_n) = X_n \beta$ | -->

<!-- `jSDM_binomial_logit`                    | presence-absence |  $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ logit(\theta_{ij}) = X_i \beta_j$ | -->

<!-- | jSDM_binomial_logit_rand_site          | presence-absence |  $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \  logit(\theta_{ij}) = X_i \beta_j \ and \ \alpha_i \sim \mathcal{N}(0,V_{\alpha})$ | -->

<!-- | jSDM_binomial_logit_lv                 | presence-absence |   $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ logit(\theta_{ij}) = X_i \beta_j + W_i \lambda_j$ | -->

<!-- | jSDM_binomial_logit_rand_site_lv       | presence-absence |  $y_{ij} \sim \mathcal{B}inomial(\theta_ij, t_i) \ where, \ logit(\theta_{ij}) = X_i \beta_j +  W_i \lambda_j + \alpha_i \ and \  \alpha_i \sim \mathcal{N}(0,V_{\alpha})$ | -->

<!-- | jSDM_binomial_probit_block_one_species | presence-absence |  $y_n \sim \mathcal{B}inomial(\theta_n, t_n) \ where, \ probit(\theta_n) = X_n \beta$ | -->

<!-- | jSDM_binomial_probit_block             | presence-absence |   $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ probit(\theta_{ij}) = X_i \beta_j$ | -->

<!-- | jSDM_binomial_probit_block_rand_site   | presence-absence | $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ probit(\theta_{ij}) = X_i \beta_j + \alpha_i \ and \  \alpha_i \sim \mathcal{N}(0,V_{\alpha})$ | -->

<!-- | jSDM_binomial_probit_block_lv          | presence-absence |  $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ probit(\theta_{ij}) = X_i \beta_j + W_i \lambda_j$ | -->

<!-- | `jSDM_binomial_probit_block_rand_site_lv` | presence-absence |  $y_{ij} \sim \mathcal{B}inomial(\theta_{ij}, t_i) \ where, \ probit(\theta_{ij}) = X_i \beta_j +  W_i \lambda_j + \alpha_i \ and \  \alpha_i \sim \mathcal{N}(0,V_{\alpha})$ | -->

## References

<div id="refs" class="references">

<div id="ref-Warton2015">

Warton, D.I., Blanchet, F.G., O’Hara, R.B., Ovaskainen, O., Taskinen,
S., Walker, S.C. & Hui, F.K. (2015) So many variables: Joint modeling in
community ecology. *Trends in Ecology & Evolution*, **30**, 766–779.

</div>

</div>
