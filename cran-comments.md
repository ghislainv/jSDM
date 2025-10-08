# jSDM v0.2.7

Jeanne Clément has notified the CRAN team that the package has a new
maintainer (email to CRAN-submissions@R-project.org on 27/07/2023).

The new package maintainer is Ghislain Vieilledent
<ghislain.vieilledent@cirad.fr>.

── R CMD check results ──────────────────────────────────────────────── jSDM 0.2.7 ────
Duration: 13m 26s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>’
  
  New maintainer:
    Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
  Old maintainer(s):
    Jeanne Clément <jeanne.clement16@laposte.net>

❯ checking installed package size ... NOTE
    installed size is 60.3Mb
    sub-directories of 1Mb or more:
      doc    2.4Mb
      libs  57.2Mb

❯ checking examples ... [40s/17s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                                    user system elapsed
  jSDM_binomial_probit             3.772  6.252   2.489
  jSDM_gaussian                    2.559  5.291   1.578
  jSDM_binomial_probit_long_format 4.034  2.815   3.577

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

# jSDM v0.2.6

## clang-UBSAN

In jSDM_poisson_log tests, Rcpp_jSDM_poisson_log_traits_rand_site_lv :
runtime error: 3.11579e+10 is outside the range of representable
values of type 'unsigned int'.

I have changed the type of the count data (Y) from unsigned int to
unsigned long int in all C++ functions called by the jSDM_poisson_log
function, in order to avoid this range problem and I've put back all
the jSDM_poisson_log function tests that were deleted for the previous
submission.

I have updated the date and removed commented lines of codes in the
examples.

Using [R-hub builder](https://builder.r-hub.io/) to build on:


## Platform: Windows Server 2022, R-devel, 64 bit

❯ checking CRAN incoming feasibility ... [26s] NOTE
  Maintainer: 'Jeanne Clément <jeanne.clement16@laposte.net>'
  
  New submission
  
  Package was archived on CRAN
  
  Possibly misspelled words in DESCRIPTION:
    Warton (25:43)
    al (25:54)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2023-03-17 as issues were not corrected
      in time.
  
  Uses the superseded package: 'snow'

❯ checking installed package size ... NOTE
    installed size is  5.2Mb
    sub-directories of 1Mb or more:
      libs   3.7Mb

❯ checking HTML version of manual ... NOTE
  

❯ checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 5 notes ✖


##  Platform: Ubuntu Linux 20.04.1 LTS, R-release, GCC

❯ checking CRAN incoming feasibility ... [7s/40s] NOTE
  Maintainer: ‘Jeanne Clément <jeanne.clement16@laposte.net>’
  
  New submission
  
  Package was archived on CRAN
  
  Possibly misspelled words in DESCRIPTION:
    al (25:54)
    Warton (25:43)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2023-03-17 as issues were not corrected
      in time.
  
  Uses the superseded package: ‘snow’

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Platform : Fedora Linux, R-devel, clang, gfortran

❯ checking package dependencies ... NOTE
  Package suggested but not available for checking: ‘boral’

❯ checking installed package size ... NOTE
    installed size is 38.3Mb
    sub-directories of 1Mb or more:
      libs  36.8Mb

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

## Platform : Fedora Linux, R-devel, clang, gfortran

❯ checking CRAN incoming feasibility ... [7s/40s] NOTE
  Maintainer: ‘Jeanne Clément <jeanne.clement16@laposte.net>’
  
  New submission
  
  Package was archived on CRAN
  
  Possibly misspelled words in DESCRIPTION:
    al (25:54)
    Warton (25:43)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2023-03-17 as issues were not corrected
      in time.
  
  Uses the superseded package: ‘snow’

❯ checking installed package size ... NOTE
    installed size is  58.3Mb
    sub-directories of 1Mb or more:
      libs   58.3Mb

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

# jSDM v0.2.5

## clang-UBSAN

In test jSDM_poisson_log, Rcpp_jSDM_poisson_log_traits_lv : 
runtime error: 2.86885e+10 is outside the range of representable values of type 'unsigned int'. 

I have fixed this issues. 

I have updated the date. 

Always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'C++'

The LICENSE file is only needed if you have additional restrictions to the license which you have not? In that case omit the file and its reference in the DESCRIPTION file.

Write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names. 

Add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)

You write information messages to the console that cannot be easily suppressed. It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object. 
Instead of print()/cat() rather use message()/warning() or if(verbose)cat(..) (or maybe stop()) if you really have to write text to the console. (except for print, summary, interactive functions)

Make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.

Always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos

I have fixed this issues.

# jSDM v0.2.3

## R CMD check results

Using [R-hub builder](https://builder.r-hub.io/) to build on :

* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit   
* Windows Server 2008 R2 SP1, R-release, 32/64 bit   
* Windows Server 2022, R-devel, 64 bit   
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC  
* Ubuntu Linux 20.04.1 LTS, R-release, GCC  
  
There were no ERRORs or WARNINGs.   

There was sometimes 1 NOTE:

* checking installed package size ... NOTE  
  installed size is  5.2Mb
  sub-directories of 1Mb or more:
    libs   3.8Mb

## clang-UBSAN

In test jSDM_poisson_log, Rcpp_jSDM_poisson_log_traits_lv : 
runtime error: 2.81295e+15 is outside the range of representable values of type 'unsigned int'. 

I have fixed this issues. 

I have also removed progress bars in non-interactive use, thank you for the advice.

# jSDM v0.2.2

## R CMD check results

Using [R-hub builder](https://builder.r-hub.io/) to build on :

* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit   
* Windows Server 2008 R2 SP1, R-release, 32/64 bit   
* Windows Server 2022, R-devel, 64 bit   
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC  
* Ubuntu Linux 20.04.1 LTS, R-release, GCC  
  
There were no ERRORs or WARNINGs.   

There was sometimes 1 NOTE:

* checking installed package size ... NOTE  
  installed size is  5.2Mb
  sub-directories of 1Mb or more:
    libs   3.8Mb

# jSDM v0.2.1

## R CMD check results

Using [R-hub builder](https://builder.r-hub.io/) to build on :

* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit   
* Windows Server 2008 R2 SP1, R-release, 32/64 bit   
* Windows Server 2022, R-devel, 64 bit   
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC  
* Ubuntu Linux 20.04.1 LTS, R-release, GCC  
  
There were no ERRORs or WARNINGs.   

There was sometimes 1 NOTE:

* checking installed package size ... NOTE  
  installed size is 4.5Mb  
  sub-directories of 1Mb or more:  
    - libs   3.2Mb 
    
 
## clang-UBSAN

jSDM-Ex.Rout:Rcpp_jSDM_poisson_log_rand_site_lv.cpp:78:17: runtime
error: nan is outside the range of representable values of type
'unsigned int'

and lots more.

I have fixed this issues. 
    
# jSDM v0.2.0

## R CMD check results

Using [R-hub builder](https://builder.r-hub.io/) to build on :

* Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit   
* Windows Server 2008 R2 SP1, R-release, 32/64 bit   
* Windows Server 2022, R-devel, 64 bit   
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC  
* Ubuntu Linux 20.04.1 LTS, R-release, GCC  
  
There were no ERRORs or WARNINGs.   

There was sometimes 2 NOTES:

* checking installed package size ... NOTE  
  installed size is 4.5Mb  
  sub-directories of 1Mb or more:  
    - libs   3.2Mb 

* checking CRAN incoming feasibility ... NOTE
  Maintainer: '
  Jeanne Clément <jeanne.clement16@laposte.net>'
  New maintainer:
  Jeanne Clément <jeanne.clement16@laposte.net>
  Old maintainer(s):
  Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>

Ghislain Vieilledent wants to transfer the maintainership of the package to Jeanne Clément. 
