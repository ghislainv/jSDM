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