##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2011-present Ghislain Vieilledent
## 
##########################################################################

.onAttach <- function(...) {
   # echo output to screen
   packageStartupMessage("##\n## jSDM R package \n",
                         "## For joint species distribution models \n",
   					     "## https://ecology.ghislainv.fr/jSDM \n",
                         "##\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("jSDM", libpath)
}





