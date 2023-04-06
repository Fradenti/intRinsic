## Version 1.0.0

* This release is accompanied by a publication in the `Journal of Statistical Software`.
  Please note that the DOI in the CITATION is for a new JSS publication that will be registered `after` the package release on CRAN.
* Also, we have fixed minor bugs and have updated the documentation adding the most recent references
  
## R CMD check results

Please note that the DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN. As a result, I currently obtain the following notes across the various platforms:

* Found the following (possibly) invalid URLs:
   URL: https://doi.org/10.18637/jss.v106.i09
     From: README.md
     Status: 404
     Message: Not Found

* Found the following (possibly) invalid DOIs:
   DOI: 10.18637/jss.v106.i09
     From: DESCRIPTION
           inst/CITATION
     Status: 404
     Message: Not Found

Aside from the notes above, running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()` and `devtools::check_win_release()` produces 
no ERRORs, WARNINGs, or additional NOTEs.  

I removed the recent note 

* checking C++ specification ...
  NOTE Specified C++11: please drop specification unless essential

following https://www.tidyverse.org/blog/2023/03/cran-checks-compiled-code/, updating the `src/Makevars` and `msrc/Makevars.win` files. 

The function `devtools::check_rhub()` produces no ERRORs, no WARNINGs, and one NOTE. Specifically, only on `Fedora Linux, R-devel, clang, gfortran`, we obtained

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

Note that only on `Debian Linux, R-devel, GCC ASAN/UBSAN` I obtained a PREPERROR message. Inspection of the log reveals that the issue is that:

* ERROR: dependency ‘RcppArmadillo’ is not available for package ‘intRinsic’

If I instead run `rhub::check(platform = "debian-gcc-devel-nold")`, the package compiles with no issues.

Finally, this package, in its current state, also passes all the standard checks performed via GitHub actions.

## Downstread dependencies

There are currently no downstream dependencies for this package.
