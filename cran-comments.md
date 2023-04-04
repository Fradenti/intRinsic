## Version 1.0.0

* This release is accompanied by a publication on the `Journal of Statistical Software` 
  Please note that the DOI in the CITATION is for a new JSS publication that will be registered `after` the package release on CRAN.
* Also, we have fixed minor bugs and have updated the documentations adding the most recent references
  
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
no ERRORs, WARNINGs, nor NOTEs.  

The functions `devtools::check_rhub()`,  produce no ERRORs, no WARNINGs, and one NOTE. Specifically, only on `Fedora Linux, R-devel, clang, gfortran`, we obtained

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

Finally, this package, in its current state, also passes all the standard 
checks performed via GitHub actions.

## Downstread dependencies

There are currently no downstream dependencies for this package.
