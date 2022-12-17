## Version 0.2.2

In this release, 

* we have addressed some minor bugs (see `NEWS.md` for more details)
* we have improved the documentation by correcting some typos
* we have added some additional input checks 
* we have updated the documentations adding the most recent references
  
## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions
`devtools::check_win_devel()` and `devtools::check_win_release()` produces 
no ERRORs, WARNINGs, nor NOTEs.  

The functions `devtools::check_rhub()`,  produce no ERRORs, no WARNINGs, and one NOTE. Specifically, only on `Fedora Linux, R-devel, clang, gfortran`, we obtained

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found


Finally, this package, in its current state, also passes all the standard 
checks performed via GitHub actions.

## Downstread dependencies

There are currently no downstream dependencies for this package.
