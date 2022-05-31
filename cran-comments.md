## Version 0.2.1

In this release, 

* we have addressed some minor bugs
* we have added basic methods functions (e.g., `print()`, `plot()`,...) to ease 
  the user's access to the results provided by the main functions.
 
  
## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions
`devtools::check_win_devel()` and `devtools::check_win_release()` produce 
no ERRORs, WARNINGs, nor NOTEs.  

The functions `devtools::check_rhub()`,  produce no ERRORs, no WARNINGs, and one NOTE.  
Specifically, only on `Platform: Windows Server 2022, R-devel, 64 bit`, we obtained

* checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

Finally, this package, in its current state, also passes all the standard 
checks performed via GitHub action.

## Downstread dependencies

There are currently no downstream dependencies for this package.
