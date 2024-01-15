## Version 1.0.2

* We removed the strong dependency of `intRinsic` from the `MCMCpack` package, which has been scheduled for archival.
* We fixed minor typos in the documentation.
  
## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()` and `devtools::check_win_release()` does not produce 
any ERRORs, WARNINGs, or NOTEs.

The function `devtools::check_rhub()` produces no ERRORs, no WARNINGs, and the following NOTEs. 

- On `Windows Server 2022, R-devel, 64 bit`, we obtained

* checking for non-standard things in the check directory ... NOTE
Found the following files/directories: ''NULL''
* checking for detritus in the temp directory ... NOTE
Found the following files/directories: 'lastMiKTeXException'

- On `Fedora Linux, R-devel, clang, gfortran`, we obtained

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found


- On `Ubuntu Linux 20.04.1 LTS, R-release, GCC`, we obtained

* checking installed package size ... NOTE
  installed size is  5.7Mb
  sub-directories of 1Mb or more:
    libs   5.2Mb

Finally, this package, in its current state, also passes all the standard checks performed via GitHub actions.

## Downstream dependencies

There are currently no downstream dependencies for this package.
