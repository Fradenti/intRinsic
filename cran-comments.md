## Version 1.0.1

* I am kindly asking for an early package update on CRAN since the 
current CRAN version of `intRinsic (v1.0.0)` contains a bug that makes the script published in the Journal of Statistical Software (`JSS`) not fully reproducible. Therefore, I would like to address this issue as soon as possible to avoid causing any disruption to potential users. In detail, we:
* Fixed a typo in `compute_mus` when computing ratios from `dist_matrix`. Now the function handles also the `dissimilarity` class
* Fixed a typo in `compute_mus` involving the class of the returned object
* Adjusted some indentations in `print` methods

  
## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()` and `devtools::check_win_release()` produces 
no ERRORs, WARNINGs, no NOTEs.

The function `devtools::check_rhub()` produces no ERRORs, no WARNINGs, and one NOTE. 
Specifically, only on `Fedora Linux, R-devel, clang, gfortran`, we obtained

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

Note that only on `Debian Linux, R-devel, GCC ASAN/UBSAN` I obtained a PREPERROR message. Inspection of the log reveals that the issue is that:

* ERROR: dependency ‘RcppArmadillo’ is not available for package ‘intRinsic’

If I instead run `rhub::check(platform = "debian-gcc-devel-nold")`, the package compiles with no issues.

Finally, this package, in its current state, also passes all the standard checks performed via GitHub actions.

## Downstream dependencies

There are currently no downstream dependencies for this package.
