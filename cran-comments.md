## Version 1.1.0

With this version, we:

* Re-organized the cpp code into multiple files
* Solved an indexing problem with `log_Zeta` and `log_corr`, sometimes causing `NA` with single-manifold data when using `Hidalgo()`
* Translated some `R` parts of the `Hidalgo()` code into `C++`
* Corrected a bug affecting the `Truncated-pointmass` approach
* Updated the README file
* Set up the new `rhub` checks workflows

## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()`, `devtools::check_win_release()`, and `devtools::check_mac_release()` does not produce any ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
