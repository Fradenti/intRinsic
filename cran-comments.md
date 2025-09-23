## Version 1.1.1

With this version, we:

- Replaced `arma::is_finite` with `std::isfinite` to ensure compatibility with the latest `RcppArmadillo` release and CRAN policies.
- Implemented small code adjustments to address and resolve compiler warnings.

## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()`, `devtools::check_win_release()`, and `devtools::check_mac_release()` does not produce any ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
