## Version 1.1.2

With this version, we:

- removed the strong dependency of `intRinsic` from the `latex2exp` package, which has been scheduled for archival.

## R CMD check results

Running `devtools::check(args = c('--as-cran','--no-manual'))` locally and the functions `devtools::check_win_devel()`, 
`devtools::check_win_release()`, and `devtools::check_mac_release()` does not produce any ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
