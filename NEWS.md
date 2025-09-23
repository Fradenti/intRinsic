# intRinsic 1.1.1

* Replaced `arma::is_finite` with `std::isfinite` to ensure compatibility with the latest `RcppArmadillo` release and CRAN policies.
* Implemented small code adjustments to address and resolve compiler warnings.

# intRinsic 1.1.0

* Re-organized the cpp code into multiple files
* Solved an indexing problem with `log_Zeta` and `log_corr`, sometimes causing `na` with single-manifold data when using `Hidalgo()`
* Translated some `R` parts of the `Hidalgo()` code into `C++`
* Corrected a bug affecting the `Truncated-pointmass` approach
* Updated the README file
* Set up the new `rhub` checks workflows

# intRinsic 1.0.2

* Removed the dependency from the package `MCMCpack` 
* Fixed and updated the documentation (minor changes)

# intRinsic 1.0.1

* Fixed a typo in `compute_mus` when computing ratios from `dist_matrix`. Now the function handles also the `dissimilarity` class
* Fixed a typo in `compute_mus` involving the class of the returned object
* Adjusted some indentations in `print` methods
* Adjusted the `autoplot` method related to `gride_mle`

# intRinsic 1.0.0

* This update marks the acceptance of the `intRinsic` vignette in the `Journal of Statistical Software`. Note: the DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN
* Fixed some spelling typos and new-line inconsistencies
* Removed the dependencies on the `as_tibble()` function (deprecated)
* The function `twonn_decimated()` is now deprecated and will be removed in future releases. The function to use is `twonn_decimation()`
* Now `gride_evolution()` and `twonn_decimation()` remove duplicated observations by default

# intRinsic 0.2.2

* Fixed bug in `Hidalgo()`, causing errors when setting `verbose = FALSE`
* Fixed bug in `Hidalgo()` initialization. Added warning if `alpha_Dirichlet` is set so small it causes underflow problems
* Fixed bug in `clustering()`. Now the returned `K` when the option `salso` is selected, is correct
* Fixed bug in `print.mus()`. Now it checks if the passed object is a list
* Corrected a few typos in the documentation
* Updated `doi` of recently published references
* Added additional checks in `compute_mus`: NA, symmetry, non-negativity

# intRinsic 0.2.1

* Fixed minor bugs involving the `gride()` family of functions
* Added several methods (`print()`, `summary()`, `plot()`,...) to facilitate the
  extraction of the results
* Now, the `autoplot()` method is directly exported from the `ggplot2` package

# intRinsic 0.2.0

* Submitted to CRAN
* Released to the public the first official version

