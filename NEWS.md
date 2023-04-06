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

