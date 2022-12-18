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

