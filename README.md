
<!-- README.md is generated from README.Rmd. Please edit that file -->

# intRinsic v1.1.1 <img src="man/figures/intLogo.png" align="right" width="120" />

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/intRinsic)](https://cran.r-project.org/package=intRinsic)
[![Last
Commit](https://img.shields.io/github/last-commit/fradenti/intRinsic)](https://github.com/fradenti/intRinsic)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/intRinsic?color=brightgreen)](https://www.r-pkg.org/pkg/intRinsic)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/intRinsic?color=brightgreen)](https://www.r-pkg.org/pkg/intRinsic)
[![JSS](https://img.shields.io/badge/JSS-10.18637%2Fjss.v040.i08-brightgreen)](https://www.jstatsoft.org/article/view/v106i09)
<!-- badges: end -->

A package with functions to estimate the intrinsic dimension of a
dataset via likelihood-based approaches. Specifically, the package
implements the `TWO-NN` and `Gride` estimators and the `Hidalgo`
Bayesian mixture model.

To install the package from CRAN, run

``` r
install.packages("intRinsic")
```

To install the package from this GitHub repository, run

``` r
# install.packages("remotes")

#Turn off warning-error-conversion regarding package versions
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

#install from github
remotes::install_github("Fradenti/intRinsic")
```

Simple example on Swissroll dataset

``` r
library(intRinsic)
X <- Swissroll(2000)
twonn(X)
```

The vignette for this package has been published in the
`Journal of Statistical Software`. The article can be found at [this
link](https://doi.org/10.18637/jss.v106.i09).

Please help me improve this package by reporting suggestions, typos, and
issues at [this link](https://github.com/Fradenti/intRinsic/issues).

Please note that the previous versions of the package (from `v0.1.0` to
`v1.0.2`) are still available as GitHub Releases at [this
page](https://github.com/Fradenti/intRinsic/releases).
