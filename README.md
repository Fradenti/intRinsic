# intRinsic <img src="man/figures/intLogo.png" align="right" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Fradenti/intRinsic_dev/workflows/R-CMD-check/badge.svg)](https://github.com/Fradenti/intRinsic_dev/actions)
[![CRAN](https://www.r-pkg.org/badges/version/intRinsic)](https://cran.r-project.org/package=intRinsic)
[![Last Commit](https://img.shields.io/github/last-commit/fradenti/intRinsic)](https://github.com/fradenti/intRinsic)
[![Downloads (monthly)](https://cranlogs.r-pkg.org/badges/intRinsic?color=brightgreen)](https://www.r-pkg.org/pkg/intRinsic)
[![Downloads (total)](https://cranlogs.r-pkg.org/badges/grand-total/intRinsic?color=brightgreen)](https://www.r-pkg.org/pkg/intRinsic)
<!-- [![JSS](https://img.shields.io/badge/JSS-10.18637%2Fjss.v040.i08-brightgreen)]() -->
<!-- badges: end -->

A package with functions to estimate the intrinsic dimension of a dataset via likelihood-based approaches. 
Specifically, the package implements the `TWO-NN` and `Gride` estimators and the `Hidalgo` Bayesian mixture model.

To install the package from CRAN, run
```r
install.packages("intRinsic")
```

To install the package from this GitHub repository, run
```r
# install.packages("remotes")

#Turn off warning-error-conversion regarding package versions
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

#install from github
remotes::install_github("Fradenti/intRinsic")
```

Simple example on Swissroll dataset

```r
library(intRinsic)
X <- Swissroll(2000)
twonn(X)
```

The updated vignette for this package is available on ArXiv, at [this link](https://arxiv.org/pdf/2102.11425.pdf).

Please note that the previous versions of the package (v0.1.0 and v0.2.0) are still available under the GitHub release at [this page](https://github.com/Fradenti/intRinsic/releases).
