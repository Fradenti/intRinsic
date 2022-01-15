# intRinsic <img src="man/figures/intLogo.png" align="right" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Fradenti/intRinsic/workflows/R-CMD-check/badge.svg)](https://github.com/Fradenti/intRinsic/actions)
<!-- badges: end -->

A package with functions to estimate the intrinsic dimension of a dataset via likelihood-based approaches. 
Specifically, the package implements the `TWO-NN` and `Gride` estimators and the `Hidalgo` Bayesian mixture model.

<!--
A `.pdf` vignette for this package can be found on Arxiv at [this link](https://arxiv.org/pdf/2102.11425.pdf).
-->

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
remotes::install_github("Fradenti/intRinsic_dev")
```

Simple example on Swissroll dataset

```r
library(intRinsic)
X <- Swissroll(2000)
twonn(X)
```


Please note that the original version of the package (v0.1.0) is still available under the GitHub release with tag `v0.1.0` at [this page](https://github.com/Fradenti/intRinsic/releases).
