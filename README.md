# intRinsic - public repository <img src="man/figures/intLogo.png" align="right" width="120" />

<!-- badges: start -->
<!--
[![R-CMD-check](https://github.com/Fradenti/intRinsic_dev/workflows/R-CMD-check/badge.svg)](https://github.com/Fradenti/intRinsic_dev/actions)
-->
<!-- badges: end -->

Provides functions to estimate the intrinsic dimension of a dataset via likelihood-based approaches. 
Specifically, the package implements the `TWO-NN` and `Gride` estimators and the `Hidalgo` Bayesian mixture model.

<!--
A `.pdf` vignette for this package can be found on Arxiv at [this link](https://arxiv.org/pdf/2102.11425.pdf).
-->
  
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
