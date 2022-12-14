#' The Generalized Ratio distribution
#'
#' Density function and random number generator for the Generalized Ratio
#' distribution with NN orders equal to \code{n1} and \code{n2}.
#' See \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}
#' for more details.
#'
#' @aliases dgera rgera
#'
#' @name generalized_ratios_distribution
#'
#' @param x vector of quantiles.
#' @param nsim integer, the number of observations to generate.
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param d value of the intrinsic dimension.
#' @param log logical, if \code{TRUE}, it returns the log-density
#'
#' @return \code{dgera} gives the density. \code{rgera} returns a vector of
#' random observations sampled from the generalized ratio distribution.
#'
#' @references
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
#'
#' @examples
#' draws   <- rgera(100,3,5,2)
#' density <- dgera(3,3,5,2)
#'
NULL

#' @name generalized_ratios_distribution
#' @export
rgera <- function(nsim,
                  n1 = 1,
                  n2 = 2,
                  d) {
  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }
  a <- n2 - n1
  b <- n1
  x <- stats::rgamma(nsim, a, 1)
  y <- stats::rgamma(nsim, b, 1)
  z <- x / y
  sample <- (z + 1) ^ (1 / d)
  attr(sample, "n1") <- n1
  attr(sample, "n2") <- n2
  attr(sample, "upper_D") <- "unknown"
  structure(sample, class = c("mus", class(sample)))
}

#' @name generalized_ratios_distribution
#' @export
dgera <- function(x,
                  n1 = 1,
                  n2 = 2,
                  d,
                  log = FALSE) {

  x <- c(x)
  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }
  d_n12    <- (n2 - n1)

  if (d_n12 == 1L) {
    lognum   <- log(d)
    logden   <- ((n2 - 1) * d + 1) * (log(x))
    log_dens <- lognum - logden + (log(x > 1))
  } else{
    logB     <- sum(log(1:(n2 - n1 - 1))) + sum(log(1:(n1-1))) -
      sum(log(1:(n2-1)))
    lognum   <- log(d) + (d_n12 - 1) * (log(x ^ d - 1))
    logden   <- ((n2 - 1) * d + 1) * (log(x))
    log_dens <- lognum - logden + sum(log(x > 1)) - logB
  }
  if (!log) {
    log_dens <- exp(log_dens)
  }

  structure(log_dens, class = c("dgera", class(log_dens)))

  return(log_dens)
}
