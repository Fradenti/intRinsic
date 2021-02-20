#' Maximum Likelihood Estimator for the TWO-NN model
#'
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param unbiased logical, if \code{TRUE} the MLE is corrected to ensure unbiasedness.
#' @param conf_lev the confidence level for the computation of the confidence interval.
#' @param trimmed logical, if \code{TRUE}, the trimming of the most extreme observations is performed.
#' @param alpha_trimmed proportion of trimmed observations.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return The MLE for the intrinsic dimension of the dataset \code{X} with its confidence interval.
#' @export
mle_twonn <-
  function(X = NULL,
           unbiased = T,
           conf_lev = .95,
           trimmed = F,
           alpha_trimmed = 0.01,
           DM = NULL) {
    mu <- generate_mus(X, dupl_remove = T, DM = DM)
    n <- base::length(mu)
    alpha_considered <- 1 - alpha_trimmed
    if (trimmed) {
      mu <- sort(mu)[1:floor(n * alpha_considered)]
      n <- base::length(mu)
    }
    alpha <- 1 - conf_lev
    logmu <- base::log(mu)
    if (unbiased) {
      Z <- (n - 1) / base::sum(logmu)
      UB <-
        Z * stats::qgamma(1 - alpha / 2, shape = n, rate = n - 1)
      LB <- Z * stats::qgamma(alpha / 2, shape = n, rate = n - 1)
      Res <-
        cbind(
          `Lower Bound` = LB,
          `Estimate` = Z,
          `Upper Bound` = UB
        )
    } else {
      Z <- (n) / base::sum(logmu)
      UB <- Z * stats::qgamma(1 - alpha / 2, shape = n, rate = n)
      LB <- Z * stats::qgamma(alpha / 2, shape = n, rate = n)
      Res <-
        cbind(
          `Lower Bound` = LB,
          `Estimate` = Z,
          `Upper Bound` = UB
        )
    }
    return(Res)
  }



#' Least Square Estimator for the Intrinsic Dimension
#'
#'
#' @param X  a n \eqn{times} D data matrix.
#' @param trimmed logical, if \code{TRUE}, the trimming of the most extreme observations is performed.
#' @param alpha_trimmed proportion of trimmed observations.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return the estimate for the ID and, if requested, the plot depicting the linear fit.
#'
#' @export
linfit_twonn <-
  function(X,
           trimmed = F,
           alpha_trimmed = .01,
           DM = NULL) {
    alpha_considered <- 1 - alpha_trimmed
    mu <- generate_mus(X, DM = DM)
    n <- base::length(mu)
    if (trimmed) {
      mu <- base::sort(mu)[1:floor(n * alpha_considered)]
      n <- base::length(mu)
    }
    F_mui <- (0:(n - 1)) / n
    y <- -base::log(1 - (F_mui))
    x <- base::sort(base::log(mu))
    modlin1 <- stats::lm(y ~ x - 1)
    Res <- numeric(3)
    Res[c(1, 3)] <- stats::confint.lm(modlin1)
    Res[2] <- stats::coef(modlin1)
    names(Res) <- c("Lower Bound", "Estimate", "Upper Bound")
    Res <- list(Estimates = Res,
                Lin_mod = modlin1)

    class(Res) <- "lnf_twonn"
    return(Res)
  }




#' Bayesian estimator for the Instrisic Dimension
#'
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param a_d shape parameter of the Gamma prior on \eqn{d}.
#' @param b_d rate parameter of the Gamma prior on \eqn{d}.
#' @param alpha the order of the quantile for the credible set.
#' @param trimmed logical, if \code{TRUE}, the trimming of the most extreme observations is performed.
#' @param alpha_trimmed proportion of trimmed observations.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return quantiles, mean, mode and median of the posterior distribution on \eqn{d}.
#'
#' @export
#'
bayesfit_twonn <- function(X,
                           a_d = .001,
                           b_d = .001,
                           alpha = .95,
                           trimmed = F,
                           DM = NULL,
                           alpha_trimmed = .01) {
  alpha_considered <- 1 - alpha_trimmed


  mu <- generate_mus(X, DM = DM)
  n <- base::length(mu)
  if (trimmed) {
    mu <- sort(mu)[1:floor(n * alpha_considered)]
    n <- base::length(mu)
  }
  logmu <- base::log(mu)
  Res <- numeric(5)
  Res[2] <- (a_d + length(logmu)) / (b_d + sum(logmu))
  Res[3] <-
    stats::qgamma(.5, shape = (a_d + length(logmu)), rate = (b_d + sum(logmu)))
  Res[4] <- (a_d + length(logmu) - 1) / (b_d + sum(logmu))
  Res[c(1, 5)] <-
    stats::qgamma(c((1 - alpha) / 2, (1 + alpha) / 2),
                  shape = (a_d + length(logmu)),
                  rate = (b_d + sum(logmu)))
  names(Res) <-
    c("Lower Bound", "Mean", "Median", "Mode", "Upper Bound")
  Res <- list(
    Estimates = Res,
    hp_prior = c(a_d, b_d),
    hp_posterior = c(a_d + length(logmu), b_d + sum(logmu))
  )
  class(Res) <- "byf_twonn"
  return(Res)
}


