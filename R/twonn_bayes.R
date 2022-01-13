#' Bayesian estimator for the \code{TWO-NN} model
#'
#' The function fits the \code{TWO-NN} model via Bayesian estimation, employing
#' a conjugate prior. The formulas can be found in
#' \href{https://arxiv.org/abs/2104.13832}{Denti et al., 2021+}.
#'
#' @param mus vector of second to first NN distance ratios.
#' @param a_d shape parameter of the Gamma prior on the parameter \code{d}.
#' @param b_d rate parameter of the Gamma prior on the parameter \code{d}.
#' @param alpha posterior probability contained in the computed credible
#' interval.
#' @param c_trimmed proportion of trimmed observations.
#'
#' @return object of class \code{twonn_bayes}, which is a list containing the
#' {(1 + \code{alpha}) / 2} and 1 - \code{\alpha} / 2 quantiles, mean, mode and
#' median of the posterior distribution of \code{d}.
#'
#' @seealso \code{\link{twonn}}, \code{\link{autoplot.twonn_bayes}}
#'
#' @keywords Internal
#' @noRd
#'
#' @references
#' Denti F, Doimo D, Laio A, Mira A (2022+). "Distributional Results for
#' Model-Based Intrinsic Dimension Estimators."
#' arXiv preprint. 2104.13832, \url{https://arxiv.org/abs/2104.13832}.
#'
twonn_bayes <- function(mus,
                        a_d = 0.001,
                        b_d = 0.001,
                        alpha = 0.95,
                        c_trimmed = 0.01) {
  n <- n_original <- base::length(mus)
  c_considered <- 1 - c_trimmed

  n  <- base::length(mus)

  if (c_trimmed) {
    mus <- sort(mus)[1:floor(n * c_considered)]
    n <- base::length(mus)
  }

  logmus  <- base::log(mus)
  Res    <- numeric(5)
  Res[2] <- (a_d + length(logmus)) / (b_d + sum(logmus))
  Res[3] <-
    stats::qgamma(.5, shape = (a_d + length(logmus)),
                  rate = (b_d + sum(logmus)))
  Res[4] <- (a_d + length(logmus) - 1) / (b_d + sum(logmus))
  Res[c(1, 5)] <-
    stats::qgamma(c((1 - alpha) / 2, (1 + alpha) / 2),
                  shape = (a_d + length(logmus)),
                  rate = (b_d + sum(logmus)))
  Res <- list(
    est = Res,
    alpha     = alpha,
    hp_prior  = c(a_d, b_d),
    hp_posterior = c(a_d + length(logmus), b_d + sum(logmus)),
    c_trimmed = c_trimmed,
    n_original = n_original,
    n = n
  )

  structure(Res, class = c("twonn_bayes", class(Res)))
}

#' Print \code{TWO-NN} Bayes object
#'
#' @param x object of class \code{twonn_bayes}, obtained from the function
#' \code{twonn_bayes()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the Bayesian TWO-NN to console.
#'
#' @export
print.twonn_bayes <- function(x, ...) {
  cat("Model: TWO-NN\n")
  cat("Method: Bayesian Estimation\n")
  cat(paste0(
    "Sample size: ",
    x[["n_original"]],
    ", Obs. used: ",
    x[["n"]],
    ". Trimming proportion: ",
    100 * x[["c_trimmed"]],
    "%\n"
  ))
  cat(paste0("Prior d ~ Gamma(",
             x[["hp_prior"]][1],
             ", ",
             x[["hp_prior"]][2],
             ")\n"))
  cat(paste0(
    "Credibile Interval quantiles: ",
    (1 - x[["alpha"]]) / 2 * 100,
    "%, ",
    (1 + x[["alpha"]]) / 2 * 100,
    "%\n"
  ))
  cat(paste0("Posterior ID estimates:"))
  y <- cbind(
    `Lower Bound` = x[["est"]][1],
    `Mean` = x[["est"]][2],
    `Median` = x[["est"]][3],
    `Mode` = x[["est"]][4],
    `Upper Bound` = x[["est"]][5]
  )
  print(knitr::kable(y))
}
