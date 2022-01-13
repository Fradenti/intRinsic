#' Maximum Likelihood Estimator for the \code{TWO-NN} model
#'
#' The function fits the \code{TWO-NN} model via maximum likelihood estimation,
#' as discussed in
#' \href{https://arxiv.org/abs/2104.13832}{Denti et al., 2021+}.
#' This function is not exported, and can be accessed with \code{twonn()},
#' specifying \code{method = "mle"}.
#'
#' @param mus vector of second to first NN distance ratios.
#' @param unbiased logical, if \code{TRUE} the MLE is corrected to ensure
#' unbiasedness.
#' @param alpha confidence level needed for the computation of the confidence
#' interval.
#' @param c_trimmed proportion of trimmed observations.
#'
#' @return list of class \code{twonn_mle} containing the maximum likelihood
#' estimation of the intrinsic dimension with its confidence interval.
#' @keywords Internal
#'
#' @noRd
#'
#' @seealso \code{\link{twonn}}
#'
#' @references
#' #' Facco E, D'Errico M, Rodriguez A, Laio A (2017). "Estimating the intrinsic
#' dimension of datasets by a minimal neighborhood information."
#' Scientific Reports, 7(1), 1-8.
#' ISSN 20452322, doi: 10.1038/s41598-017-11873-y.
#'
#' Denti F, Doimo D, Laio A, Mira A (2022+). "Distributional Results for
#' Model-Based Intrinsic Dimension Estimators."
#' arXiv preprint. 2104.13832, \url{https://arxiv.org/abs/2104.13832}.
#'
twonn_mle <-
  function(mus,
           unbiased = TRUE,
           alpha = .95,
           c_trimmed = 0.01) {
    n <- n_original <- base::length(mus)
    c_considered <- 1 - c_trimmed

    if (c_trimmed > 0) {
      mus <- sort(mus)[1:floor(n * c_considered)]
      n  <- base::length(mus)
    }

    alpha  <- 1 - alpha
    logmus <- base::log(mus)

    if (unbiased) {
      Z  <- (n - 1) / base::sum(logmus)
      UB <-
        Z * stats::qgamma(1 - alpha / 2, shape = n, rate = n - 1)
      LB <- Z * stats::qgamma(alpha / 2, shape = n, rate = n - 1)
    } else {
      Z <- (n) / base::sum(logmus)
      UB <- Z * stats::qgamma(1 - alpha / 2, shape = n, rate = n)
      LB <- Z * stats::qgamma(alpha / 2, shape = n, rate = n)
    }
    results <- c(LB, Z, UB)

    Res <-
      list(
        est = results,
        cl = 1 - alpha,
        c_trimmed = c_trimmed,
        n_original = n_original,
        n = n
      )
    structure(Res, class = c("twonn_mle", class(Res)))
  }


#' Print \code{TWO-NN} MLE output
#'
#' @param x object of class \code{twonn_mle}, obtained from the function
#' \code{twonn_mle()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the TWO-NN estimated via
#' maximum likelihood to console.
#'
#' @export
print.twonn_mle <- function(x, ...) {
  cat("Model: TWO-NN\n")
  cat("Method: MLE\n")
  cat(paste0(
    "Sample size: ",
    x[["n_original"]],
    ", Obs. used: ",
    x[["n"]],
    ". Trimming proportion: ",
    100 * x[["c_trimmed"]],
    "%\n"
  ))
  cat(paste0("ID estimates (confidence level: ", x[["cl"]], ")"))
  y <- cbind(
    `Lower Bound` = x[["est"]][1],
    `Estimate` = x[["est"]][2],
    `Upper Bound` = x[["est"]][3]
  )
  print(knitr::kable(y))
}
