#' Least Squares Estimator for the \code{TWO-NN} model
#'
#' The function fits the \code{TWO-NN} model via least squares estimation, as
#' originally proposed in
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al., 2017}.
#'
#' @param mus vector of second to first NN distance ratios.
#' @param alpha confidence level needed for the computation of the confidence
#' interval.
#' @param c_trimmed proportion of trimmed observations.
#'
#' @return object of class \code{twonn_linfit}, which is a list containing the
#' least squares estimate of the intrinsic dimension, along with the \code{lm()}
#' output used for the computation.
#'
#' @seealso \code{\link{twonn}}, \code{\link{autoplot.twonn_linfit}}
#' @keywords Internal
#'
#' @noRd
#'
#' @references
#' Facco E, D'Errico M, Rodriguez A, Laio A (2017). "Estimating the intrinsic
#' dimension of datasets by a minimal neighborhood information."
#' Scientific Reports, 7(1), 1-8.
#' ISSN 20452322, doi: 10.1038/s41598-017-11873-y.
#'
twonn_linfit <- function(mus,
                         alpha = 0.95,
                         c_trimmed = .01) {
  c_considered    <- 1 - c_trimmed
  n <- n_original <- base::length(mus)

  if (c_trimmed > 0) {
    mus <- base::sort(mus)[1:floor(n * c_considered)]
    n  <- base::length(mus)
  }

  F_musi <- (0:(n - 1)) / n
  y <- -base::log(1 - (F_musi))
  x <- base::sort(base::log(mus))

  modlin1 <- stats::lm(y ~ x - 1)

  Res          <- numeric(3)
  Res[c(1, 3)] <- stats::confint.lm(modlin1, level = alpha)
  Res[2]       <- stats::coef(modlin1)
  Res <- list(
    est = Res,
    lm_output = modlin1,
    conflev = alpha,
    c_trimmed = c_trimmed,
    n_original = n_original,
    n = n
  )

  structure(Res, class = c("twonn_linfit", class(Res)))
}

#' Print TWO-NN Least Squares output
#'
#' @param x object of class \code{twonn_linfit}, obtained from the function
#' \code{twonn_linfit()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the TWO-NN estimated via least
#' squares to console.
#'
#' @export
print.twonn_linfit <- function(x, ...) {
  cat("Model: TWO-NN\n")
  cat("Method: Least Square Estimation\n")
  cat(paste0(
    "Sample size: ",
    x[["n_original"]],
    ", Obs. used: ",
    x[["n"]],
    ". Trimming proportion: ",
    100 * x[["c_trimmed"]],
    "%\n"
  ))
  cat(paste0("ID estimates (confidence level: ", x[["conflev"]], ")"))
  y <- cbind(
    `Lower Bound` = x[["est"]][1],
    `Estimate` = x[["est"]][2],
    `Upper Bound` = x[["est"]][3]
  )
  print(knitr::kable(y))
}
