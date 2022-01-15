#' Generalized Ratios ID Estimation via MLE
#'
#' The function fits the frequentist Gride model. To run this function, use the
#' high-level \code{gride()} and specify \code{method = "mle"}.
#' The function finds the maximum likelihood estimates for \code{d}, and
#' subsequently simulate parametric bootstrap samples for uncertainty
#' quantification.
#' See \href{https://arxiv.org/abs/2104.13832}{Denti et al., 2021+}
#' for more details.
#'
#'
#' @param mus_n1_n2 vector of generalized order NN distance ratios.
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param nsim number of bootstrap simulations to consider.
#' @param alpha confidence level for the computation of the confidence interval.
#' @param upper_D nominal dimension of the dataset (upper bound for the
#' maximization routine).
#'
#' @return MLE estimate obtained via numeric optimization along with the
#' bootstrap confidence interval.
#' @keywords internal
#' @noRd
#'
#' @references
#' Denti F, Doimo D, Laio A, Mira A (2022+). "Distributional Results for
#' Model-Based Intrinsic Dimension Estimators."
#' arXiv preprint. 2104.13832, \url{https://arxiv.org/abs/2104.13832}.
#'
#' @seealso \code{\link{gride}}
#'
gride_mle <- function(mus_n1_n2 = NULL,
                      n1 = 1,
                      n2 = 2,
                      nsim = 2000,
                      alpha = .95,
                      upper_D = NULL) {
  if (class(mus_n1_n2)[1] == "mus") {
    n1 <- attr(mus_n1_n2, which = "n1")
    n2 <- attr(mus_n1_n2, which = "n2")
  }

  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }

  one_m_alpha <- 1 - alpha

  bs <- gride_bootstrap(
    mus_n1_n2 = mus_n1_n2,
    n1 = n1,
    n2 = n2,
    nsim = nsim,
    upper_D = upper_D
  )

  qq  <- base::unname(stats::quantile(bs$boot_sample,
                                      probs = c(one_m_alpha / 2,
                                                1 - one_m_alpha / 2)))

  Res <- list(
    est = c(
      lb = qq[1],
      mle = bs$mle,
      ub = qq[2]
    ),
    cl = alpha,
    boot_sample = bs$boot_sample,
    n1 = n1,
    n2 = n2,
    nsim = nsim
  )
  structure(Res, class = c("gride_mle", class(Res)))

}


#' Print \code{Gride} MLE object
#'
#' @param x object of class \code{gride_mle}, obtained from the function
#' \code{gride_mle()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the Gride estimated via maximum
#' likelihood to console.
#' @export
print.gride_mle <- function(x, ...) {
  cat(paste0("Model: Gride(", x[["n1"]], ",", x[["n2"]], ")\n"))
  cat("Method: MLE\n")
  cat(paste0("CI obtained with a parametric bootstrap sample of size ",
             x[["nsim"]], "\n"))
  cat(paste0("ID estimates (confidence level: ", x[["cl"]], ")"))
  y <- cbind(
    `Lower Bound` = x[["est"]][1],
    `Estimate` = x[["est"]][2],
    `Upper Bound` = x[['est']][3]
  )
  print(knitr::kable(y))
}
