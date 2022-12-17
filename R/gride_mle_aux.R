#' Bootstrap sample generator for the \code{Gride} MLE
#'
#' @param mus_n1_n2 vector of generalized order NN distance ratios.
#' @param n1 order of the first NN considered.
#' @param n2 order of the second NN considered.
#' @param nsim integer, the number of bootstrap replications to consider.
#' @param upper_D nominal dimension of the dataset (upper bound for the
#' maximization routine).
#'
#' @keywords Internal
#' @noRd
#'
#' @return list containing the MLE, confidence interval, and the
#' bootstrap sample used for the estimation.
#'
gride_bootstrap <- function(mus_n1_n2,
                            n1 = 1,
                            n2 = 2,
                            nsim = 2000,
                            upper_D = NULL) {

  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }

  if (is.null(upper_D)) {
    stop("Please provide the nominal dimension of the dataset D in upper_D",
         call. = FALSE)
  }

  n       <- length(mus_n1_n2)

  mle_est <- gride_mle_point(
    mus_n1_n2 = mus_n1_n2,
    n1 = n1,
    n2 = n2,
    upper_D = upper_D
  )

  boot_mus_samples <- replicate(nsim, rgera(
    nsim = n,
    n1 = n1,
    n2 = n2,
    d = mle_est
  ))

  bootstrap_sample <- apply(boot_mus_samples,
                            2,
                            function(x)
                              stats::optimize(
                                gride_log_likelihood,
                                interval = c(1.01, upper_D),
                                n1 = n1,
                                n2 = n2,
                                mus_n1_n2 = x,
                                maximum = TRUE
                              )$max)

  Res <- list(mle = mle_est, boot_sample = bootstrap_sample)

  structure(Res, class = c("bootstrap_sample", "gride", class(Res)))

}

#' Maximum Likelihood Estimator (MLE) for the id using generic ratio statistics
#'
#' @param mus_n1_n2 vector of generalized order NN distance ratios statistics.
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param upper_D nominal dimension of the dataset (upper bound for the
#' maximization routine).
#'
#' @keywords Internal
#' @noRd
#'
#' @return Gride MLE point estimate obtained via numeric optimization.
#'
gride_mle_point <- function(mus_n1_n2,
                            n1 = 1,
                            n2 = 2,
                            upper_D = NULL) {
  if (class(mus_n1_n2)[1] == "mus") {
    n1 <- attr(mus_n1_n2, which = "n1")
    n2 <- attr(mus_n1_n2, which = "n2")
  }

  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }

  if (is.null(upper_D)) {
    stop("Please provide the nominal dimension of the dataset D in upper_D",
         call. = FALSE)
  }

  O       <- stats::optimize(
    gride_log_likelihood,
    interval = c(1.01, upper_D),
    n1 = n1,
    n2 = n2,
    mus_n1_n2 = mus_n1_n2,
    maximum = TRUE
  )


  structure(O$max, class = c("mle_point", "gride", class(O$max)))

}
