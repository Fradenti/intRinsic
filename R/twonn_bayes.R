#' Bayesian estimator for the \code{TWO-NN} model
#'
#' The function fits the \code{TWO-NN} model via Bayesian estimation, employing
#' a conjugate prior. The formulas can be found in
#' \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}.
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
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic
#' dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
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


#' @name twonn
#'
#' @param x object of class \code{twonn_bayes}, obtained from the function
#' \code{twonn_bayes()}.
#' @param ... ignored.
#'
#'
#' @export
print.twonn_bayes <- function(x, ...) {
  y <- c("TWONN - Bayes - Posterior Mean" = x[["est"]][2])
  cat(y)
  cat("\n")
  invisible(x)
}

#' @name twonn
#'
#' @param object object of class \code{twonn_bayes}, obtained from the function
#' \code{twonn_bayes()}.
#' @param ... ignored.
#'
#' @export
summary.twonn_bayes <- function(object, ...) {
  y <- c(
    `Original sample size` = object[["n_original"]],
    `Used sample size` = object[["n"]],
    `Trimming proportion` = object[["c_trimmed"]],

    `Prior shape` = object[["hp_prior"]][1],
    `Prior scale` = object[["hp_prior"]][2],

    `Credible interval level` = object[["alpha"]],

    `Lower bound` = object[["est"]][1],
    `Posterior mean` = object[["est"]][2],
    `Posterior median` = object[["est"]][3],
    `Posterior mode` = object[["est"]][4],
    `Upper bound` = object[["est"]][5]
  )
  structure(y, class = c("summary.twonn_bayes","matrix"))
}


#' @name twonn
#'
#' @param x object of class \code{twonn_bayes}, obtained from the function
#' \code{twonn_bayes()}.
#' @param ... ignored.
#'
#' @export
print.summary.twonn_bayes <- function(x, ...) {
  cat("Model: TWO-NN\n")
  cat("Method: Bayesian Estimation\n")
  cat(paste0(
    "Sample size: ",
    x[1],
    ", Obs. used: ",
    x[2],
    ". Trimming proportion: ",
    100 * x[3],
    "%\n"
  ))
  cat(paste0("Prior d ~ Gamma(",
             x[4],
             ", ",
             x[5],
             ")\n"))
  cat(paste0(
    "Credible Interval quantiles: ",
    (1 - x[6]) / 2 * 100,
    "%, ",
    (1 + x[6]) / 2 * 100,
    "%\n"
  ))
  cat(paste0("Posterior ID estimates:"))
  y <- cbind(
    `Lower Bound` = x[7],
    `Mean` = x[8],
    `Median` = x[9],
    `Mode` = x[10],
    `Upper Bound` = x[11]
  )
  rownames(y) <- NULL
  print(knitr::kable(y))
  cat("\n")
  invisible(x)
}



#' @name twonn
#'
#' @param x object of class \code{twonn_bayes}, the output of the
#' \code{twonn} function when \code{method = "bayes"}.
#' @param plot_low lower bound of the interval on which the posterior density
#' is plotted.
#' @param plot_upp upper bound of the interval on which the posterior density
#' is plotted.
#' @param by step-size at which the sequence spanning the interval is
#' incremented.
#' @param ... other arguments passed to specific methods.
#'
#' @export
#'
plot.twonn_bayes <-
  function(x,
           plot_low = 0.001,
           plot_upp = NULL,
           by = .05,
           ...) {

    if (is.null(plot_upp)) {
      plot_upp <- x$est[5] + 3
    }

    xx <- seq(plot_low, plot_upp, by = by)
    y0 <- stats::dgamma(xx,
                        shape = x$hp_prior[1],
                        rate =  x$hp_prior[2])
    y <- stats::dgamma(xx,
                       shape = x$hp_posterior[1],
                       rate =  x$hp_posterior[2])

    plot(y0~xx, type="l", col = 4,
         xlab=("Intrinsic Dimension"),
         ylab=("Posterior density"), ylim = c(0, max(c(y,y0)) ))
    lines(y~xx, type="l", col = 1)
    abline(v = x$est,
           lty = 2,
           col = 2)
    graphics::title("Bayesian TWO-NN")
    invisible()
  }


