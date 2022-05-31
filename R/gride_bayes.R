#' Generalized Ratios ID Estimation via Bayesian approach
#'
#' The function fits the Bayesian Gride model. To run this function, use the
#' high-level \code{gride()} and specify \code{method = "bayes"}. The function
#' runs a Metropolis-Hasting (MH) algorithm over \code{log(d)}, adopting a
#' Normal distribution with pre-specified standard deviation \code{sigma} as
#' proposal distribution.
#' See \href{https://arxiv.org/abs/2104.13832}{Denti et al., 2022+}
#' for more details.
#'
#' @param mus_n1_n2 vector of generalized order NN distance ratios.
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param nsim number of MCMC iterations to collect in the sample.
#' @param burn_in number of iterations to discard from the MCMC sample.
#' @param sigma standard deviation of the Gaussian proposal used in the MH step.
#' @param start_d initial value for the MCMC chain. If \code{NULL}, the MLE is
#' computed.
#' @param a_d shape parameter of the Gamma prior distribution for \code{d}.
#' @param b_d rate parameter of the Gamma prior distribution for \code{d}.
#' @param alpha  posterior probability contained in computed credible interval.
#' @param upper_D upper bound for the id parameter, needed if \code{start_d}
#' is not provided to initiate the chain with the MLE estimate.
#'
#' @return a list containing the MCMC sample and the summary of the
#' specified arguments.
#'
#' @keywords Internal
#' @noRd
#'
#' @references
#' Denti F, Doimo D, Laio A, Mira A (2022+). "Distributional Results for
#' Model-Based Intrinsic Dimension Estimators."
#' arXiv preprint. 2104.13832, \url{https://arxiv.org/abs/2104.13832}.
#'
#' @seealso \code{\link{gride}}
#'
gride_bayes <- function(mus_n1_n2 = NULL,
                        n1 = 1,
                        n2 = 2,
                        nsim = 5000,
                        burn_in = 2000,
                        sigma = 0.5,
                        start_d = NULL,
                        a_d = 1,
                        b_d = 1,
                        alpha = 0.95,
                        upper_D = NULL) {


  if (n2 < n1) {
    stop("n2 should be greater than n1", call. = FALSE)
  }


  if (is.null(start_d)) {
    start_d <- gride_mle_point(
      mus_n1_n2 = mus_n1_n2,
      n1 = n1,
      n2 = n2,
      upper_D = upper_D
    )
    start_d <- log(start_d - 1)
  } else{
    if (start_d <= 1) {
      stop("Invalid starting point: it has to be > 1", call. = FALSE)
    }
    start_d  <- log(start_d - 1)
  }

  post_sample <- gride_mh_sampler(
    start_d = start_d,
    n1 = n1,
    n2 = n2,
    a_d = a_d,
    b_d = b_d,
    mus_n1_n2 = mus_n1_n2,
    nsim = nsim,
    burn_in = burn_in,
    sigma = sigma
  )

  sam         <- (exp(post_sample) + 1)
  est         <- numeric(5)
  est[2]      <- mean(sam)
  est[3]      <- stats::median(sam)
  dens        <- stats::density(sam)
  est[4]      <- dens$x[which.max(dens$y)]
  est[c(1, 5)] <- stats::quantile(sam, c((1 - alpha) / 2, (1 + alpha) /
                                           2))
  Res <- list(
    est = est,
    post_sample = sam,
    hp_prior = c(a_d, b_d),
    alpha = alpha,
    n1 = n1,
    n2 = n2,
    nsim = nsim,
    sigma = sigma
  )
  structure(Res, class = c("gride_bayes", class(Res)))
}



#' @name gride
#'
#' @param x object of class \code{gride_bayes}, obtained from the function
#' \code{gride_bayes()}.
#' @param ... ignored.
#'
#'
#' @export
print.gride_bayes <- function(x, ...) {
  cat(paste0("Gride(",x[["n1"]],",",x[["n2"]],") - Bayes - Posterior Mean\n"))
  cat(x[["est"]][2])
  invisible(x)
}

#' @name gride
#'
#' @param object object of class \code{gride_bayes}, obtained from the function
#' \code{gride_bayes()}.
#' @param ... ignored.
#'
#' @export
summary.gride_bayes <- function(object, ...) {
  y <- c(
    `MCMC Iterations` = object[["nsim"]],


    `Prior shape` = object[["hp_prior"]][1],
    `Prior scale` = object[["hp_prior"]][2],

    `Credible interval level` = object[["alpha"]],

    `n1` = object[["n1"]],
    `n2` = object[["n2"]],

    `Lower bound` = object[["est"]][1],
    `Posterior mean` = object[["est"]][2],
    `Posterior median` = object[["est"]][3],
    `Posterior mode` = object[["est"]][4],
    `Upper bound` = object[["est"]][5]
  )
  structure(y, class = c("summary.gride_bayes","matrix"))
}


#' @name gride
#'
#' @param x object of class \code{gride_bayes()}, obtained from the function
#' \code{gride_bayes()}.
#' @param ... ignored.
#'
#'
#' @export
print.summary.gride_bayes <- function(x, ...) {
  cat(paste0("Model: Gride(", x[5], ",", x[6], ")\n"))
  cat("Method: Bayesian Estimation\n")
  cat(paste0("Prior d ~ Gamma(",
             x[2],
             ", ",
             x[3],
             ")\n"))
  cat(paste0("MCMC posterior sample size: ", x[["nsim"]], "\n"))
  cat(paste0(
    "Credibile Interval quantiles: ",
    (1 - x[4]) / 2 * 100,
    "%, ",
    (1 + x[4]) / 2 * 100,
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
  print(knitr::kable(y))
  invisible(x)
}



#' @name gride
#'
#' @param x object of class \code{gride_bayes}.
#' It is obtained using the output of the \code{gride} function when
#' \code{method = "bayes"}.
#'
#' @export
#'
plot.gride_bayes <- function(x,
                             ...) {

    on.exit({par(my_par)}, add = TRUE, after = TRUE)
    my_par <- par(mfrow = c(2, 1))

    ID <- c(x$post_sample)
    cmm <- cumsum(ID) / seq_along(ID)

    plot(ID ,col="gray",type="l", xlab = "MCMC Iteration")
    lines(stats::ts(cmm), col = "darkblue",lwd=1.3)
    graphics::title("Bayesian Gride: Traceplot")


    dx <- density(ID)
    plot(dx, xlab = "Intrinsic Dimension" , ylab = "Simulated Posterior Density",
         col="darkblue",lwd=1.3, main = "")
    polygon(c(dx$x), c(dx$y),
            col = "lightgray", border = "darkblue", main = "")
    graphics::title("Bayesian Gride: Posterior Density")

    invisible()
      }

