#' Maximum Likelihood Estimator for the \code{TWO-NN} model
#'
#' The function fits the \code{TWO-NN} model via maximum likelihood estimation,
#' as discussed in
#' \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}.
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
#' Scientific Reports, 7(1).
#' ISSN 20452322, \doi{10.1038/s41598-017-11873-y}.
#'
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic
#' dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
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


#' @name twonn
#'
#' @param x object of class \code{twonn_mle}, obtained from the function
#' \code{twonn_mle()}.
#' @param ... ignored.
#'
#'
#' @export
print.twonn_mle <- function(x, ...) {
  y <- c("TWONN - MLE" = x[["est"]][2])
  print((y))
  invisible(x)
}

#' @name twonn
#'
#' @param object object of class \code{twonn_mle}, obtained from the function
#' \code{twonn_mle()}.
#' @param ... ignored.
#'
#' @export
summary.twonn_mle <- function(object, ...) {
  y <- cbind(
    `Original sample size` = object[["n_original"]],
    `Used sample size` = object[["n"]],
    `Trimming proportion` = object[["c_trimmed"]],
    `Confidence level` = object[["cl"]],
    `Lower Bound` = object[["est"]][1],
    `Estimate` = object[["est"]][2],
    `Upper Bound` = object[["est"]][3]
  )
  structure(y, class = c("summary.twonn_mle","matrix"))
}


#' @name twonn
#'
#' @param x object of class \code{twonn_mle}, obtained from the function
#' \code{twonn_mle()}.
#' @param ... ignored.
#'
#' @export
print.summary.twonn_mle <- function(x, ...) {
  cat("Model: TWO-NN\n")
  cat("Method: MLE\n")
  cat(paste0(
    "Sample size: ",
    x[1],
    ", Obs. used: ",
    x[2],
    ". Trimming proportion: ",
    100 * x[3],
    "%\n"
  ))
  cat(paste0("ID estimates (confidence level: ", x[4], ")"))
  y <- cbind(
    `Lower Bound` = x[5],
    `Estimate` = x[6],
    `Upper Bound` = x[7]
  )
  print(knitr::kable(y))
  invisible(x)
}


#' @name twonn
#' @param x object of class \code{twonn_mle}, the output of the
#' \code{twonn} function when \code{method = "mle"}.
#'
#'
#' @importFrom graphics abline legend lines matplot par points polygon
#' @importFrom stats density ts
#'
#' @export
#'
plot.twonn_mle <-
  function(x,
           ...) {

    plot(rep(0,3) ~ x$est,type="n",xlab = ("Maximum Likelihood Estimation"),ylab="")
    graphics::abline(v = x$est,col = "gray",lwd=1,lty=2)
    graphics::points(0~x$est[2],cex = 2,pch=21,bg=1)
    graphics::arrows(y0 = 0,code = 3,y1 = 0,x0 = x$est[1],x1 = x$est[3],angle=90,lwd=2)
    graphics::title("MLE TWO-NN")
    invisible()
  }


