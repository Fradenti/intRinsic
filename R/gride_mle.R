#' Generalized Ratios ID Estimation via MLE
#'
#' The function fits the frequentist Gride model. To run this function, use the
#' high-level \code{gride()} and specify \code{method = "mle"}.
#' The function finds the maximum likelihood estimates for \code{d}, and
#' subsequently simulate parametric bootstrap samples for uncertainty
#' quantification.
#' See \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}
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
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic
#' dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
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




#' @name gride
#'
#' @param object object of class \code{gride_mle}, obtained from the function
#' \code{gride_mle()}.
#' @param ... ignored.
#'
#'
#' @export
print.gride_mle <- function(x, ...) {
  y <- c("Gride - MLE" = unname(x[["est"]][2]))
  print((y))
  invisible(x)
}

#' @name gride
#'
#' @param object object of class \code{gride_mle}, obtained from the function
#' \code{gride_mle()}.
#' @param ... ignored.
#'
#' @export
summary.gride_mle <- function(object, ...) {
  y <- cbind(
    `NN order 1` = object[["n1"]],
    `NN order 2` = object[["n2"]],
    `Bootstrap simulations` = object[["nsim"]],
    `Confidence level` = object[["cl"]],

    `Lower Bound` = object[["est"]][1],
    `Estimate` = object[["est"]][2],
    `Upper Bound` = object[["est"]][3]
  )
  structure(y, class = c("summary.gride_mle","matrix"))
}


#' @name gride
#'
#' @param x object of class \code{twonn_mle}, obtained from the function
#' \code{twonn_mle()}.
#' @param ... ignored.
#'
#' @export
print.summary.gride_mle <- function(x, ...) {
  cat(paste0("Model: Gride(", x[1], ",", x[2], ")\n"))
  cat("Method: MLE\n")
  cat(paste0("CI obtained with a parametric bootstrap sample of size ",
             x[3], "\n"))
  cat(paste0("ID estimates (confidence level: ", x[4], ")"))
  y <- cbind(
    `Lower Bound` = x[5],
    `Estimate` = x[6],
    `Upper Bound` = x[7]
  )
  print(knitr::kable(y))
  invisible(x)
  }


#' @name gride
#'
#' @param x object of class \code{gride_mle}.
#' It is obtained using the output of the \code{gride} function when
#' \code{method = "mle"}.
#'
#' @param ... other arguments passed to specific methods.
#'
#' @export
#'
plot.gride_mle <- function(x,
                           ...) {
  ID <- x$boot_sample
  dx <- density(ID)
  plot(dx, xlab = "Intrinsic Dimension" , ylab = "Bootstrap Density",
       col="darkblue",lwd=1.3, main = "")
  polygon(c(dx$x), c(dx$y),
          col = "lightgray", border = "darkblue", main = "")
  abline(v =  c(x$est, x$lb, x$ub),
         lty = 2,
         col = 2)
  graphics::title("MLE Gride: Bootstrap sample")
  invisible()
}

