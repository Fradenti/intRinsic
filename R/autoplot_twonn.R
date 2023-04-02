#' Plot the output of the \code{TWO-NN} model estimated via least squares
#'
#' Use this method without the \code{.twonn_linfit} suffix.
#' The function returns the representation of the linear
#' regression that is fitted with the \code{linfit} method.
#'
#' @param object object of class \code{twonn_linfit}, the output of the
#' \code{twonn} function when \code{method = "linfit"}.
#' @param title string used as title of the plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.twonn_linfit
#'
#' @seealso \code{\link{twonn}}
#'
#' @return a \code{\link[ggplot2]{ggplot2}} object displaying the goodness of
#' the linear fit of the TWO-NN model.
#'
#' @family autoplot methods
#'
#' @export
autoplot.twonn_linfit <- function(object,
                                  title = "TWO-NN Linear Fit",
                                  ...) {
  lmod <- object$lm_output
  x    <- lmod$model$x
  y    <- lmod$model$y
  Res  <- object$est

  p1 <- ggplot2::ggplot()+
    ggplot2::geom_point(ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(
      intercept = 0,
      slope = lmod$coefficients,
      col = I("red")
    ) +
    ggplot2::ylab("-log(1-(i/N))") +
    ggplot2::xlab(expression(log(mu))) +
    ggplot2::annotate(
      "label",
      -Inf,
      Inf,
      label = paste("ID:", round(Res[2], 3)),
      hjust = -0.05,
      vjust = 1.1
    ) +
    ggplot2::annotate(
      "label",
      Inf,
      -Inf,
      label = paste("R^2:", round(summary(lmod)$r.squared, 3)),
      parse = TRUE,
      hjust = 1.05,
      vjust = -0.1
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15),
      title = ggplot2::element_text(size = 15)
    )

  p1 + ggplot2::ggtitle(title)
}


#' Plot the output of the \code{TWO-NN} model estimated via the Bayesian
#' approach
#'
#' Use this method without the \code{.twonn_bayes} suffix.
#' The function returns the density plot of the
#' posterior distribution computed with the \code{bayes} method.
#'
#' @param object object of class \code{twonn_bayes}, the output of the
#' \code{twonn} function when \code{method = "bayes"}.
#' @param title character string used as title of the plot.
#' @param plot_low lower bound of the interval on which the posterior density
#' is plotted.
#' @param plot_upp upper bound of the interval on which the posterior density
#' is plotted.
#' @param by step-size at which the sequence spanning the interval is
#' incremented.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.twonn_bayes
#'
#' @seealso \code{\link{twonn}}
#'
#' @return \code{\link[ggplot2]{ggplot2}} object displaying the posterior
#' distribution of the intrinsic dimension parameter.
#'
#' @family autoplot methods
#'
#' @export
autoplot.twonn_bayes <-
  function(object,
           plot_low = 0,
           plot_upp = NULL,
           by = .05,
           title = "Bayesian TWO-NN",
           ...) {
    if (is.null(plot_upp)) {
      plot_upp <- object$est[5] + 3
    }

    x <- seq(plot_low, plot_upp, by = by)
    y0 <- stats::dgamma(x,
                        shape = object$hp_prior[1],
                        rate =  object$hp_prior[2])
    y <- stats::dgamma(x,
                       shape = object$hp_posterior[1],
                       rate =  object$hp_posterior[2])
    G1 <- ggplot2::ggplot(dplyr::tibble(x, y)) +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y0), col = 4) +
      ggplot2::geom_line(ggplot2::aes(x = x, y = y)) +
      ggplot2::xlab("Intrinsic Dimension") +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(xintercept = object$est,
                          lty = 2,
                          col = 2) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(size = 20),
        title = ggplot2::element_text(size = 20)
      )
    G1 + ggplot2::ggtitle(title)
  }


#' Plot the output of the \code{TWO-NN} model estimated via the Maximum
#' Likelihood approach
#'
#' Use this method without the \code{.twonn_mle} suffix.
#' The function returns the point estimate along with the confidence bands
#' computed via the \code{mle} method.
#'
#' @param object object of class \code{twonn_mle}, the output of the
#' \code{twonn} function when \code{method = "mle"}.
#' @param title character string used as title of the plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.twonn_mle
#'
#' @seealso \code{\link{twonn}}
#'
#' @return \code{\link[ggplot2]{ggplot2}} object displaying the point estimate
#' and confidence interval obtained via the maximum likelihood approach of the
#' \code{id} parameter.
#'
#' @family autoplot methods
#'
#' @export
autoplot.twonn_mle <-
  function(object,
           title = "MLE TWO-NN",
           ...) {
    D  <- data.frame(x = object$est, y = 0)
    G1 <- ggplot2::ggplot() +
      ggplot2::geom_segment(ggplot2::aes(
        x = D[1, 1],
        xend = D[3, 1],
        y = D[1, 2],
        yend = D[3, 2]
      )) +
      ggplot2::geom_vline(
        xintercept = D[, 1],
        col = "gray",
        lty = 2,
        size = .5
      ) +
      ggplot2::geom_point(ggplot2::aes(x = D[c(1), 1], y = D[c(1), 2]), size =
                            10, pch = "[") +
      ggplot2::geom_point(ggplot2::aes(x = D[c(3), 1], y = D[c(3), 2]), size =
                            10, pch = "]") +
      ggplot2::geom_point(ggplot2::aes(x = D[c(2), 1], y = D[c(2), 2]), size =
                            10) +
      ggplot2::xlab("Maximum Likelihood Estimation") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 20),
        title = ggplot2::element_text(size = 20)
      )
    G1 + ggplot2::ggtitle(title)
  }




