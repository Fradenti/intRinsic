#' Plot a summary of the Linear Fit TWO-NN model
#'
#' @param object an object of class \code{lnf_twonn}, the output of the \code{linfit_twonn} function.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.lnf_twonn
#'
#' @import  ggplot2
#'
#' @return a plot that shows the goodness of the linear fit
#' @export
#'
#' @examples
#' \dontrun{
#' m <- linfit_twonn(data)
#' plot(m)
#' }
autoplot.lnf_twonn <- function(object, ...) {
  lmod <- object$Lin_mod
  x <- lmod$model$x
  y <- lmod$model$y
  Res <- object$Estimates

  p1 <- qplot(x, y) +
    theme_bw() +
    geom_abline(
      intercept = 0,
      slope = lmod$coefficients,
      col = I("red")
    ) +
    ylab("-log(1-(i/N))") +
    xlab(expression(log(mu))) +
    annotate(
      "label",
      -Inf,
      Inf,
      label = paste("ID:", round(Res[2], 3)),
      hjust = -0.05,
      vjust = 1.1
    ) +
    annotate(
      "label",
      Inf,
      -Inf,
      label = paste("R^2:", round(summary(lmod)$r.squared, 3)),
      parse = T,
      hjust = 1.05,
      vjust = -0.1
    )

  p1
}

#' Plot the TWO-NN posterior distribution
#'
#' @param object an object of class \code{byf_twonn}, the output of the \code{bayesfit_twonn} function.
#' @param plot_low the lower bound of the interval on which the posterior density is plotted.
#' @param plot_upp the upper bound of the interval on which the posterior density is plotted.
#' @param by the step-size at which the sequence spanning the interval is incremented
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.byf_twonn
#'
#' @import ggplot2
#'
#' @return a plot of the posterior distribution for the TWO-NN model
#' @export
#'
#' @examples
#' \dontrun{
#' m <- bayesfit_twonn(data)
#' plot(m)
#' }
autoplot.byf_twonn <-
  function(object,
           plot_low = 0,
           plot_upp = NULL,
           by = .05,
           ...) {
    if (is.null(plot_upp)) {
      plot_upp <- object$Estimates[5] + 3
    }

    x <- seq(plot_low, plot_upp, by = by)
    y0 <- stats::dgamma(x,
                        shape = object$hp_prior[1],
                        rate =  object$hp_prior[2])
    y <- stats::dgamma(x,
                       shape = object$hp_posterior[1],
                       rate =  object$hp_posterior[2])
    G1 <- ggplot(dplyr::as_tibble(x, y)) +
      geom_line(aes(x = x, y = y0), col = 4) +
      geom_line(aes(x = x, y = y)) +
      xlab("Intrinsic Dimension") +
      ylab("Posterior density") +
      theme_bw() +
      geom_vline(xintercept = object$Estimates,
                 lty = 2,
                 col = 2)
    G1
  }
