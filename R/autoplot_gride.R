#' Plot the simulated MCMC chains for the Bayesian \code{Gride}
#'
#' Use this method without the \code{.gride_bayes} suffix.
#' It displays the traceplot of the chain generated
#' with Metropolis-Hasting updates to visually assess mixing and convergence.
#' Alternatively, it is possible to plot the posterior density.
#'
#' @param object object of class \code{gride_bayes}.
#' It is obtained using the output of the \code{gride} function when
#' \code{method = "bayes"}.
#' @param traceplot logical. If \code{FALSE}, the function returns a plot of the
#' posterior density. If \code{TRUE}, the function returns the traceplots of the
#' MCMC used to simulate from the posterior distribution.
#' @param title optional string to display as title.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.gride_bayes
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}.
#' It could represent the traceplot of the posterior simulations for the
#' Bayesian \code{Gride} model (\code{traceplot = TRUE}) or a density plot
#' of the simulated posterior distribution (\code{traceplot = FALSE}).
#'
#' @seealso \code{\link{gride}}
#'
#' @family autoplot methods
#'
#' @export
autoplot.gride_bayes <- function(object,
                                 traceplot = FALSE,
                              title = "Bayesian Gride - Posterior distribution",
                                 ...) {
  if (traceplot) {
    sam <- c(object$post_sample)
    cmm <- cumsum(sam) / seq_along(sam)
    G1 <-
      ggplot2::ggplot(dplyr::tibble(sam, ind = seq_along(sam), cmm = cmm)) +
      ggplot2::geom_line(ggplot2::aes(x = .data$ind, y = .data$sam), col =
                           "lightgray") +
      ggplot2::geom_line(ggplot2::aes(x = .data$ind, y = .data$cmm), col =
                           "darkblue") +
      ggplot2::theme_bw() +
      ggplot2::ylab("Intrinsic Dimension") +
      ggplot2::xlab("MCMC Iteration") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 15),
        axis.title.y = ggplot2::element_text(size = 15),
        title = ggplot2::element_text(size = 15)
      ) +
      ggplot2::ggtitle(title,
        subtitle = bquote(
n[1] == .(object$n1) ~ "," ~ n[2] == .(object$n2) ~ "," ~ sigma == .(object$sigma)
      ))


  } else{
    sam <- object$post_sample
    G1 <- ggplot2::ggplot(dplyr::tibble(c(sam))) +
      ggplot2::geom_density(ggplot2::aes(x = .data$value),
                            col = "darkblue",
                            fill = "lightgray") +
      ggplot2::xlab("Intrinsic Dimension") +
      ggplot2::ylab("Simulated Posterior Density") +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(xintercept = object$est,
                          lty = 2,
                          col = 2) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 15),
        axis.title.y = ggplot2::element_text(size = 15),
        title = ggplot2::element_text(size = 15)
      ) +
      ggplot2::ggtitle(title,
                       subtitle = bquote(
n[1] == .(object$n1) ~ "," ~ n[2] == .(object$n2) ~ "," ~ sigma == .(object$sigma)
      ))


  }
  return(G1)
}




#' Plot the simulated bootstrap sample for the MLE \code{Gride}
#'
#' Use this method without the \code{.gride_mle} suffix.
#' It displays the density plot of sample obtained via
#' parametric bootstrap for the \code{Gride} model.
#'
#' @param object object of class \code{gride_mle}.
#' It is obtained using the output of the \code{gride} function when
#' \code{method = "mle"}.
#' @param title title for the plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.gride_mle
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}. It displays the
#' density plot of the sample generated via parametric bootstrap to help the
#' visual assessment of the uncertainty of the \code{id} estimates.
#'
#'
#' @seealso \code{\link{gride}}
#'
#' @export
autoplot.gride_mle <- function(object,
                               title = "MLE Gride - Bootstrap sample",
                               ...) {
  G1 <- ggplot2::ggplot(dplyr::tibble(object$boot_sample)) +
    ggplot2::geom_density(ggplot2::aes(x = .data$value),
                          col = "darkblue",
                          fill = "lightgray") +
    ggplot2::xlab("Intrinsic Dimension") +
    ggplot2::ylab("Bootstrap Density") +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
      xintercept = c(object$est, object$lb, object$ub),
      lty = 2,
      col = 2
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15),
      title = ggplot2::element_text(size = 15)
    ) +
    ggplot2::ggtitle(title,
          subtitle = bquote(n[1] == .(object$n1) ~ "," ~ n[2] == .(object$n2)))

  return(G1)
}


