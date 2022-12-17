#' \code{Gride} evolution based on Maximum Likelihood Estimation
#'
#' The function allows the study of the evolution of the \code{id} estimates
#' as a function of the scale of a dataset. A scale-dependent analysis
#' is essential to identify the correct number of relevant directions in noisy
#' data. To increase the average distance from the second NN (and thus the
#' average neighborhood size) involved in the estimation, the function computes
#' a sequence of \code{Gride} models with increasing NN orders, \code{n1} and
#' \code{n2}.
#' See also \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}
#' for more details.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param vec_n1 vector of integers, containing the smaller NN orders considered
#' in the evolution.
#' @param vec_n2 vector of integers, containing the larger NN orders considered
#' in the evolution.
#' @param upp_bound upper bound for the interval used in the numerical
#' optimization (via \code{optimize}). Default is set to 50.
#'
#' @return list containing the Gride evolution, the corresponding NN distance
#' ratios, the average n2-th NN order distances, and the NN orders considered.
#'
#' @name gride_evolution
#'
#' @export
#'
#' @references
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic
#' dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
#'
#' @examples
#' \donttest{
#' X       <-  replicate(5,rnorm(10000,0,.1))
#' gride_evolution(X = X,vec_n1 = 2^(0:5),vec_n2 = 2^(1:6))
#'}
#'
gride_evolution <- function(X, vec_n1, vec_n2, upp_bound = 50) {
  if (any(vec_n2 < vec_n1)) {
    stop("at least one n2 is lower than one n1", call. = FALSE)
  }
  if (length(vec_n1) != length(vec_n2)) {
    stop("the considered NN orders in n1 and n2 do not match", call. = FALSE)
  }

  D        <- ncol(X)
  n        <- nrow(X)
  W        <- length(vec_n1)
  MUdots   <- matrix(NA, n, W)
  K        <- FNN::get.knn(X, k = max(vec_n2))
  path     <- numeric(W)
  avg_distance_n2 <- numeric(W)

  for (w in 1:W) {
    MUdots[, w] <- (K$nn.dist[, vec_n2[w]]) /
      (K$nn.dist[, vec_n1[w]])

    avg_distance_n2[w] <- mean(K$nn.dist[, vec_n2[w]])
    path[w] <- stats::optimize(
      gride_log_likelihood,
      interval = c(0.01, min(D, upp_bound) + 1),
      n1 = vec_n1[w],
      n2 = vec_n2[w],
      mus_n1_n2 = MUdots[, w],
      maximum = TRUE
    )$max
  }

  res <- list(
    path     = path,
    MUdots   = MUdots,
    NNorders = rbind(vec_n1, vec_n2),
    avg_distance_n2 = avg_distance_n2
  )

  structure(res, class = c("gride_evolution", class(res)))

}


#' @name gride_evolution
#'
#' @param x object of class \code{gride_evolution}, obtained from the function
#' \code{gride_evolution()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the Gride evolution to
#' console.
#'
#' @export
print.gride_evolution <- function(x, ...) {
  cat(paste0("Model: Gride evolution\n"))
  cat(paste0("Smaller NN order ranging from ",
             min(x[["NNorders"]][1, ]),
             " to ",
             max(x[["NNorders"]][1, ]), "\n"))
  cat(paste0("Larger NN order ranging from ",
             min(x[["NNorders"]][2, ]),
             " to ",
             max(x[["NNorders"]][2, ]), "\n"))
  cat(paste0(
    "Average distance from the n2-th NN ranging from ",
    round(min(x[["avg_distance_n2"]]), 4),
    " to ",
    round(max(x[["avg_distance_n2"]]), 4),
    "\n"
  ))

  invisible(x)
}



#' @name gride_evolution
#'
#' @param x an object of class \code{gride_evolution}.
#'
#' @param ... other arguments passed to specific methods.
#'
#' @export
#'
plot.gride_evolution <- function(x,
                                     ...) {
  id     <- x$path
  avg_n2 <- x$avg_distance_n2

  plot(id~avg_n2,type = "b", col = "darkblue",log = "x",
       ylab = ("Intrinsic dimension"),
       xlab = Log[10] ~ average ~ n[2] ~ distance)
  graphics::title("Gride Evolution")
  invisible()
}

#' Plot the evolution of \code{Gride} estimates
#'
#' Use this method without the \code{.gride_evolution} suffix.
#' It plots the evolution of the \code{id}
#' estimates as a function of the average distance from the furthest NN of
#' each point.
#'
#' @param object an object of class \code{gride_evolution}.
#' @param title an optional string to customize the title of the plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.gride_evolution
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}. It displays the
#' the evolution of the Gride maximum likelihood estimates as a function
#' of the average distance from \code{n2}.
#'
#' @export
#'
autoplot.gride_evolution <- function(object,
                                     title = "Gride Evolution",
                                     ...) {
  D <- data.frame(id = object$path,
                  n2 = object$avg_distance_n2)
  G1 <- ggplot2::ggplot(D) +
    ggplot2::geom_path(ggplot2::aes(x = .data$n2,
                                    y = .data$id),
                       col = "darkblue") +
    ggplot2::geom_point(ggplot2::aes(x = .data$n2,
                                     y = .data$id),
                        col = "darkblue") +
    ggplot2::ylab("Intrinsic dimension") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::xlab(Log[10] ~ average ~ n[2] ~ distance) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      title = ggplot2::element_text(size = 20)
    ) +
    ggplot2::ggtitle(title)

  return(G1)
}
