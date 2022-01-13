#' Decimated \code{TWO-NN} evolution with proportions
#'
#' The estimation of the intrinsic dimension is related to the scale of the
#' dataset. To escape the local reach of the \code{TWO-NN} estimator,
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al. (2017)}
#' proposed to subsample the original dataset in order to induce greater
#' distances between the data points. By investigating the estimates' evolution
#' as a function of the size of the neighborhood, it is possible to obtain
#' information about the validity of the modeling assumptions and the robustness
#' of the model in the presence of noise.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param proportions vector containing the fractions of the dataset to be
#' considered.
#' @param seed random seed controlling the sequence of subsampling.
#'
#' @return list containing the \code{TWO-NN} evolution (maximum likelihood
#' estimation and confidence intervals), the average distance from the second
#' NN, and the vector of proportions that were considered.
#'
#' @keywords internal
#' @noRd
#'
twonn_dec_prop <- function(X,
                           proportions = 1,
                           seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (any(proportions  > 1) | any(proportions  < 0)) {
    stop(
      "The vector proportions must contains values between 0 (excluded) and 1 (included)",
      call. = FALSE
    )
  }

  n               <- nrow(X)
  if (floor(min(proportions) * n)  <= 2) {
    stop(
      "Proportions are too low, no observations left. Please lower the number of steps considered"
    )
  }

  W               <- length(proportions)
  twonns          <- matrix(NA, W, 3)
  avg_distance_n2 <- numeric(W)
  sub_ind        <-  1:n
  n_new           <- n

  # Classic TWO-NN
  K                  <- FNN::get.knn(X, k = 2)
  mudots             <- K$nn.dist[, 2] / K$nn.dist[, 1]
  avg_distance_n2[1] <- mean(K$nn.dist[, 2])
  ests               <- twonn_mle(mudots)
  twonns[1, ]         <- ests$est

  for (w in 2:(W)) {
    prop       <- proportions[(w)] / proportions[w - 1]

    sub_ind    <- sample(x = sub_ind,
                         size = floor(n_new * prop))

    n_new      <- length(sub_ind)

    subX       <- X[sub_ind, ]
    K          <- FNN::get.knn(subX, k = 2)
    mudots     <- K$nn.dist[, 2] / K$nn.dist[, 1]
    avg_distance_n2[w] <- mean(K$nn.dist[, 2])
    ests       <- twonn_mle(mudots)
    twonns[w, ] <- ests$est

  }

  res <- list(
    decimated_twonn = twonns,
    avg_distance_n2 = avg_distance_n2,
    proportions = proportions
  )
  structure(res, class = c("twonn_dec_prop", class(res)))

}



#' Print \code{TWO-NN} evolution object decimated via vector of proportions
#'
#' @param x object of class \code{twonn_dec_prop}, obtained from the function
#' \code{twonn_dec_prop()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the decimated TWO-NN to console.
#'
#' @export
print.twonn_dec_prop <- function(x, ...) {
  cat(paste0("Model: decimated TWO-NN\n"))
  if (length(x[["proportions"]]) > 10) {
    cat(paste0(
      "Decimating proportions ranging from ",
      round(max(x[["proportions"]]), 5),
      " to ",
      round(min(x[["proportions"]]), 5)
    ),
    ".\n")
  } else{
    cat("Decimating proportions: ", round(x[["proportions"]], 5), ".\n")
  }
  cat(paste0(
    "Average distance from the n2-th NN ranging from ",
    round(min(x[["avg_distance_n2"]]), 4),
    " to ",
    round(max(x[["avg_distance_n2"]]), 4),
    ".\n"
  ))
  invisible(x)
}
