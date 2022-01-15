#' Decimated TWO-NN evolution with halving steps
#'
#' The estimation of the \code{id} is related to the scale of the
#' dataset. To escape the local reach of the \code{TWO-NN} estimator,
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al. (2017)}
#' propose to subsample the original dataset in order to induce greater
#' distances between the data points. By investigating the estimates' evolution
#' as a function of the size of the neighborhood, it is possible to obtain
#' information about the validity of the modeling assumptions and the robustness
#' of the model in the presence of noise.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param steps number of times the dataset is halved.
#' @param seed random seed controlling the sequence of subsampling.
#'
#' @return object of class \code{twonn_decimation_steps}, which is a list
#' containing the \code{TWO-NN} evolution (maximum likelihood estimation and
#' confidence intervals), the average distance from the second NN, and the
#' number of steps.
#'
#' @keywords internal
#' @noRd
#'
twonn_dec_by <- function(X,
                         steps = 2,
                         seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n               <- nrow(X)
  if (floor(2 ^ (-steps) * n)  <= 2) {
    stop(
"Too many steps, no observations left. Please lower the number of steps considered",
      call. = FALSE
    )
  }

  W               <- steps + 1
  twonns          <- matrix(NA, W, 3)
  avg_distance_n2 <- numeric(W)
  inds            <- 1:n

  # Classic TWO-NN
  K                  <- FNN::get.knn(X, k = 2)
  mudots             <- K$nn.dist[, 2] / K$nn.dist[, 1]
  avg_distance_n2[1] <- mean(K$nn.dist[, 2])
  ests               <- twonn_mle(mudots)
  twonns[1, ]         <- ests$est


  for (w in 2:(W)) {
    sub_ind    <- sample(x = inds,
                         size = floor(.5 * length(inds)))

    subX       <- X[sub_ind, ]
    K          <- FNN::get.knn(subX, k = 2)
    mudots     <- K$nn.dist[, 2] / K$nn.dist[, 1]
    avg_distance_n2[w] <- mean(K$nn.dist[, 2])
    ests       <- twonn_mle(mudots)
    twonns[w, ] <- ests$est
    inds <- sub_ind

  }

  res <- list(
    decimated_twonn = twonns,
    avg_distance_n2 = avg_distance_n2,
    steps = steps
  )
  structure(res, class = c("twonn_dec_by", class(res)))

}




#' Print \code{TWO-NN} evolution object decimated via halving steps
#'
#' @param x object of class \code{twonn_dec_prop}, obtained from the function
#' \code{twonn_dec_prop()}.
#' @param ... ignored.
#'
#' @return the function prints a summary of the decimated TWO-NN to console.
#'
#' @export
print.twonn_dec_by <- function(x, ...) {
  cat(paste0("Model: decimated TWO-NN\n"))
  cat(paste0("Dataset halved ", x[["steps"]], " times \n"))
  cat(paste0(
    "Average distance from the n2-th NN ranging from ",
    round(min(x[["avg_distance_n2"]]), 4),
    " to ",
    round(max(x[["avg_distance_n2"]]), 4),
    "\n"
  ))
  invisible(x)
}

