#' Generates a noise-free Swiss Roll dataset
#'
#' @param N number of observations contained in the output dataset
#'
#' @return a three-dimensional Swiss Roll dataset
#' @export
#'
#' @examples Swissroll_maker(1000)
Swissroll_maker <- function(N) {
  X <- stats::runif(N, 0, 10)
  Y <- stats::runif(N, 0, 10)
  return(SR = cbind(
    x = X * cos(X),
    y = Y,
    z = X * sin(X)
  ))
}
