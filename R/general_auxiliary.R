#' Generates a noise-free Swiss roll dataset
#'
#' The function creates a three-dimensional dataset with coordinates
#' following the Swiss roll mapping, transforming random uniform data points
#' sampled on the interval \code{(0,10)}.
#'
#' @param n number of observations contained in the output dataset.
#'
#' @return a three-dimensional \code{data.frame} containing the coordinates of
#' the points generated via the Swiss roll mapping.
#'
#' @export
#'
#' @examples
#' Data <- Swissroll(1000)
#'
Swissroll <- function(n) {
  X <- stats::runif(n, 0, 10)
  Y <- stats::runif(n, 0, 10)
  return(SR = data.frame(
    x = X * cos(X),
    y = Y,
    z = X * sin(X)
  ))
}

