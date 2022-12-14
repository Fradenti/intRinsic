#' Compute the ratio statistics needed for the intrinsic dimension estimation
#'
#' The function \code{compute_mus} computes the ratios of distances between
#' nearest neighbors (NNs) of generic order, denoted as
#' \code{mu(n_1,n_2)}.
#' This quantity is at the core of all the likelihood-based methods contained
#' in the package.
#'
#' @param X a dataset with \code{n} observations and \code{D} variables.
#' @param dist_mat a distance matrix computed between \code{n} observations.
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param Nq logical indicator. If \code{TRUE}, it provides the \code{N^q}
#' matrix needed for fitting the Hidalgo model.
#' @param q integer, number of NN considered to build \code{N^q}.
#'
#' @references
#' Facco E, D'Errico M, Rodriguez A, Laio A (2017). "Estimating the intrinsic
#' dimension of datasets by a minimal neighborhood information."
#' Scientific Reports, 7(1).
#' ISSN 20452322, \doi{10.1038/s41598-017-11873-y}.
#'
#' Denti F, Doimo D, Laio A, Mira A (2022). "The generalized ratios intrinsic
#' dimension estimator."
#' Scientific Reports, 12(20005).
#' ISSN  20452322, \doi{10.1038/s41598-022-20991-1}.
#'
#' @return a vector containing the ratio statistics, an object of class
#' \code{mus}. The length of the vector is equal to the number of observations
#' considered, unless ties are present in the dataset. In that case, the
#' duplicates are removed. Optionally, if \code{Nq} is \code{TRUE}, the function
#' returns a list containing both the ratio statistics and the adjacency
#' matrix \code{N^q}.
#
#' @export
#'
#'
#' @examples
#' X           <- replicate(2,rnorm(1000))
#' mu          <- compute_mus(X, n1 = 1, n2 = 2)
#' mudots      <- compute_mus(X, n1 = 4, n2 = 8)
#' pre_hidalgo <- compute_mus(X, n1 = 4, n2 = 8, Nq = TRUE, q = 3)
compute_mus <- function(X = NULL,
                        dist_mat = NULL,
                        n1 = 1,
                        n2 = 2,
                        Nq = FALSE,
                        q  = 3) {
  if (n2 < n1) {
    stop("n2 is lower than n1", call. = FALSE)
  }

  if (is.null(dist_mat) & is.null(X)) {
    stop("Please provide either a dataset X or a distance matrix",
         call. = FALSE)
  }

  D <- NULL

  if (is.null(dist_mat)) {
    if(any(is.na(X))){
      stop("There are missing values in the provided dataset.
Please remove the problematic observations and try again.",
           call. = FALSE)
    }

    X     <- as.matrix(X)
    D     <- ncol(X)
    n     <- n0 <- nrow(X)
    check <- 0
    check <- duplicated(X)

    if (sum(check) > 0) {
      X <- X[-which(check),]
      n <- nrow(X)
      warning(
        paste0(
          "\n  Duplicates are present and will be removed.
        \n  Original sample size: ",
          n0,
          ". New sample size: ",
          n,
          "."
        ),
        call. = FALSE
      )
    }

    if (!Nq) {
      K   <- FNN::get.knn(X, k = n2)
      mus <- K$nn.dist[, n2] / K$nn.dist[, n1]

    } else{
      K   <- FNN::get.knn(X, k = max(n2, q))
      mus <- K$nn.dist[, n2] / K$nn.dist[, n1]
      NQ  <- matrix(0, n, n)

      for (h in 1:n) {
        NQ[h, K$nn.index[h,]] <- 1
      }
      mus <- list(mus = mus, NQ = NQ)

    }

    attr(mus, "upper_D") <- D

  } else {
    # checks on dist_mat

    # if it is of class dist, then transform it into a matrix
    if(inherits(dist_mat,"dist")){
      dist_mat <- as.matrix(dist_mat)
    }

    ## does it contain non-negative distances?
    if(!all(dist_mat>=0)){
      stop("Negative distances detected in dist_mat. Please provide a valid distance matrix",
           call. = FALSE)
    }
    # is it symmetric?
    if(!isSymmetric(dist_mat)){
      stop("The provided distance matrix is not symmetric.
Please provide a valid distance matrix",
           call. = FALSE)
    }
    # NA?
    if(any(is.na(X))){
      stop("There are missing values in the provided distance matrix.
Please remove the problematic observations and try again.",
           call. = FALSE)
    }
    n0    <- nrow(dist_mat)
    dummy <- dist_mat
    dummy[lower.tri(dummy, diag = TRUE)] <- -1
    inds     <- unique(which(dummy == 0, arr.ind = TRUE)[, 2])

    if (length(inds) > 0) {
      dist_mat <- dist_mat[-inds, -inds]
      n        <- nrow(dist_mat)

      warning(
        paste0(
          "\n  Duplicates are present and will be removed.
             Original sample size: ",
          n0,
          ". New sample size: ",
          n,
          "."
        )
      )
    }

    if (!Nq) {
      sDistMat <- apply(dist_mat, 1, function(x)
        sort(x, index = TRUE))
      mus      <- unlist(lapply(sDistMat,
                                function(z)
                                  z$x[n2 + 1] / z$x[n1 + 1]))

    } else{
      NQ       <- matrix(0, n, n)
      for (h in 1:n) {
        NQ[h, (sDistMat[[h]]$ix)[2:q]] <- 1
      }
      mus <- list(mus = mus, Nq = NQ)

    }
  }

  attr(mus, "n1") <- n1
  attr(mus, "n2") <- n2
  structure(mus, class = c("mus", class(mus)))

}



#' @name compute_mus
#'
#' @param x object of class \code{mus}, obtained from the
#' function \code{compute_mus()}.
#' @param ... ignored.
#'
#' @export
print.mus <- function(x, ...) {

  if(is.list(x)){ #check if the mus object contains Nq
    nn <- length(x[[1]])
  }else{
    nn <- length(x)
  }

  cat("Ratio statistics mu's:\n")
  cat(paste0("NN orders: n1 = ", attr(x, "n1"), ", n2 = ",
             attr(x, "n2"), ".\n"))
  cat(paste0("Sample size: ", nn, "."))
  if (!is.null(attr(x, "upper_D"))) {
    cat(paste0("\nNominal Dimension: ", attr(x, "upper_D"), "."))
  }

  invisible(x)
  }




#' @name compute_mus
#'
#' @importFrom graphics hist curve
#'
#' @param x object of class \code{mus}, obtained from the
#' function \code{compute_mus()}.
#' @param range_d a sequence of values for which the generalized ratios density
#' is superimposed to the histogram of \code{mus}.
#' @param ... ignored.
#'
#' @export
plot.mus <- function(x, range_d = NULL, ...) {


  n1 <- attr(x, "n1")
  n2 <- attr(x, "n2")
  hist(x,
       breaks = 30,
       freq = F,
       col = "white",
       main = "", xlab = latex2exp::TeX("$\\mu$"))

  if(!is.null(range_d)){

  for (i in seq_along(range_d)) {
    y <- dgera(sort(x),
               n1 = n1,
               n2 = n2,
               d = range_d[i])
    lines(y ~ sort(x),
          col = range_d[i] - min(range_d) + 1,
          lwd = 2,
          lty = 3)
  }

  legend(
    "topright",
    legend = range_d,
    col = (range_d) - min(range_d) + 1,
    lty = 3,
    lwd = 2
  )
  }
  invisible()
}


