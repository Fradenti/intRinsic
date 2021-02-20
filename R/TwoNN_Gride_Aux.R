#' Computes the ratio statistics \eqn{\mu} needed for the Intrinsic Dimension estimation
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param dupl_remove logical, should duplicates in the data be removed automatically?
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return a vector of length \eqn{n}.
#' @export
generate_mus <- function(X, dupl_remove = T, DM = NULL) {
  check <- (duplicated(X))
  if (sum(check) > 0) {
    if (dupl_remove == F) {
      cat(
        "\nDuplicates are present in the dataset.\nI cannot proceed!\nThe problematic indexes are:\n"
      )
      return(which(check))
    } else {
      cat("\nDuplicates are present. Removing them.\n")
      X <- X[-which(check),]
    }
  }
  if (is.null(DM)) {
    K <- FNN::get.knn(X, k = 2)
    res <- K$nn.dist[, 2] / K$nn.dist[, 1]
  } else {
    sDistMat <- apply(DM, 1, function(x)
      sort(x, index = T))
    res <- purrr::map_dbl(sDistMat, ~ (.x$x)[3] / (.x$x)[2])
  }
  return(res)
}

#' Random sampler for observations generated from a Beta prime distribution
#'
#' @param nsim integer, the number of observations to generate.
#' @param a the first shape parameter of the distribution.
#' @param b the second shape parameter of the distribution.
#'
#' @return a vector of length \eqn{nsim} sampled from a Beta prime distribution.
#' @export
#'
#' @examples
#' Y <- beta_prime(100,3,2)
#'
beta_prime <- function(nsim,a,b){

  x <- stats::rgamma(nsim,a,1)
  y <- stats::rgamma(nsim,b,1)

  x/y
}

#' Random sampler for the generalized ratio distribution of \eqn{\dot{\mu}}
#'
#' @param nsim integer, the number of observations to generate.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param d the value of the intrinsic dimension.
#'
#' @return a vector of random observations sampled from the generalized ratio distribution.
#' @export
#'
#' @examples
#' simulate_mudot(100,3,5,2)
simulate_mudot <- function(nsim,n1=1,n2=2,d){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  a <- n2-n1
  b <- n1
  z <- beta_prime(nsim,a,b)
  (z+1)^(1/d)
}





#' Ratio statistics \eqn{\dot{\mu}} needed for the Intrinsic Dimension estimation
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param dupl_remove logical, should duplicates in the data be removed automatically?
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return a vector of length \eqn{n}.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_mudots(dataset, n1 = 1, n2 = 2)
#' }
generate_mudots <- function(X, n1, n2, dupl_remove = T, DM = NULL){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  check      <- (duplicated(X))
  if(sum(check)>0){
    if(dupl_remove==F){
      stop("\nDuplicates are present in the dataset.")
    }else{
      cat("Duplicates are present. Removing them.")
      X <- X[ - which(check) , ]
    }
  }
  if(is.null(DM)){
    K   <- FNN::get.knn(X,k=n2)
    res <- K$nn.dist[,n2]/K$nn.dist[,n1]
  }else{
    sDistMat <- apply(DM,1,function(x) sort(x,index=T))
    res <- purrr::map_dbl(sDistMat, ~ c((.x$x)[n2]/(.x$x)[n1]))
  }
  return(res)
}
