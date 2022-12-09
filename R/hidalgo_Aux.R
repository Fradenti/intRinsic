#' Auxiliary functions for the \code{Hidalgo} model
#'
#' Collection of functions used to extract meaningful information from the object returned
#' by the function \code{Hidalgo}
#'
#' @name auxHidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#'
#' @return \code{posterior_mean} returns the observation-specific \code{id} posterior means estimated with \code{Hidalgo}.
#'
#' @export
posterior_means <- function(x){
  if (class(x)[1] != "Hidalgo") {
    stop("object is not of class 'Hidalgo'", call. = FALSE)
  }
  return(x$id_summary$MEAN)

}

#' @name auxHidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#'
#' @return \code{initial_values} returns a list with the parameter specification
#' passed to the model.
#'
#' @export
initial_values <- function(x){
  if (class(x)[1] != "Hidalgo") {
    stop("object is not of class 'Hidalgo'", call. = FALSE)
  }
  return(x$recap)

}


#' @name auxHidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#'
#' @return \code{posterior_median} returns the observation-specific \code{id} posterior medians estimated with  \code{Hidalgo}.
#'
#' @export
posterior_medians <- function(x){
  if (class(x)[1] != "Hidalgo") {
    stop("object is not of class 'Hidalgo'", call. = FALSE)
  }
  out <- x$id_summary$MEDIAN
  return(out)

}


#' @name auxHidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param alpha posterior probability contained in the computed credible
#' interval.
#'
#' @return \code{credible_interval} returns the observation-specific credible intervals for a specific
#' probability \code{alpha}.
#'
#' @export
credible_intervals <- function(x, alpha = .95){
  if (class(x)[1] != "Hidalgo") {
    stop("object is not of class 'Hidalgo'", call. = FALSE)
  }
  out <- t(apply(x$id_postpr,2,function(x) stats::quantile(x, probs = c((1 - alpha) / 2, (1 + alpha) / 2))))
  return(out)
}
