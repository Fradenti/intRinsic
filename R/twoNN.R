#' \code{TWO-NN} estimator
#'
#' The function can fit the two-nearest neighbor estimator within the maximum
#' likelihood and the Bayesian frameworks. Also, one can obtain the estimates
#' using least squares estimation, depending on the specification of the
#' argument \code{method}. This model has been originally presented in
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al., 2017}
#' . See also \href{https://www.nature.com/articles/s41598-022-20991-1}{Denti et al., 2022}
#' for more details.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param dist_mat distance matrix computed between the \code{n} observations.
#' @param mus vector of second to first NN distance ratios.
#' @param method chosen estimation method. It can be
#'  \describe{
#'        \item{\code{"mle"}}{for maximum likelihood estimator;}
#'        \item{\code{"linfit"}}{for estimation via the least squares approach;}
#'        \item{\code{"bayes"}}{for estimation with the Bayesian approach.}
#'               }
#' @param alpha the confidence level (for \code{mle} and least squares fit) or
#' posterior probability in the credible interval (\code{bayes}).
#' @param c_trimmed the proportion of trimmed observations.
#' @param unbiased logical, applicable when \code{method = "mle"}.
#' If \code{TRUE}, the MLE is corrected to ensure unbiasedness.
#' @param a_d shape parameter of the Gamma prior on the parameter \code{d},
#' applicable when \code{method = "bayes"}.
#' @param b_d rate parameter of the Gamma prior on the parameter \code{d},
#' applicable when \code{method = "bayes"}.
#' @param ... additional arguments for the different methods.
#'
#'
#' @return list characterized by a class type that depends on the \code{method}
#' chosen. Regardless of the \code{method}, the output list always contains the
#' object \code{est}, which provides the estimated intrinsic dimension along
#' with uncertainty quantification. The remaining objects vary with the
#' estimation method. In particular, if
#' \describe{
#' \item{\code{method = "mle"}}{the output reports the MLE and the relative
#' confidence interval;}
#' \item{\code{method = "linfit"}}{the output includes the \code{lm()} object used for the computation;}
#' \item{\code{method = "bayes"}}{the output contains the (1 + \code{alpha}) / 2 and (1 - \code{alpha}) / 2 quantiles, mean, mode, and median of the posterior distribution of \code{d}.}
#' }
#'
#' @export
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
#' @examples
#' # dataset with 1000 observations and id = 2
#' X <- replicate(2,rnorm(1000))
#' twonn(X)
#' # dataset with 1000 observations and id = 3
#' Y <- replicate(3,runif(1000))
#' #  Bayesian and least squares estimate from distance matrix
#' dm <- as.matrix(dist(Y,method = "manhattan"))
#' twonn(dist_mat = dm,method = "bayes")
#' twonn(dist_mat = dm,method = "linfit")
#'
twonn <- function(X = NULL,
                  dist_mat = NULL,
                  mus = NULL,
                  method = c("mle", "linfit", "bayes"),
                  alpha = 0.95,
                  c_trimmed = .01,
                  unbiased = TRUE,
                  a_d = 0.001,
                  b_d = 0.001,
                  ...) {

  if (!is.null(mus)) {
    if( !all(mus >= 1) ){
      stop("Detected some values in mus below 1.
           Please provide a proper vector of ratios.",
           call. = FALSE)
    }
    }

  if (is.null(mus)) {
    if (is.null(X) & is.null(dist_mat)) {
      stop("Please provide either a dataset X or a distance matrix",
           call. = FALSE)
    }

    mus          <- compute_mus(
      X,
      dist_mat = dist_mat,
      n1 = 1,
      n2 = 2,
      Nq = FALSE
    )
  }

  method <- match.arg(method)


  switch(
    method,
    mle    = twonn_mle(
      mus = mus,
      alpha = alpha,
      c_trimmed = c_trimmed,
      unbiased = unbiased,
      ...
    ),
    linfit = twonn_linfit(
      mus = mus,
      alpha = alpha,
      c_trimmed = c_trimmed,
      ...
    ),
    bayes  = twonn_bayes(
      mus = mus,
      alpha = alpha,
      c_trimmed = c_trimmed,
      a_d = a_d,
      b_d = b_d,
      ...
    )
  )


}
