#' \code{Gride}: the Generalized Ratios ID Estimator
#'
#' The function can fit the Generalized ratios ID estimator under both the
#' frequentist and the Bayesian frameworks, depending on the specification of
#' the argument \code{method}. The model is the direct extension of the
#' \code{TWO-NN} method presented in
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al., 2017}
#' . See also \href{https://arxiv.org/abs/2104.13832}{Denti et al., 2021+} \
#' for more details.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param dist_mat distance matrix computed between the \code{n} observations.
#' @param mus_n1_n2 vector of generalized order NN distance ratios.
#' @param method the chosen estimation method. It can be
#' \describe{
#' \item{\code{"mle"}}{maximum likelihood estimation;}
#' \item{\code{"bayes"}}{estimation with the Bayesian approach.}
#' }
#' @param n1 order of the first NN considered. Default is 1.
#' @param n2 order of the second NN considered. Default is 2.
#' @param alpha confidence level (for \code{mle}) or posterior probability in
#' the credible interval (\code{bayes}).
#' @param upper_D nominal dimension of the dataset (upper bound for the
#' maximization routine).
#' @param nsim number of bootstrap samples or posterior simulation to consider.
#' @param burn_in number of iterations to discard from the MCMC sample.
#' Applicable if \code{method = "bayes"}.
#' @param sigma standard deviation of the Gaussian proposal used in the MH step.
#' Applicable if \code{method = "bayes"}.
#' @param start_d initial value for the MCMC chain. If \code{NULL},
#' the MLE is used. Applicable if \code{method = "bayes"}.
#' @param a_d shape parameter of the Gamma prior distribution for \code{d}.
#' Applicable if \code{method = "bayes"}.
#' @param b_d rate parameter of the Gamma prior distribution for \code{d}.
#' Applicable if \code{method = "bayes"}.
#' @param ... additional arguments for the different methods.
#'
#'
#' @return a list containing the \code{id} estimate obtained with the Gride
#' method, along with the relative confidence or credible interval
#' (object \code{est}). The class of the output object changes according to the
#' chosen \code{method}. Similarly,
#' the remaining elements stored in the list reports a summary of the key
#' quantities involved in the estimation process, e.g.,
#' the NN orders \code{n1} and \code{n2}.
#'
#' @export
#'
#' @references
#' Facco E, D'Errico M, Rodriguez A, Laio A (2017). "Estimating the intrinsic
#' dimension of datasets by a minimal neighborhood information."
#' Scientific Reports, 7(1), 1-8.
#' ISSN 20452322, doi: 10.1038/s41598-017-11873-y.
#'
#' Denti F, Doimo D, Laio A, Mira A (2022+). "Distributional Results for
#' Model-Based Intrinsic Dimension Estimators."
#' arXiv preprint. 2104.13832, \url{https://arxiv.org/abs/2104.13832}.
#'
#' @examples
#' \donttest{
#'  X  <- replicate(2,rnorm(500))
#'  dm <- as.matrix(dist(X,method = "manhattan"))
#'  res <- gride(X, nsim = 500)
#'  res
#'  plot(res)
#'  gride(dist_mat = dm, method = "bayes", upper_D =10,
#'  nsim = 500, burn_in = 100)
#' }
gride <- function(X = NULL,
                  dist_mat = NULL,
                  mus_n1_n2 = NULL,
                  method = c("mle", "bayes"),
                  n1 = 1,
                  n2 = 2,
                  alpha = 0.95,
                  nsim = 5000,
                  upper_D = 50,
                  burn_in = 2000,
                  sigma = 0.5,
                  start_d = NULL,
                  a_d = 1,
                  b_d = 1,
                  ...) {


  if (is.null(mus_n1_n2)) {
    if (is.null(X) & is.null(dist_mat)) {
      stop("Please provide either a dataset X or a distance matrix",
           call. = FALSE)}

      mus_n1_n2  <- compute_mus(
        X = X,
        dist_mat = dist_mat,
        n1 = n1,
        n2 = n2
        )

  }else{
    if (class(mus_n1_n2)[1] == "mus") {
      n1 <- attr(mus_n1_n2, which = "n1")
      n2 <- attr(mus_n1_n2, which = "n2")
      if(attr(mus_n1_n2, which = "upper_D") != "unknown"){
        upper_D <- attr(mus_n1_n2, which = "upper_D") + 5
      }

    }


    }


  if (!is.null(X)) {
      upper_D <- ncol(X) + 1
  }
  if (is.null(X) & is.null(upper_D)) {
    stop("Please provide the nominal dimension of the dataset D in upper_D",
         call. = FALSE)
  }

  method <- match.arg(method)
  switch(
    method,
    mle    = gride_mle(
      mus_n1_n2 = mus_n1_n2,
      n1 = n1,
      n2 = n2,
      alpha = alpha,
      upper_D = upper_D,
      nsim = nsim,
      ...
    ),
    bayes  = gride_bayes(
      mus_n1_n2 = mus_n1_n2,
      alpha = alpha,
      n1 = n1,
      n2 = n2,
      upper_D = upper_D,
      nsim = nsim,
      burn_in = burn_in,
      sigma = sigma,
      start_d = start_d,
      a_d = a_d, b_d = b_d,
      ...
    )
  )

}
