#' Estimate the decimated \code{TWO-NN} evolution with halving steps or vector of
#' proportions
#'
#' The estimation of the \code{id} is related to the scale of the
#' dataset. To escape the local reach of the \code{TWO-NN} estimator,
#' \href{https://www.nature.com/articles/s41598-017-11873-y}{Facco et al. (2017)}
#' proposed to subsample the original dataset in order to induce greater
#' distances between the data points. By investigating the estimates' evolution
#' as a function of the size of the neighborhood, it is possible to obtain
#' information about the validity of the modeling assumptions and the robustness
#' of the model in the presence of noise.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param method method to use for decimation:
#' \describe{
#'   \item{\code{"steps"}}{the number of times the dataset is halved;}
#'   \item{\code{"proportion"}}{the dataset is subsampled according to a vector
#'   of proportions.}
#' }
#' @param proportions vector containing the fractions of the dataset to be
#' considered.
#' @param steps number of times the dataset is halved.
#' @param seed random seed controlling the sequence of sub-sampled observations.
#'
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
#' @return list containing the \code{TWO-NN} evolution
#' (maximum likelihood estimation and confidence intervals), the average
#' distance from the second NN, and the vector of proportions that were
#' considered. According to the chosen estimation method, it is accompanied with
#' the vector of proportions or halving steps considered.
#'
#' @export
#'
#' @seealso \code{\link{twonn}}
#'
#' @examples
#' X <- replicate(4,rnorm(1000))
#' twonn_decimated(X,,method = "proportions",
#'                 proportions = c(1,.5,.2,.1,.01))
#'
twonn_decimated <- function(X,
                            method = c("steps", "proportions"),
                            steps = 0,
                            proportions = 1,
                            seed = NULL) {
  method <- match.arg(method)

  if (steps == 0 & length(proportions) == 1) {
    if (proportions == 1)
      method <- "mle"
  }

  res <- switch(
    method,
    steps = twonn_dec_by(X = X,
                         steps = steps,
                         seed = seed),
    proportions = twonn_dec_prop(
      X = X,
      proportions = proportions,
      seed = seed
    ),
    mle = twonn(X = X,
                method = "mle")
  )
  return(res)
}
