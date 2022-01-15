#' Plot the output of the \code{Hidalgo} function
#'
#' Use this method without the \code{.Hidalgo} suffix and after loading the
#' \code{ggplot2} package. It produces several plots to explore the output of
#' the \code{Hidalgo} model.
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param type character that indicates the type of plot that is requested.
#' It can be:
#' \describe{
#'  \item{\code{"raw_chains"}}{plot the MCMC and the ergodic means NOT corrected
#'  for label switching;}
#'  \item{\code{"point_estimates"}}{plot the posterior mean and median of the id
#'  for each observation, after the chains are processed for label switching;}
#'  \item{\code{"class_plot"}}{plot the estimated id distributions stratified by
#'  the groups specified in the class vector;}
#'  \item{\code{"clustering"}}{plot the posterior coclustering matrix. Rows and
#'  columns can be stratified by and external class and/or a clustering
#'  solution.}
#'  }
#' @param class_plot_type if \code{type} is chosen to be \code{"class_plot"},
#' one can plot the stratified id estimates with a \code{"density"} plot or a
#' \code{"histogram"}, or using \code{"boxplots"} or \code{"violin"} plots.
#' @param class factor variable used to stratify observations according to
#' their the \code{id} estimates.
#' @param psm posterior similarity matrix containing the posterior probability
#' of coclustering.
#' @param clust vector containing the cluster membership labels.
#' @param title character string used as title of the plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.Hidalgo
#'
#' @return a \code{\link[ggplot2]{ggplot2}} object produced by the function
#' according to the \code{type} chosen.
#' More precisely, if
#' \describe{
#'  \item{\code{method = "raw_chains"}}{The functions produces the traceplots
#'  of the parameters \code{d_k}, for \code{k=1...K}.
#'  The ergodic means for all the chains are superimposed. The \code{K} chains
#'  that are plotted are not post-processed.
#'  Ergo, they are subjected to label switching;}
#'  \item{\code{method = "point_estimates"}}{The function returns two
#'  scatterplots displaying
#' the posterior mean and median \code{id} for each observation, after that the
#' MCMC has been postprocessed to handle label switching;}
#'  \item{\code{method = "class_plot"}}{The function returns a plot that can be
#'  used to visually assess the relation between the posterior \code{id}
#'  estimates and an external, categorical variable. The type of plot varies
#'  according to the specification of \code{class_plot_type}, and it can be
#'  either a set of boxplots or violin plots, or a collection of overlapping
#'  densities or histograms;}
#'  \item{\code{method = "clustering"}}{The function displays the posterior
#'  similarity matrix, to allow the study of the clustering structure present in
#'  the data estimated via the mixture model. Rows and columns can be stratified
#'  by and external class and/or a clustering structure.}
#'  }
#'
#'
#' @seealso \code{\link{Hidalgo}}
#'
#' @importFrom rlang .data
#'
#'
autoplot.Hidalgo <- function(object,
                             type = c("raw_chains",
                                      "point_estimates",
                                      "class_plot",
                                      "clustering"),
                             class_plot_type = c("histogram", "density",
                                                 "boxplot", "violin"),
                             class = NULL,
                             psm = NULL,
                             clust = NULL,
                             title = NULL,
                             ...) {
  type <- match.arg(type)

  if (type == "class_plot") {
    class_plot_type <- match.arg(class_plot_type)
    if (is.null(class))
      stop("Please provide the factor variable to stratify the id estimates")
  }

  G1 <- switch(
    type,
    raw_chains      = ggHid_chains(object),
    point_estimates = ggHid_mean_median(object),
    class_plot      = ggHid_class(object, class, class_plot_type),
    clustering      = ggHid_psm(object, psm, class, clust, ...)
  )

  if (!is.null(title)) {
    G1 <- G1 + ggplot2::ggtitle(title)
  }

  G1
}
