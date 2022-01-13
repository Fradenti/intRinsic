#' Plot the raw traceplots of the \code{id} parameters
#'
#' The functions produces the traceplots of the parameters
#' \code{d_k}, for \code{k=1...K}. The ergodic means for all the chains
#' are superimposed.
#' The \code{K} chains that are plotted are not post-processed.
#' Ergo, they are subjected to label switching.
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param ... other arguments passed to specific methods.
#'
#' @importFrom rlang .data
#'
#' @seealso \code{\link{autoplot.Hidalgo}}
#' @rdname autoplot.gride_bayes
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}, which displays the
#' chains of the \code{id} parameters sampled from the mixture model.
#' @keywords Internal
#' @noRd
#'
ggHid_chains <- function(object, ...) {
  ID  <- object$id_raw
  cmm <- (apply(ID, 2, function(x)
    cumsum(x) / seq_along(x)))
  D   <- reshape2::melt(ID)
  D1  <- reshape2::melt(cmm)

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = D,
      ggplot2::aes(
        x = .data$Var1,
        y = .data$value,
        group = .data$Var2
      ),
      col = "gray",
      alpha = .2
    ) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Raw MCMC - Intrinsic Dimension") +
    ggplot2::xlab("MCMC Iteration") +
    ggplot2::geom_line(
      data = D1,
      ggplot2::aes(
        x = .data$Var1,
        y = .data$value,
        group = .data$Var2,
        col = factor(.data$Var2)
      ),
      alpha = 1,
      lwd = 1
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      title = ggplot2::element_text(size = 20),
      legend.position = "none"
    )
}


#' Plot a summary of the distributions of re-arranged chains
#'
#' The functions produces two panels, reporting the means (left) and the medians
#' (right) of the processed chains. Each observation is mapped to its own
#' intrinsic dimension value assumed at each iteration \code{t} of the MCMC,
#' denoted as \code{d(t,z_i)}. The 90% credible intervals are also depicted
#' with gray lines.
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param ... other arguments passed to specific methods.
#' @seealso \code{\link{autoplot.Hidalgo}}
#'
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}. It displays two
#' scatterplots containing the posterior mean and median \code{id} for each
#' observation, after that the MCMC has been postprocessed to handle label
#' switching.
#'
#' @keywords internal
#' @noRd
#'
ggHid_mean_median <- function(object, ...) {
  a    <- object$id_summary
  data <-  rbind(
    data.frame(
      x = a$OBS,
      low = a$Q.05,
      est = a$MEAN,
      upp = a$Q.95,
      type = "Mean"
    ),
    data.frame(
      x = a$OBS,
      low = a$Q.05,
      est = a$MEDIAN,
      upp = a$Q.95,
      type = "Median"
    )
  )
  ggplot2::ggplot(data = data) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$x,
        xend = .data$x,
        y = .data$low,
        yend = .data$upp
      ),
      col = "gray",
      alpha = .4
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.position = "none") +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$x, y = .data$est),
      col = "darkblue",
      alpha = 1,
      pch = 21
    ) +
    ggplot2::facet_wrap(~ type) +
    ggplot2::xlab("Observation") +
    ggplot2::ylab("Intrinsic Dimension")
}

#' Plot posterior \code{id} for each observation stratified by external factor
#'
#' The function produces different plots to investigate the relationship between
#' the posterior estimates of the \code{id} and an external, categorical
#' variable \code{class}.
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param class factor variable used to stratify observations according to their
#' \code{id} estimates.
#' @param class_plot_type a string indicating the preferred type of plot used
#' for the comparison. It can be:
#' \describe{
#' \item{\code{"histogram"} or \code{"density"}}{which produces overlapping
#' plots of the stratified distributions}
#' \item{\code{"boxplot"} or \code{"violin"}}{which produces side-to-side
#' boxplots or violin plots}
#' }
#' @param ... other arguments passed to specific methods.
#'
#' @seealso \code{\link{autoplot.Hidalgo}}
#' @keywords internal
#' @noRd
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}. It can be use to
#' visually study the relation between the posterior \code{id} estimates and an
#' external, categorical variable.
#' The type of plot varies according to the specification of
#' \code{class_plot_type}, and it can be either a set of boxplot or
#' violin plots, or a collection of overlapping densities or histograms.
#'
ggHid_class <- function(object,
                        class,
                        class_plot_type = c("histogram", "density",
                                            "boxplot", "violin"),
                        ...) {
  class_plot_type <- match.arg(class_plot_type)
  D <- object$id_summary
  D <- rbind(
    data.frame(
      Class = as.factor(class),
      est = D$MEAN,
      type = "Mean"
    ),
    data.frame(
      Class = as.factor(class),
      est = D$MEDIAN,
      type = "Median"
    )
  )

  if (class_plot_type == "histogram") {
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = D,
        ggplot2::aes(x = .data$est,
                     fill = .data$Class),
        position = "identity",
        col = 1,
        alpha = .5,
        bins = 25
      )
  } else if (class_plot_type == "density") {
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_density(
        data = D,
        ggplot2::aes(x = .data$est,
                     fill = .data$Class),
        position = "identity",
        col = 1,
        alpha = .5
      )
  } else if (class_plot_type == "boxplot") {
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_boxplot(
        data = D,
        ggplot2::aes(
          x = .data$est,
          y = .data$Class,
          fill = .data$Class
        ),
        col = 1
      )
  } else if (class_plot_type == "violin") {
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_violin(
        data = D,
        ggplot2::aes(
          x = .data$est,
          y = .data$Class,
          fill = .data$Class
        ),
        col = 1
      )
  }

  p1 +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20),
                   legend.position = "none") +
    ggplot2::facet_wrap(~ type) +
    ggplot2::xlab("ID posterior estimate") +
    ggplot2::ylab("")

}

#' Plot the posterior similarity matrix
#'
#' The function produces a heatmap of the posterior similarity (coclustering)
#' matrix (psm) computed from the MCMC output of the function \code{Hidalgo()}.
#' Rows and columns can be organized according to a clustering solution or to an
#' external categorical variable. To plot the coclustering matrix, the user
#' needs to load the \code{pheatmatp} package.
#'
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param psm posterior similarity matrix that can be provided directly to the
#' function.
#' @param clust clustering solution used to stratify the psm plot.
#' @param class factor variable used to stratify observations according to their
#' the \code{id} estimates.
#' @param id_names vector of label identifying each observation.
#' It is optional, and if the number of observations is above 75, it is ignored
#' for better graphical representation.
#' @param ... other arguments passed to specific methods.
#'
#' @keywords internal
#' @noRd
#'
#' @return plot produced by the function \code{\link[pheatmap]{pheatmap}}.
#' The plotted psm allows to study the clustering structure present in the data
#' estimated via the mixture model.
#'
#' @seealso \code{\link{autoplot.Hidalgo}}
#'
ggHid_psm <- function(object,
                      psm = NULL,
                      clust = NULL,
                      class = NULL,
                      id_names = NULL,
                      ...) {
  if (is.null(psm)) {
    psm <- psm_and_cluster(object)$psm
  }

  n   <- nrow(psm)

  if (is.null(id_names)) {
    id_names <- paste("Observation", 1:n)
  }


  if (!is.null(clust)) {
    cl     <- clust
    cl_col <- data.frame(Cluster = cl)
    rownames(cl_col) <- id_names
    cl_col_ord <- cl_col[order(cl_col), , drop = FALSE]
  }

  if (!is.null(class)) {
    cl_row <- data.frame(Class = class)
    rownames(cl_row) <- id_names
    cl_row_ord <- cl_row[order(cl_row), , drop = FALSE]
  }

  estimate_clustering <- !is.null(clust)
  lab <- n <= 75

  if (!is.null(class) & estimate_clustering) {
    psm_ord <- psm[order(cl_row), order(cl_col)]
    rownames(psm_ord) <- id_names[order(cl_row)]
    colnames(psm_ord) <- id_names[order(cl_col)]

    Q <- pheatmap::pheatmap(
      psm_ord,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = cl_row_ord,
      annotation_col = cl_col_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab
    )
  } else if (!is.null(class) & !estimate_clustering) {
    psm_ord <- psm[order(cl_row), ]
    rownames(psm_ord) <- id_names[order(cl_row)]

    Q <- pheatmap::pheatmap(
      psm_ord,
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      annotation_row = cl_row_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab
    )
  } else if (is.null(class) & estimate_clustering) {
    psm_ord <- psm[, order(cl_col)]
    colnames(psm_ord) <- id_names[order(cl_col)]

    Q <- pheatmap::pheatmap(
      psm_ord,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = cl_col_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab
    )
  } else {
    Q <- pheatmap::pheatmap(
      psm,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = lab,
      show_colnames = lab
    )
  }
  return(Q)
}
