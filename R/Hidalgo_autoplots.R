#' Plot the intrinsic dimension chains
#'
#' @param object an object of class \code{Hidalgo}, the output of the \code{Hidalgo} function.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.Hidalgo
#'
#' @return Plots the MCMC of the different intrinsic dimensions.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @export
autoplot.Hidalgo <- function(object, ...) {
  ID  <- object$intrinsic_dimension
  cmm <- (apply(ID, 2, function(x) cumsum(x) / seq_along(x)))
  D   <- reshape2::melt(ID)
  D1  <- reshape2::melt(cmm)

  mycolors <- colorRampPalette(brewer.pal(9, "Blues")[6:9])(ncol(ID))

  ggplot() +
    geom_line(data = D, aes(
      x = .data$Var1,
      y = .data$value,
      group = .data$Var2
    ), col = "gray", alpha = .5) +
    theme_bw() +
    ylab("Raw MCMC - Intrinsic Dimension") +
    xlab("MCMC Iteration") +
    geom_line(
      data = D1, aes(
        x = .data$Var1,
        y = .data$value,
        group = .data$Var2,
        col = factor(.data$Var2)
      ),
      alpha = 1, lwd = 1
    ) +
    scale_color_manual(values = mycolors) +
    theme(legend.position = "none")
}

#' Plot the postprocessed MCMC trace plots
#'
#' @param object an object of class \code{hid_all_mcmc}, the output of the \code{Hidalgo_postpr_chains} function when all_chains is \code{TRUE}.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.hid_all_mcmc
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @return the trace plots of the intrinsic dimension  each observation after the postprocessing
#' @export
#'
autoplot.hid_all_mcmc <- function(object, ...) {
  n <- nrow(object$all_MCMC)
  nd <- ncol(object$all_MCMC)
  Iteration <- 1:n

  Melted <- reshape2::melt(object$all_MCMC)


  mycolors <- colorRampPalette(brewer.pal(9, "Blues")[2:9])(nd)

  qplot(
    data = Melted, x = rep(Iteration, nd),
    y = .data$value,
    group = .data$Var2, lwd=I(1),
    geom = "line", col = factor(.data$Var2)
  ) +
    theme_bw() + theme(legend.position = "none") +
    xlab("MCMC Iteration") + ylab("Intrinsic Dimension") +
    scale_color_manual(values = mycolors) # scale_color_brewer(palette = "Blues")
}

#' Plot mean and median of the \code{id} for each observation after postprocessing
#'
#' @param object an object of class \code{hid_sum_mcmc}, the output of the \code{Hidalgo_postpr_chains} function when all_chains is \code{FALSE}.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.hid_sum_mcmc
#'
#' @import ggplot2
#'
#' @return a plot displaying median and mean for each observations computed over the postprocessed \code{id} chains.
#' @export
#'
autoplot.hid_sum_mcmc <- function(object, ...) {
  a <- object$ID_summary
  aa <- ggplot() +
    geom_segment(
      data = a,
      aes(
        x = .data$OBS, xend = .data$OBS,
        y = .data$Q.05, yend = .data$Q.95
      ),
      col = "gray", alpha = .4
    ) +
    geom_point(data = a, aes(x = .data$OBS, y = .data$MEDIAN), col = "darkblue", alpha = 1, shape = 2) +
    geom_point(data = a, aes(x = .data$OBS, y = .data$MEAN), col = "steelblue", alpha = 1, shape = 4) +
    theme_bw() +
    xlab("Observation") +
    ylab("Intrinsic Dimension")
  aa
}

#' Plot posterior \code{id} for each observation stratified by an external factor
#'
#' @param object an object of class \code{hid_class}, the output of the \code{Hidalgo_ID_class} function.
#' @param class a factor according to the observations can be stratified by.
#' @param type character, provides the type of the plot. It must be one of the following options: histogram, density, boxplot, or violin plot.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.hid_class
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @return a plot of the \code{id} estimates stratified by class.
#' @export
autoplot.hid_class <- function(object, class, type = "histogram", ...) {
  D <- object$ID_perClass

  if (type == "histogram") {
    p1 <- ggplot() +
      geom_histogram(
        data = D,
        aes(
          x = .data$X,
          fill = .data$G
        ),
        position = "identity",
        col = 1,
        alpha = .5,
        bins = 25
      )
  } else if (type == "density") {
    p1 <- ggplot() +
      geom_density(
        data = D,
        aes(
          x = .data$X,
          fill = .data$G
        ),
        position = "identity",
        col = 1,
        alpha = .5
      )
  } else if (type == "boxplot") {
    p1 <- ggplot() +
      geom_boxplot(
        data = D,
        aes(
          x = .data$X,
          y = .data$G,
          fill = .data$G
        ),
        col = 1
      )
  } else if (type == "violin") {
    p1 <- ggplot() +
      geom_violin(
        data = D,
        aes(
          x = .data$X,
          y = .data$G,
          fill = .data$G
        ),
        col = 1
      )
  } else {
    return(cat("Invalid plot type: try again with one of these: histogram, density, boxplot, violin"))
  }

  nd <- length(unique(class))
  mycolors <- colorRampPalette(brewer.pal(9, "Blues")[2:9])(nd)

  p1 <- p1 + scale_fill_manual("Class",values = mycolors) +
    xlab("ID posterior estimate") +
    ylab("") + theme_bw()

  return(p1)
}

#' Plot the Posterior Coclustering Matrix
#'
#' @param object an object of class \code{hid_psm}, the output of the function \code{Hidalgo_coclustering_matrix}.
#' @param class factor according to the observations can be stratified by.
#' @param id_names a vector of label identifying each observation.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.hid_psm
#'
#' @return Heatmap displaying the Posterior Coclustering Matrix. Columns and rows can be ordered by an external class and/or the optimal partition
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Hidalgo_coclustering_matrix(output, class = group_variable, VI = F, plot = F)
#' }
autoplot.hid_psm <- function(object, id_names = NULL, class = NULL, ...) {
  n <- nrow(object$psm)

  if (is.null(id_names)) {
    id_names <- paste("Observation", 1:n)
  }


  if (!is.null(object$optimalCL)) {
    cl <- object$optimalCL
    cl_col <- data.frame(Cluster = cl)
    rownames(cl_col) <- id_names
    cl_col_ord <- cl_col[order(cl_col), , drop = F]
  }

  if (!is.null(class)) {
    cl_row <- data.frame(Class = class)
    rownames(cl_row) <- id_names
    cl_row_ord <- cl_row[order(cl_row), , drop = F]
  }

  estimate_clustering <- !is.null(object$optimalCL)
  lab <- n < 75

  nd <- 100
  mycolors <- colorRampPalette(brewer.pal(9, "Blues")[2:9])(nd)


  if (!is.null(class) & estimate_clustering) {
    psm_ord <- object$psm[order(cl_row), order(cl_col)]
    rownames(psm_ord) <- id_names[order(cl_row)]
    colnames(psm_ord) <- id_names[order(cl_col)]

    pheatmap::pheatmap(
      psm_ord,
      cluster_rows = F,
      cluster_cols = F,
      annotation_row = cl_row_ord,
      annotation_col = cl_col_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab,
      color = mycolors
    )
  } else if (!is.null(class) & !estimate_clustering) {
    psm_ord <- object$psm[order(cl_row), ]
    rownames(psm_ord) <- id_names[order(cl_row)]

    pheatmap::pheatmap(
      psm_ord,
      cluster_rows = F,
      cluster_cols = T,
      annotation_row = cl_row_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab,
      color = mycolors
    )
  } else if (is.null(class) & estimate_clustering) {
    psm_ord <- object$psm[, order(cl_col)]
    colnames(psm_ord) <- id_names[order(cl_col)]

    pheatmap::pheatmap(
      psm_ord,
      cluster_rows = T,
      cluster_cols = F,
      annotation_col = cl_col_ord,
      angle_col = "90",
      show_rownames = lab,
      show_colnames = lab,
      color = mycolors
    )
  } else {
    pheatmap::pheatmap(object$psm,
                       cluster_rows = T,
                       cluster_cols = T,
                       show_rownames = lab,
                       show_colnames = lab,
                       color = mycolors
    )
  }
}
