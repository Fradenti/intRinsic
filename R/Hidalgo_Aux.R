#' Preprocess the data for the Hidalgo model
#'
#' @param X  a \code{n} \eqn{`times} \code{D} data matrix.
#' @param q  integer, the local homogeneity parameter.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#'
#' @return A list containing
#'       \itemize{
#'       \item \code{mu} the vector of the ratios NN/1NN.
#'       \item \code{Nq} the matrix \eqn{N^q}.
#'       \item \code{q} the provided local homogeneity parameter.
#'       }
#'
#' @export
Hidalgo_data_preprocessing <- function(X = NULL, q, DM = NULL) {
  if (is.null(DM)) {
    check <- (duplicated(X))
    if (sum(check) > 0) {
      cat(paste(
        "Duplicates are present. Removing the following observations:\n",
        which(check)
      ), "\n")
      X <- X[-which(check), ]
    }
    K <- FNN::get.knn(X, k = q)
    mu <- K$nn.dist[, 2] / K$nn.dist[, 1]
    n <- length(mu)
    NQ <- matrix(0, n, n)
    for (x in 1:n) {
      NQ[x, K$nn.index[x, ]] <- 1
    }
  } else {
    sDistMat <- apply(DM, 1, function(x) sort(x, index = T))
    mu <- purrr::map_dbl(sDistMat, ~ (.x$x)[3] / (.x$x)[2])
    for (x in 1:n) {
      NQ[x, (sDistMat[[x]]$ix)[2:q]] <- 1
    }
  }
  return(list(mu = mu, Nq = NQ, q = q))
}

#' Postprocess the \code{id} chains to disentangle the label switching
#'
#' @param output the output of the \code{Hidalgo} function.
#' @param all_chains : logical, if TRUE all the chains, after the disentanglement for label switching, are reported.  If \code{FALSE}, a dataset with summary statistics is returned instead.
#'
#' @return a dataset \code{nsim} \eqn{times} \code{n}, with the observation specific chain of the id.
#' @export
#'
Hidalgo_postpr_chains <- function(output, all_chains = F) {
  nsim <- dim(output$membership_labels)[1]
  n <- dim(output$membership_labels)[2]
  Res <- t(sapply(1:nsim, function(i) output$intrinsic_dimension[i, output$membership_labels[i, ]]))
  if (all_chains == T) {
    Res <- list(all_MCMC = Res)
    class(Res) <- "hid_all_mcmc"
  } else {
    avg <- colMeans(Res)
    med <- robustbase::colMedians(Res)
    sd1 <- apply((Res), 2, function(z) stats::quantile(z, .05))
    sd2 <- apply((Res), 2, function(z) stats::quantile(z, .25))
    sd3 <- apply((Res), 2, function(z) stats::quantile(z, .75))
    sd4 <- apply((Res), 2, function(z) stats::quantile(z, .95))
    Res <- dplyr::tibble(Q.05 = sd1, Q.25 = sd2,
                         MEAN = avg, MEDIAN = med,
                         Q.75 = sd3, Q.95 = sd4, OBS = 1:n)
    Res <- list(ID_summary = Res)
    class(Res) <- "hid_sum_mcmc"
  }
  return(Res)
}

#' Compute mean and median ID of the postprocessed chains stratified by an external factor
#'
#' @param output the output of the \code{Hidalgo} function.
#' @param class a factor according to the observations can be stratified by.
#' @param avg logical, if \code{TRUE} the posterior average is computed for each observation's chain. Otherwise, the median is used.
#'
#' @return a tibble containing the \code{id} estimates stratified by class.
#' @export
Hidalgo_ID_class <- function(output, class, avg = T) {
  class <- factor(class)

  REV <-
    Hidalgo_postpr_chains(
      output = output,
      all_chains = F
    )$ID_summary
  if (avg) {
    D <- dplyr::mutate(dplyr::select(REV, .data$MEAN), X = .data$MEAN, G = class)
  } else {
    D <- dplyr::mutate(dplyr::select(REV, .data$MEDIAN), X = .data$MEDIAN, G = class)
  }
  Res <- list(ID_perClass = D)
  class(Res) <- "hid_class"
  return(Res)
}

#' Compute the Posterior Coclustering Matrix
#'
#' @param object the output of the \code{Hidalgo} function.
#' @param estimate_clustering logical, if \code{TRUE}, the optimal clustering based on the minimization of the Variation of Information is computed. Can be slow with more than 1000 data points.
#' @param greed logical, if \code{TRUE} the greedy version of the VI algorithm is used. Use this option when the dataset is in the order of the hundreds
#'
#'
#' @return a list containing the posterior coclustering matrix and, if requested, the optimal partition.
#' @export
#'
#' @examples
#' \dontrun{
#' Hidalgo_coclustering_matrix(output, class = group_variable, estimate_clustering = F, plot = F)
#' }
Hidalgo_coclustering_matrix <- function(object,
                                        estimate_clustering = TRUE,
                                        greed = F) {
  psm <- PSM(object$membership_labels)
  cl <- NULL
  n <- nrow(psm)
  if (estimate_clustering & greed) {
    cl <- factor(mcclust.ext::minVI(psm, method = "greedy")$cl)
  } else if (estimate_clustering) {
    cl <- factor(mcclust.ext::minVI(psm)$cl)
  }
  Res <- list(psm = psm, optimalCL = cl)
  class(Res) <- "hid_psm"
  return(Res)
}


