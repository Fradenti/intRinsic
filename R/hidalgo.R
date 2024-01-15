#' Fit the \code{Hidalgo} model
#'
#' The function fits the Heterogeneous intrinsic dimension algorithm, developed
#' in Allegra et al., 2020. The model is a Bayesian mixture of Pareto
#' distribution with modified likelihood to induce homogeneity across
#' neighboring observations. The model can segment the observations into
#' multiple clusters characterized by different
#' intrinsic dimensions. This permits to capture hidden patterns in the data.
#' For more details on the algorithm, refer to
#' \href{https://www.nature.com/articles/s41598-020-72222-0}{Allegra et al., 2020}.
#' For an example of application to basketball data, see
#' \href{https://imstat.org/journals-and-publications/annals-of-applied-statistics/annals-of-applied-statistics-next-issues/}{Santos-Fernandez et al., 2021}.
#'
#' @param X data matrix with \code{n} observations and \code{D} variables.
#' @param dist_mat distance matrix computed between the \code{n} observations.
#' @param K integer, number of mixture components.
#' @param nsim number of MCMC iterations to run.
#' @param burn_in number of MCMC iterations to discard as burn-in period.
#' @param thinning integer indicating the thinning interval.
#' @param verbose logical, should the progress of the sampler be printed?
#' @param q integer, first local homogeneity parameter. Default is 3.
#' @param xi real number between 0 and 1, second local homogeneity parameter.
#' Default is 0.75.
#' @param alpha_Dirichlet parameter of the symmetric Dirichlet prior
#' on the mixture weights. Default is 0.05, inducing a sparse mixture.
#' Values that are too small (i.e., lower than 0.005) may cause underflow.
#' @param a0_d shape parameter of the Gamma prior on \code{d}.
#' @param b0_d rate parameter of the Gamma prior on \code{d}.
#' @param prior_type character, type of Gamma prior on \code{d}, can be
#' \describe{
#'    \item{\code{"Conjugate"}}{a conjugate Gamma distribution is elicited;}
#'    \item{\code{"Truncated"}}{the conjugate Gamma prior is truncated over the
#'    interval \code{(0,D)};}
#'    \item{\code{"Truncated_PointMass"}}{same as \code{"Truncated"}, but a
#'    point mass is placed on \code{D}, to allow the \code{id} to be
#'    identically equal to the nominal dimension.}
#' }
#' @param D integer, the maximal dimension of the dataset.
#' @param pi_mass probability placed a priori on \code{D} when
#' \code{Truncated_PointMass} is chosen.
#'
#'
#' @name Hidalgo
#'
#' @seealso \code{\link{id_by_class}} and \code{\link{clustering}}
#' to understand how to further postprocess the results.
#'
#' @return object of class \code{Hidalgo}, which is a list containing
#'\describe{
#'     \item{\code{cluster_prob}}{chains of the posterior mixture weights;}
#'     \item{\code{membership_labels}}{chains of the membership labels for all
#'     the observations;}
#'     \item{\code{id_raw}}{chains of the \code{K} intrinsic dimensions
#'     parameters, one per mixture component;}
#'     \item{\code{id_postpr}}{a chain for each observation, corrected for
#'     label switching;}
#'     \item{\code{id_summary}}{a matrix containing, for each observation, the
#'     value of posterior mean and the 5\%, 25\%, 50\%, 75\%, 95\% quantiles;}
#'     \item{\code{recap}}{a list with the objects and specifications passed to
#'     the function used in the estimation.}
#'  }
#' @export
#'
#' @references
#' Allegra M, Facco E, Denti F, Laio A, Mira A (2020).
#' “Data segmentation based on the local intrinsic dimension.”
#' Scientific Reports, 10(1), 1–27. ISSN 20452322,
#' \doi{10.1038/s41598-020-72222-0},
#'
#' Santos-Fernandez E, Denti F, Mengersen K, Mira A (2021).
#' “The role of intrinsic dimension in high-resolution player tracking data –
#' Insights in basketball.” Annals of Applied Statistics - Forthcoming, –
#' ISSN 2331-8422, 2002.04148, \doi{10.1038/s41598-022-20991-1}
#'
#' @examples
#' \donttest{
#' set.seed(1234)
#' X            <- replicate(5,rnorm(500))
#' X[1:250,1:2] <- 0
#' X[1:250,]    <- X[1:250,] + 4
#' oracle       <- rep(1:2,rep(250,2))
#' # this is just a short example
#' # increase the number of iterations to improve mixing and convergence
#' h_out        <- Hidalgo(X, nsim = 500, burn_in = 500)
#' plot(h_out, type =  "B")
#' id_by_class(h_out, oracle)
#' }
#'
#'
Hidalgo <- function(X  = NULL,
                    dist_mat = NULL,
                    K  = 10,
                    nsim     = 5000,
                    burn_in  = 5000,
                    thinning = 1,
                    verbose  = TRUE,
                    q   = 3,
                    xi = .75,
                    alpha_Dirichlet = .05,
                    a0_d = 1,
                    b0_d = 1,
                    prior_type = c("Conjugate",
                                   "Truncated",
                                   "Truncated_PointMass"),
                    D = NULL,
                    pi_mass = .5) {

  prior_type <- match.arg(prior_type)

  if (verbose)
    cat("Computing ratios and Nq...\n")

  if (!is.null(X)) {
    D         <- ncol(X)
    InputList <-
      compute_mus(
        X = X,
        q = q,
        dist_mat = NULL,
        Nq = TRUE
      )
    mus       <- InputList$mus
    Nq        <- InputList$NQ

  } else if (!is.null(dist_mat)) {
    if (is.null(D) & prior_type != "Conjugate") {
      cat("Please provide the maximal dimension of the dataset")
    }

    InputList <-
      compute_mus(
        X = NULL,
        q = q,
        dist_mat = dist_mat,
        Nq = TRUE
      )
    mus       <- InputList$mus
    Nq        <- InputList$NQ

  } else if (is.null(X) & is.null(dist_mat)) {
    stop("Please provide either a dataset or a distance matrix", call. = FALSE)
  }


  if (verbose)
    cat("Done! \n")
  if (verbose)
    cat("Initializing the algorithm... \n")
  n  <-  length(mus)
  QQ <- (log(xi) - log(1 - xi))
  #  Initialization of Output:
  ALL_D  <- matrix(NA_real_, K, nsim)
  ALL_Pi <- matrix(NA_real_, K, nsim)
  ALL_Ci <- matrix(NA_real_, n, nsim)
  lev    <- seq_len(K)

  IRC.list <-  index_row_col(Nq, q, n)
  rm(Nq)

  pl  <- c(.mcmc_pack_rdirichlet(1, rep(alpha_Dirichlet,K)))
  # in case alpha_Dirichlet is too low, initialize as equiprobable
  if( any(is.na(pl)) ){

    warning(paste("alpha_Dirichlet is too low and caused underflows when initializing.
Consider increasing its value in the next run."),
            call. = FALSE)

    pl <-  rep(alpha_Dirichlet,K)
  }

  Ci    <- sample(K, n, TRUE, pl)
  d     <- stats::rgamma(K, a0_d, 1 / b0_d)

  if (verbose)
    cat("Done! \n")


  if (verbose)
    cat("Computing recurring quantities to save you time... \n")
  log_Zeta <-
    log(sapply(seq_len(n), function(zz)
      Norm_Constant_Z_l2(
        Nzi_l = zz,
        N = n,
        xi = xi,
        q = q
      )))
  log_corr <-
    (log_Zeta[seq_len(n - 1)] - log_Zeta[seq(2, n)]) *  (seq_len(n - 1) - 1)
  if (verbose)
    cat("Done! \n")


  if (verbose) {
    cat("MCMC progress:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 1,
                                  max = nsim * thinning + burn_in,
                                  style = 3)
    ipbar <- 0
  }

  time1 <-  Sys.time()

  # Gibbs sampler loop

  for (sim in seq_len(nsim * thinning + burn_in)) {

    # STEP 1 - sample  indicators

    Ci <- Update_memberships_faster(
      mu_obser = mus,
      dl =  d,
      pl = pl,
      K  = K,
      N  = n,
      q  = q,
      Ci = Ci,
      QQ = QQ,
      possible_label = lev,
      index_row =  IRC.list$IR,
      index_col =  IRC.list$IC,
      log_Precomp_Z = log_Zeta,
      log_Precomp_ratios = log_corr
    )

    # STEP 2, updating Pi^*

    Ci     <- factor(Ci, levels = seq_len(K))
    N_Slog <- Groups_quantities(mu_obser = mus,
                                Ci = Ci,
                                K = K)
    nl     <- N_Slog[, 1]
    sLog   <- N_Slog[, 2]
    pl     <- .mcmc_pack_rdirichlet(1, alpha_Dirichlet + nl)

    # Step 3  --  Update d

    a.star <- nl   + a0_d
    b.star <- sLog + b0_d


    if (prior_type == "Conjugate") {
      d <- stats::rgamma(K, a.star, b.star)
    } else if (prior_type == "Truncated") {
      p_upper <-  stats::pgamma(D, a.star, b.star)
      p_unif  <-  stats::runif(K, 0, p_upper)
      d       <-  stats::qgamma(p_unif, a.star, b.star)
    } else if (prior_type == "Truncated_PointMass") {
      log1 <- log(1 - pi_mass) +
        a0_d * log(b0_d) -
        lgamma(a0_d) +
        lgamma(a.star) -
        a.star * log(b.star) +
        stats::pgamma(D, a.star, b.star, log.p = TRUE) -
        stats::pgamma(D, a0_d, b0_d, log.p = TRUE)
      log0 <- log(pi_mass)   + nl * log(D) - D * sLog
      LOGP <- cbind(log0, log1)
      RES  <-
        apply(LOGP, 1, function(x)
          sample(c(0, 1), 1, TRUE, exp(x - max(x))))
      for (ll in seq_len(K)) {
        if (RES[ll] == 1L) {
          p_upper <-  stats::pgamma(D, a.star[ll], b.star[ll])
          p_unif  <-  stats::runif(1, 0, p_upper)
          d[ll]   <-  stats::qgamma(p_unif, a.star[ll], b.star[ll])
        } else{
          d[ll] <- D
        }
      }
    }

    # Step 4  --  Store simulation post burn-in

    if (sim > burn_in && ((sim - burn_in) %% thinning == 0L)) {
      rr            <- floor((sim - burn_in) / thinning)
      ALL_Ci[, rr]   <-  Ci
      ALL_Pi[, rr]   <-  pl
      ALL_D[, rr]    <-  d
    }
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  } # End of Gibbs loop
  time2 <-  Sys.time()

  Recap <- list(
    X = X,
    dist_mat = dist_mat,
    mus = mus,
    K = K,
    nsim = nsim,
    burn_in = burn_in,
    thinning = thinning,
    q = q,
    xi = xi,
    alpha_Dirichlet = alpha_Dirichlet,
    a0_d = a0_d,
    b0_d = b0_d,
    prior_type = prior_type,
    D = D,
    pi_mass = pi_mass,
    elapsed = time2 - time1
  )

  # Step 5  --  Postprocessing of the MCMC

  if (verbose)
    cat("\nPosterior sampling complete. Postprocessing the chains...\n")
  intrinsic_dimension  <- t(ALL_D)
  cluster_prob         <- t(ALL_Pi)
  membership_labels    <- t(ALL_Ci)

  postpro_chains <-
    t(sapply(seq_len(nsim), function(i)
      intrinsic_dimension[i, membership_labels[i, ]]))

  avg     <- colMeans(postpro_chains)
  med     <- apply(postpro_chains, 2, stats::median)
  sd1     <-
    apply(postpro_chains, 2, function(z)
      stats::quantile(z, .05))
  sd2     <-
    apply(postpro_chains, 2, function(z)
      stats::quantile(z, .25))
  sd3     <-
    apply(postpro_chains, 2, function(z)
      stats::quantile(z, .75))
  sd4     <-
    apply(postpro_chains, 2, function(z)
      stats::quantile(z, .95))
  summa   <- data.frame(
    Q.05 = sd1,
    Q.25 = sd2,
    MEAN = avg,
    MEDIAN = med,
    Q.75 = sd3,
    Q.95 = sd4,
    OBS = seq_len(n)
  )
  if (verbose){
    cat("Done! \n")
    close(pbar)
  }
  output <- list(
    cluster_prob         = cluster_prob,
    membership_labels    = membership_labels,
    id_raw               = intrinsic_dimension,
    id_postpr            = postpro_chains,
    id_summary           = summa,
    recap                = Recap
  )
  structure(output, class = c("Hidalgo", class(output)))

}



#' @name Hidalgo
#'
#' @param x an object of class \code{Hidalgo}, obtained from the function
#' \code{Hidalgo()}.
#' @param ... other arguments passed to specific methods.
#'
#'
#' @export
print.Hidalgo <- function(x, ...) {
  cat("Model: Hidalgo\n")
  cat("Method: Bayesian Estimation\n")
  cat(paste0(
    "Prior d ~ Gamma(",
    x$recap[["a0_d"]],
    ", ",
    x$recap[["b0_d"]],
    "), type = ",
    x$recap[["prior_type"]],
    "\n"
  ))
  cat(
    paste0(
      "Prior on mixture weights: Dirichlet(",
      x$recap[["alpha_Dirichlet"]],
      ") with ",
      x$recap[["K"]],
      " mixture components\n"
    )
  )
  cat(
    paste0(
      "MCMC details:\nTotal iterations: ",
      x$recap[["nsim"]] + x$recap[["burn_in"]],
      ", Burn in: ",
      x$recap[["burn_in"]],
      ", Thinning: ",
      x$recap[["thinning"]],
      "\nUsed iterations: ",
      x$recap[["nsim"]],
      "\nElapsed time: ",
      round(x$recap[["elapsed"]], 4),
      " ",
      attr(x$recap[["elapsed"]], "units"),"\n"
    )
  )
  invisible(x)
}


#' @name Hidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param type character that indicates the type of plot that is requested.
#' It can be:
#' \describe{
#'  \item{\code{"A"}}{plot the MCMC and the ergodic means NOT corrected
#'  for label switching;}
#'  \item{\code{"B"}}{plot the posterior mean and median of the id
#'  for each observation, after the chains are processed for label switching;}
#'  \item{\code{"C"}}{plot the estimated id distributions stratified by
#'  the groups specified in the class vector;}
#'  }
#' @param class factor variable used to stratify observations according to
#' their the \code{id} estimates.
#' @param ... other arguments passed to specific methods.
#'
#' @importFrom graphics boxplot plot segments matplot
#'
#'
#' @export
plot.Hidalgo <- function(x,
                         type = c("A","B","C"),
                         class = NULL,
                         ...) {
  type <- match.arg(type)

  if (type == "C" & is.null(class)){
      stop("Please provide the factor variable to stratify the id estimates")
  }


  if (type == "A") {
    ID  <- x$id_raw
    cmm <- (apply(ID, 2, function(x)
      cumsum(x) / seq_along(x)))
    matplot(
      ts(x$id_raw),
      type = "l",
      col = "lightgray",
      lty = 1,
      xlab = "MCMC Iteration",
      ylab = "Raw MCMC - Intrinsic Dimension"
    )
    matplot(
      cmm,
      type = "l",
      lty = 1,
      lwd = 2,
      add = T
    )
  } else if (type == "B") {

    on.exit({par(my_par)}, add = TRUE, after = TRUE)
    my_par <- par(mfrow = c(2, 1))

    ID  <- x$id_summary
    plot(
      ID$MEAN,
      ylim = c(min(ID$Q.05), max(ID$Q.95)),
      type = "n",
      main = "Posterior Means",
      ylab = "Intrinsic Dimension",
      xlab = "Observation"
    )
    segments(
      x0 = ID$OBS,
      x1 = ID$OBS,
      y0 = ID$Q.05,
      ID$Q.95,
      col = "lightgray"
    )
    points(ID$MEAN, ylim = c(min(ID$Q.05), max(ID$Q.95)), col = "darkblue")



    plot(
      ID$MEDIAN,
      ylim = c(min(ID$Q.05), max(ID$Q.95)),
      type = "n",
      main = "Posterior Medians",
      ylab = "Intrinsic Dimension",
      xlab = "Observation"
    )
    segments(
      x0 = ID$OBS,
      x1 = ID$OBS,
      y0 = ID$Q.05,
      ID$Q.95,
      col = "lightgray"
    )
    points(ID$MEDIAN, ylim = c(min(ID$Q.05), max(ID$Q.95)), col = "darkblue")

  }else if(type == "C"){

    on.exit({par(my_par)}, add = TRUE, after = TRUE)
    my_par <- par(mfrow = c(2, 1))

    ID  <- x$id_summary
    boxplot(
      ID$MEAN ~ class,
      col  = as.numeric(class)+1,
      ylab = ("ID posterior estimate") ,
      xlab = "Class",
      main = "Posterior Means"
    )
    boxplot(
      ID$MEDIAN ~ class,
      col = as.numeric(class)+1,
      ylab = ("ID posterior estimate"),
      xlab = "Class",
      main = "Posterior Medians"
    )
  }
  invisible()

}

#' @name Hidalgo
#'
#' @param object object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param ... other arguments passed to specific methods.
#'
#' @export
summary.Hidalgo <- function(object, ... ){

  out <- object$id_summary[,1:6]
  structure(out, class = c("summary.Hidalgo",class(out)))

}

#' @name Hidalgo
#'
#' @param x object of class \code{Hidalgo}, the output of the
#' \code{Hidalgo()} function.
#' @param ... other arguments passed to specific methods.
#'
#' @export
print.summary.Hidalgo <- function(x, ...){
  cat("Hidalgo - Posterior id estimates\n")
  print(summary(as.data.frame(x)))
  cat("\n")
  invisible(x)
}


