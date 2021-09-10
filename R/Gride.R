#' Density for the generic ratio of two NNs' distances
#'
#' @param mu.dot real positive value, the quantile to be considered.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param d real value, the ID parameter.
#' @param log logical, if \cite{TRUE}, the log-density is returned.
#'
#' @return the value of the (log-)density in the specified quantile.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' d_mudot(x=4, n1 = 1, n2 = 2, d=4)
#' }
#'
d_mudot <- function(mu.dot, n1, n2, d, log=F){

  if(n2<n1){
    stop("n2 is lower than n1")
  }

  d_n12    <- (n2-n1)
  if(d_n12==1){
    ch       <- n1
    lognum   <- log(d)
    logden   <- ((n2-1)*d+1) * log(mu.dot)
    log_dens <- log(ch) +  lognum - logden + log(mu.dot>1)
  }else{
    logch    <- sum(log(1:n2))-sum(log(1:n1))-sum(log(1:(n2-n1)))
    lognum   <- log(d) + (d_n12 - 1) * log(mu.dot ^ d - 1)
    logden   <- ((n2 - 1) * d + 1) * log(mu.dot)
    log_dens <- log(d_n12) + logch +  lognum - logden + log(mu.dot > 1)
  }
  if(log){
    log_dens
  }else{
    exp(log_dens)
  }
}


#' Log-likelihood generated from the generic ratio of two NN's distances
#'
#' @param par the ID, variable of the objective function,
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param mu_dots the ratios of distance of generic order.
#'
#'
#' @return the value of the log-likelihood evaluated in \code{par}.
#'
#' @examples
#' \dontrun{
#' loglikelihood_gride(d_value,n1=1,n2=2,mudotsample)
#' }
#'
loglikelihood_gride_R <- function(par, n1, n2, mu_dots){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  nn       <- length(mu_dots)
  d_n12    <- (n2-n1)

  if (d_n12 == 1) {
    lognum   <- nn * log(par)
    logden   <- ((n2 - 1) * par + 1) * sum(log(mu_dots))
    log_dens <- lognum - logden + sum(log(mu_dots > 1))
  } else{
    logch    <- sum(log(1:n2))-sum(log(1:n1))-sum(log(1:(n2-n1)))
    lognum   <- nn * log(par) + (d_n12 - 1) * sum(log(mu_dots ^ par - 1))
    logden   <- ((n2 - 1) * par + 1) * sum(log(mu_dots))
    log_dens <- logch * nn + nn * log(d_n12)  +
                lognum - logden + sum(log(mu_dots > 1))
  }

  log_dens
}









#' Unnormalized log-posterior for Grimes estimator on the real scale
#'
#' @param z the transformed ID value, according to \eqn{z=\log(d-1)}.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param mu_dots the sample of distance ratios \eqn{\dot{\mu}}.
#' @param a_d shape parameter of the Gamma prior distribution for \eqn{d}.
#' @param b_d rate parameter of the Gamma prior distribution for \eqn{d}.
#'
#' @return the value of the unnormalized log-posterior distribution evaluated in \eqn{z}.
#'
#' @examples
#' \dontrun{
#' unn_posterior_gride_real(z=3, n1=1, n2=2, mu_dots=mu_sample)
#' }
#'
unn_posterior_gride_real <- function(z, n1, n2, mu_dots, a_d=1, b_d=1){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  par      <- exp(z)+1
  nn       <- length(mu_dots)
  d_n12    <- (n2-n1)

  if (d_n12 == 1) {
    ch       <- n1
    lognum   <- nn * log(par)
    logden   <- ((n2 - 1) * par + 1) * sum(log(mu_dots))
    log_dens <- lognum - logden + sum(log(mu_dots > 1))
  } else{
    logch    <- sum(log(1:n2))-sum(log(1:n1))-sum(log(1:(n2-n1)))
    lognum   <- nn * log(par) + (d_n12 - 1) * sum(log(mu_dots ^ par - 1))
    logden   <- ((n2 - 1) * par + 1) * sum(log(mu_dots))
    log_dens <- logch * nn + nn * log(d_n12)  +
      lognum - logden + sum(log(mu_dots > 1))
  }

  log_dens + stats::dgamma(par,shape = a_d,rate = b_d,log = T) + z
}



#' Posterior Metropolis Hasting sampler for the Grimes model
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param nsim the number of MCMC iterations to collect in the sample.
#' @param burn_in the number of iterations to discard from the MCMC sample.
#' @param sigma standard deviation of the Gaussian proposal used in the MH step.
#' @param start_d initial value for the MCMC chain. If \eqn{NULL}, the MLE is used.
#' @param a_d shape parameter of the Gamma prior distribution for \eqn{d}.
#' @param b_d rate parameter of the Gamma prior distribution for \eqn{d}.
#'
#' @return a list with the coda MCMC sample and the traceplot.
#' @export
#'
#' @import ggplot2
#' @importFrom coda mcmc
#'
#' @examples
#' \dontrun{
#' posterior_sampler_gride(dataset, n1 = 1, n2 = 2)
#' }
posterior_sampler_gride <- function(X,
                                   DM = NULL,
                                   n1 = 1,
                                   n2 = 2,
                                   nsim = 5000,
                                   burn_in = 2000,
                                   sigma = .5,
                                   start_d = NULL,
                                   a_d = 1,
                                   b_d = 1) {

  if(n2<n1){
    stop("n2 is lower than n1")
  }

  mu    <- generate_mudots(
    X = X,
    n1 = n1,
    n2 = n2,
    DM = DM,
    dupl_remove = T
  )

  PAR <- numeric(nsim + burn_in)

  if (is.null(start_d)) {
    start_d <- mle_gride_opt(X = X, n1 = n1, n2 = n2)
    start_d <- log(start_d-1)
    PAR[1]  <- start_d
  } else{
    if(start_d<=1){stop("Invalid starting point: it has to be > 1")}
    start_d  <- log(start_d-1)
    PAR[1]   <- start_d
  }

  for (i in 2:(nsim + burn_in)) {
    oldpar <- PAR[i - 1]
    newpar <- stats::rnorm(1, oldpar, sigma)
    alpha  <- unn_posterior_gride_real(
      z = newpar,
      n1 = n1,
      n2 = n2,
      mu_dots = mu,
      a_d = a_d,
      b_d = b_d
    ) -
      unn_posterior_gride_real(
        z = oldpar,
        n1 = n1,
        n2 = n2,
        mu_dots = mu,
        a_d = a_d,
        b_d = b_d
      )
    if (log(stats::runif(1)) < alpha) {
      PAR[i] <- newpar
    } else{
      PAR[i] <- oldpar
    }
  }

  PAR        <- PAR[-c(1:burn_in)]
  MCMCSample <- coda::mcmc(exp(PAR) + 1)

  class(MCMCSample) <- "mcmc_gride"
  return(MCMCSample)
}


#' Gride estimator for ID under the Bayesian framework
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param DM a distance matrix describing the nearest neighbors structure. Default is \code{NULL}.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param sigma standard deviation of the Gaussian proposal used in the MH step.
#' @param a_d shape parameter of the Gamma prior distribution for \eqn{d}.
#' @param b_d rate parameter of the Gamma prior distribution for \eqn{d}.
#' @param alpha the amount of probability contained in the posterior credible interval.
#' @param nsim the number of MCMC iterations to collect in the sample.
#' @param burn_in the number of iterations to discard from the MCMC sample.
#'
#' @return the vector of summary statistics from the posterior MCMC sample.
#' @export
#'
#' @importFrom hBayesDM estimate_mode
#'
#'
#' @examples
#' \dontrun{
#' bayesfit_gride(X=dataset, n1 = 4, n2 = 8)
#' }
bayesfit_gride <- function(X,
                           DM=NULL,
                           n1,n2,
                           sigma=.5,
                           a_d = 1,
                           b_d = 1,
                           alpha=.95,
                           nsim = 5000,
                           burn_in = 2000){

  if(n2<n1){
    stop("n2 is lower than n1")
  }


  sam <- posterior_sampler_gride(X = X,DM = DM,
                                n1 = n1,n2 = n2,sigma = sigma,
                                nsim = nsim,burn_in = burn_in,
                                start_d = NULL)

  Res       <- numeric(5)
  Res[2]    <- mean(sam)
  Res[3]    <- stats::median(sam)
  Res[4]    <- hBayesDM::estimate_mode(sam)
  Res[c(1,5)] <- stats::quantile(sam,c((1-alpha)/2,(1+alpha)/2))
  names(Res) <- c("Lower Bound", "Mean", "Median", "Mode","Upper Bound")
  Res <- list(Estimates = Res, mcmc = sam)
  class(Res) <- "byf_gride"
  return(Res)
}







#' Bootstrap sample generator for the Gride maximum likelihood estimate of the intrinsic dimension
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param nsim integer,the number of bootstrap replications to consider.
#'
#' @return a matrix of dimension \code{n} \eqn{\times} \code{nsim} with the bootstrap samples organized by column.
#' @export
#'
#' @examples
#' \dontrun{
#' Bootstrap_gride(Data,2,5)
#' }
bootstrap_gride <- function(X,
                            n1,n2,
                            nsim=2000){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  n  <- nrow(X)
  mle.est <- mle_gride_opt(X,n1 = n1,n2 = n2)
  SAMPL <- replicate(nsim,simulate_mudot(nsim = n,n2 = n2,n1 = n1,mle.est))
  bs <- apply(SAMPL,
              2,
              function(x)  stats::optimize( loglikelihood_gride,
                                            interval = c(0.01, ncol(X) + 1),
                                            n1 = n1,
                                            n2 = n2,
                                            mu_dots = x,
                                            maximum = T)$max)
  return(bs)
}


#' Maximum Likelihood Estimator (MLE) for the ID using generic ratios
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#'
#' @return the Gride MLE estimate obtained via numeric optimization.
#' @export
#'
#' @examples
#' \dontrun{
#' mle_gride_opt(dataset,n1=1,n2=2)
#' }
#'
mle_gride_opt <- function(X,n1=2,n2=2){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  mudot <- generate_mudots(X,n1 = n1,n2 = n2)
  O <- stats::optimize( loglikelihood_gride,
                        interval = c(0.01, ncol(X) + 1),
                        n1 = n1,
                        n2 = n2,
                        mu_dots = mudot,
                        maximum = T
  )
  return(c(mle=O$max))

}

#' Maximum Likelihood Estimator (MLE) for the ID using generic ratios
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param n1 order of the NN according to the first distance is computed.
#' @param n2 order of the NN according to the second distance is computed.
#' @param nsim the number of bootstrap samples to consider.
#' @param conf_lev the confidence level for the computation of the confidence interval.
#'
#' @return the MLE estimate obtained via numeric optimization along with the bootstrap confidence interval.
#' @export
#' @examples
#' \dontrun{
#' mle_gride(dataset,n1=1,n2=2)
#' }
#'
mle_gride <- function(X,n1=1,
                      n2=2, nsim=2000, conf_lev=.95){
  if(n2<n1){
    stop("n2 is lower than n1")
  }
  n  <- nrow(X)
  mle.est <- mle_gride_opt(X,n1,n2)
  alpha <- 1-conf_lev
  bs <- bootstrap_gride(X = X,
                  n1 = n1,n2 = n2,
                  nsim=nsim)
  qq  <- stats::quantile(bs-mle.est,probs = c(alpha/2,1-alpha/2))
  Res <- c(mle-qq[2],mle.est,mle-qq[1])
  names(Res) <- c("Lower Bound", "Estimate","Upper Bound")
  return(Res)
}




#' Plateau function
#'
#' @param X a \code{n} \eqn{times} \code{D} data matrix.
#' @param n.max integer,
#' @param by integer, step
#'
#' @return a matrix, with id estimates and average distances
#' @export
#'
#' @examples
#' \dontrun{
#' plateau(dataset,n.max=100,by=2)
#' }
plateau <- function(X,n.max,by=2){

  K    <- FNN::get.knn(X,k=n.max)
  inds <- round(seq(by,n.max,by = by))
  nind <- length(inds)
  IND  <- cbind(inds/by,inds)
  res  <- matrix(NA, nind, 3)
  res  <- cbind(res, IND)
  for(i in 1:nind){

    mudots <- K$nn.dist[,IND[i,2]]/K$nn.dist[,IND[i,1]]
    res[i,3] <- mean(K$nn.dist[,IND[i,2]])
    res[i,2] <- mean(K$nn.dist[,IND[i,1]])
    res[i,1] <- stats::optimize( loglikelihood_gride,
                                 interval = c(0.01, ncol(X)+5 ),
                                 n1 = IND[i,1],
                                 n2 = IND[i,2],
                                 mu_dots = mudots,
                                 maximum = T)$max
  }
  return(res)

}


