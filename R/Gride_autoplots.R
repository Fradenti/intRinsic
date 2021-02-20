#' Plot the simulated MCMC chains
#'
#' @param object an object of class \code{autoplot.mcmc_gride}, the output of the \code{posterior_sampler_gride} function.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.mcmc_gride
#'
#' @return the trace plots of the posterior simulations for the Gride model
#' @export
#'
autoplot.mcmc_gride <- function(object, ...){
  mcmc <- c(object)
  cmm <- cumsum(mcmc)/seq_along(mcmc)
  G2 <- ggplot(dplyr::tibble(mcmc,ind=1:length(mcmc),cmm=cmm))+
    geom_line(aes(x=.data$ind,y=.data$mcmc),col="lightgray")+
    geom_line(aes(x=.data$ind,y=.data$cmm),col="darkblue")+
    theme_bw()+
    ylab("Intrinsic Dimension")+
    xlab("MCMC Iteration")
  G2
}

#' Plot the posterior distribution for the Gride model
#'
#' @param object an object of class \code{byf_gride}, the output of the function \code{bayesfit_gride}.
#' @param chain logical. If \code{TRUE}, the MCMC chain is reported. Otherwise, it displays the posterior density.
#' @param ... other arguments passed to specific methods.
#'
#' @rdname autoplot.byf_gride
#'
#' @import ggplot2
#'
#' @export
#'
autoplot.byf_gride <- function(object, chain = F, ...){

  sam <- object$mcmc
  Res <- object$Estimates

  if(chain){

    sam <- c(sam)
    cmm <- cumsum(sam)/seq_along(sam)
    G1 <- ggplot(dplyr::tibble(sam,ind=1:length(sam),cmm=cmm))+
      geom_line(aes(x=.data$ind,y=.data$sam),col="lightgray")+
      geom_line(aes(x=.data$ind,y=.data$cmm),col="darkblue")+
      theme_bw()+
      ylab("Intrinsic Dimension")+
      xlab("MCMC Iteration")

  }else{

    G1 <- ggplot(dplyr::as_tibble(c(sam))) +
      geom_density(aes(x=.data$value),col="darkblue",fill="lightgray")+
      xlab("Intrinsic Dimension") +
      ylab("Simulated Posterior Density") + theme_bw() +
      geom_vline(xintercept = Res,
                 lty = 2,
                 col = 2)

  }

  G1
}
