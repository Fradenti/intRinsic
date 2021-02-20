#' Gibbs sampler for the Hidalgo model
#'
#' @param X \code{n} \eqn{times} \code{D} data matrix.
#' @param DM \code{n} \eqn{times} \code{n} distance matrix (can be alternative to \code{X}).
#' @param K integer, number of mixture components.
#' @param nsim number of MCMC iterations to run.
#' @param burn_in number of MCMC iterations to discard as burn-in period.
#' @param thinning thinning interval.
#' @param verbose logical, should the progress of the sampler be printed?
#' @param q integer, first local homogeneity parameter.
#' @param csi real between 0 and 1, second local homogeneity parameter.
#' @param alpha_Dirichlet parameter of the Dirichlet prior on the mixture weights.
#' @param a0_d shape parameter of the Gamma prior on \code{d}.
#' @param b0_d rate parameter of the Gamma prior on \code{d}.
#' @param prior_type character, type of Gamma prior on \code{d}, can be \code{Conjugate}, \code{Truncated}, or \code{Truncated_PointMass}.
#' @param D integer, the maximal dimension of the dataset.
#' @param pi_mass probability placed a priori on \code{D} when \code{Truncated_PointMass} is chosen.
#' @param recap logical, includes a list of the chosen parameters in the output.
#' @param if_coda logical, if \code{TRUE} the MCMC is returned as a coda object.
#'
#' @importFrom coda mcmc
#'
#' @return a list containing
#'     \itemize{
#'     \item \code{membership_labels} chains of the membership labels for all the observations;
#'     \item \code{cluster_prob} chains of the posterior mixture weights;
#'     \item \code{intrinsic_dimension} chains of the K intrinsic dimensions, one per mixture component.
#'     }
#' @export
Hidalgo <- function(X  = NULL,
                    DM = NULL,
                    K  = 10,
                    nsim     = 1000,
                    burn_in  = 1000,
                    thinning = 1,
                    verbose  =T,
                    q   = 3,
                    csi = .75,
                    alpha_Dirichlet = .05,
                    a0_d = 1,
                    b0_d = 1,
                    prior_type = "Conjugate",
                    D = NULL,
                    pi_mass = .5,
                    recap = T,
                    if_coda = T){

  if(!(prior_type %in% c("Conjugate",
                         "Truncated",
                         "Truncated_PointMass") )){
    stop("Please Provide a proper prior type")
  }
  if(verbose) cat("Computing ratios and Nq...\n")


  if(!is.null(X)){
    D         <- ncol(X)
    InputList <- Hidalgo_data_preprocessing(X = X, q = q, DM = NULL)
    mus       <- InputList$mu
    Nq        <- InputList$Nq
  }else if(!is.null(DM)){
    InputList <- Hidalgo_data_preprocessing(X = NULL, q = q, DM = DM)
    mus       <- InputList$mu
    Nq        <- InputList$Nq
  }else if(is.null(X) & is.null(DM)){
    stop("Dude, you need to give me something to work with!")
  }


  if(verbose) cat("Done! \n")
  if(verbose) cat("Initializing the algorithm... \n")
  n  <-  length(mus)
  QQ <- (log(csi)-log(1-csi))
  #  Initialization of Output:
  ALL_D  <- matrix(NA_real_,K,nsim)
  ALL_Pi <- matrix(NA_real_,K,nsim)
  ALL_Ci <- matrix(NA_real_,n,nsim)
  lev    <- 1:K
  # ////////////////////////////////////////////////////////////////////////////////////////////
  IRC.list = index_row_col(Nq, q, n);
  rm(Nq)
  pl    = MCMCpack::rdirichlet( K, alpha_Dirichlet)
  Ci    = sample(K,n,TRUE, pl)
  d     = stats::rgamma(K,a0_d, 1/b0_d)
  u_rep = stats::runif(1,0,1)
  if(verbose) cat("Done! \n")


  if(verbose) cat("Computing recurring quantities to save you time... \n")
  log_Zeta <- log(sapply( 1:n, function(zz) Norm_Constant_Z_l2(Nzi_l = zz,N = n,csi = csi,q = q) ))
  log_corr <- (log_Zeta[1:(n-1)]-log_Zeta[2:(n)]) * ( (1:(n-1)) - 1 )
  if(verbose) cat("You are welcome! \n")

  ############################## GIBBS SAMPLER LOOP ################################################
  if (verbose) {
    cat("MCMC progress:\n")
    utils::flush.console()
    pbar <- utils::txtProgressBar(min = 1,
                                  max = nsim*thinning + burn_in,
                                  style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }


  for(sim in 1:(nsim*thinning + burn_in)){

    ####################################################################################################################################################################
    # STEP 1 - sample  indicators
    ####################################################################################################################################################################for(j in 1:J){
    Ci = Update_memberships_faster(mu_obser = mus,
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
                                   log_Precomp_ratios = log_corr);

    ####################################################################################################################################################################
    # STEP 2, uptading Pi^*
    ####################################################################################################################################################################
    Ci     <- factor(Ci,levels = 1:K)
    N_Slog <- Groups_quantities(mu_obser = mus,Ci = Ci,K = K)
    nl     <- N_Slog[,1]
    sLog   <- N_Slog[,2]
    pl     <- MCMCpack::rdirichlet(1,alpha_Dirichlet+nl)

    ####################################################################################################################################################################
    # Step 3  --  Update d
    ####################################################################################################################################################################
    a.star <- nl   + a0_d
    b.star <- sLog + b0_d

    if(prior_type=="Conjugate"){
      d = stats::rgamma(K,a.star,b.star)
    }else if(prior_type=="Truncated"){
      p_upper <-  stats::pgamma(D,a.star,b.star)
      p_unif  <-  stats::runif(K,0,p_upper)
      d       <-  stats::qgamma(p_unif,a.star,b.star)
    }else if(prior_type=="Truncated_PointMass"){
      log1 <- log(1-pi_mass) +
        a0_d*log(b0_d)-
        lgamma(a0_d)+
        lgamma(a.star)-
        a.star*log(b.star)+
        stats::pgamma(D,a.star,b.star,log.p = T)-
        stats::pgamma(D,a0_d,b0_d,log.p = T)
      log0 <- log(pi_mass)   + nl*log(D)-D*sLog
      LOGP <- cbind(log0,log1)
      RES  <- apply(LOGP,1, function(x) sample(0:1,1,T,exp(x-max(x))))
      for(ll in 1:K){
        if(RES[ll]==1){
          p_upper <-  stats::pgamma(D,a.star[ll],b.star[ll])
          p_unif  <-  stats::runif(1,0,p_upper)
          d[ll]   <-  stats::qgamma(p_unif,a.star[ll],b.star[ll])
        }else{
          d[ll] <- D
        }}
    }


    if (sim > burn_in && ((sim - burn_in) %% thinning == 0)) {
      rr            <- floor((sim - burn_in)/thinning);
      ALL_Ci[,rr]   <-  Ci;
      ALL_Pi[,rr]   <-  pl;
      ALL_D[,rr]    <-  d;
    }
    if (verbose) {
      ipbar <- ipbar + 1
      utils::setTxtProgressBar(pbar, ipbar)
    }
  }

  Recap <- list(X=X, DM=DM,mus=mus,
                K=K,
                nsim=nsim,
                burn_in=burn_in,
                thinning=thinning,
                verbose=verbose,
                q=q, csi=csi,
                alpha_Dirichlet=alpha_Dirichlet,
                a0_d = a0_d,
                b0_d = b0_d,
                prior_type= prior_type,
                D=D,
                pi_mass=pi_mass)

  if(if_coda){
  output <- list(
    cluster_prob         = t(ALL_Pi),
    membership_labels    = t(ALL_Ci),
    intrinsic_dimension  = t(ALL_D))
  }else{
    output <- list(
      cluster_prob         = coda::mcmc(t(ALL_Pi)),
      membership_labels    = coda::mcmc(t(ALL_Ci)),
      intrinsic_dimension  = coda::mcmc(t(ALL_D)))
  }



  if(recap){
    output[["Recap"]] <- Recap
  }

  class(output) <- "Hidalgo"
  return(output)


}
