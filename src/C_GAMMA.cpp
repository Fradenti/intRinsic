#include "C_GAMMA.h"

// [[Rcpp::export]]
arma::colvec gam_trunc(int D, int K,
                       double a0_d, double b0_d,
                       arma::colvec n_l,
                       arma::colvec sLog){
  arma::colvec d(K);
  arma::colvec p_upper(K);
  arma::colvec p_unif(K);

  arma::colvec a_star = n_l + a0_d;
  arma::colvec b_star = sLog + b0_d;

  for(int k=0; k<K; k++){
    p_upper[k] = R::pgamma( D, a_star[k], 1.0/b_star[k], true, false);
    p_unif[k]  = R::runif( 0, p_upper[k]);
    d[k]       = R::qgamma(p_unif[k], a_star[k], 1.0/b_star[k], true, false);
  }
  return(d);
}

// [[Rcpp::export]]
arma::colvec gam_trunc_pmass(int D, int K,
                             double a0_d, double b0_d,
                             arma::colvec n_l,
                             arma::colvec sLog,
                             double pi_mass){

  arma::colvec d(K);
  arma::colvec p_upper(K);
  arma::colvec p_unif(K);
  arma::colvec log_prob(2), prob(2);
  arma::colvec a_star = n_l + a0_d;
  arma::colvec b_star = sLog + b0_d;

  arma::vec dtrunc = gam_trunc(D, K, a0_d, b0_d,
                               n_l, sLog);

  for(int ll=0; ll<K; ll++) {
    log_prob[0] =log(pi_mass) + n_l[ll] * log(D) - D * sLog[ll];
    log_prob[1] =
        log(1 - pi_mass) +
        a0_d * log(b0_d) -
        lgamma(a0_d) +
        lgamma(a_star[ll]) -
        a_star[ll] * log(b_star[ll]); //+
      // R::pgamma(D, a_star[ll], b_star[ll], true, true) -
      // R::pgamma(D, a0_d, b0_d, true, true);
      // not needed, as I am sampling from a truncated gamma
      // (this ratio is part of the distribution!)

    prob = exp(log_prob) / arma::accu(exp(log_prob));
  if(R::runif(0,1)<prob[0]){
    d[ll] = D;
  }else{
    d[ll] = dtrunc[ll];
  }
  }
  return(d);
  }
