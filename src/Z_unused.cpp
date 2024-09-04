/*
//////////////////////////////////////////////////////////////////////////////////////////////
  // [[Rcpp::export]]
double rtgamma_once(double shape, double rate, double lb, double ub) {
  // Draws from a truncated Gamma(shape, rate) distribution on (lb, ub)
  // weprovide the shape and rate parameters of the gamma
  // these are converted internally to scale parameters used by Rf_{dpq}gamma
  double lub = Rf_pgamma(lb, shape, 1.0/rate, 1,0);
  double uub = Rf_pgamma(ub, shape, 1.0/rate, 1,0);
  double u = Rf_runif(lub, uub);
  double draw = Rf_qgamma(u, shape, 1.0/rate, 1,0);
  return draw;
}
// [[Rcpp::export]]
arma::colvec UPD_d_TRUNC( List AIFD,  double a0, double b0, int K, double D){

  arma::colvec D_new(K);

  arma::colvec LOGSUM = as<arma::colvec>(AIFD["SL0"]);
  arma::colvec n10    = as<arma::colvec>(AIFD["nl0"]);
  arma::colvec a_star = a0 + n10;
  arma::colvec b_star = b0 + LOGSUM;
  for( int l=0; l<K; l++){
    D_new[l] = rtgamma_once(a_star[l], b_star[l], 0.0, D);
  }
  return(D_new);
}
// [[Rcpp::export]]
arma::colvec UPD_d_TRUNC_MASS( List AIFD,  double a0, double b0, int K, double D, double piMass){

  arma::colvec D_new(K);
  arma::colvec LOGSUM = as<arma::colvec>(AIFD["SL0"]);
  arma::colvec n10    = as<arma::colvec>(AIFD["nl0"]);
  arma::colvec a_star = a0 + n10;
  arma::colvec b_star = b0 + LOGSUM;
  arma::colvec lp(2);
  for( int l=0; l<K; l++){
    //non tengo conto che è troncata nel peso!
      double log_mx  = (::lgamma(a_star[l])) - (::lgamma(a0)) +
        a0*log(b0) - (a_star[l]) * log(b_star[l]) - LOGSUM[l];
      lp(0) = (log_mx)      + log(piMass)  ;
      lp(1) = log(1-piMass) + n10[l] * log(D) - (D)*LOGSUM[l]   ;
      lp    = exp(lp-max(lp));
      double  U = R::runif(0,1);
      double totpeso = lp(0)/accu(lp);

      // Rcout << totpeso << "..." <<    "\n";
      if( U < totpeso){
        D_new[l] = rtgamma_once(a_star[l], b_star[l], 0.0, D);
      }else{
        D_new[l] = D;
      }
  }
  return(D_new);
}

*/
