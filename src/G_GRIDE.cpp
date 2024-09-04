#include <RcppArmadillo.h>

// [[Rcpp::export]]
double gride_log_likelihood(double d,
                            int n1,
                            int n2,
                            arma::colvec mus_n1_n2){


  if(n2<n1){
    Rcpp::stop("n2 is lower than n1");
  }

  int d_n12 = (n2-n1);
  int nn    = mus_n1_n2.n_elem;
  double log_dens, lognum, logden, logB;

  arma::colvec vn1  = arma::linspace(1, n1-1, n1-1);
  arma::colvec vn2  = arma::linspace(1, n2-1, n2-1);
  arma::colvec vn12 = arma::linspace(1, d_n12-1, d_n12-1);


  if (d_n12 == 1) {
    lognum   = nn * log(d);
    logden   = ((n2 - 1) * d + 1) * arma::accu( log( mus_n1_n2 ) );
    log_dens = lognum - logden + sum( log( mus_n1_n2 > 1 ) );
  } else{
    logB     = arma::accu(log(vn1)) + arma::accu(log(vn12))-arma::accu(log(vn2));
    lognum   = nn * log(d) +
      (d_n12 - 1) * arma::accu( log( pow( mus_n1_n2, d)  - 1 ) );
    logden   = ((n2 - 1) * d + 1) * sum(log(mus_n1_n2));
    log_dens =  lognum - logden + sum(log(mus_n1_n2 > 1.)) - logB * nn;
  }

  return(log_dens);
}


// [[Rcpp::export]]
double gride_log_posterior(double z,
                           int n1,
                           int n2,
                           double a_d,
                           double b_d,
                           arma::colvec mus_n1_n2){


  if(n2 < n1){
    Rcpp::stop("n2 should be greater than n1");
  }

  double d = exp(z) + 1;
  double log_post = gride_log_likelihood(d, n1, n2, mus_n1_n2) +
    z +
    R::dgamma( d, a_d, 1/b_d, 1 );

  return(log_post);
}

// [[Rcpp::export]]
arma::colvec gride_mh_sampler(double start_d,
                              int n1,
                              int n2,
                              double a_d,
                              double b_d,
                              arma::colvec mus_n1_n2,
                              int nsim,
                              int burn_in,
                              double sigma){

  arma::colvec mh_sample(nsim);

  double oldpar, newpar, alpha;
  oldpar = start_d;
  for(int i=0; i<(nsim + burn_in); i++) {
    newpar = arma::randn(1)[0] * sigma + oldpar;
    alpha  = gride_log_posterior(
      newpar,
      n1,
      n2,
      a_d,
      b_d,
      mus_n1_n2) -
        gride_log_posterior(
          oldpar,
          n1,
          n2,
          a_d,
          b_d,
          mus_n1_n2);

    if (log(arma::randu(1)[0]) < alpha) {
      oldpar = newpar;
    }
    if(i >= burn_in){
      mh_sample(i-burn_in) = oldpar;
    }
  }
  return(mh_sample);
}

