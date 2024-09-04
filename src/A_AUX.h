#ifndef A_AUX
#define A_AUX

#include <RcppArmadillo.h>

double log_Likelihood_double(double mu_obs, double d);
arma::rowvec Stratified_operations_0(arma::vec x,
                                     arma::vec col1,
                                     int val1);
arma::mat Groups_quantities(arma::colvec mu_obser,
                            arma::vec  Ci,
                            int K);
Rcpp::List index_row_col(arma::mat Nq,
                   int q,
                   int N);

double Norm_Constant_Z_l2( int Nzi_l, int N, double xi, int q);

arma::colvec rdir_cpp(arma::colvec alpha);

  #endif

