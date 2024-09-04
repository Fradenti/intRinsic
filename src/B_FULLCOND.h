#ifndef B_FULLCOND
#define B_FULLCOND

#include <RcppArmadillo.h>
#include "A_AUX.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

arma::colvec  Update_memberships_faster(arma::colvec mu_obser,
                                        arma::colvec dl,
                                        arma::colvec pl,
                                        int K, int N,int q,
                                        arma::colvec possible_label,
                                        arma::colvec Ci, double QQ,
                                        arma::umat index_row,
                                        Rcpp::List index_col,
                                        arma::colvec log_Precomp_Z,
                                        arma::colvec log_Precomp_ratios);
#endif
