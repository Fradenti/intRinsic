#ifndef C_GAMMA
#define C_GAMMA

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::colvec gam_trunc(int D, int K,
                       double a0_d, double b0_d,
                       arma::colvec n_l,
                       arma::colvec sLog);

arma::colvec gam_trunc_pmass(int D, int K,
                             double a0_d, double b0_d,
                             arma::colvec n_l,
                             arma::colvec sLog,
                             double pi_mass);
#endif
