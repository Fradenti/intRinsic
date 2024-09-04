#include "B_FULLCOND.h"

// [[Rcpp::export]]
arma::colvec  Update_memberships_faster(arma::colvec mu_obser,
                                        arma::colvec dl,
                                        arma::colvec pl,
                                        int K, int N,int q,
                                        arma::colvec possible_label,
                                        arma::colvec Ci, double QQ,
                                        arma::umat index_row,
                                        Rcpp::List index_col,
                                        arma::colvec log_Precomp_Z,
                                        arma::colvec log_Precomp_ratios) {
  Rcpp::checkUserInterrupt();
  arma::colvec n_i_in_l(K), m_i_in_l(K), log_Zeta(K), log_corr(K),
  log_Nq_given_z_i(K), pp(K), pp2(K), logLik_mu_i(K), Nzi_lAUX(K);
  for(int i=0; i<N; i++){  // for every label
    arma::colvec fakeC = Ci;
    logLik_mu_i.fill(0);
    n_i_in_l.fill(0);
    m_i_in_l.fill(0);
    unsigned int Nzi_l = 0 ;
    arma::uvec IR=index_row.col(i), IC=index_col(i);
    for(int l=0; l<(K); l++){ // for every possible value taken by Ci
      logLik_mu_i(l) = log_Likelihood_double( mu_obser[i], dl[l]  );
      fakeC[i] = l+1;
      Nzi_l = accu(fakeC==(l+1));
      //////////////////////////////////////////////////////////////////////////////////
      n_i_in_l[l] = accu( fakeC.elem(IR) == (l+1) );
      m_i_in_l[l] = accu( fakeC.elem(IC) == (l+1) );
      log_Zeta[l] = log_Precomp_Z[Nzi_l-1];     //Norm_Constant_Z_l2(Nzi_l,  N, xi, q);
      log_corr[l] = log_Precomp_ratios[Nzi_l - 1];//Norm_Constant_Z_l2(Nzi_l-1,  N, xi, q)/Zeta[l];
      Nzi_lAUX[l] = Nzi_l-1;
      //////////////////////////////////////////////////////////////////////////////////
    }


    log_Nq_given_z_i =  (n_i_in_l+m_i_in_l) * QQ - log_Zeta + log_corr;
    pp = logLik_mu_i + log_Nq_given_z_i + log(pl);
    if(pp.has_nan()){
      Rcpp::Rcout << log_Nq_given_z_i <<"C\n";
      Rcpp::Rcout << log_Zeta <<"A\n";
      Rcpp::Rcout << log_corr <<"B\n";
      Rcpp::Rcout << Nzi_lAUX <<"D\n";
      Rcpp::Rcout << log_Precomp_ratios[Nzi_l-1];
    }
    // closing loop on labels to compute probabilities
    if(arma::is_finite(max(pp))){
      pp2 = exp( pp - max(pp));
      Ci(i) = RcppArmadillo::sample(
        possible_label, 1, TRUE, pp2)[0];
    }else{
      Rcpp:: Rcout << "***!!!***";
      Ci(i) = RcppArmadillo::sample(possible_label, 1, TRUE)[0];
    }
  } // closing loop on observations
  return Ci;
}
