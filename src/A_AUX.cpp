#include "A_AUX.h"

// [[Rcpp::export]]
double log_Likelihood_double(double mu_obs, double d){
  return log(d) - (d+1.0) * log(mu_obs) ;
}

arma::rowvec Stratified_operations_0(arma::vec x,
                                     arma::vec col1,
                                     int val1){

  arma::uvec inds1  = find( col1 == val1 ) ;
  arma::vec  x10    = x.elem(inds1);
  arma::rowvec Res(2); Res.zeros();
  Res(0) = x10.size();
  if(Res(0) > 0){
    Res(1) = accu(log(x10));
  }
  return Res;
}

// [[Rcpp::export]]
arma::mat Groups_quantities(arma::colvec mu_obser,
                            arma::vec  Ci,
                            int K){ // performs stratified operations for every k=1,..,K
  arma::mat A(K,2);
  for(int ls = 0; ls<K; ls++){
    A.row(ls) = Stratified_operations_0(mu_obser, Ci, ls+1);;
  }
  return(A);
}

// [[Rcpp::export]]
double Norm_Constant_Z_l2( int Nzi_l, int N, double xi, int q){
  int Nb = Nzi_l-1;
  int Nw = N-Nzi_l;
  double res=0; ;
  for( int nb=0; nb<(q+1); nb++ ){
    res +=  Rf_choose(Nb, nb) * Rf_choose(Nw, q-nb) * pow(xi, nb) * pow(1-xi,q-nb);
  }
  return (res);
}

// [[Rcpp::export]]
arma::vec log_Zeta_maker(int N, double xi, int q){

  arma::vec NCZ(N);
  for(int i=0; i<N; i++){
    NCZ[i] = Norm_Constant_Z_l2(i+1, N, xi, q);
  }
  return(log(NCZ));
}


arma::mat Postprocessing_Chains( arma::mat labels, arma::mat IDs){
  int N    = labels.n_rows;
  int NSIM = labels.n_cols;
  arma::mat trackedID(N,NSIM);
  int ind;
  for(int i =0; i<N; i++){
    for(int j=0; j<NSIM; j++){
      ind = labels(i,j);
      trackedID(i,j) =  IDs(ind-1,j);
    }}
  return trackedID;
}


// [[Rcpp::export]]
Rcpp::List index_row_col(arma::mat Nq, int q, int N){

  arma::umat index_row(q,N);
  arma::field<arma::uvec> index_col(N);

  for(int i=0; i<N; i++){
    index_row.col(i) = find(Nq.row(i)==1);
    arma::uvec index_col2 = find(Nq.col(i)==1);
    index_col(i) = index_col2;
  }
  return(
      Rcpp::List::create(
      Rcpp::_["IR"] =  index_row,
      Rcpp::_["IC"] =  index_col ) );
}


// [[Rcpp::export]]
arma::colvec rdir_cpp(arma::colvec alpha){

  unsigned int K = alpha.n_elem;
  arma::colvec dirvec(K);
  double su = 0.0;
  for(int i=0; i<K; i++){
    dirvec[i] = R::rgamma(alpha[i],1.0);
    su += dirvec[i];
  }
  return(dirvec/su);
}
