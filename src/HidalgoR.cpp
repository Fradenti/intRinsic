#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double log_Likelihood_double(double mu_obs, double d){
  return log(d)-(d+1)*log(mu_obs) ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
arma::colvec Stratified_operations_0(arma::vec x, arma::vec col1, int val1){

  arma::uvec inds1      = find( col1 == val1 ) ;
  arma::vec  x10        = x.elem(inds1);
  arma::colvec Res(2);
  Res(0) = x10.size();
  if(Res(0) > 0){
    Res(1) = accu(log(x10));
  }else{
    Res(1) = 0.0;
  }
  return Res;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat Groups_quantities(arma::colvec mu_obser, arma::vec  Ci, int K){ // performs stratified operations for every k=1,..,K
  arma::colvec Res;
  arma::mat A(K,2);
  for(int ls = 0; ls<K; ls++){
    Res = Stratified_operations_0(mu_obser, Ci, ls+1);
    A(ls,0) = Res(0);
    A(ls,1) = Res(1);
  }
  return(A);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//' Compute the normalizing constant for the
//'
//' @param Nzi_l the number of observations assigned to each manifold
//' @param N total number of
//' @param csi the strenght of the homogenity enforcement
//' @param q the number of NN considered for each observation
//'
//'
//' @export
// [[Rcpp::export]]
double Norm_Constant_Z_l2( int Nzi_l, int N, double csi, int q){
  int Nb = Nzi_l-1;
  int Nw = N-Nzi_l;
  double res=0; ;
  for( int nb=0; nb<(q+1); nb++ ){
    res +=  Rf_choose(Nb, nb) * Rf_choose(Nw, q-nb) * pow(csi, nb) * pow(1-csi,q-nb);
  }
  return (res);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//' Summarize the Nq matrix into indexes of rows and columns
//'
//' @param Nq the binary adjacency matrix
//' @param q number of neighbors considered
//' @param N number of observations
//'
//' @export
// [[Rcpp::export]]
List index_row_col(arma::mat Nq, int q, int N){
  arma::umat index_row(q,N);
  arma::field<arma::uvec> index_col(N);
  for(int i=0; i<N; i++){
    index_row.col(i) = find(Nq.row(i)==1);
    arma::uvec index_col2 = find(Nq.col(i)==1);
    index_col(i) = index_col2;
  }
  return(List::create(
      _["IR"]    =  index_row,
      _["IC"]    =  index_col
  ));
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::colvec  Update_memberships(arma::colvec mu_obser,
  arma::colvec dl,
  arma::colvec pl,
  int K, int N,int q,
  arma::colvec possible_label,
  arma::colvec Ci, double CSI,
  arma::umat index_row,
  List index_col) {
  R_CheckUserInterrupt();
  arma::colvec n_i_in_l(K), m_i_in_l(K), Zeta(K), corr(K),
  log_Nq_given_z_i(K), pp(K), pp2(K), logLik_mu_i(K);
  for(int i=0; i<N; i++){  // for every label
    arma::colvec fakeC=Ci;
    logLik_mu_i.fill(0);
    n_i_in_l.fill(0);
    m_i_in_l.fill(0);
    arma:: uvec IR=index_row.col(i),
      IC=index_col(i);
    for(int l=0; l<(K); l++){ // for every possible value taken by Ci
      logLik_mu_i(l) = log_Likelihood_double( mu_obser[i], dl[l]  );
      fakeC[i] = l+1;
      int Nzi_l = accu(fakeC==(l+1));
      //////////////////////////////////////////////////////////////////////////////////
      n_i_in_l[l] = accu( fakeC.elem(IR) == (l+1) );
      m_i_in_l[l] = accu( fakeC.elem(IC) == (l+1) );
      Zeta[l]       = Norm_Constant_Z_l2(Nzi_l,  N, CSI, q);
      corr[l]       = Norm_Constant_Z_l2(Nzi_l-1,  N, CSI, q)/Zeta[l];
      corr[l]       = (Nzi_l-1)*log(corr[l]);
      //////////////////////////////////////////////////////////////////////////////////
    }
    //Rcout << n_i_in_l << "..," << m_i_in_l;
    log_Nq_given_z_i =  (n_i_in_l+m_i_in_l) * (log(CSI)-log(1-CSI)) - log(Zeta) + corr;
    pp = logLik_mu_i + log_Nq_given_z_i + log(pl);
    // closing loop on labels to compute probabilities
    if(arma::is_finite(max(pp))){
      pp2 = exp( pp - max(pp));
      Ci(i) = RcppArmadillo::sample(
        possible_label, 1, TRUE, pp2)[0];
    }else{
      Ci(i) = RcppArmadillo::sample(possible_label, 1, TRUE)[0];
    }
  } // closing loop on observations
  return Ci;
}

// [[Rcpp::export]]
arma::colvec  Update_memberships_faster(arma::colvec mu_obser,
                                 arma::colvec dl,
                                 arma::colvec pl,
                                 int K, int N,int q,
                                 arma::colvec possible_label,
                                 arma::colvec Ci, double QQ,
                                 arma::umat index_row,
                                 List index_col,
                                 arma::colvec log_Precomp_Z,
                                 arma::colvec log_Precomp_ratios) {
  R_CheckUserInterrupt();
  arma::colvec n_i_in_l(K), m_i_in_l(K), log_Zeta(K), log_corr(K),
  log_Nq_given_z_i(K), pp(K), pp2(K), logLik_mu_i(K);
  for(int i=0; i<N; i++){  // for every label
    arma::colvec fakeC=Ci;
    logLik_mu_i.fill(0);
    n_i_in_l.fill(0);
    m_i_in_l.fill(0);
    arma:: uvec IR=index_row.col(i),
      IC=index_col(i);
    for(int l=0; l<(K); l++){ // for every possible value taken by Ci
      logLik_mu_i(l) = log_Likelihood_double( mu_obser[i], dl[l]  );
      fakeC[i] = l+1;
      int Nzi_l = accu(fakeC==(l+1));
      //////////////////////////////////////////////////////////////////////////////////
      n_i_in_l[l] = accu( fakeC.elem(IR) == (l+1) );
      m_i_in_l[l] = accu( fakeC.elem(IC) == (l+1) );
      log_Zeta[l]     = log_Precomp_Z[Nzi_l-1];     //Norm_Constant_Z_l2(Nzi_l,  N, CSI, q);
      log_corr[l]     = log_Precomp_ratios[Nzi_l-1];//Norm_Constant_Z_l2(Nzi_l-1,  N, CSI, q)/Zeta[l];
      //////////////////////////////////////////////////////////////////////////////////
    }
    log_Nq_given_z_i =  (n_i_in_l+m_i_in_l) * QQ - log_Zeta + log_corr;
    pp = logLik_mu_i + log_Nq_given_z_i + log(pl);
    // closing loop on labels to compute probabilities
    if(arma::is_finite(max(pp))){
      pp2 = exp( pp - max(pp));
      Ci(i) = RcppArmadillo::sample(
        possible_label, 1, TRUE, pp2)[0];
    }else{
      Ci(i) = RcppArmadillo::sample(possible_label, 1, TRUE)[0];
    }
  } // closing loop on observations
  return Ci;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
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

//' Efficient computation of the Posterior Coclustering Matrix
//'
//' @param inds a matrix containing the memebership labels for each observation (on the columns) in each MCMC iteration (on the rows)
//' @export
// [[Rcpp::export]]
arma::mat PSM(arma::mat inds){
  int nsim = inds.n_rows;
  int n= inds.n_cols;
  arma::mat PSM(n,n), D(n,n);
  D.eye(); PSM.zeros();

  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      arma::colvec Z = (inds.col(i)-inds.col(j));
      arma::uvec success = find(Z==0);
      PSM(i,j) = success.n_elem;
      PSM(j,i) = PSM(i,j);
    }
  }

  return(PSM/nsim + D);

}

//' Log-likelihood generated from the generic ratio of two NN's distances
//'
//' @param par the ID, variable of the objective function,
//' @param n1 order of the NN according to the first distance is computed.
//' @param n2 order of the NN according to the second distance is computed.
//' @param mu_dots the ratios of distance of generic order.
//'
//' @export
// [[Rcpp::export]]
double loglikelihood_gride(double par, int n1, int n2, arma::colvec mu_dots){


  if(n2<n1){
    stop("n2 is lower than n1");
  }
  int d_n12 = (n2-n1);
  int nn    = mu_dots.n_elem;
  double log_dens, lognum, logden, logch;

  arma::colvec vn1  = arma::linspace(1, n1, n1);
  arma::colvec vn2  = arma::linspace(1, n2, n2);
  arma::colvec vn12 = arma::linspace(1, d_n12, d_n12);


    if (d_n12 == 1) {
      lognum   = nn * log(par);
      logden   = ((n2 - 1) * par + 1) * arma::accu( log( mu_dots ) );
      log_dens = lognum - logden + sum( log( mu_dots > 1 ) );
    } else{
      logch    = arma::accu(log(vn2))-arma::accu(log(vn1))-arma::accu(log(vn12));
      // Rcout << logch << "--\n"  ;
      lognum   = nn * log(par) +
                 (d_n12 - 1) * arma::accu( log( pow( mu_dots, par)  - 1 ) );
      // Rcout << lognum << "--\n" ;
      logden   = ((n2 - 1) * par + 1) * sum(log(mu_dots));
      // Rcout << logden << "--\n"  ;
      log_dens = logch * nn + nn * log(d_n12)  +
                 lognum - logden + sum(log(mu_dots > 1.));
    }

    return(log_dens);
}
