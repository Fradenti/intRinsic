#ifndef G_GRIDE
#define G_GRIDE


double gride_log_likelihood(double d,
                            int n1,
                            int n2,
                            arma::colvec mus_n1_n2);

double gride_log_posterior(double z,
                           int n1,
                           int n2,
                           double a_d,
                           double b_d,
                           arma::colvec mus_n1_n2);

arma::colvec gride_mh_sampler(double start_d,
                              int n1,
                              int n2,
                              double a_d,
                              double b_d,
                              arma::colvec mus_n1_n2,
                              int nsim,
                              int burn_in,
                              double sigma);

#endif
