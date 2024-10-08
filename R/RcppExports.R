# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

log_Likelihood_double <- function(mu_obs, d) {
    .Call(`_intRinsic_log_Likelihood_double`, mu_obs, d)
}

Groups_quantities <- function(mu_obser, Ci, K) {
    .Call(`_intRinsic_Groups_quantities`, mu_obser, Ci, K)
}

Norm_Constant_Z_l2 <- function(Nzi_l, N, xi, q) {
    .Call(`_intRinsic_Norm_Constant_Z_l2`, Nzi_l, N, xi, q)
}

log_Zeta_maker <- function(N, xi, q) {
    .Call(`_intRinsic_log_Zeta_maker`, N, xi, q)
}

index_row_col <- function(Nq, q, N) {
    .Call(`_intRinsic_index_row_col`, Nq, q, N)
}

rdir_cpp <- function(alpha) {
    .Call(`_intRinsic_rdir_cpp`, alpha)
}

Update_memberships_faster <- function(mu_obser, dl, pl, K, N, q, possible_label, Ci, QQ, index_row, index_col, log_Precomp_Z, log_Precomp_ratios) {
    .Call(`_intRinsic_Update_memberships_faster`, mu_obser, dl, pl, K, N, q, possible_label, Ci, QQ, index_row, index_col, log_Precomp_Z, log_Precomp_ratios)
}

gam_trunc <- function(D, K, a0_d, b0_d, n_l, sLog) {
    .Call(`_intRinsic_gam_trunc`, D, K, a0_d, b0_d, n_l, sLog)
}

gam_trunc_pmass <- function(D, K, a0_d, b0_d, n_l, sLog, pi_mass) {
    .Call(`_intRinsic_gam_trunc_pmass`, D, K, a0_d, b0_d, n_l, sLog, pi_mass)
}

gride_log_likelihood <- function(d, n1, n2, mus_n1_n2) {
    .Call(`_intRinsic_gride_log_likelihood`, d, n1, n2, mus_n1_n2)
}

gride_log_posterior <- function(z, n1, n2, a_d, b_d, mus_n1_n2) {
    .Call(`_intRinsic_gride_log_posterior`, z, n1, n2, a_d, b_d, mus_n1_n2)
}

gride_mh_sampler <- function(start_d, n1, n2, a_d, b_d, mus_n1_n2, nsim, burn_in, sigma) {
    .Call(`_intRinsic_gride_mh_sampler`, start_d, n1, n2, a_d, b_d, mus_n1_n2, nsim, burn_in, sigma)
}

