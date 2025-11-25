functions {
#include utils/age_struct.stan
#include utils/lpdfs.stan
}
data {
  //--- survey data  ---
  int N; // n_patches * n_time
  int n_ages; // number of ages
  int n_patches; // number of patches
  int n_time; // years for training
  //--- toggles ---
  vector[N] y;
  int<lower = 0, upper = 1> cloglog;
  int<lower = 0, upper = 1> rho_mu;
  int<lower = 0, upper = 4> likelihood;
  //--- for vectorizing zero-inflation ----
  int N_nz;
  int N_z;
  array[N_nz] int id_nz;
  array[N_z] int id_z;
  //--- environmental data ----
  //--- * for mortality ----
  int<lower = 1> K_t;
  matrix[N, K_t] X_t;
}
transformed data {
  int N_st = N_z + N_nz;
}
parameters {
  vector[N] mu;
  vector[K_t] beta_t;
  array[rho_mu] real lxi;
  array[likelihood == 0 ? 1 : 0] real<lower = 0> sigma_obs;
  array[likelihood > 0 ? 1 : 0] real<lower = 0> phi;
}
transformed parameters {
}
generated quantities {
  vector[N_st] log_lik;
  {
    array[rho_mu] real xi;
    if (rho_mu)
      xi = - exp(lxi);
    vector[N] rho;
    // rho now hasa "regression like" type
    if (cloglog) {
      if (rho_mu) {
        rho =
          inv_cloglog(X_t * beta_t + xi[1] .* log(mu));
      } else {
        rho =
          inv_cloglog(X_t * beta_t);
      }
    } else {
      if (rho_mu) {
        rho =
          inv_logit(X_t * beta_t + xi[1] .* log(mu));
      } else {
        rho =
          inv_logit(X_t * beta_t);
      }
    }
    // zeros
    for (n in 1:N_z) {
      log_lik[n] = log(rho[id_z[n]]);
    }
    // non-zeros
    for (n in 1:N_nz) {
      if (likelihood == 0) {
        real loc_par;
        loc_par = log(mu[id_nz[n]]) + square(sigma_obs[1]) / 2;
        log_lik[n + N_z] = log1m(rho[id_nz[n]]) +
          lognormal_lpdf(y[id_nz[n]] | loc_par, sigma_obs[1]);
      } else if (likelihood == 1) {
        real mu_ln;
        real sigma_ln;
        sigma_ln = sqrt(log1p(phi[1] * inv_square(mu[id_nz[n]])));
        mu_ln = log(square(mu[id_nz[n]]) * inv_sqrt(square(mu[id_nz[n]]) + phi[1]));
        log_lik[n + N_z] = log1m(rho[id_nz[n]]) +
          ln_mu_lpdf(y[id_nz[n]] | mu[id_nz[n]], phi[1]);
      } else if (likelihood == 2) {
        real gamma_beta;
        gamma_beta = phi[1] / mu[id_nz[n]];
        log_lik[n + N_z] = log1m(rho[id_nz[n]]) +
          gamma_lpdf(y[id_nz[n]] | phi[1], gamma_beta);
      } else if (likelihood == 3) {
        real a_ll;
        real b_ll;
        b_ll = phi[1] + 1;
        a_ll = sin(pi() / b_ll) * mu[id_nz[n]] * b_ll * inv(pi());
        log_lik[n + N_z] = log1m(rho[id_nz[n]]) +
          loglogistic_lpdf(y[id_nz[n]] | a_ll, b_ll);
      } else {
        array[2] real aux_tn = rep_array(0.0, 2);
        aux_tn[2] = normal_rng(mu[id_nz[n]], phi[1]);
        log_lik[n + N_z] = log1m(rho[id_nz[n]]) +
          normal_lpdf(y[id_nz[n]] | mu[id_nz[n]], phi[1]) +
          normal_lccdf(0.0 | mu[id_nz[n]], phi[1]);
      }
    }
  }
}
