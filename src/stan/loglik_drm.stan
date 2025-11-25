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
  int<lower = 0, upper = 1> movement;
  int<lower = 0, upper = 1> est_surv; // estimate mortality?
  int<lower = 0, upper = 1> est_init; // estimate "initial cohort"
  int<lower = 0, upper = 1> cloglog;
  int<lower = 0, upper = 1> minit;
  int<lower = 0, upper = 1> rho_mu;
  int<lower = 0, upper = 4> likelihood;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  array[est_surv ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  array[movement ? n_ages : 0] int ages_movement;
  vector[n_ages] selectivity_at_age;
  //--- initial cohort (if not estimated) ----
  array[est_init ? 0 : n_ages - 1] real init_data;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
  int<lower = 1> K_t;
  matrix[N, K_t] X_t;
}
transformed data {
  matrix[est_surv ? 0 : n_time, est_surv ? 0 : n_patches] fixed_m;
  if (!est_surv)
    fixed_m = rep_matrix(- m[1], n_time, n_patches);
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
  vector[N] log_lik;
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
    for (n in 1:N) {
      if (likelihood == 0) {
        real loc_par;
        loc_par = log(mu[n]) + square(sigma_obs[1]) / 2;
        if (y[n] == 0) {
          // only evaluate density if there are length comps to evaluate
          log_lik[n] = log(rho[n]);
        } else {
          log_lik[n] = log1m(rho[n]) +
            lognormal_lpdf(y[n] | loc_par, sigma_obs[1]);
        }
      } else if (likelihood == 1) {
        real mu_ln;
        real sigma_ln;
        sigma_ln = sqrt(log1p(phi[1] * inv_square(mu[n])));
        mu_ln = log(square(mu[n]) * inv_sqrt(square(mu[n]) + phi[1]));
        if (y[n] == 0) {
          log_lik[n] = log(rho[n]);
        } else {
          log_lik[n] = log1m(rho[n]) +
            ln_mu_lpdf(y[n] | mu[n], phi[1]);
        }
      } else if (likelihood == 2) {
        real gamma_beta;
        gamma_beta = phi[1] / mu[n];
        if (y[n] == 0) {
          log_lik[n] = log(rho[n]);
        } else {
          log_lik[n] = log1m(rho[n]) +
            gamma_lpdf(y[n] | phi[1], gamma_beta);
        }
      } else if (likelihood == 3) {
        real a_ll;
        real b_ll;
        b_ll = phi[1] + 1;
        a_ll = sin(pi() / b_ll) * mu[n] * b_ll * inv(pi());
        if (y[n] == 0) {
          log_lik[n] = log(rho[n]);
        } else {
          log_lik[n] = log1m(rho[n]) +
            loglogistic_lpdf(y[n] | a_ll, b_ll);
        }
      } else {
        array[2] real aux_tn = rep_array(0.0, 2);
        aux_tn[2] = normal_rng(mu[n], phi[1]);
        if (y[n] == 0) {
          log_lik[n] = log(rho[n]) +
            normal_lpdf(0.0 | mu[n], phi[1]) +
            normal_lccdf(0.0 | mu[n], phi[1]);
        } else {
          log_lik[n] = log1m(rho[n]) +
            normal_lpdf(y[n] | mu[n], phi[1]) +
            normal_lccdf(0.0 | mu[n], phi[1]);
        }
      }
    }
  }
}
