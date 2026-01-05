functions {
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
  int<lower = 0, upper = 1> cloglog;
  int<lower = 0, upper = 1> rho_mu;
  int<lower = 0, upper = 4> likelihood;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
  int<lower = 1> K_t;
  matrix[N, K_t] X_t;
}
transformed data {
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
      int is_zero = 0;
      if (y[n] == 0) {
        is_zero += 1;
      }
      log_lik[n] = ptziloglik_lpdf(y[n] | likelihood, is_zero,
                                   mu[n], rho[n],
                                   likelihood == 0 ? sigma_obs[1] : phi[1]);
    }
  }
}
