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
    for (n in 1:N_z) {
      log_lik[id_z[n]] =
        ptziloglik_lpdf(y[id_z[n]] | likelihood, 1,
                        mu[id_z[n]], rho[id_z[n]],
                        likelihood == 0 ? sigma_obs[1] : phi[1]);
    }
    for (n in 1:N_nz) {
      log_lik[id_nz[n]] =
        ptziloglik_lpdf(y[id_nz[n]] | likelihood, 0,
                        mu[id_nz[n]], rho[id_nz[n]],
                        likelihood == 0 ? sigma_obs[1] : phi[1]);
    }
  }
}
