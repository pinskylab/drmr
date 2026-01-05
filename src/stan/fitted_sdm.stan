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
  int<lower = 0, upper = 1> cloglog;
  int<lower = 0, upper = 1> rho_mu;
  int<lower = 0, upper = 4> likelihood;
  //--- environmental data ----
  int<lower = 1> K_z;
  matrix[N, K_z] Z;
}
transformed data {
}
parameters {
  vector[N] mu;
  vector[K_z] beta_t;
  array[rho_mu] real lxi;
  array[likelihood == 0 ? 1 : 0] real<lower = 0> sigma_obs;
  array[likelihood > 0 ? 1 : 0] real<lower = 0> phi;
}
transformed parameters {
}
generated quantities {
  vector[N] y_pp;
  {
    array[rho_mu] real xi;
    if (rho_mu)
      xi = - exp(lxi);
    vector[N] rho;
    // rho now hasa "regression like" type
    if (cloglog) {
      if (rho_mu) {
        rho =
          inv_cloglog(Z * beta_t + xi[1] .* log(mu));
      } else {
        rho =
          inv_cloglog(Z * beta_t);
      }
    } else {
      if (rho_mu) {
        rho =
          inv_logit(Z * beta_t + xi[1] .* log(mu));
      } else {
        rho =
          inv_logit(Z * beta_t);
      }
    }
    for (n in 1:N) {
      y_pp[n] = drmsdm_rng(mu[n], rho[n],
                           likelihood == 0 ? sigma_obs[1] : phi[1],
                           likelihood);
    }
  }
}
