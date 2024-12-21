functions {
}
data {
  //--- survey data  ---
  int N;
  int n_ages; // number of ages
  int n_patches; // number of patches
  int n_time; // years for forecasts
  int n_time_train; // years for training
  array[N] int time;
  //--- toggles ---
  int<lower = 0, upper = 1> p_error;
  int<lower = 0, upper = 1> cloglog; // use cloglog instead of logit for theta
  int<lower = 0, upper = 3> likelihood; // (0 = Original LN, 1 = repar LN, 2 =
                                        // Gamma, 3 = log-Logistic)
  //--- suitability (for theta) ----
  int<lower = 1> K_z;
  matrix[N, K_z] Z;
  //--- environmental data ----
  int<lower = 1> K_x;
  matrix[N, K_x] X;
}
parameters {
  //--- "regression" coefficients ----
  vector[K_z] coef_t;
  vector[K_x] coef_r;
  //--- parameters from AR process ----
  vector[p_error ? n_time_train : 0] rec_dev;
  array[p_error] real rho;
  array[p_error] real sigma_r;
  //--- additional parameter for different lik functions ----
  array[likelihood > 0 ? 1 : 0] real phi;
  array[likelihood == 0 ? 1 : 0] real<lower = 0> sigma_obs;
}
generated quantities {
  //--- projected expected density ----
  vector[n_time * n_patches] mu_proj;
  //--- projected total density ----
  vector[n_time * n_patches] y_proj;
  //--- projected absence probability ----
  vector[n_time * n_patches] theta_proj;
  //--- AR term ----
  vector[p_error ? n_time : 0] rec_proj;
  if (p_error) {
    {
      vector[n_time] raw;
      vector[p_error ? n_time - 1 : 0] lagged_rec;
      raw[1] = std_normal_rng();
      rec_proj[1] = rho[1] * rec_dev[n_time_train] +
        sigma_r[1] * raw[1];
      for (tp in 2:n_time) {
        raw[tp] = std_normal_rng();
        lagged_rec[tp - 1] = rec_proj[tp - 1];
        rec_proj[tp] = rho[1] * lagged_rec[tp - 1] +
          sigma_r[1] * raw[tp];
      }
    }
  }
  //--- mu_proj calculations ----
  {
    vector[N] lmu;
    lmu = X * coef_r;
    if (p_error) {
      for (n in 1:N)
        lmu[n] += rec_dev[time[n]];
    }
    mu_proj = exp(lmu);
  }
  //--- absence probabilities ----
  if (cloglog) {
    theta_proj =
      inv_cloglog(Z * coef_t);
  } else {
    theta_proj =
      inv_logit(Z * coef_t);
  }
  //--- y_proj calculations ----
  for (n in 1:N) {
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu_proj[n]) + square(sigma_obs[1]) / 2;
      y_proj[n] = (1 - bernoulli_rng(theta_proj[n])) *
        lognormal_rng(loc_par, sigma_obs[1]);
    } else if (likelihood == 1) {
      real mu_ln;
      real sigma_ln;
      sigma_ln = sqrt(log1p(phi[1] * inv_square(mu_proj[n])));
      mu_ln = log(square(mu_proj[n]) * inv_sqrt(square(mu_proj[n]) + phi[1]));
      y_proj[n] = (1 - bernoulli_rng(theta_proj[n])) *
        lognormal_rng(mu_ln, sigma_ln);
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi[1] / mu_proj[n];
      y_proj[n] = (1 - bernoulli_rng(theta_proj[n])) *
        gamma_rng(phi[1], gamma_beta);
    } else if (likelihood == 3) {
      real a_ll;
      real b_ll;
      b_ll = phi[1] + 1;
      a_ll = sin(pi() / b_ll) * mu_proj[n] * inv(pi() * b_ll);
      y_proj[n] = (1 - bernoulli_rng(theta_proj[n])) *
        loglogistic_rng(a_ll, b_ll);
    }
  }
}
