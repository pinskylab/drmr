functions {
  // counting number of zeros
  int num_non_zero_fun(vector y) {
    int A = 0;
    int N = size(y);
    
    for (n in 1 : N) {
      if (y[n] != 0) {
        A += 1;
      }
    }
    return A;
  }
  
  array[] int non_zero_index_fun(vector y, int A) {
    int N = size(y);
    array[A] int non_zero_index;
    int counter = 0;
    for (n in 1 : N) {
      if (y[n] != 0) {
        counter += 1;
        non_zero_index[counter] = n;
      }
    }
    return non_zero_index;
  }
  
  array[] int zero_index_fun(vector y, int Z) {
    int N = size(y);
    array[Z] int zero_index;
    int counter = 0;
    for (n in 1 : N) {
      if (y[n] == 0) {
        counter += 1;
        zero_index[counter] = n;
      }
    }
    return zero_index;
  }
  /**
   * Penalized complexity prior for the AR(1) correlation. This prior assumes:
   * P(|x| > u) = alpha A reasonable strategy here is to set u as a possibly high
   * autocorrelation and alpha as a small probability.
   *
   * Reference: https://doi.org/10.1111/jtsa.12242
   * 
   * @param x autocorrelation parameter
   * @param u upper bound for |x|
   * @param alpha prior probability of exceeding the upper bound u.
   * 
   * @return an array of numbers by age, year and patch
   */
  real pcp_ar0_lpdf(real x, real alpha, real u) {
    real aux = - log(1 - square(x));
    real rho = - log(alpha) / sqrt(- log(1 - square(u)));
    return log(rho) - log(2) - rho * sqrt(aux) +
      log(abs(x)) + aux - 0.5 * log(aux);
  }
  /**
   * Penalized complexity prior for the log-transformed standard deviation of a
   * random effect. It assumes: P(exp(x) > u) = alpha.
   * 
   * @param x log-sd parameter
   * @param u upper bound for conditional SD
   * @param alpha prior probability of exceeding the upper bound u.
   * 
   * @return an array of numbers by age, year and patch
   */
  real pcp_logsd_lpdf(real x, real alpha, real u) {
    real lambda = - log(alpha) / u;
    return log(lambda) + x - lambda * exp(x);
  }

  /**
   * Gamma lpdf reparametrized.
   *
   * Reference: https://doi.org/10.1111/jtsa.12242
   * 
   * @param x non-negative random variable
   * @param mu theoretical mean of X
   * @param phi inverse scale parameter.
   * 
   * @return a log-pdf
   */
  real gamma_mu_lpdf(real x, real mu, real phi) {
    return - lgamma(phi) + lmultiply(phi, phi) - lmultiply(phi, mu) +
      lmultiply(phi, x) - phi * x * inv(mu);
  }

  /**
   * Log-Normal lpdf reparametrized.
   *
   * Reference: https://doi.org/10.1111/jtsa.12242
   * 
   * @param x non-negative random variable
   * @param mu theoretical mean of X
   * @param phi weird parameter parameter.
   * 
   * @return a log-pdf
   */
  real ln_mu_lpdf(real x, real mu, real phi) {
    real f_mu_phi;
    f_mu_phi = log1p(phi * inv_square(mu));
    real output = 0;
    output += log2() + log(pi()) + log(f_mu_phi) + log(x) +
      square(log(x) + log2() - 0.5 * f_mu_phi) * inv(f_mu_phi);
    return - 0.5 * output;
  }

  /**
   * Inverse Gaussian lpdf reparametrized.
   *
   * Reference: https://doi.org/10.1111/jtsa.12242
   * 
   * @param x non-negative random variable
   * @param mu theoretical mean of X
   * @param phi weird parameter parameter.
   * 
   * @return a log-pdf
   */
  real igaus_mu_lpdf(real x, real mu, real phi) {
    real f_mu_phi;
    real output = 0;
    output += log2() + log(pi()) + 3 * log(x) - log(phi) +
      log(phi) * square(x - mu) * inv(mu * mu * x);
    return - 0.5 * output;
  }

  real igaus_mu_rng(real mu, real phi) {
    real nu;
    real output;
    real z;
    nu = std_normal_rng();
    z = uniform_rng(0, 1);
    output = square(nu);
    output = mu + square(mu) * output * inv(2 * phi) -
      mu * inv(2 * phi) *
      sqrt(4 * mu * phi * output + square(mu * output));
    if (z <= mu * inv(mu + output)) {
      return output;
    } else {
      return square(mu) * inv(output);
    }
  }
}
data {
  //--- survey data  ---
  int N; // n_patches * n_time
  int n_patches; // number of patches
  int n_time; // years for training
  array[N] int<lower = 1, upper = n_time> time;
  vector[N] y;
  //--- toggles ---
  int<lower = 0, upper = 1> ar_re;
  int<lower = 0, upper = 1> cloglog; // use cloglog instead of logit for rho
  int<lower = 0, upper = 3> likelihood; // (0 = Original LN, 1 = repar LN, 2 =
                                        // Gamma, 3 = log-Logistic)
  //--- suitability (for rho) ----
  int<lower = 1> K_z;
  matrix[N, K_z] Z;
  //--- environmental data ----
  //--- * for counts ----
  int<lower = 1> K_x;
  matrix[N, K_x] X;
  //--- priors hyperparameters ----
  real pr_sigma_obs_mu; 
  real pr_sigma_obs_sd; // formerly sigma_obs_cv
  real pr_phi_a; // revise this. It is for phi in gamma, loglogistic, and
                 // inverse-gaussian
  real pr_phi_b;
  // * now AR SD parameters have pcpriors
  real pr_lsigma_t_mu; 
  real pr_lsigma_t_sd;
  real pr_alpha_a; 
  real pr_alpha_b;
  vector[K_z] pr_beta_t_mu;
  vector[K_z] pr_beta_t_sd;
  vector[K_x] pr_beta_r_mu;
  vector[K_x] pr_beta_r_sd;
}
transformed data {
  //--- Vectorizing zero-inflation ----
  int N_nz;
  N_nz = num_non_zero_fun(y);
  int N_z = N - N_nz;
  array[N_nz] int id_nz;
  array[N_z] int id_z;
  id_nz = non_zero_index_fun(y, N_nz);
  id_z = zero_index_fun(y, N_z);
}
parameters {
  // parameter associated to the likelihood 
  array[likelihood == 0 ? 1 : 0] real<lower = 0> sigma_obs;
  array[likelihood > 0 ? 1 : 0] real<lower = 0> phi;
  // coefficients for recruitment (it is a log-linear model)
  vector[K_x] beta_r;
  // parameter associated with "encounter probability"
  vector[K_z] beta_t;
  //--- * AR process parameters ----
  // conditional SD
  array[ar_re] real log_sigma_t;
  // autocorrelation
  array[ar_re] real<lower = 0, upper = 1> alpha;
  // aux latent variable
  vector[ar_re ? n_time : 0] w_t;
}
transformed parameters {
  //--- AR process ----
  array[ar_re] real sigma_t;
  vector[ar_re ? n_time : 0] z_t;
  {
    if (ar_re) {
      sigma_t[1] = exp(log_sigma_t[1]);
      z_t[1] = sigma_t[1] * w_t[1];
      for (tp in 2:n_time) {
        z_t[tp] += alpha[1] * z_t[tp - 1] +
          sigma_t[1] * w_t[tp];
      }
    }
  }
  //--- quantities used in the likelihood ----
  // probability of encounter at specific time/patch combinations
  vector[N] rho;
  // Expected density at specific time/patch combinations
  vector[N] mu =
    rep_vector(0.0, N);
  // rho now hasa "regression like" type
  if (cloglog) {
    rho =
      inv_cloglog(Z * beta_t);
  } else {
    rho =
      inv_logit(Z * beta_t);
  }
  {
    vector[N] lmu;
    lmu = X * beta_r;
    if (ar_re) {
      for (n in 1:N)
        lmu[n] += z_t[time[n]];
    }
    mu = exp(lmu);
  }
}
// close transformed parameters block
model {
  //--- AR process ----
  if (ar_re) {
    target += std_normal_lpdf(w_t);
    target += normal_lpdf(log_sigma_t[1] | pr_lsigma_t_mu, pr_lsigma_t_sd);
    target += beta_lpdf(alpha[1] | pr_alpha_a, pr_alpha_b); 
  }
  //--- Counts ----
  target += normal_lpdf(beta_r | pr_beta_r_mu, pr_beta_r_sd);
  //--- Suitability ----
  target += normal_lpdf(beta_t | pr_beta_t_mu, pr_beta_t_sd);
  //--- Likelihood ----
  if (likelihood == 0) {
    target += normal_lpdf(sigma_obs[1] | pr_sigma_obs_mu, pr_sigma_obs_sd) -
      1.0 * normal_lccdf(0 | pr_sigma_obs_mu, pr_sigma_obs_sd);
  } else {
    // change these parameters (PC prior for exponential?)
    target += gamma_lpdf(phi[1] | pr_phi_a, pr_phi_b);
  }
  // only evaluate density if there are length comps to evaluate
  target += sum(log(rho[id_z]));
  if (likelihood == 0) {
    vector[N_nz] loc_par;
    loc_par = log(mu[id_nz]) + square(sigma_obs[1]) / 2;
    target += log1m(rho[id_nz]);
    target += lognormal_lpdf(y[id_nz] | loc_par, sigma_obs[1]);
  } else if (likelihood == 1) {
    vector[N_nz] mu_ln;
    vector[N_nz] sigma_ln;
    sigma_ln = sqrt(log1p(phi[1] * inv_square(mu[id_nz])));
    mu_ln = log(square(mu[id_nz]) .*
                inv_sqrt(square(mu[id_nz]) + phi[1]));
    target += log1m(rho[id_nz]);
    target += lognormal_lpdf(y[id_nz] | mu_ln, sigma_ln);
  } else if (likelihood == 2) {
    vector[N_nz] b_g;
    b_g = phi[1] / mu[id_nz];
    target += log1m(rho[id_nz]);
    target += gamma_lpdf(y[id_nz] | phi[1], b_g);
  } else {
    vector[N_nz] a_ll;
    real b_ll;
    b_ll = phi[1] + 1;
    a_ll = sin(pi() / b_ll) * mu[id_nz] * b_ll * inv(pi());
    target += log1m(rho[id_nz]);
    target += loglogistic_lpdf(y[id_nz] | a_ll, b_ll);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_pp;
  for (n in 1:N) {
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu[n]) + square(sigma_obs[1]) / 2;
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        lognormal_rng(loc_par, sigma_obs[1]);
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
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        lognormal_rng(mu_ln, sigma_ln);
      if (y[n] == 0) {
        log_lik[n] = log(rho[n]);
      } else {
        log_lik[n] = log1m(rho[n]) +
          ln_mu_lpdf(y[n] | mu[n], phi[1]);
      }
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi[1] / mu[n];
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        gamma_rng(phi[1], gamma_beta);
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
      a_ll = sin(pi() / b_ll) * mu[n]  * b_ll * inv(pi());
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        loglogistic_rng(a_ll, b_ll);
      if (y[n] == 0) {
        log_lik[n] = log(rho[n]);
      } else {
        log_lik[n] = log1m(rho[n]) +
          loglogistic_lpdf(y[n] | a_ll, b_ll);
      }
    }
  }
}
