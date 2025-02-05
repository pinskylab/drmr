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
    real theta = - log(alpha) / sqrt(- log(1 - square(u)));
    return log(theta) - log(2) - theta * sqrt(aux) +
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
  int<lower = 0, upper = 1> time_ar;
  int<lower = 0, upper = 1> cloglog; // use cloglog instead of logit for theta
  int<lower = 0, upper = 1> qr_z; // use qr parametrization for theta?
  int<lower = 0, upper = 1> qr_x; // use qr parametrization for logrec?
  int<lower = 0, upper = 3> likelihood; // (0 = Original LN, 1 = repar LN, 2 =
                                        // Gamma, 3 = log-Logistic)
  //--- suitability (for theta) ----
  int<lower = 1> K_z;
  matrix[N, K_z] Z;
  //--- environmental data ----
  //--- * for counts ----
  int<lower = 1> K_x;
  matrix[N, K_x] X;
  //--- priors hyperparameters ----
  real pr_sigma_obs_mu; 
  real pr_sigma_obs_sd; // formerly sigma_obs_cv
  real pr_phi_mu; // revise this. It is for phi in gamma, loglogistic, and
                 // inverse-gaussian
  real pr_phi_sd;
  // * now AR SD parameters have pcpriors
  real pr_logsd_r_mu; 
  real pr_logsd_r_sd;
  real pr_alpha_a; 
  real pr_alpha_b;
  vector[K_z] pr_coef_t_mu;
  vector[K_z] pr_coef_t_sd;
  vector[K_x] pr_coef_r_mu;
  vector[K_x] pr_coef_r_sd;
}
transformed data {
  //--- thin-QR parametrization (improves sampler at almost no cost) ----
  // * detection prob
  matrix[qr_z ? N : 0,
         qr_z ? K_z : 0] Q_z;
  matrix[qr_z ? K_z : 0,
         qr_z ? K_z : 0] R_z;
  matrix[qr_z ? K_z : 0,
         qr_z ? K_z : 0] R_z_inv;
  if (qr_z) {
    Q_z = qr_thin_Q(Z) * sqrt(N - 1);
    R_z = qr_thin_R(Z) / sqrt(N - 1);
    R_z_inv = inverse(R_z);
  }
  // * Count
  matrix[qr_x ? N : 0, qr_x ? K_x : 0] Q_x;
  matrix[qr_x ? K_x : 0, qr_x ? K_x : 0] R_x;
  matrix[qr_x ? K_x : 0, qr_x ? K_x : 0] R_x_inv;
  if (qr_x) {
    Q_x = qr_thin_Q(X) * sqrt(N - 1);
    R_x = qr_thin_R(X) / sqrt(N - 1);
    R_x_inv = inverse(R_x);
  }
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
  array[likelihood > 0 ? 1 : 0] real log_phi;
  // coefficients for recruitment (it is a log-linear model)
  vector[K_x] coef_r0;
  // parameter associated with "encounter probability"
  vector[K_z] coef_t0;
  //--- * AR process parameters ----
  // conditional SD
  array[time_ar] real log_tau;
  // autocorrelation
  array[time_ar] real<lower = 0, upper = 1> shift_alpha;
  // aux latent variable
  vector[time_ar ? n_time : 0] raw;
}
transformed parameters {
  array[likelihood > 0 ? 1 : 0] real phi;
  if (likelihood > 0)
    phi = exp(log_phi);
  //--- AR process ----
  array[time_ar] real tau;
  array[time_ar] real alpha;
  vector[time_ar ? n_time : 0] z_t;
  {
    vector[time_ar ? n_time : 0] lagged_z_t;
    if (time_ar) {
      tau[1] = exp(log_tau[1]);
      alpha[1] = 2.0 * shift_alpha[1] - 1.0;
      z_t = tau[1] * raw;
      for (tp in 2:n_time) {
        lagged_z_t[tp] = z_t[tp - 1];
        z_t[tp] += alpha[1] * lagged_z_t[tp];
      }
    }
  }
  //--- quantities used in the likelihood ----
  // probability of encounter at specific time/patch combinations
  vector[N] theta;
  // Expected density at specific time/patch combinations
  vector[N] mu =
    rep_vector(0.0, N);
  // theta now hasa "regression like" type
  if (cloglog) {
    matrix[N, K_z] Z_aux;
    Z_aux = qr_z ? Q_z : Z;
    theta =
      inv_cloglog(Z_aux * coef_t0);
  } else {
    matrix[N, K_z] Z_aux;
    Z_aux = qr_z ? Q_z : Z;
    theta =
      inv_logit(Z_aux * coef_t0);
  }
  {
    matrix[N, K_x] X_aux;
    vector[N] lmu;
    X_aux = qr_x ? Q_x : X;
    lmu = X_aux * coef_r0;
    if (time_ar) {
      for (n in 1:N)
        lmu[n] += z_t[time[n]];
    }
    mu = exp(lmu);
  }
}
// close transformed parameters block
model {
  //--- AR process ----
  if (time_ar) {
    target += std_normal_lpdf(raw);
    target += normal_lpdf(log_tau[1] | pr_logsd_r_mu, pr_logsd_r_sd);
    target += beta_lpdf(shift_alpha[1] | pr_alpha_a, pr_alpha_b); 
  }
  //--- Counts ----
  target += normal_lpdf(coef_r0 | pr_coef_r_mu, pr_coef_r_sd);
  //--- Suitability ----
  target += normal_lpdf(coef_t0 | pr_coef_t_mu, pr_coef_t_sd);
  //--- Likelihood ----
  if (likelihood == 0) {
    target += normal_lpdf(sigma_obs[1] | pr_sigma_obs_mu, pr_sigma_obs_sd) -
      1.0 * normal_lccdf(0 | pr_sigma_obs_mu, pr_sigma_obs_sd);
  } else {
    // change these parameters (PC prior for exponential?)
    target += student_t_lpdf(log_phi[1] | 3, pr_phi_mu, pr_phi_sd);
  }
  // only evaluate density if there are length comps to evaluate
  target += sum(log(theta[id_z]));
  if (likelihood == 0) {
    vector[N_nz] loc_par;
    loc_par = log(mu[id_nz]) + square(sigma_obs[1]) / 2;
    target += log1m(theta[id_nz]);
    target += lognormal_lpdf(y[id_nz] | loc_par, sigma_obs[1]);
  } else if (likelihood == 1) {
    vector[N_nz] mu_ln;
    vector[N_nz] sigma_ln;
    sigma_ln = sqrt(log1p(phi[1] * inv_square(mu[id_nz])));
    mu_ln = log(square(mu[id_nz]) .*
                inv_sqrt(square(mu[id_nz]) + phi[1]));
    target += log1m(theta[id_nz]);
    target += lognormal_lpdf(y[id_nz] | mu_ln, sigma_ln);
  } else if (likelihood == 2) {
    vector[N_nz] b_g;
    b_g = phi[1] / mu[id_nz];
    target += log1m(theta[id_nz]);
    target += gamma_lpdf(y[id_nz] | phi[1], b_g);
  } else {
    vector[N_nz] a_ll;
    real b_ll;
    b_ll = phi[1] + 1;
    a_ll = sin(pi() / b_ll) * mu[id_nz] * inv(pi() * b_ll);
    target += log1m(theta[id_nz]);
    target += loglogistic_lpdf(y[id_nz] | a_ll, b_ll);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_pp;
  vector[K_z] coef_t;
  vector[K_x] coef_r;
  for (n in 1:N) {
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu[n]) + square(sigma_obs[1]) / 2;
      y_pp[n] = (1 - bernoulli_rng(theta[n])) *
        lognormal_rng(loc_par, sigma_obs[1]);
      if (y[n] == 0) {
        // only evaluate density if there are length comps to evaluate
        log_lik[n] = log(theta[n]);
      } else {
        log_lik[n] = log1m(theta[n]) +
          lognormal_lpdf(y[n] | loc_par, sigma_obs[1]);
      }
    } else if (likelihood == 1) {
      real mu_ln;
      real sigma_ln;
      sigma_ln = sqrt(log1p(phi[1] * inv_square(mu[n])));
      mu_ln = log(square(mu[n]) * inv_sqrt(square(mu[n]) + phi[1]));
      y_pp[n] = (1 - bernoulli_rng(theta[n])) *
        lognormal_rng(mu_ln, sigma_ln);
      if (y[n] == 0) {
        log_lik[n] = log(theta[n]);
      } else {
        log_lik[n] = log1m(theta[n]) +
          ln_mu_lpdf(y[n] | mu[n], phi[1]);
      }
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi[1] / mu[n];
      y_pp[n] = (1 - bernoulli_rng(theta[n])) *
        gamma_rng(phi[1], gamma_beta);
      if (y[n] == 0) {
        log_lik[n] = log(theta[n]);
      } else {
        log_lik[n] = log1m(theta[n]) +
          gamma_lpdf(y[n] | phi[1], gamma_beta);
      }
    } else if (likelihood == 3) {
      real a_ll;
      real b_ll;
      b_ll = phi[1] + 1;
      a_ll = sin(pi() / b_ll) * mu[n] * inv(pi() * b_ll);
      y_pp[n] = (1 - bernoulli_rng(theta[n])) *
        loglogistic_rng(a_ll, b_ll);
      if (y[n] == 0) {
        log_lik[n] = log(theta[n]);
      } else {
        log_lik[n] = log1m(theta[n]) +
          loglogistic_lpdf(y[n] | a_ll, b_ll);
      }
    }
  }
  if (qr_z) {
    coef_t = R_z_inv * coef_t0;
  } else coef_t = coef_t0;
  if (qr_x) {
    coef_r = R_x_inv * coef_r0;
  } else coef_r = coef_r0;
}
