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
   * @title Generate theoretical mean according to the simplest model possible
   *
   * @description
   * 
   * @param n_patches number of patches
   * @param n_time number of years of training data
   * @param n_ages number of age classes
   * @param f_a_t f_{at}?
   * @param neg_mort minus natural mortality (instantaneous) rate
   * @param init a n_ages - 1 array.
   * 
   * @return an array of numbers by age, year and patch
   */
  array[] matrix simplest(int n_patches,
                          int n_time,
                          int n_ages,
                          // Mortality parameter
                          matrix f_a_t,
                          matrix neg_mort,
                          // initialization
                          array[] real init,
                          int init_type,
                          matrix recruitment) {
    // initializing output with zeros
    array[n_ages] matrix[n_time, n_patches] output
      = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
    // initializing recruitment elements
    for (p in 1:n_patches) {
      for (i in 1:n_time) {
        output[1, i, p] = recruitment[i, p];
      }
    }
    // Think about the indexes below for init_type 1 & 2
    if (init_type == 1) {
      // initializes every element at the recruitment
      output = rep_array(recruitment, n_ages);
    } else if (init_type == 2) {
      for (a in 2 : n_ages) {
        // makes more sense
        output[a, 1, ] = rep_row_vector(init[a - 1], n_patches);
      }
    } else if (init_type == 4) {
      for (a in 2 : n_ages) {
        // makes more sense
        output[a, 1, ] = rep_row_vector(init[a - 1], n_patches);
      }
    } else if (init_type == 3) {
      for (a in 2 : n_ages) {
        output[a, 1, ] = 
          rep_row_vector(init[a - 1], n_patches) .*
          recruitment[1, ];
      }
    } else if (init_type == 5) {
      for (a in 2 : n_ages) {
        output[a, 1, ] = 
          rep_row_vector(init[a - 1], n_patches) .*
          recruitment[1, ];
      }
    }
    // given the initialization, does the order of the for loops matter?
    for (p in 1 : n_patches) {
      for (i in 2 : n_time) {
        for (a in 2 : n_ages) {
          output[a, i, p] = output[a - 1, i - 1, p]
            * exp(neg_mort[i - 1, p] - f_a_t[a - 1, i - 1]);
        }
      }
    }
    return output;
  }
  /**
   * @title Adding process error
   *
   * @description
   * 
   * @param lrec log recruitment for each site
   * @param z "process error"
   * 
   * @return an array of numbers by age, year and patch
   */
  matrix add_pe(vector lrec, vector z) {
    array[3] int dimensions;
    int nt = num_elements(z);
    int np = num_elements(lrec) %/% nt; // %/% is integer division!
    matrix[nt, np] output = to_matrix(lrec, nt, np);
    for (p in 1:np) {
      output[, p] += z;
    }
    return output;
  }

  /**
   * @title Applying movement
   *
   * @description
   * 
   * @param lambda array of number of individuals per age, time, and patch
   * @param M movement matrix
   * @param mov_age ages at which movement starts (this can be generalized)
   * 
   * @return an array of numbers by age, year and patch
   */
  array[] matrix apply_movement(array[] matrix lambda, matrix M,
                                array[] int mov_age) {
    array[3] int dimensions;
    dimensions = dims(lambda);
    array[dimensions[1]] matrix[dimensions[2], dimensions[3]] output;
    output = lambda;
    for (a in 1:dimensions[1]) {
      if (mov_age[a]) {
        for (time in 1:dimensions[2]) {
          output[a, time, 1:dimensions[3]] =
            lambda[a, time, 1:dimensions[3]] * M';
        }
      }
    }
    return output;
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
  int n_ages; // number of ages
  int n_patches; // number of patches
  int n_time; // years for training
  vector[N] y;
  //--- toggles ---
  int<lower = 0, upper = 1> time_ar;
  int<lower = 0, upper = 1> movement;
  int<lower = 0, upper = 1> est_mort; // estimate mortality?
  int<lower = 0, upper = 1> est_init; // estimate "initial cohort"
  int<lower = 0, upper = 1> cloglog; // use cloglog instead of logit for rho
  int<lower = 0, upper = 1> qr_t; // use qr parametrization for rho?
  int<lower = 0, upper = 1> qr_r; // use qr parametrization for logrec?
  int<lower = 0, upper = 1> qr_m; // use qr parametrization for mortality?
  int<lower = 0, upper = 4> likelihood; // (0 = Original LN, 1 = repar LN, 2 =
                                        // Gamma, 3 = log-Logistic. 4 = truncated normal)
  int<lower = 0, upper = 5> init_type;
  //--- suitability (for rho) ----
  int<lower = 1> K_t;
  matrix[N, K_t] X_t;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  array[est_mort ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  matrix[movement ? n_patches: 1, movement ? n_patches : 1] adj_mat;
  array[movement ? n_ages : 0] int ages_movement;
  vector[n_ages] selectivity_at_age;
  //--- initial cohort (if not estimated) ----
  array[est_init ? 0 : n_ages - 1] real init_data;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_mort ? 1 : 0] int<lower = 1> K_m;
  matrix[est_mort ? N : 1, est_mort ? K_m[1] : 1] X_m;
  //--- * for recruitment ----
  int<lower = 1> K_r;
  matrix[N, K_r] X_r;
  //--- priors hyperparameters ----
  real pr_sigma_obs_mu; 
  real pr_sigma_obs_sd; // formerly sigma_obs_cv
  real pr_phi_mu; // revise this. It is for phi in gamma, loglogistic, and
  // inverse-gaussian
  real pr_phi_sd;
  real pr_logsd_r_mu; 
  real pr_logsd_r_sd;
  real pr_alpha_a; 
  real pr_alpha_b;
  real pr_zeta_a; 
  real pr_zeta_b;
  vector[K_t] pr_coef_t_mu;
  vector[K_t] pr_coef_t_sd;
  vector[est_mort ? K_m[1] : 0] pr_coef_m_mu;
  vector[est_mort ? K_m[1] : 0] pr_coef_m_sd;
  vector[K_r] pr_coef_r_mu;
  vector[K_r] pr_coef_r_sd;
}
transformed data {
  int na2 = n_ages - 1;
  //--- Movement ----
  // `identity_mat` is a "patch X patch" identity matrix used to apply the
  // "reflexive" movement.
  matrix[movement ? n_patches : 0, movement ? n_patches : 0] identity_mat;
  if (movement)
    identity_mat = identity_matrix(n_patches);
  //--- Mortality ----
  // `fixed_m` is a "time X patch" constant matrix. All of its positions are
  // equal to the fishing mortality `m`. The purpose of this matrix is the make
  // the changes in the code when not estimating `m` minimal.
  matrix[est_mort ? 0 : n_time, est_mort ? 0 : n_patches] fixed_m;
  if (!est_mort)
    fixed_m = rep_matrix(- m[1], n_time, n_patches);
  //--- thin-QR parametrization (improves sampler at almost no cost) ----
  // * detection prob
  matrix[qr_t ? N : 0,
         qr_t ? K_t : 0] Q_t;
  matrix[qr_t ? K_t : 0,
         qr_t ? K_t : 0] R_t;
  matrix[qr_t ? K_t : 0,
         qr_t ? K_t : 0] R_t_inv;
  if (qr_t) {
    Q_t = qr_thin_Q(X_t) * sqrt(N - 1);
    R_t = qr_thin_R(X_t) / sqrt(N - 1);
    R_t_inv = inverse(R_t);
  }
  // * recruitment
  matrix[qr_r ? N : 0, qr_r ? K_r : 0] Q_r;
  matrix[qr_r ? K_r : 0, qr_r ? K_r : 0] R_r;
  matrix[qr_r ? K_r : 0, qr_r ? K_r : 0] R_r_inv;
  if (qr_r) {
    Q_r = qr_thin_Q(X_r) * sqrt(N - 1);
    R_r = qr_thin_R(X_r) / sqrt(N - 1);
    R_r_inv = inverse(R_r);
  }
  // * mortality
  matrix[qr_m ? N : 0, qr_m ? K_m[1] : 0] Q_m;
  matrix[qr_m ? K_m[1] : 0, qr_m ? K_m[1] : 0] R_m;
  matrix[qr_m ? K_m[1] : 0, qr_m ? K_m[1] : 0] R_m_inv;
  if (qr_m) {
    Q_m = qr_thin_Q(X_m) * sqrt(N - 1);
    R_m = qr_thin_R(X_m) / sqrt(N - 1);
    R_m_inv = inverse(R_m);
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
  vector[K_r] coef_r0;
  // parameter associated with "encounter probability"
  vector[K_t] coef_t0;
  // coefficients for mortality/survival (it is a log-linear model)
  vector[est_mort ? K_m[1] : 0] coef_m0;
  //--- * initialization parameter ----
  array[est_init ? na2 : 0] real log_init;
  //--- * AR process parameters ----
  // conditional SD
  array[time_ar] real log_tau;
  // autocorrelation
  array[time_ar] real<lower = 0, upper = 1> shift_alpha;
  // aux latent variable
  vector[time_ar ? n_time : 0] raw;
  //--- * Movement parameter ----
  // logit of the probability of staying in the same patch
  array[movement] real<lower = 0, upper = 1> zeta;
}
transformed parameters {
  array[likelihood > 0 ? 1 : 0] real phi;
  array[likelihood > 0 ? 1 : 0] real lprior_phi;
  if (likelihood > 0) {
    phi = exp(log_phi);
    lprior_phi[1] = student_t_lpdf(log_phi[1] | 3, pr_phi_mu, pr_phi_sd);
  }
  //--- Recruitment ----
  // we are working with "log recruitment" here
  vector[N] log_rec;
  if (qr_r) {
    log_rec = Q_r * coef_r0;
  } else {
    log_rec = X_r * coef_r0;
  }
  //--- Initialization ----
  array[est_init ? na2 : 0] real init_par;
  array[est_init && init_type > 3 ? na2 : 0] real init_ord;
  if (est_init) {
    if (init_type > 3) {
      init_ord[1] = log_init[1];
      for (i in 2 : na2) {
        init_ord[i] = log_init[i - 1] + exp(log_init[i]);
      }
      init_par = exp(init_ord);
    } else {
      init_par = exp(log_init);
    }
  }
  //--- AR process ----
  array[time_ar] real tau;
  array[time_ar] real lprior_tau;
  array[time_ar] real alpha;
  array[time_ar] real lprior_alpha;
  vector[time_ar ? n_time : 0] z_t;
  vector[time_ar ? n_time : 0] lagged_z_t;
  if (time_ar) {
    tau[1] = exp(log_tau[1]);
    lprior_tau[1] = normal_lpdf(log_tau[1] | pr_logsd_r_mu, pr_logsd_r_sd);
    alpha[1] = 2.0 * shift_alpha[1] - 1.0;
    lprior_alpha[1] = beta_lpdf(shift_alpha[1] | pr_alpha_a, pr_alpha_b);
    z_t = tau[1] * raw;
    for (tp in 2:n_time) {
      lagged_z_t[tp] = z_t[tp - 1];
      z_t[tp] += alpha[1] * lagged_z_t[tp];
    }
  }
  //--- Mortality ----
  vector[est_mort ? N : 0] mortality;
  if (est_mort)
    mortality = qr_m ? (Q_m * coef_m0) : (X_m * coef_m0);
  // Expected density at specific time/patch combinations by age
  array[n_ages] matrix[n_time, n_patches] lambda;
  // filling lambda according to our "simplest model"
  lambda =
    simplest(n_patches, n_time, n_ages,
             f,
             est_mort ? to_matrix(mortality, n_time, n_patches) : fixed_m,
             est_init ? init_par : init_data,
             init_type,
             time_ar ?
             exp(add_pe(log_rec, z_t)) :
             to_matrix(exp(log_rec), n_time, n_patches));
  //--- Movement ----
  // probability of staying in the current patch
  // movement matrix
  matrix[movement ? n_patches : 0, movement ? n_patches : 0] mov_mat;
  array[movement] real lprior_zeta;
  if (movement) {
    // probability of movement is evenly distributed across neighbors
    real d = (1 - zeta[1]);
    mov_mat = zeta[1] * identity_mat;
    lprior_zeta[1] = beta_lpdf(zeta[1] | pr_zeta_a, pr_zeta_b);
    // It is super important that the adj_mat is setup correctly. Note that,
    // this matrix is "row standardized". That is, its rows add up to 1.
    mov_mat += d * adj_mat;
    lambda = apply_movement(lambda, mov_mat, ages_movement);
  }
  //--- quantities used in the likelihood ----
  // probability of encounter at specific time/patch combinations
  vector[N] rho;
  // Expected density at specific time/patch combinations
  vector[N] mu =
    rep_vector(0.0, N);
  // rho now has a "regression like" type
  if (cloglog) {
    matrix[N, K_t] X_aux;
    X_aux = qr_t ? Q_t : X_t;
    rho =
      inv_cloglog(X_aux * coef_t0);
  } else {
    matrix[N, K_t] X_aux;
    X_aux = qr_t ? Q_t : X_t;
    rho =
      inv_logit(X_aux * coef_t0);
  }
  {
    matrix[n_time, n_patches] mu_aux =
      rep_matrix(0.0, n_time, n_patches);
    //--- filling mu ----
    for (time in 1 : n_time) {
      for (p in 1 : n_patches) {
        real mu_aux2;
        mu_aux2 = dot_product(to_vector(lambda[1:n_ages, time, p]),
                              selectivity_at_age);
        if (!is_nan(mu_aux2)) {
          mu_aux[time, p] = mu_aux2;
        }
      } // close patches
    }
    mu = to_vector(mu_aux);
  }
}
// close transformed parameters block
model {
  //--- initialization parameters ----
  if (est_init) {
    target += std_normal_lpdf(log_init);
    /*
      The link below provides info about the transformation of an unordered into
  an ordered vector.
  https://mc-stan.org/docs/reference-manual/transforms.html#absolute-jacobian-determinant-of-the-ordered-inverse-transform
    */
    if (init_type > 3)
      target += sum(init_ord[2 : na2]);
  }
  //--- AR process ----
  if (time_ar) {
    target += std_normal_lpdf(raw);
    target += lprior_tau[1];
    target += lprior_alpha[1];
  }
  //--- Movement ----
  if (movement) {
    target += lprior_zeta[1];
  }
  //--- Mortality ----
  if (est_mort)
    target += normal_lpdf(coef_m0 | pr_coef_m_mu, pr_coef_m_sd);
  //--- Recruitment ----
  target += normal_lpdf(coef_r0 | pr_coef_r_mu, pr_coef_r_sd);
  //--- suitability ----
  target += normal_lpdf(coef_t0 | pr_coef_t_mu, pr_coef_t_sd);
  //--- Likelihood ----
  if (likelihood == 0) {
    target += normal_lpdf(sigma_obs[1] | pr_sigma_obs_mu, pr_sigma_obs_sd) -
      1.0 * normal_lccdf(0 | pr_sigma_obs_mu, pr_sigma_obs_sd);
  } else {
    // change these parameters (PC prior for exponential?)
    target += lprior_phi[1];
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
  } else if (likelihood == 3) {
    vector[N_nz] a_ll;
    real b_ll;
    b_ll = phi[1] + 1;
    a_ll = sin(pi() / b_ll) * mu[id_nz] * inv(pi() * b_ll);
    target += log1m(rho[id_nz]);
    target += loglogistic_lpdf(y[id_nz] | a_ll, b_ll);
  } else {
    target += log1m(rho[id_nz]);
    target += normal_lpdf(y | mu, phi[1]) -
      normal_lccdf(rep_vector(0.0, N) | mu, phi[1]);
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_pp;
  vector[K_t] coef_t;
  vector[K_r] coef_r;
  vector[est_mort ? K_m[1] : 0] coef_m;
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
      a_ll = sin(pi() / b_ll) * mu[n] * inv(pi() * b_ll);
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        loglogistic_rng(a_ll, b_ll);
      if (y[n] == 0) {
        log_lik[n] = log(rho[n]);
      } else {
        log_lik[n] = log1m(rho[n]) +
          loglogistic_lpdf(y[n] | a_ll, b_ll);
      }
    } else {
      array[2] real aux_tn = rep_array(0.0, 2);
      aux_tn[2] = normal_rng(mu[n], phi[1]);
      y_pp[n] = (1 - bernoulli_rng(rho[n])) *
        max(aux_tn);
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
  if (qr_t) {
    coef_t = R_t_inv * coef_t0;
  } else coef_t = coef_t0;
  if (qr_r) {
    coef_r = R_r_inv * coef_r0;
  } else coef_r = coef_r0;
  if (qr_m) {
    coef_m = R_m_inv * coef_m0;
  } else if (est_mort) {
    coef_m = coef_m0;
  }
}
