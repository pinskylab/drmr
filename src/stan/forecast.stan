functions {
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
   * @param init a n_time by n_patches matrix of initialization for age-group
   *   "1". In the simplest case, it can be the recruitment.
   * 
   * @return an array of numbers by age, year and patch
   */
  array[] matrix simplest(int n_patches,
                          int n_time,
                          int n_ages,
                          // Mortality parameter
                          matrix f_a_t,
                          matrix neg_mort,
                          // initialization (currently with recruitment)
                          matrix init) {
    // initializing output with zeros
    array[n_ages] matrix[n_time, n_patches] output
      = rep_array(init, n_ages);
    /* output[1] += init; */
    for (i in 2 : n_time) {
      for (p in 1 : n_patches) {
        for (a in 2 : n_ages) {
          output[a, i, p] = output[a - 1, i - 1, p] +
            neg_mort[i - 1, p] - f_a_t[a - 1, i - 1];
        }
      }
    }
    return exp(output);
  }

  array[] matrix forecast_simplest(int n_patches,
                                   int n_time,
                                   int n_ages,
                                   // Mortality parameter
                                   matrix f_a_t,
                                   matrix neg_mort,
                                   // initialization (currently with recruitment)
                                   matrix init,
                                   // from past
                                   matrix lambda_past,
                                   matrix f_past,
                                   vector neg_mort_past) {
    // initializing output with zeros
    array[n_ages] matrix[n_time, n_patches] output
      = rep_array(init, n_ages);
    int past_last_time;
    past_last_time = cols(f_past[1]);
    /* output[1] += init; */
    for (i in 1 : n_time) {
      for (p in 1 : n_patches) {
        for (a in 2 : n_ages) {
          if (i == 1) {
            output[a, i, p] = lambda_past[a - 1, p]
              + neg_mort_past[p] - f_past[a - 1, past_last_time];
          } else {
            output[a, i, p] = output[a - 1, i - 1, p]
              + neg_mort[i - 1, p] - f_a_t[a - 1, i - 1];
          }
        }
      }
    }
    return exp(output);
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
   * @param mov_age ages at which movement happens (this can be generalized)
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
  int<lower = 0, upper = 1> time_ar;
  int<lower = 0, upper = 1> movement;
  int<lower = 0, upper = 1> est_surv; // estimate mortality?
  int<lower = 0, upper = 1> cloglog; // use cloglog instead of logit for rho
  int<lower = 0, upper = 4> likelihood; // (0 = Original LN, 1 = repar LN, 2 =
                                        // Gamma, 3 = log-Logistic, 4 = truncated normal)
  //--- suitability (for rho) ----
  int<lower = 1> K_t;
  matrix[N, K_t] X_t;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  matrix[n_ages, n_time_train] f_past;
  array[est_surv ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  array[movement ? n_ages : 0] int ages_movement;
  vector[n_ages] selectivity_at_age;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
  matrix[est_surv ? n_patches : 1, est_surv ? K_m[1] : 1] X_m_past;
  //--- * for recruitment ----
  int<lower = 1> K_r;
  matrix[N, K_r] X_r;
}
parameters {
  //--- "regression" coefficients for absence and recr ----
  vector[K_t] beta_t;
  vector[K_r] beta_r;
  //--- pop dyn parameters ----
  matrix[n_ages, n_patches] lambda;
  //--- additional parameter for different lik functions ----
  array[likelihood > 0 ? 1 : 0] real phi;
  array[likelihood == 0 ? 1 : 0] real<lower = 0> sigma_obs;
  //--- movement matrix ---
  matrix[movement ? n_patches : 0, movement ? n_patches : 0] mov_mat;
  //--- reg for mortality ---
  vector[est_surv ? K_m[1] : 0] beta_s;
  //--- parameters from AR process ----
  vector[time_ar ? n_time_train : 0] z_t;
  array[time_ar] real alpha;
  array[time_ar ? 1 : 0] real tau;
}
generated quantities {
  //--- projected expected density by age ---
  array[n_ages] matrix[n_time, n_patches] lambda_proj;
  //--- projected expected density ----
  vector[n_time * n_patches] mu_proj;
  //--- projected total density ----
  vector[n_time * n_patches] y_proj;
  //--- projected absence probability ----
  vector[n_time * n_patches] rho_proj;
  //--- AR term ----
  vector[time_ar ? n_time : 0] z_tp;
  if (time_ar) {
    {
      vector[n_time] raw;
      vector[time_ar ? n_time - 1 : 0] lagged_rec;
      raw[1] = std_normal_rng();
      z_tp[1] = alpha[1] * z_t[n_time_train] +
        tau[1] * raw[1];
      for (tp in 2:n_time) {
        raw[tp] = std_normal_rng();
        lagged_rec[tp - 1] = z_tp[tp - 1];
        z_tp[tp] = alpha[1] * lagged_rec[tp - 1] +
          tau[1] * raw[tp];
      }
    }
  }
  //--- lambda_proj calculations ----
  {
    vector[N] log_rec;
    log_rec = X_r * beta_r;
    if (time_ar) {
      for (n in 1:N)
        log_rec[n] += z_tp[time[n]];
    }
    vector[n_patches] past_m;
    matrix[n_time, n_patches] current_m;
    if (!est_surv) {
      current_m = rep_matrix(- m[1], n_time, n_patches);
      past_m = rep_vector(- m[1], n_patches);
    } else {
      current_m = to_matrix(X_m * beta_s, n_time, n_patches);
      past_m = X_m_past * beta_s;
    }
    // forecast_simplest is a function in the utils/theoretical_mean.stan file
    lambda_proj = forecast_simplest(n_patches,
                                    n_time,
                                    n_ages,
                                    f,
                                    current_m,
                                    to_matrix(log_rec,
                                              n_time,
                                              n_patches),
                                    lambda,
                                    f_past,
                                    past_m);
    if (movement) {
      lambda_proj =
        apply_movement(lambda_proj, mov_mat, ages_movement);
    }
    //--- mu_proj calculations ----
    matrix[n_time, n_patches] mu_aux =
      rep_matrix(0.0, n_time, n_patches);
    //--- filling mu ----
    for (tp in 1 : n_time) {
      for (p in 1 : n_patches) {
        real mu_aux2;
        mu_aux2 = dot_product(to_vector(lambda_proj[1:n_ages, tp, p]),
                              selectivity_at_age);
        if (!is_nan(mu_aux2)) {
          mu_aux[tp, p] = mu_aux2;
        }
      } // close patches
    }
    mu_proj = to_vector(mu_aux);
  }
  //--- absence probabilities ----
  if (cloglog) {
    rho_proj =
      inv_cloglog(X_t * beta_t);
  } else {
    rho_proj =
      inv_logit(X_t * beta_t);
  }
  //--- y_proj calculations ----
  for (n in 1:N) {
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu_proj[n]) + square(sigma_obs[1]) / 2;
      y_proj[n] = (1 - bernoulli_rng(rho_proj[n])) *
        lognormal_rng(loc_par, sigma_obs[1]);
    } else if (likelihood == 1) {
      real mu_ln;
      real sigma_ln;
      sigma_ln = sqrt(log1p(phi[1] * inv_square(mu_proj[n])));
      mu_ln = log(square(mu_proj[n]) * inv_sqrt(square(mu_proj[n]) + phi[1]));
      y_proj[n] = (1 - bernoulli_rng(rho_proj[n])) *
        lognormal_rng(mu_ln, sigma_ln);
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi[1] / mu_proj[n];
      y_proj[n] = (1 - bernoulli_rng(rho_proj[n])) *
        gamma_rng(phi[1], gamma_beta);
    } else if (likelihood == 3) {
      real a_ll;
      real b_ll;
      b_ll = phi[1] + 1;
      a_ll = sin(pi() / b_ll) * mu_proj[n] * inv(pi() * b_ll);
      y_proj[n] = (1 - bernoulli_rng(rho_proj[n])) *
        loglogistic_rng(a_ll, b_ll);
    } else {
      array[2] real aux_tn = rep_array(0.0, 2);
      aux_tn[2] = normal_rng(mu_proj[n], phi[1]);
      y_proj[n] = (1 - bernoulli_rng(rho_proj[n])) *
        max(aux_tn);
    }
  }
}
