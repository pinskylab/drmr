functions {
#include utils/age_struct.stan
}
data {
  //--- survey data  ---
  int N; // n_patches * n_time
  int n_ages; // number of ages
  int n_patches; // number of patches
  int n_time; // years for training
  //--- toggles ---
  int<lower = 0, upper = 1> movement;
  int<lower = 0, upper = 1> est_surv; // estimate mortality?
  int<lower = 0, upper = 1> est_init; // estimate "initial cohort"
  int<lower = 0, upper = 1> minit;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  array[est_surv ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  matrix[movement ? n_patches: 1, movement ? n_patches : 1] adj_mat;
  array[movement ? n_ages : 0] int ages_movement;
  vector[n_ages] selectivity_at_age;
  //--- initial cohort (if not estimated) ----
  array[est_init ? 0 : n_ages - 1] real init_data;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
}
transformed data {
  matrix[est_surv ? 0 : n_time, est_surv ? 0 : n_patches] fixed_m;
  if (!est_surv)
    fixed_m = rep_matrix(- m[1], n_time, n_patches);
  matrix[movement ? n_patches : 0, movement ? n_patches : 0] identity_mat;
  if (movement)
    identity_mat = identity_matrix(n_patches);
}
parameters {
  // coefficients for recruitment (it is a log-linear model)
  vector[N] log_rec;
  // coefficients for mortality/survival (it is a log-linear model)
  array[est_surv] vector[est_surv ? K_m[1] : 0] beta_s;
  //--- * movement ----
  array[movement] real<lower = 0, upper = 1> zeta;
  //--- * initialization parameter ----
  array[est_init ? n_ages - 1 : 0] real log_init;
}
transformed parameters {
}
generated quantities {
  // Expected density at specific time/patch combinations
  array[n_ages] matrix[n_time, n_patches] lambda;
  {
    //--- Initialization ----
    array[est_init ? n_ages - 1 : 0] real init_par;
    if (est_init)
      init_par = exp(log_init);
    //--- Mortality ----
    vector[est_surv ? N : 0] mortality;
    if (est_surv)
      mortality = X_m * beta_s[1];

    // filling lambda according to our "simplest model"
    lambda =
      simplest(n_patches, n_time, n_ages,
               f,
               est_surv ? to_matrix(mortality, n_time, n_patches) : fixed_m,
               est_init ? init_par : init_data,
               to_matrix(log_rec, n_time, n_patches),
               minit);
  }
  //--- Movement ----
  if (movement) {
    matrix[movement ? n_patches : 0, movement ? n_patches : 0] mov_mat;
    real d = (1 - zeta[1]);
    mov_mat = zeta[1] * identity_mat;
    mov_mat += d * adj_mat;
    lambda =
      apply_movement(lambda, mov_mat, ages_movement);
  }
}
