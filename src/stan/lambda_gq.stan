functions {
#include utils/age_struct.stan
#include utils/age_st_forecast.stan
}
data {
  //--- survey data  ---
  int N; // n_sites * n_time
  int n_ages; // number of ages
  int n_sites; // number of sites
  int n_time; // years for training
  array[N] int time;
  array[N] int site;
  //--- toggles ---
  int<lower = 0, upper = 1> movement;
  int<lower = 0, upper = 1> est_surv; // estimate mortality?
  int<lower = 0, upper = 1> est_init; // estimate "initial cohort"
  int<lower = 0, upper = 1> minit;
  int<lower = 0, upper = 3> ar_re;
  int<lower = 0, upper = 3> iid_re;
  int<lower = 0, upper = 3> sp_re;
  int<lower = 0, upper = 1> proj;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  array[est_surv ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  matrix[movement ? n_sites: 1, movement ? n_sites : 1] adj_mat;
  array[movement ? n_ages : 0] int ages_movement;
  vector[n_ages] selectivity_at_age;
  //--- initial cohort (if not estimated) ----
  array[est_init ? 0 : n_ages - 1] real init_data;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
  //--- variables for projection ----
  array[proj] int n_proj;
  array[proj ? n_sites * n_proj[1] : 0] int time_proj;
  array[proj ? n_sites * n_proj[1] : 0] int site_proj;
  matrix[proj ? n_ages : 0, proj ? n_proj[1] : 0] f_proj;
  int<lower = 1> K_r;
  matrix[proj ? n_proj[1] * n_sites : 0, proj ? K_r : 0] X_r;
  matrix[est_surv * proj ? n_proj[1] * n_sites : 1, est_surv * proj ? K_m[1] : 1] X_mproj;
}
transformed data {
  matrix[est_surv ? 0 : n_time, est_surv ? 0 : n_sites] fixed_m;
  if (!est_surv)
    fixed_m = rep_matrix(- m[1], n_time, n_sites);
  matrix[movement ? n_sites : 0, movement ? n_sites : 0] identity_mat;
  if (movement)
    identity_mat = identity_matrix(n_sites);
  array[proj] int N_proj;
  if (proj)
    N_proj[1] = n_proj[1] * n_sites;
}
parameters {
  // coefficients for recruitment (it is a log-linear model)
  vector[N] log_rec;
  vector[proj ? K_r : 0] beta_r;
  vector[ar_re > 0 ? n_time : 0] z_t;
  array[proj * ar_re > 0 ? 1 : 0] real alpha;
  array[proj * ar_re > 0 ? 1 : 0] real sigma_t;
  array[iid_re > 0 ? 1 : 0] vector[n_sites] z_i;
  vector[sp_re > 0 ? n_sites : 0] z_s;
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
  // Expected density at specific time/site combinations
  array[n_ages] matrix[n_time, n_sites] lambda;
  array[proj ? n_ages : 0] matrix[proj ? n_proj[1] : 0, proj ? n_sites : 0] lambda_proj;
  {
    //--- Initialization ----
    array[est_init ? n_ages - 1 : 0] real init_par;
    if (est_init)
      init_par = log_init;
    //--- Mortality ----
    matrix[n_time, n_sites] mortality;
    if (est_surv) {
      vector[N] m_aux;
      m_aux = X_m * beta_s[1];
      if (ar_re == 2) {
        for (n in 1:N)
          m_aux[n] += z_t[time[n]];
      }
      if (iid_re == 2) {
        for (n in 1:N)
          m_aux[n] += z_i[1][site[n]];
      }
      if (sp_re == 2) {
        for (n in 1:N)
          m_aux[n] += z_s[site[n]];
      }
      mortality = to_matrix(-log1p(exp(-m_aux)), n_time, n_sites);
    } else {
      mortality = fixed_m;
    }
    // filling lambda according to our "simplest model"
    lambda =
      simplest(n_sites, n_time, n_ages,
               f,
               mortality,
               est_init ? init_par : init_data,
               to_matrix(log_rec, n_time, n_sites),
               minit);
    if (proj) {
      vector[ar_re > 0 ? n_proj[1] : 0] z_tp;
      if (ar_re > 0) {
        {
          vector[n_proj[1]] w_t;
          for (tp in 1:n_proj[1]) w_t[tp] = std_normal_rng();
          z_tp[1] = alpha[1] * z_t[n_time] +
            sigma_t[1] * w_t[1];
          for (tp in 2:n_proj[1]) {
            z_tp[tp] = alpha[1] * z_tp[tp - 1] +
              sigma_t[1] * w_t[tp];
          }
        }
      }
      vector[N_proj[1]] log_rec_proj;
      log_rec_proj = X_r * beta_r;
      if (ar_re == 1) {
        for (n in 1:N_proj[1])
          log_rec_proj[n] += z_tp[time_proj[n]];
      }
      if (iid_re == 1) {
        for (n in 1:N_proj[1])
          log_rec_proj[n] += z_i[1][site_proj[n]];
      }
      if (sp_re == 1) {
        for (n in 1:N_proj[1])
          log_rec_proj[n] += z_s[site_proj[n]];
      }
      matrix[n_proj[1], n_sites] current_m;
      vector[n_sites] past_m;
      if (!est_surv) {
        current_m = rep_matrix(- m[1], n_proj[1], n_sites);
        past_m = rep_vector(- m[1], n_sites);
      } else {
        vector[N_proj[1]] m_aux_proj;
        m_aux_proj = X_mproj * beta_s[1];
        if (ar_re == 2) {
          for (n in 1:N_proj[1])
            m_aux_proj[n] += z_tp[time_proj[n]];
        }
        if (iid_re == 2) {
          for (n in 1:N_proj[1])
            m_aux_proj[n] += z_i[1][site_proj[n]];
        }
        if (sp_re == 2) {
          for (n in 1:N_proj[1])
            m_aux_proj[n] += z_s[site_proj[n]];
        }
        current_m = to_matrix(-log1p(exp(-m_aux_proj)), n_proj[1], n_sites);
        // Reconstruction of past_m (last year of training)
        past_m = mortality[n_time]';
      }
      matrix[n_ages, n_sites] lambda_last;
      for (a in 1:n_ages) {
        lambda_last[a] = lambda[a, n_time]; 
      }
      lambda_proj = forecast_simplest(n_sites,
                                      n_proj[1],
                                      n_ages,
                                      f_proj,
                                      current_m,
                                      to_matrix(log_rec_proj,
                                                n_proj[1],
                                                n_sites),
                                      lambda_last,
                                      f,
                                      past_m);
    }
  }
  //--- Movement ----
  if (movement) {
    matrix[movement ? n_sites : 0, movement ? n_sites : 0] mov_mat;
    real d = (1 - zeta[1]);
    mov_mat = zeta[1] * identity_mat;
    mov_mat += d * adj_mat;
    lambda =
      apply_movement(lambda, mov_mat, ages_movement);
    if (proj)
      lambda_proj =
        apply_movement(lambda_proj, mov_mat, ages_movement);
  }
}
