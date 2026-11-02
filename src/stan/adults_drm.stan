functions {
#include utils/age_struct.stanfunctions
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
  int<lower = 0, upper = 2> rec_dd;   // 0 for Ricker, 1 for Beverton-Holt,
                                      // 2 for none
  int<lower = 0, upper = 1> acc_dd;
  int<lower = 0, upper = 3> ar_re;
  int<lower = 0, upper = 3> iid_re;
  int<lower = 0, upper = 3> sp_re;
  //--- fish mortality data ----
  matrix[n_ages, n_time] f;
  array[est_surv ? 0 : 1] real m; // total mortality
  //--- movement related quantities ----
  matrix[movement ? n_sites: 1, movement ? n_sites : 1] adj_mat;
  int<lower = 0> n_edges_adj;
  array[movement ? n_ages : 0] int ages_movement;
  //--- maturity and weight at age ----
  vector[n_ages] amat;
  vector[n_ages] weight;
  //--- initial cohort (if not estimated) ----
  array[est_init ? 0 : n_ages - 1] real init_data;
  //--- environmental data ----
  //--- * for mortality ----
  array[est_surv ? 1 : 0] int<lower = 1> K_m;
  matrix[est_surv ? N : 1, est_surv ? K_m[1] : 1] X_m;
}
transformed data {
  matrix[est_surv ? 0 : n_time, est_surv ? 0 : n_sites] fixed_m;
  if (!est_surv)
    fixed_m = rep_matrix(- m[1], n_time, n_sites);
  vector[movement ? n_edges_adj : 0] w_adj;
  array[movement ? n_edges_adj : 0] int v_adj;
  array[movement ? n_sites + 1 : 0] int u_adj;
  matrix[movement ? n_sites : 0, movement ? n_sites : 0] adj_mat_std;
  if (movement) {
    for (i in 1:n_sites) {
      real row_s = sum(adj_mat[i]);
      if (row_s > 0) {
        adj_mat_std[i] = adj_mat[i] / row_s;
      } else {
        adj_mat_std[i] = adj_mat[i];
      }
    }
    w_adj = csr_extract_w(adj_mat_std);
    v_adj = csr_extract_v(adj_mat_std);
    u_adj = csr_extract_u(adj_mat_std);
  }
}
parameters {
  // log recruitment (including its random effects)
  vector[N] log_rec;
  vector[ar_re > 0 ? n_time : 0] z_t;
  array[iid_re > 0 ? 1 : 0] vector[n_sites] z_i;
  vector[sp_re > 0 ? n_sites : 0] z_s;
  // coefficients for mortality/survival (it is a log-linear model)
  array[est_surv] vector[est_surv ? K_m[1] : 0] beta_s;
  //--- * movement ----
  array[movement] real<lower = 0, upper = 1> zeta;
  //--- * density-dependence ----
  array[rec_dd < 2 ? 1 : 0] real kappa;
  //--- * initialization parameter ----
  array[est_init ? n_ages - 1 : 0] real log_init;
}
generated quantities {
  // Density (or biomass, when weight at age is informed) of adults in
  // reproductive age. It follows the same ordering as the data used for
  // fitting the model (i.e., by site and then time).
  vector[N] adults;
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
    //--- Expected density by age ----
    array[n_ages] matrix[n_time, n_sites] lambda;
    if (rec_dd < 2) {
      if (movement) {
        lambda =
          pop_rec_dd_movement(n_sites, n_time, n_ages,
                              f,
                              mortality,
                              est_init ? init_par : init_data,
                              to_matrix(log_rec, n_time, n_sites),
                              amat, weight,
                              kappa[1], rec_dd,
                              zeta[1], w_adj, v_adj, u_adj,
                              ages_movement,
                              acc_dd);
      } else {
        lambda =
          pop_rec_dd(n_sites, n_time, n_ages,
                     f,
                     mortality,
                     est_init ? init_par : init_data,
                     to_matrix(log_rec, n_time, n_sites),
                     amat, weight,
                     kappa[1], rec_dd,
                     acc_dd);
      }
    } else {
      if (movement) {
        lambda =
          simplest_movement(n_sites, n_time, n_ages,
                            f,
                            mortality,
                            est_init ? init_par : init_data,
                            to_matrix(log_rec, n_time, n_sites),
                            minit,
                            zeta[1], w_adj, v_adj, u_adj,
                            ages_movement);
      } else {
        lambda =
          simplest(n_sites, n_time, n_ages,
                   f,
                   mortality,
                   est_init ? init_par : init_data,
                   to_matrix(log_rec, n_time, n_sites),
                   minit);
      }
    }
    //--- Adults (same quantity driving the density-dependence) ----
    matrix[n_time, n_sites] adults_aux =
      rep_matrix(0.0, n_time, n_sites);
    for (a in 1:n_ages) {
      adults_aux += lambda[a] * (amat[a] * weight[a]);
    }
    adults = to_vector(adults_aux);
  }
}
