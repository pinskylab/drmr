/**
 * @title Forecast theoretical mean according to the simplest model possible
 *
 * @description
 * 
 * @param n_patches number of patches
 * @param n_time number of years of training data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param init initialization (currently with recruitment)
 * @param lambda_past a n_ages by n_patches matrix containing the expected
 * density at the last training year.
 * @param f_past a n_ages by n_time_train matrix
 * @param neg_mort_past a n_patches vector
 * 
 * @return an array of numbers by age, year and patch
 */
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
  past_last_time = cols(f_past);
  /* output[1] += init; */
  for (i in 1 : n_time) {
    for (p in 1 : n_patches) {
      for (a in 2 : n_ages) {
        if (i == 1) {
          output[a, i, p] = log(lambda_past[a - 1, p])
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
 * @title Forecast theoretical mean according to the simplest model possible with movement
 *
 * @description Mechanistic version
 * 
 * @param n_patches number of patches
 * @param n_time number of years of training data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param recruitment a n_time by n_patches matrix;
 * @param lambda_past a n_ages by n_patches matrix containing the expected
 * density at the last training year.
 * @param f_past a n_ages by n_time_train matrix
 * @param neg_mort_past a n_patches vector
 * @param zeta probability of staying in the current site
 * @param w_adj sparse CSR vector of non-zero entries of adjacency matrix
 * @param v_adj sparse CSR array of column indices
 * @param u_adj sparse CSR array of row starting indices
 * @param mov_age ages at which movement starts
 * 
 * @return an array of numbers by age, year and patch
 */
array[] matrix forecast_simplest_movement(int n_patches,
                                          int n_time,
                                          int n_ages,
                                          matrix f_a_t,
                                          matrix neg_mort,
                                          matrix recruitment,
                                          matrix lambda_past,
                                          matrix f_past,
                                          vector neg_mort_past,
                                          real zeta,
                                          vector w_adj,
                                          array[] int v_adj,
                                          array[] int u_adj,
                                          array[] int mov_age) {
  // initializing output with zeros
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
  int past_last_time = cols(f_past);
  
  for (i in 1 : n_time) {
    // Recruitment
    output[1, i] = exp(recruitment[i]);
    
    for (a in 2 : n_ages) {
      row_vector[n_patches] lambda_prev;
      row_vector[n_patches] surv;
      if (i == 1) {
        lambda_prev = lambda_past[a - 1];
        surv = exp(to_row_vector(neg_mort_past) - f_past[a - 1, past_last_time]);
      } else {
        lambda_prev = output[a - 1, i - 1];
        surv = exp(neg_mort[i - 1] - f_a_t[a - 1, i - 1]);
      }
      
      row_vector[n_patches] lambda_surv = lambda_prev .* surv;
      
      if (mov_age[a]) {
        vector[n_patches] adj_x = csr_matrix_times_vector(n_patches, n_patches, w_adj, v_adj, u_adj, lambda_surv');
        output[a, i] = (zeta * lambda_surv' + (1 - zeta) * adj_x)';
      } else {
        output[a, i] = lambda_surv;
      }
    }
  }
  return output;
}
