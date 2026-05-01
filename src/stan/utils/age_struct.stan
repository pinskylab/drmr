/**
 * @title Generate theoretical mean according to the simplest model possible
 *
 * @description
 * 
 * @param n_patches number of patches
 * @param n_time number of years of training data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param init a n_ages - 1 array.
 * @param recruitment a n_time by n_patches matrix;
 * @param init_type a n_time by n_patches matrix;
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
                        matrix recruitment,
                        int minit) {
  /*
    This function could be updated to have patches and time as the input as
    opposed to n_patches and n_time. The advantage would be that the ordering
    of these variables in the input dataset would not be important. The would
    come at the cost of increasing the complexity of the function.
   */
  // initializing output with zeros
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
  output[1] = recruitment;
  if (minit) {
    for (p in 1:n_patches) {
      for (a in 2:n_ages) {
        output[a, 1, p] = output[a - 1, 1, p] +
          neg_mort[1, p] - f_a_t[a - 1, 1];
      }
    }
  } else {
    for (a in 1 : (n_ages - 1)) {
      output[a + 1, 1, ] = rep_row_vector(init[a], n_patches);
    }
  }
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
  int n_ages = size(lambda);
  array[n_ages] matrix[rows(lambda[1]), cols(lambda[1])] output = lambda;
  for (a in 1:n_ages) {
    if (mov_age[a]) {
      output[a] = lambda[a] * M';
    }
  }
  return output;
}

/**
 * @title Generate theoretical mean according to the simplest model possible with movement
 *
 * @description This version includes mechanistic movement.
 * 
 * @param n_patches number of patches
 * @param n_time number of years of training data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param init a n_ages - 1 array.
 * @param recruitment a n_time by n_patches matrix;
 * @param minit 1 if mortality is stable at the beginning
 * @param zeta probability of staying in the current site
 * @param w_adj sparse CSR vector of non-zero entries of adjacency matrix
 * @param v_adj sparse CSR array of column indices
 * @param u_adj sparse CSR array of row starting indices
 * @param mov_age ages at which movement starts
 * 
 * @return an array of numbers by age, year and patch
 */
array[] matrix simplest_movement(int n_patches,
                                 int n_time,
                                 int n_ages,
                                 matrix f_a_t,
                                 matrix neg_mort,
                                 array[] real init,
                                 matrix recruitment,
                                 int minit,
                                 real zeta,
                                 vector w_adj,
                                 array[] int v_adj,
                                 array[] int u_adj,
                                 array[] int mov_age) {
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
    
  // Time t=1
  output[1, 1] = exp(recruitment[1]);
  if (minit) {
    for (p in 1:n_patches) {
      for (a in 2:n_ages) {
        output[a, 1, p] = output[a - 1, 1, p] *
          exp(neg_mort[1, p] - f_a_t[a - 1, 1]);
      }
    }
  } else {
    for (a in 1 : (n_ages - 1)) {
      output[a + 1, 1] = rep_row_vector(exp(init[a]), n_patches);
    }
  }
  
  for (i in 2 : n_time) {
    // Recruitment at time i
    output[1, i] = exp(recruitment[i]);
    
    // Survival and Movement from time i-1 to i
    for (a in 2 : n_ages) {
      // Survival first (at each patch)
      row_vector[n_patches] surv = exp(neg_mort[i - 1] - f_a_t[a - 1, i - 1]);
      row_vector[n_patches] lambda_surv = output[a - 1, i - 1] .* surv;
      
      if (mov_age[a]) {
        // Mechanistic movement: survivors move
        vector[n_patches] adj_x =
          csr_matrix_times_vector(n_patches, n_patches, w_adj,
                                  v_adj, u_adj, lambda_surv');
        output[a, i] = (zeta * lambda_surv' + (1 - zeta) * adj_x)';
      } else {
        output[a, i] = lambda_surv;
      }
    }
  }
  return output;
}

/**
 * @title Generate theoretical mean according to a density-dependent recruitment model
 *
 * @description Mechanistic version with endogenous recruitment (Ricker or Beverton-Holt).
 * 
 * @param n_patches number of patches
 * @param n_time number of years of data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param init initialization for ages 2 to n_ages at time t=1
 * @param recruitment_env log-productivity (alpha) matrix [n_time, n_patches]
 * @param mat maturity-at-age vector [n_ages]
 * @param weight weight-at-age vector [n_ages]
 * @param beta density-dependence coefficient
 * @param rec_type 0 for Ricker, 1 for Beverton-Holt
 * 
 * @return an array of numbers by age, year and patch
 */
array[] matrix pop_rec_dd(int n_patches,
                          int n_time,
                          int n_ages,
                          matrix f_a_t,
                          matrix neg_mort,
                          array[] real init,
                          matrix recruitment_env,
                          vector mat,
                          vector weight,
                          real beta,
                          int rec_type) {
  // Initializing output with zeros
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
    
  // Time t=1
  // Recruitment at t=1 is purely environmental (no previous stock known)
  output[1, 1] = exp(recruitment_env[1]);
  
  // Initialization of other ages at t=1
  for (a in 1 : (n_ages - 1)) {
    output[a + 1, 1] = rep_row_vector(exp(init[a]), n_patches);
  }
  
  for (i in 2 : n_time) {
    for (p in 1 : n_patches) {
      // 1. Calculate Spawning Stock Biomass (SSB) at time i-1
      real ssb_prev = 0;
      for (a in 1 : n_ages) {
        ssb_prev += output[a, i - 1, p] * mat[a] * weight[a];
      }
      
      // 2. Density-dependent Recruitment
      // recruitment_env[i, p] is log(alpha)
      if (ssb_prev > 1e-10) {
        real log_S = log(ssb_prev);
        if (rec_type == 0) {
          // Ricker: R = alpha * S * exp(-beta * S)
          output[1, i, p] = exp(recruitment_env[i, p] + log_S - beta * ssb_prev);
        } else {
          // Beverton-Holt: R = (alpha * S) / (1 + beta * S)
          output[1, i, p] = exp(recruitment_env[i, p] + log_S - log1p(beta * ssb_prev));
        }
      } else {
        output[1, i, p] = 0.0;
      }
      
      // 3. Survival transition
      for (a in 2 : n_ages) {
        output[a, i, p] = output[a - 1, i - 1, p] *
          exp(neg_mort[i - 1, p] - f_a_t[a - 1, i - 1]);
      }
    }
  }
  
  return output;
}

/**
 * @title Generate theoretical mean according to a density-dependent recruitment model with movement
 *
 * @description Mechanistic version with endogenous recruitment and movement.
 * 
 * @param n_patches number of patches
 * @param n_time number of years of data
 * @param n_ages number of age classes
 * @param f_a_t fishing mortality at age "a" and time "t"
 * @param neg_mort minus natural mortality (instantaneous) rate
 * @param init initialization for ages 2 to n_ages at time t=1
 * @param recruitment_env log-productivity (alpha) matrix [n_time, n_patches]
 * @param mat maturity-at-age vector [n_ages]
 * @param weight weight-at-age vector [n_ages]
 * @param beta density-dependence coefficient
 * @param rec_type 0 for Ricker, 1 for Beverton-Holt
 * @param zeta probability of staying in the current site
 * @param w_adj sparse CSR vector of non-zero entries of adjacency matrix
 * @param v_adj sparse CSR array of column indices
 * @param u_adj sparse CSR array of row starting indices
 * @param mov_age ages at which movement starts
 * 
 * @return an array of numbers by age, year and patch
 */
array[] matrix pop_rec_dd_movement(int n_patches,
                                   int n_time,
                                   int n_ages,
                                   matrix f_a_t,
                                   matrix neg_mort,
                                   array[] real init,
                                   matrix recruitment_env,
                                   vector mat,
                                   vector weight,
                                   real beta,
                                   int rec_type,
                                   real zeta,
                                   vector w_adj,
                                   array[] int v_adj,
                                   array[] int u_adj,
                                   array[] int mov_age) {
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
    
  // Time t=1
  output[1, 1] = exp(recruitment_env[1]);
  for (a in 1 : (n_ages - 1)) {
    output[a + 1, 1] = rep_row_vector(exp(init[a]), n_patches);
  }
  
  for (i in 2 : n_time) {
    // 1. Endogenous Recruitment (happens locally in each patch)
    row_vector[n_patches] ssb_prev = rep_row_vector(0.0, n_patches);
    for (a in 1 : n_ages) {
      ssb_prev += output[a, i - 1] .* (mat[a] * weight[a]);
    }
    
    for (p in 1 : n_patches) {
      if (ssb_prev[p] > 1e-10) {
        real log_S = log(ssb_prev[p]);
        if (rec_type == 0) {
          output[1, i, p] = exp(recruitment_env[i, p] + log_S - beta * ssb_prev[p]);
        } else {
          output[1, i, p] = exp(recruitment_env[i, p] + log_S - log1p(beta * ssb_prev[p]));
        }
      } else {
        output[1, i, p] = 0.0;
      }
    }
    
    // 2. Survival and Movement (transition from i-1 to i)
    for (a in 2 : n_ages) {
      row_vector[n_patches] surv = exp(neg_mort[i - 1] - f_a_t[a - 1, i - 1]);
      row_vector[n_patches] lambda_surv = output[a - 1, i - 1] .* surv;
      
      if (mov_age[a]) {
        vector[n_patches] adj_x =
          csr_matrix_times_vector(n_patches, n_patches, w_adj,
                                  v_adj, u_adj, lambda_surv');
        output[a, i] = (zeta * lambda_surv' + (1 - zeta) * adj_x)';
      } else {
        output[a, i] = lambda_surv;
      }
    }
  }
  
  return output;
}
