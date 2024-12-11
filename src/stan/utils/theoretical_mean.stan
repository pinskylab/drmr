/**
 * @title Generate theoretical mean according to the Ricker model
 * @return an array of numbers by age, year and patch
 */
// array[] matrix ricker_model() {
// }

/**
 * @title Generate theoretical mean according to the BH model
 *
 * @description
 * @param r0 recruitment if population were not fished
 * @param maturity_at_age ?? (used for SSB)
 * @param wt_at_age ?? (used for SSB)
 * @param h steepness parameter of the Beverton-Holt model; it controls how
 *    sensitive recruitment was to spawning stock biomass
 * @param ssb0 theoretical max of spawning stock biomass
 
 * @return an array of numbers by age, year and patch
 */
// array[] matrix bh_model(int n_patches,
//                         int n_time,
//                         int n_ages,
//                         // Mortality parameter
//                         matrix f_a_y,
//                         real mort,
//                         // Initializing N
//                         int use_init_n_at_age,
//                         matrix init_n_at_age,
//                         matrix log_rec,
//                         real r0,
//                         vector maturity_at_age,
//                         vector wt_at_age,
//                         real h,
//                         real ssb0) {
// }

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
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
  output[1] += init;
  for (i in 2 : n_time) {
    for (p in 1 : n_patches) {
      for (a in 2 : n_ages) {
        output[a, i, p] = output[a - 1, i - 1, p]
          * exp(neg_mort[i - 1, p] - f_a_t[a - 1, i - 1]);
      }
    }
  }
  return output;
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
                                 array[] matrix lambda_past,
                                 matrix f_past,
                                 vector neg_mort_past) {
  // initializing output with zeros
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
  int past_last_time;
  past_last_time = rows(lambda_past[1]);
  output[1] += init;
  for (i in 1 : n_time) {
    for (p in 1 : n_patches) {
      for (a in 2 : n_ages) {
        if (i == 1) {
          output[a, i, p] = lambda_past[a - 1, past_last_time, p]
            * exp(neg_mort_past[p] - f_past[a - 1, past_last_time]);
        } else {
          output[a, i, p] = output[a - 1, i - 1, p]
            * exp(neg_mort[i - 1, p] - f_a_t[a - 1, i - 1]);
        }
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
 * @param mov_age age at which movement starts (this can be generalized)
 * 
 * @return an array of numbers by age, year and patch
 */
array[] matrix apply_movement(array[] matrix lambda, matrix M,
                              int mov_age) {
  array[3] int dimensions;
  dimensions = dims(lambda);
  array[dimensions[1]] matrix[dimensions[2], dimensions[3]] output;
  output = lambda;
  for (a in mov_age:dimensions[1]) {
    for (time in 1:dimensions[2]) {
      output[a, time, 1:dimensions[3]] =
        lambda[a, time, 1:dimensions[3]] * M';
    }
  }
  return output;
}
