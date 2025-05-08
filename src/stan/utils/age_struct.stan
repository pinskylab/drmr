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
                        matrix recruitment) {
  // initializing output with zeros
  array[n_ages] matrix[n_time, n_patches] output
    = rep_array(rep_matrix(0.0, n_time, n_patches), n_ages);
  for (p in 1:n_patches) {
    for (i in 1:n_time) {
      output[1, i, p] = recruitment[i, p];
    }
  }
  for (a in 1 : (n_ages - 1)) {
    output[a, 1 : a, ] = rep_matrix(to_vector(init[1 : a]), n_patches);
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
