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

// * selectivity at length *//
vector sel_length(vector lengths, real par1,
                  real length_50_sel, real sel_delta) {
  int n = num_elements(lengths);
  vector[n] out;
  out = par1 * ((lengths - length_50_sel) / sel_delta);
  out = inv_logit(out);
  return out;
}

//* Von Bertalanffy curve *//
vector vb_curve(vector age, real l_inf,
                real k, real t0) {
  int n = num_elements(age);
  vector[n] out;
  out = - l_inf * expm1(-k * (age - t0));
  return out;
}

//* selectivity at age *//
vector sel_age(vector age, real par1,
               real length_50_sel,
               real sel_delta, real l_inf,
               real k, real t0) {
  int n = num_elements(age);
  vector[n] lengths;
  lengths = vb_curve(age, l_inf, k, t0);
  vector[n] out;
  out = par1 * ((lengths - length_50_sel) / sel_delta);
  out = inv_logit(out);
  return out;
}

// finally: function to calculate range quantiles
// real calculate_range_quantile(int np, vector patches,
//                               array[] real dens_by_patch,
//                               real quantile_out) {
//   vector[np] csum_dens;
//   vector[np] csum_to_edge;
//   real quant_position;
//   real cutoff;
//   int cutoff_id;
//   cutoff = sum(dens_by_patch[1 : np]) * quantile_out; // where along the
//   // vector of counts does
//   // the quantile cutoff
//   // fall? (then need to
//   // translate that into a
//   // patch position,
//   // below)
//   for (i in 1 : np) {
//     csum_dens[i] = csum(dens_by_patch[1 : np], i); // calculate cumulative sum
//     // of density along each
//     // patch
//     if (csum_dens[i] <= cutoff) {
//       csum_to_edge[i] = csum_dens[i]; // keep only the patches below the edge
//     } else {
//       csum_to_edge[i] = 0;
//     }
//   }
//   cutoff_id = min(which_equal(csum_to_edge, 0)); // get lowest patch with a 0
//   // (that's the patch where
//   // the edge falls)
//   quant_position = ((cutoff - max(csum_to_edge)) / dens_by_patch[cutoff_id])
//     + cutoff_id - 1; // calculate what proportion of the edge-containing patch
//   // is "filled in" by the actual fish up to the weighted
//   // quantile (that's the decimal) and then add in the other
//   // patches
//   return quant_position;
// } // close quantile function

  // selectivity at length
  // par1 = log(19)
