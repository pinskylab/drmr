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

real ziloglik_lpdf(vector y,
                   data int likelihood,
                   data int N_nz,
                   data int N,
                   data array[] int id_nz,
                   data array[] int id_z,
                   vector mu,
                   vector rho,
                   real phi) {
  real out = 0.0;
  out += sum(log(rho[id_z]));
  out += sum(log1m(rho[id_nz]));
  if (likelihood == 0) {
    vector[N_nz] loc_par;
    loc_par = log(mu[id_nz]) + square(phi) / 2;
    out += lognormal_lpdf(y[id_nz] | loc_par, phi);
  } else if (likelihood == 1) {
    vector[N_nz] mu_ln;
    vector[N_nz] sigma_ln;
    sigma_ln = sqrt(log1p(phi * inv_square(mu[id_nz])));
    mu_ln = log(square(mu[id_nz]) .*
                inv_sqrt(square(mu[id_nz]) + phi));
    out += lognormal_lpdf(y[id_nz] | mu_ln, sigma_ln);
  } else if (likelihood == 2) {
    vector[N_nz] b_g;
    b_g = phi / mu[id_nz];
    out += gamma_lpdf(y[id_nz] | phi, b_g);
  } else if (likelihood == 3) {
    vector[N_nz] a_ll;
    real b_ll;
    b_ll = phi + 1;
    a_ll = sin(pi() / b_ll) * mu[id_nz] * b_ll / inv(pi());
    out += loglogistic_lpdf(y[id_nz] | a_ll, b_ll);
  } else {
    out += normal_lpdf(y | mu, phi) -
      normal_lccdf(rep_vector(0.0, N) | mu, phi);
  }
  return out;
}

real ptziloglik_lpdf(real y,
                     data int likelihood,
                     data int is_zero,
                     real mu,
                     real rho,
                     real phi) {
  real out = 0.0;
  if (is_zero) {
    out += log(rho);
  } else {
    out += log1m(rho);
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu) + square(phi) / 2;
      out += lognormal_lpdf(y | loc_par, phi);
    } else if (likelihood == 1) {
      out += ln_mu_lpdf(y | mu, phi);
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi / mu;
      out += gamma_lpdf(y | phi, gamma_beta);
    } else if (likelihood == 3) {
      real a_ll;
      real b_ll;
      b_ll = phi + 1;
      a_ll = sin(pi() / b_ll) * mu * b_ll * inv(pi());
      out += loglogistic_lpdf(y | a_ll, b_ll);
    } else {
      out += normal_lpdf(y | mu, phi) -
        normal_lccdf(0.0 | mu, phi);
    }
  }
  return out;
}

real drmsdm_rng(real mu,
                real rho,
                real phi,
                data int likelihood) {
  real out;
  out = (1 - bernoulli_rng(rho));
  if (out > 0) {
    if (likelihood == 0) {
      real loc_par;
      loc_par = log(mu) + square(phi) / 2;
      out *= lognormal_rng(loc_par, phi);
    } else if (likelihood == 1) {
      real mu_ln;
      real sigma_ln;
      sigma_ln = sqrt(log1p(phi * inv_square(mu)));
      mu_ln = log(square(mu) * inv_sqrt(square(mu) + phi));
      out *= lognormal_rng(mu_ln, sigma_ln);
    } else if (likelihood == 2) {
      real gamma_beta;
      gamma_beta = phi / mu;
      out *= gamma_rng(phi, gamma_beta);
    } else if (likelihood == 3) {
      real a_ll;
      real b_ll;
      b_ll = phi + 1;
      a_ll = sin(pi() / b_ll) * mu * b_ll * inv(pi());
      out *= loglogistic_rng(a_ll, b_ll);
    } else {
      array[2] real aux_tn = rep_array(0.0, 2);
      aux_tn[2] = normal_rng(mu, phi);
      out *= max(aux_tn);
    }
  }
  return out;
}
