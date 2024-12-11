/**
 * Penalized complexity prior for the AR(1) correlation. This prior assumes:
 * P(|x| > u) = alpha A reasonable strategy here is to set u as a possibly high
 * autocorrelation and alpha as a small probability.
 *
 * Reference: https://doi.org/10.1111/jtsa.12242
 * 
 * @param x autocorrelation parameter
 * @param u upper bound for |x|
 * @param alpha prior probability of exceeding the upper bound u.
 * 
 * @return an array of numbers by age, year and patch
 */
real pcp_ar0_lpdf(real x, real alpha, real u) {
  real aux = - log(1 - square(x));
  real theta = - log(alpha) / sqrt(- log(1 - square(u)));
  return log(theta) - log(2) - theta * sqrt(aux) +
    log(abs(x)) + aux - 0.5 * log(aux);
}

/**
 * Penalized complexity prior for the log-transformed standard deviation of a
 * random effect. It assumes: P(exp(x) > u) = alpha.
 * 
 * @param x log-sd parameter
 * @param u upper bound for conditional SD
 * @param alpha prior probability of exceeding the upper bound u.
 * 
 * @return an array of numbers by age, year and patch
 */
real pcp_logsd_lpdf(real x, real alpha, real u) {
  real lambda = - log(alpha) / u;
  return log(lambda) + x - lambda * exp(x);
}
