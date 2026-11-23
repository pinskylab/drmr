##' @title Random number generation from a truncated Student's t distribution
##' @description Generates random numbers from a truncated Student's t
##'   distribution.
##' @details For details on the method used, see:
##'   [https://stats.stackexchange.com/questions/567944/how-can-i-sample-from-a-shifted-and-scaled-student-t-distribution-with-a-specifi](https://stats.stackexchange.com/questions/567944/how-can-i-sample-from-a-shifted-and-scaled-student-t-distribution-with-a-specifi)
##' @param n An \code{integer} specifying the number of random numbers to
##'   generate.
##' @param mean A \code{numeric} scalar representing the mean of the underlying t
##'   distribution. Defaults to 0.
##' @param sd A \code{numeric} scalar representing the standard deviation of the
##'   underlying t distribution. Defaults to 1.
##' @param df A \code{numeric} scalar representing the degrees of freedom of the
##'   underlying t distribution. Defaults to 3.
##' @param range A \code{numeric vector} of length 2 specifying the lower and
##'   upper truncation bounds. Defaults to `c(0, Inf)`, indicating truncation
##'   from below at 0.
##' @return A \code{numeric vector} of length \code{n} containing random numbers
##'   drawn from the specified truncated t distribution.
##' @seealso [stats::rt()], [stats::qt()]
##' @author Anonymous
rtt <- function(n, mean = 0, sd = 1,
                df = 3, range = c(0, Inf)) {
  ll <- stats::pt(q = min(range - mean) / sd, df = df)
  hh <- stats::pt(q = max(range - mean) / sd, df = df)
  u  <- stats::runif(n, min = ll, max = hh)
  mean + stats::qt(p = u, df = df) * sd
}

##' @title Random number generation from a truncated Normal distribution
##' @description Generates random numbers from a truncated Normal distribution.
##' @inheritParams rtt n range
##' @param mean A \code{numeric} scalar representing the mean of the underlying
##'   normal distribution. Defaults to 0.
##' @param sd A \code{numeric} scalar representing the standard deviation of the
##'   underlying normal distribution. Defaults to 1.
##' @return A \code{numeric vector} of length \code{n} containing random numbers
##'   drawn from the specified truncated normal distribution.
##' @seealso [stats::rnorm()], [stats::qnorm()]
##' @author Anonymous
rtn <- function(n, mean = 0, sd = 1, range = c(0, Inf)) {
  lo <- min(range)
  hi <- max(range)
  if (lo > mean) {
    ## the interval lies above the mean: upper-tail probabilities do not
    ## round to 1, so the draws stay finite
    pl <- stats::pnorm(q = lo, mean = mean, sd = sd, lower.tail = FALSE)
    ph <- stats::pnorm(q = hi, mean = mean, sd = sd, lower.tail = FALSE)
    u  <- stats::runif(n, min = ph, max = pl)
    return(stats::qnorm(p = u, mean = mean, sd = sd, lower.tail = FALSE))
  }
  ll <- stats::pnorm(q = lo, mean = mean, sd = sd)
  hh <- stats::pnorm(q = hi, mean = mean, sd = sd)
  u  <- stats::runif(n, min = ll, max = hh)
  stats::qnorm(p = u, mean = mean, sd = sd)
}
