##' @title Density of a truncated Student's t distribution
##' @description Calculates the density of a truncated Student's t distribution.
##' @param x A \code{numeric vector} at which to evaluate the density.
##' @param mean A \code{numeric} scalar representing the mean of the underlying t
##'   distribution. Defaults to 0.
##' @param sd A \code{numeric} scalar representing the standard deviation of the
##'   underlying t distribution. Defaults to 1.
##' @param df A \code{numeric} scalar representing the degrees of freedom of the
##'   underlying t distribution. Defaults to 3.
##' @param range A \code{numeric} vector of length 2 specifying the lower and
##'   upper truncation bounds. Defaults to \code{c(0, Inf)}, indicating
##'   truncation from below at 0.
##' @param log \code{Logical}; if \code{TRUE}, the natural logarithm of the
##'   density is returned.  Defaults to \code{FALSE}.
##' @return A numeric vector of the same length as \code{x}, containing the
##'   (log) density values of the truncated t distribution.
##' @seealso [stats::dt()], [stats::pt()]
##' @author lcgodoy
dtt <- function(x, mean = 0, sd = 1,
                df = 3, range = c(0, Inf),
                log = FALSE) {
  ll <- stats::pt(q = min(range - mean) / sd, df = df)
  hh <- stats::pt(q = max(range - mean) / sd, df = df)
  ## 1/sd is the Jacobian of the transformation.
  out <- stats::dt((x - mean) / sd, df = df, log = TRUE) -
    log(hh - ll) +
    log(sd)
  if (!log)
    out <- exp(out)
  return(out)
}

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
##' @author lcgodoy
rtt <- function(n, mean = 0, sd = 1,
                df = 3, range = c(0, Inf)) {
  ll <- stats::pt(q = min(range - mean) / sd, df = df)
  hh <- stats::pt(q = max(range - mean) / sd, df = df)
  u  <- stats::runif(n, min = ll, max = hh)
  mean + stats::qt(p = u, df = df) * sd
}

##' @title Density of a truncated Normal distribution
##' @description Calculates the density of a truncated Normal distribution.
##' @param x A \code{numeric vector} at which to evaluate the density.
##' @param mean A \code{numeric} scalar representing the mean of the underlying
##'   normal distribution. Defaults to 0.
##' @param sd A \code{numeric} scalar representing the standard deviation of the
##'   underlying normal distribution. Defaults to 1.
##' @param range A \code{numeric} vector of length 2 specifying the lower and
##'   upper truncation bounds. Defaults to `c(0, Inf)`, indicating truncation
##'   from below at 0.
##' @param log \code{Logical}; if \code{TRUE}, the natural logarithm of the
##'   density is returned. Defaults to \code{FALSE}.
##' @return A numeric vector of the same length as \code{x}, containing the
##'   (log) density values of the truncated normal distribution.
##' @seealso [stats::dnorm()], [stats::pnorm()]
##' @author lcgodoy
dtn <- function(x, mean = 0, sd = 1, range = c(0, Inf),
                log = FALSE) {
  ll <- stats::pnorm(q = min(range), mean = mean, sd = sd)
  hh <- stats::pnorm(q = max(range), mean = mean, sd = sd)
  out <- stats::dnorm(x, mean = mean, sd = sd,
                      log = TRUE) - log(hh - ll)
  if (!log)
    out <- exp(out)
  return(out)
}

##' @title Random number generation from a truncated Normal distribution
##' @description Generates random numbers from a truncated Normal distribution.
##' @param n An \code{integer} specifying the number of random samples to
##'   generate.
##' @param mean A \code{numeric} scalar representing the mean of the underlying
##'   normal distribution. Defaults to 0.
##' @param sd A \code{numeric} scalar representing the standard deviation of the
##'   underlying normal distribution. Defaults to 1.
##' @param range A \code{numeric vector} of length 2 specifying the lower and
##'   upper truncation bounds. Defaults to \code{c(0, Inf)}, indicating
##'   truncation from below at 0.
##' @return A \code{numeric vector} of length \code{n} containing random numbers
##'   drawn from the specified truncated normal distribution.
##' @seealso [stats::rnorm()], [stats::qnorm()]
##' @author lcgodoy
rtn <- function(n, mean = 0, sd = 1, range = c(0, Inf)) {
  ll <- stats::pnorm(q = min(range), mean = mean, sd = sd)
  hh <- stats::pnorm(q = max(range), mean = mean, sd = sd)
  u  <- stats::runif(n, min = ll, max = hh)
  stats::qnorm(p = u, mean = mean, sd = sd)
}
