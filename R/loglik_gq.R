##' @inherit get_fitted_pars
fitted_pars_ll <- function(data_list) {
  output <- c("mu", "beta_t")
  if (data_list$rho_mu == 1)
    output <- c(output, "lxi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  return(output)
}

##' @title Computing the log-likelihood function fot `adrm` and `sdm` objets.
##' @param x an object of class \code{adrm} or \code{sdm}. These objects are
##'   generated as the output of the [fit_drm()] and [fit_sdm()] functions,
##'   respectively,
##' @param cores number of threads used to calculate the log-likelihood
##'   function. If four chains were used in the \code{fit_drm} (of
##'   \code{fit_sdm}), then four (or less) threads are recommended.
##' @return an object of class \code{"CmdStanGQ"} containing the log-likelihood
##'   function evaluated at each data point given each sample from the
##'   posterior.
##' @name loglik
##' @author lcgodoy
##' @export
log_lik <- function(x, cores = 1) UseMethod("log_lik", x)

##' @rdname loglik
##' @export
log_lik.adrm <- function(x,
                         cores = 1) {
  stopifnot(inherits(x$stanfit,
                     c("CmdStanFit", "CmdStanLaplace",
                       "CmdStanPathfinder", "CmdStanVB")))
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(x$data)
  fitted_params <-
    x$stanfit$draws(variables = pars)
  ## load compiled model
  ll_out <-
    instantiate::stan_package_model(name = "loglik_drm",
                                    package = "drmr")
  ## computing loglik
  output <- ll_out$
    generate_quantities(fitted_params = fitted_params,
                        data = x$data,
                        parallel_chains = cores)
  return(output)
}

##' @rdname loglik
##' @export
log_lik.sdm <- function(x,
                        cores = 1) {
  stopifnot(inherits(x$stanfit,
                     c("CmdStanFit", "CmdStanLaplace",
                       "CmdStanPathfinder", "CmdStanVB")))
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(x$data)
  fitted_params <-
    x$stanfit$draws(variables = pars)
  ## load compiled model
  ll_out <-
    instantiate::stan_package_model(name = "loglik_sdm",
                                    package = "drmr")
  ## computing log_lik
  output <- ll_out$
    generate_quantities(fitted_params = fitted_params,
                        data = x$data,
                        parallel_chains = cores)
  return(output)
}
