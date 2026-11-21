##' @title Computing the log-likelihood function fot `adrm` and `sdm` objets.
##' @inheritParams plot.drmrmodels x
##' @inheritParams preddrm
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
  check_gq(output)
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
  check_gq(output)
  return(output)
}
