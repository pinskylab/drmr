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

##' @title DRM log-likelihood (for model comparisons' purpose)
##'
##' @description Calculates the log-likelihood associated to a given \code{drm}
##'   model.
##'
##' @param drm A \code{list} object containing the output from the [fit_drm]
##'   function.
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{drm}, then four (or less) threads are recommended.
##'
##' @author lcgodoy
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for forecasting.
##'
##' @export
loglik_drm <- function(drm,
                       cores = 1) {
  stopifnot(inherits(drm$stanfit, "CmdStanFit"))
  ## number time points for forecasting
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(drm$data)
  fitted_params <-
    drm$stanfit$draws(variables = pars)
  ## load compiled lambda
  ll_out <-
    instantiate::stan_package_model(name = "loglik_drm",
                                    package = "drmr")
  ## computing forecast
  output <- ll_out$
    generate_quantities(fitted_params = fitted_params,
                        data = drm$data,
                        parallel_chains = cores)
  return(output)
}
