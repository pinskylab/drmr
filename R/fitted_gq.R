##' @title DRM fitted values
##'
##' @description Calculates the "fitted values" associated to an \code{adrm}
##'   object. These fitted values are useful for "hindcasts" and posterior
##'   predictive checks.
##'
##' @param object \code{adrm} object containing the output from the [fit_drm]
##'   function.
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{fit_drm}, then four (or less) threads are recommended.
##' @param ... additional parameters to be passed to \code{$generated_quantities}
##'
##' @author lcgodoy
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for the observed data.
##'
##' @export
fitted.adrm <- function(object,
                        cores = 1,
                        ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(object$data)
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ## load compiled model
  ft_out <-
    instantiate::stan_package_model(name = "fitted_drm",
                                    package = "drmr")
  ## computing fitted values
  output <- ft_out$
    generate_quantities(fitted_params = fitted_params,
                        data = object$data,
                        parallel_chains = cores,
                        ...)
  return(output)
}

##' @title SDM fitted values
##'
##' @description Calculates the log-likelihood associated to a \code{sdm}
##'   obkect.
##'
##' @param object A \code{sdm} object containing the output from the [fit_sdm]
##'   function.
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{sdm}, then four (or less) threads are recommended.
##' @param ... additional parameters to be passed to \code{$generated_quantities}
##'
##' @author lcgodoy
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for the observed data.
##'
##' @export
fitted.sdm <- function(object,
                       cores = 1,
                       ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(object$data)
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ## load compiled lambda
  ft_out <-
    instantiate::stan_package_model(name = "fitted_sdm",
                                    package = "drmr")
  ## computing forecast
  output <- ft_out$
    generate_quantities(fitted_params = fitted_params,
                        data = object$data,
                        parallel_chains = cores,
                        ...)
  return(output)
}
