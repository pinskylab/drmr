##' @title DRM fitted values
##'
##' @description Calculates the "fitted values" associated to an \code{adrm}
##'   object. These fitted values are useful for "hindcasts" and posterior
##'   predictive checks.
##'
##' @param object \code{adrm} object containing the output from the [fit_drm]
##'   function.
##' @param type type of predictions to be computed. Admitted values are
##' \itemize{
##' \item `"predictive"` (default): posterior predictive distribution;
##' \item `"expected"`: theoretical mean of the posterior predictive distribution;
##' \item `"latent"`: latent density (i.e., disconsidering the observation error);
##' }
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
                        type = "predictive",
                        cores = 1,
                        ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(type %in% c("predictive", "expected", "latent"))
  type_pred <- switch(type,
                      predictive = 0,
                      expected   = 1,
                      latent     = 2)
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(object$data)
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ## load compiled model
  ft_out <-
    instantiate::stan_package_model(name = "fitted_drm",
                                    package = "drmr")
  ## computing fitted values
  gq <- ft_out$
    generate_quantities(fitted_params = fitted_params,
                        data = c(object$data, list(type = type_pred)),
                        parallel_chains = cores,
                        ...)
  spt <- data.frame(v1 = object$cols$site_levels[object$data$site],
                    v2 = object$data$time + object$data$time_init - 1)
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
}

##' @title SDM fitted values
##'
##' @description Calculates the log-likelihood associated to a \code{sdm}
##'   obkect.
##'
##' @param object A \code{sdm} object containing the output from the [fit_sdm]
##'   function.
##' @param type type of predictions to be computed. Admitted values are
##' \itemize{
##' \item `"predictive"` (default): posterior predictive distribution;
##' \item `"expected"`: theoretical mean of the posterior predictive distribution;
##' \item `"latent"`: latent density (i.e., disconsidering the observation error);
##' }
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
                       type = "predictive",
                       cores = 1,
                       ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(type %in% c("predictive", "expected", "latent"))
  type_pred <- switch(type,
                      predictive = 0,
                      expected   = 1,
                      latent     = 2)
  ##--- pars from model fitted ----
  pars <- fitted_pars_ll(object$data)
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ## load compiled lambda
  ft_out <-
    instantiate::stan_package_model(name = "fitted_sdm",
                                    package = "drmr")
  ## computing forecast
  gq <- ft_out$
    generate_quantities(fitted_params = fitted_params,
                        data = c(object$data, list(type = type_pred)),
                        parallel_chains = cores,
                        ...)
  spt <- data.frame(v1 = object$cols$site_levels[object$data$site],
                    v2 = object$data$time + object$data$time_init - 1)
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
}
