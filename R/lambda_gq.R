##' @title Age-specific expected densities based on DRM.
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
##' @name age_dens
##' 
##' @export
ages_edens <- function(drm,
                       cores = 1) {
  stopifnot(inherits(drm$stanfit,
                     c("CmdStanFit", "CmdStanLaplace",
                       "CmdStanPathfinder", "CmdStanVB")))
  ## number time points for forecasting
  ##--- pars from model fitted ----
  pars <- fitted_pars_lambda(drm$data)
  fitted_params <-
    drm$stanfit$draws(variables = pars)
  ## load compiled lambda
  lambda_comp <-
    instantiate::stan_package_model(name = "lambda_gq",
                                    package = "drmr")
  output <- lambda_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = drm$data,
                        parallel_chains = cores)
  spt <- data.frame(v1 = drm$cols$site_levels[drm$data$site],
                    v2 = drm$data$time + drm$data$time_init - 1,
                    site_id = drm$data$site)
  colnames(spt)[1:2] <- rev(unname(unlist(drm$cols[2:3])))
  output <- list("gq" = output,
                 "spt" = spt)
  return(new_aesd(output))
}

##' @rdname age_dens
##' @param ... params to be passed to \code{age_edens}
##' @export
lambda_drm <- function(...) {
  lifecycle::deprecate_warn("0.3.1", "lambda_drm()", "ages_edens()")
  ages_edens(...)
}

##' Extract and Summarize Age-Specific Densities
##'
##' Parses the output of \code{lambda_drm} to return a data frame of
##' age-specific densities, indexed by age, time, and patch.
##'
##' @param lambda_obj A \code{CmdStanGQ} object returned by \code{lambda_drm}.
##' @param ages An optional vector of integers. If provided, the output is filtered
##'   to include only these specific ages.
##' @param probs A numeric vector of quantiles to calculate in the summary.
##'   Defaults to c(0.05, 0.5, 0.95).
##'
##' @return A \code{data.frame} containing the summary statistics and parsed indices
##'   (age, time, patch).
##' @export
summarise_adens <- function(lambda_obj, ages = NULL,
                            probs = c(0.05, 0.5, 0.95)) {
  summ <- lambda_obj$summary(NULL, "mean", "sd",
                             ~quantile(., probs = probs))
  summ <- as.data.frame(summ)
  summ <- summ[grep("^lambda\\[", summ$variable), ]
  vars <- summ$variable
  age  <- as.numeric(sub(".*\\[([0-9]+),.*", "\\1", vars))
  time <- as.integer(sub(".*,([0-9]+),.*", "\\1", vars))
  patch <- as.integer(sub(".*,([0-9]+)].*", "\\1", vars))
  out <- cbind(
    data.frame(age = age, time = time, patch = patch),
    summ
  )
  if (!is.null(ages)) {
    out <- out[out$age %in% ages, ]
  }
  rownames(out) <- NULL 
  return(out)
}
