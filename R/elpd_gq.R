##' @title (out-of-sample) Expected Log-posterior Density (ELPD) based on
##'   \code{adrm} and \code{sdm} objects.
##' @description Considering a new dataset (across the same sites), computes
##'   the out-of-sample ELPD based on the DRM passed as \code{drm}.
##'
##' @param x A \code{list} object containing the output from the [fit_drm()] (or
##'   [fit_sdm()]) function.
##' @param new_data a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions. Note that, this \code{data.frame} must contain the
##'   response variable used when fitting the DRM as well.
##' @param past_data a \code{data.frame} with the dataset last year used in
##'   model fitting. Only needed when \code{f_test} is not missing or when
##'   estimating survival.
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (columns) and timepoint (rows).
##' @param seed a seed used for the forecasts. Forecasts are obtained through
##'   Monte Carlo samples from the posterior predictive distribution. Therefore,
##'   a \code{seed} is needed to ensure the results' reproducibility.
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{drm}, then four (or less) threads are recommended.
##' @param ... additional parameters to be passed to \code{elpd}
##' 
##' @details The current version of the code assumes the data where forecasts
##'   are needed are ordered by "site" and "site" and, in addition, the sites
##'   MUST be the same as the ones used to obtain the parameters' estimates from
##'   the the \code{drm} object.
##'
##' @author lcgodoy

##' @return an object of class \code{"CmdStanGQ"} containing the ELPD
##'   function evaluated at each data point given each sample from the
##'   posterior.
##' @name elpd
##' @author lcgodoy
##' @export
elpd <- function(x, ...) UseMethod("elpd", x)

##' @rdname elpd
##' @export
elpd.adrm <- function(x,
                      new_data,
                      past_data,
                      f_test,
                      seed = 1,
                      cores = 1,
                      ...) {
  stopifnot(inherits(x$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                  "CmdStanPathfinder", "CmdStanVB")))
  ## number time points for forecasting
  ntime_for <- length(unique(new_data[[x$cols$time_col]]))
  time_for <- new_data[[x$cols$time_col]] -
    min(new_data[[x$cols$time_col]]) + 1
  if (!missing(f_test)) {
    stopifnot(NROW(f_test) == x$data[["n_ages"]])
  }
  x_tt <- stats::model.matrix(x[["formulas"]][["formula_zero"]],
                              data = new_data)
  x_rt <- stats::model.matrix(x[["formulas"]][["formula_rec"]],
                              data = new_data)
  if (missing(f_test)) {
    f_test <- matrix(0, ncol = ntime_for, nrow = x$data$n_ages)
  }
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(x$data, "drm")
  fitted_params <-
    x$stanfit$draws(variables = pars)
  ##--- list for forecast object ----
  forecast_data <-
    list(n_sites = x$data$n_sites,
         n_ages = x$data$n_ages,
         n_time = ntime_for,
         n_time_train = x$data$n_time,
         time = time_for,
         site = new_data[[x$cols$site_col]],
         rho_mu = x$data$rho_mu,
         ar_re = x$data$ar_re,
         iid_re = x$data$iid_re,
         sp_re = x$data$sp_re,
         movement = x$data$movement,
         est_surv = x$data$est_surv,
         cloglog = x$data$cloglog,
         likelihood = x$data$likelihood,
         new_y = new_data[[x$cols$y_col]],
         K_t = x$data$K_t,
         X_t = x_tt,
         f = f_test,
         f_past = x$data$f,
         m = x$data$m,
         adj_mat = x$data$adj_mat,
         ages_movement = x$data$ages_movement,
         selectivity_at_age = x$data$selectivity_at_age,
         K_r = x$data$K_r,
         X_r = x_rt)
  forecast_data$N <- forecast_data$n_sites * forecast_data$n_time
  if (length(x$data$K_m) > 0) {
    stopifnot(!missing(past_data))
    x_mpast <-
      stats::model.matrix(x[["formulas"]][["formula_surv"]],
                          data = past_data)
    x_m <- stats::model.matrix(x[["formulas"]][["formula_surv"]],
                               data = new_data)
    forecast_data$K_m <- array(NCOL(x_m), dim = 1)
    forecast_data$X_m <- x_m
    forecast_data$X_m_past <- x_mpast
  } else {
    forecast_data$K_m <- integer(0)
    forecast_data$X_m <- matrix(0)
    forecast_data$X_m_past <- matrix(0)
  }
  ## load compiled forecast
  elpd_comp <-
    instantiate::stan_package_model(name = "elpd_drm",
                                    package = "drmr")
  ## computing forecast
  output <- elpd_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = forecast_data,
                        seed = seed,
                        parallel_chains = cores)
  return(output)
}

##' @rdname elpd
##' @export
elpd.sdm <- function(x,
                     new_data,
                     seed = 1,
                     cores = 1,
                      ...) {
  stopifnot(inherits(x$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                  "CmdStanPathfinder", "CmdStanVB")))
  ## time points for forecasting
  ntime_for <- length(unique(new_data[[x$cols$time_col]]))
  time_for <- new_data[[x$cols$time_col]] -
    min(new_data[[x$cols$time_col]]) + 1
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(x$data, "sdm")
  fitted_params <-
    x$stanfit$draws(variables = pars)
  ##--- list for forecast object ----
  forecast_data <-
    list(n_sites = x$data$n_sites,
         n_time = ntime_for,
         n_time_train = x$data$n_time,
         time = time_for,
         site = new_data[[x$cols$site_col]],
         rho_mu = x$data$rho_mu,
         ar_re = x$data$ar_re,
         iid_re = x$data$iid_re,
         sp_re = x$data$sp_re,
         cloglog = x$data$cloglog,
         likelihood = x$data$likelihood,
         new_y = new_data[[x$cols$y_col]],
         Z = stats::model.matrix(x[["formulas"]][["formula_zero"]],
                                 data = new_data),
         X = stats::model.matrix(x[["formulas"]][["formula_dens"]],
                                 data = new_data))
  forecast_data$N <- forecast_data$n_sites * forecast_data$n_time
  forecast_data$K_z <- NCOL(forecast_data$Z)
  forecast_data$K_x <- NCOL(forecast_data$X)
  ## load compiled forecast
  elpd_comp <-
    instantiate::stan_package_model(name = "elpd_sdm",
                                    package = "drmr")
  ## computing forecast
  output <- elpd_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = forecast_data,
                        seed = seed,
                        parallel_chains = cores)
  return(output)
}
