##' @inherit get_fitted_pars
fitted_pars_drm <- function(data_list) {
  output <- c("coef_t", "coef_r",
              "lambda")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$movement == 1)
    output <- c(output, "mov_mat")
  if (data_list$est_mort == 1)
    output <- c(output, "coef_m")
  if (data_list$time_ar == 1)
    output <- c(output, "z_t", "alpha", "tau")
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_sdm <- function(data_list) {
  output <- c("coef_t", "coef_r")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$time_ar == 1)
    output <- c(output, "z_t", "alpha", "tau")
  return(output)
}

##' @title Retrieve parameters needed for forecasting
##'
##' @description This function identifies the parameters necessary for carrying
##'   out forecasts based on the data used to fit a DRM (or SDM).
##' 
##' @param data_list a \code{list} used as input for model fitting. Typically,
##'   the output from the [make_data] function.
##' @param model a \code{character} indicating which model forecasts are sought
##'   for. This input admits two possible entries: "drm" (default) or "sdm".
##' @return a \code{character} vector of labels indicating the parameters
##'   necessary for the forecast.
##' @author Lucas Godoy
get_fitted_pars <- function(data_list, model = "drm") {
  stopifnot(inherits(data_list, "list"))
  stopifnot(length(model) == 1)
  if (model == "drm") {
    output <- fitted_pars_drm(data_list)
  } else
    output <- fitted_pars_sdm(data_list)
  return(output)
}

##' @title Value of a covariate that maximizes the response variable in a
##'   quadratic model.
##' 
##' @description Consider a linear predictor having linear and square terms
##'   associated with a variable \eqn{x}. Assume this variable was centered
##'   before being included in the linear predictor. This functions returns the
##'   value of \eqn{x} (on its original scale) such that the linear predictor is
##'   maximized (or minimized).
##' 
##' @param drm A \code{CmdStanFit} object containing samples from the posterior
##'   distribution.
##' @param drm_data a \code{list} used as input for model fitting. Typically,
##'   the output from the [make_data] function.
##' @param ntime_for an \code{integer} denoting the number of timepoints for the
##'   forecast.
##' @param x_tt a design \code{matrix} of variables associated to the
##'   probability of absence at each site/time.
##' @param x_mt a design \code{matrix} of variables associated to survival.
##' @param x_rt a design \code{matrix} of variables associated to recruitment.
##' @param x_mpast a design \code{matrix} of variables associated to survival
##'   from the last year of the "training" perior (i.e., the period using to fit
##'   the model and obtain the parameters' estimates).
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (columns) and timepoint (rows).
##' @param seed a seed used for the forecasts. Forecasts are obtained through
##'   Monte Carlo samples from the posterior predictive distribution. Therefore,
##'   a \code{seed} is needed to ensure the results' reproducibility.
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{drm}, then four (or less) threads are recommended.
##' 
##' 
##' @details It is important that the rows of design matrices \code{x_tt},
##'   \code{x_rt} and \code{x_mt} are associated to the same patch/site and
##'   timepoints. In addition, the current version of the code assumes the data
##'   where forecasts are needed is ordered by "patch" and "site" and, in
##'   addition, its patches MUST be the same as the ones used to obtain the
##'   parameters' estimates from the the \code{drm} object.
##' 
##' @author Lucas Godoy
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for forecasting.
##' 
##' @export
predict_drm <- function(drm,
                        drm_data,
                        ntime_for,
                        x_tt,
                        x_rt,
                        x_mt = matrix(1, ncol = 1),
                        x_mpast = matrix(1, ncol = 1),
                        f_test,
                        seed = 1,
                        cores = 1) {
  ## time points for forecasting
  stopifnot(inherits(drm, "CmdStanFit"))
  stopifnot(NCOL(x_tt) == drm_data$K_t)
  stopifnot(NCOL(x_rt) == drm_data$K_r)
  if (length(drm_data$K_m) > 0) {
    stopifnot(!missing(x_mt))
    stopifnot(!missing(x_mpast))
    stopifnot(NCOL(x_mt) == drm_data$K_m)
    stopifnot(NCOL(x_mpast) == drm_data$K_m)
  }
  if (!missing(f_test)) {
    stopifnot(NROW(f_test) == drm_data$n_ages)
  }
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(drm_data, "drm")
  fitted_params <-
    drm$draws(variables = pars)
  ##--- list for forecast object ----
  forecast_data <-
    list(n_patches = drm_data$n_patches,
         n_ages = drm_data$n_ages,
         n_time = ntime_for,
         n_time_train = drm_data$n_time,
         time_ar = drm_data$time_ar,
         movement = drm_data$movement,
         est_mort = drm_data$est_mort,
         cloglog = drm_data$cloglog,
         likelihood = drm_data$likelihood,
         K_t = drm_data$K_t,
         X_t = x_tt,
         f = f_test,
         f_past = drm_data$f,
         m = drm_data$m,
         age_at_maturity = drm_data$age_at_maturity,
         selectivity_at_age = drm_data$selectivity_at_age,
         K_m = drm_data$K_m,
         X_m = x_mt,
         X_m_past = x_mpast,
         K_r = drm_data$K_r,
         X_r = x_rt)
  forecast_data$N <- forecast_data$n_patches * forecast_data$n_time
  ## load compiled forecast
  forecast_comp <-
    instantiate::stan_package_model(name = "forecast",
                                    package = "drmr")
  ## computing forecast
  output <- forecast_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = forecast_data,
                        seed = seed,
                        parallel_chains = cores)
  return(output)
}
