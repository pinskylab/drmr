##' @inherit get_fitted_pars
fitted_pars_drm <- function(data_list) {
  output <- c("beta_t", "beta_r",
              "lambda")
  if (data_list$rho_mu == 1)
    output <- c(output, "xi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$movement == 1)
    output <- c(output, "zeta")
  if (data_list$est_surv == 1)
    output <- c(output, "beta_s")
  if (data_list$ar_re > 0)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (data_list$iid_re > 0)
    output <- c(output, "z_i", "sigma_i")
  if (data_list$sp_re > 0)
    output <- c(output, "z_s")
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_sdm <- function(data_list) {
  output <- c("beta_t", "beta_r")
  if (data_list$rho_mu == 1)
    output <- c(output, "xi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$ar_re == 1)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (data_list$iid_re == 1)
    output <- c(output, "z_i")
  if (data_list$sp_re == 1)
    output <- c(output, "z_s")
  return(output)
}

##' @title Retrieve parameters needed for predictions
##'
##' @description This function identifies the parameters necessary for carrying
##'   out predictions based on the data used to fit a DRM (or SDM).
##'
##' @param data_list a \code{list} used as input for model fitting. Typically,
##'   the output from the [make_data] function.
##' @param model a \code{character} indicating which model predictions are
##'   sought for. This input admits two possible entries: "drm" (default) or
##'   "sdm".
##' @return a \code{character} vector of labels indicating the parameters
##'   necessary for the predictions.
##' @author lcgodoy
get_fitted_pars <- function(data_list, model = "drm") {
  stopifnot(inherits(data_list, "list"))
  stopifnot(length(model) == 1)
  if (model == "drm") {
    output <- fitted_pars_drm(data_list)
  } else
    output <- fitted_pars_sdm(data_list)
  return(output)
}

##' @title Predictions based on DRM.
##'
##' @description Considering a new dataset (across the same patches), computes
##'   predictions based on the DRM passed as \code{drm}.
##'
##' @param object A \code{adrm} (or \code{sdm}) object containing the output from the
##'   [fit_drm()] (or [fit_sdm()]) function.
##' @param new_data a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions.
##' @param past_data a \code{data.frame} with the dataset last year used in
##'   model fitting. Only needed when \code{f_test} is not missing or when
##'   estimating survival.
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (columns) and timepoint (rows).
##' @param type type of predictions to be computed. Admitted values are
##' \itemize{
##' \item `"predictive"` (default): posterior predictive distribution;
##' \item `"expected"`: theoretical mean of the posterior predictive distribution;
##' \item `"latent"`: latent density (i.e., disconsidering the observation error);
##' }
##' @param seed a seed used for the predictions. predictiosn are obtained
##'   through Monte Carlo samples from the posterior predictive
##'   distribution. Therefore, a \code{seed} is needed to ensure the results'
##'   reproducibility.
##' @param cores number of threads used to compute the predictions. If four
##'   chains were used in the \code{drm}, then four (or less) threads are
##'   recommended.
##' @param ... additional parameters to be passed to
##'   \code{$generated_quantities}
##'
##'
##' @details The current version of the code assumes the data where predictions
##'   are needed is ordered by "patch" and "site" and, in addition, its patches
##'   MUST be the same as the ones used to obtain the parameters' estimates from
##'   the the \code{drm} object.
##'
##' @author lcgodoy
##'
##' @name preddrm
##' 
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for predictions.
##'
##' @export
predict.adrm <- function(object,
                         new_data,
                         past_data,
                         f_test,
                         type = "predictive",
                         seed = 1,
                         cores = 1,
                         ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(type %in% c("predictive", "expected", "latent"))
  type_prep <- switch(type,
                      predictive = 0,
                      expected   = 1,
                      latent     = 2)
  ## number time points for predictions
  ntime_for <- length(unique(new_data[[object$cols$time_col]]))
  time_for <- new_data[[object$cols$time_col]] -
    min(new_data[[object$cols$time_col]]) + 1
  if (!missing(f_test)) {
    stopifnot(NROW(f_test) == object$data[["n_ages"]])
  }
  x_tt <- stats::model.matrix(object[["formulas"]][["formula_zero"]],
                              data = new_data)
  x_rt <- stats::model.matrix(object[["formulas"]][["formula_rec"]],
                              data = new_data)
  if (missing(f_test)) {
    f_test <- matrix(0, ncol = ntime_for, nrow = object$data$n_ages)
  }
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(object$data, "drm")
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ##--- list for predictions gq ----
  pred_data <-
    list(n_patches = object$data$n_patches,
         n_ages = object$data$n_ages,
         n_time = ntime_for,
         n_time_train = object$data$n_time,
         time = time_for,
         patch = new_data[[object$cols$site_col]],
         rho_mu = object$data$rho_mu,
         ar_re = object$data$ar_re,
         iid_re = object$data$iid_re,
         sp_re = object$data$sp_re,
         movement = object$data$movement,
         est_surv = object$data$est_surv,
         cloglog = object$data$cloglog,
         likelihood = object$data$likelihood,
         type = type_prep,
         K_t = object$data$K_t,
         X_t = x_tt,
         f = f_test,
         f_past = object$data$f,
         m = object$data$m,
         adj_mat = object$data$adj_mat,
         ages_movement = object$data$ages_movement,
         selectivity_at_age = object$data$selectivity_at_age,
         K_r = object$data$K_r,
         X_r = x_rt)
  pred_data$N <- pred_data$n_patches * pred_data$n_time
  if (length(object$data$K_m) > 0) {
    stopifnot(!missing(past_data))
    x_mpast <-
      stats::model.matrix(object[["formulas"]][["formula_surv"]],
                          data = past_data)
    x_m <- stats::model.matrix(object[["formulas"]][["formula_surv"]],
                               data = new_data)
    pred_data$K_m <- array(NCOL(x_m), dim = 1)
    pred_data$X_m <- x_m
    pred_data$X_m_past <- x_mpast
  } else {
    pred_data$K_m <- integer(0)
    pred_data$X_m <- matrix(0)
    pred_data$X_m_past <- matrix(0)
  }
  ## load compiled "predict"
  pred_comp <-
    instantiate::stan_package_model(name = "predict",
                                    package = "drmr")
  ## computing predictions
  output <- pred_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = pred_data,
                        seed = seed,
                        parallel_chains = cores,
                        ...)
  return(output)
}

##' @title Predictions based on SDM.
##'
##' @description Considering a new dataset (across the same patches), computes
##'   predictions based on the SDM passed as \code{sdm}.
##'
##' @param object A \code{list} object containing the output of a [fit_sdm] call.
##'
##' @param new_data a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions.
##' @param seed a seed used for the predictions. predictions are obtained
##'   through Monte Carlo samples from the posterior predictive
##'   distribution. Therefore, a \code{seed} is needed to ensure the results'
##'   reproducibility.
##' @param type type of predictions to be computed. Admitted values are
##' \itemize{
##' \item `"predictive"` (default): posterior predictive distribution;
##' \item `"expected"`: theoretical mean of the posterior predictive distribution;
##' \item `"latent"`: latent density (i.e., disconsidering the observation error);
##' }
##' @param cores number of threads used to compute the predictions. If four
##'   chains were used in the \code{drm}, then four (or less) threads are
##'   recommended.
##' @param ... additional parameters to be passed to
##'   \code{$generated_quantities}
##'
##' @details The current version of the code assumes the data where predictions
##'   are needed is ordered by "patch" and "site" and, in addition, its patches
##'   MUST be the same as the ones used to obtain the parameters' estimates from
##'   the the \code{sdm} object.
##'
##' @author lcgodoy
##'
##' @return An object of class \code{CmdStanGQ} (from the \code{instantiate}
##'   package) containing samples for the posterior predictive distribution for
##'   predictions.
##'
##' @rdname predsdm
predict.sdm <- function(object,
                        new_data,
                        type = "predictive",
                        seed = 1,
                        cores = 1,
                        ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                  "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(type %in% c("predictive", "expected", "latent"))
  type_prep <- switch(type,
                      predictive = 0,
                      expected   = 1,
                      latent     = 2)
  ## time points for predictions
  ntime_for <- length(unique(new_data[[object$cols$time_col]]))
  time_for <- new_data[[object$cols$time_col]] -
    min(new_data[[object$cols$time_col]]) + 1
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(object$data, "sdm")
  fitted_params <-
    object$stanfit$draws(variables = pars)
  pred_data <-
    list(n_patches = object$data$n_patches,
         n_time = ntime_for,
         n_time_train = object$data$n_time,
         time = time_for,
         patch = new_data[[object$cols$site_col]],
         rho_mu = object$data$rho_mu,
         ar_re = object$data$ar_re,
         iid_re = object$data$iid_re,
         sp_re = object$data$sp_re,
         cloglog = object$data$cloglog,
         likelihood = object$data$likelihood,
         type = type_prep,
         Z = stats::model.matrix(object[["formulas"]][["formula_zero"]],
                                 data = new_data),
         X = stats::model.matrix(object[["formulas"]][["formula_dens"]],
                                 data = new_data))
  pred_data$N <- pred_data$n_patches * pred_data$n_time
  pred_data$K_z <- NCOL(pred_data$Z)
  pred_data$K_x <- NCOL(pred_data$X)
  ## load compiled predict_sdm
  pred_comp <-
    instantiate::stan_package_model(name = "predict_sdm",
                                    package = "drmr")
  ## computing predictions
  output <- pred_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = pred_data,
                        seed = seed,
                        parallel_chains = cores,
                        ...)
  return(output)
}

##' @rdname preddrm
##' @param ... params to be passed to \code{predict.adrm}
##' @export
predict_drm <- function(...) {
  lifecycle::deprecate_warn("0.3.1", "predict_drm()", "predict()")
  predict.adrm(...)
}

##' @rdname predsdm
##' @param ... params to be passed to \code{predict.sdm}
##' @export
predict_sdm <- function(...) {
  lifecycle::deprecate_warn("0.3.1", "predict_sdm()", "predict()")
  predict.sdm(...)
}
