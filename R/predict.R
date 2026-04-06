##' @title Check if x is between lb and ub
##'
##' @param x A numeric vector.
##' @param lb A numeric vector of lower bounds.
##' @param ub A numeric vector of upper bounds.
##'
##' @return Returns \code{NULL} invisibly. This function stops execution if any of the
##'  following conditions are met:
##'  * lb is greater than or equal to ub.
##'  * The lengths of lb and ub are not equal.
##'  * The lengths of lb, ub, and x are not equal.
##'
check_between <- function(x, lb, ub) {
  stopifnot(all(lb <= ub))
  if (length(lb) > 1 | length(ub) > 1) {
    stopifnot(length(lb) == length(ub) &
              length(lb) == length(x))
  }
}

##' @inherit between
between_scalar <- function(x, lb, ub) {
  x >= lb & x <= ub
}

##' @title Check if elements in x are between corresponding elements in lb and
##'   ub
##'
##' @inheritParams check_between
##'
##' @return A logical vector of the same length as x, indicating whether each
##'   element of x is between the corresponding elements of lb and ub.
##'
##' @examples
##' between(1:5, 1, 5)
##' between(1:5, 2, 4)
##'
##' @export
between <- function(x, lb, ub) {
  check_between(x, lb, ub)
  if (length(lb) == 1) {
    sapply(x, between_scalar, lb, ub)
  } else {
    mapply(between_scalar, x, lb, ub)
  }
}

##' @title Calculate the interval score
##'
##' @description This function calculates the interval score for a given set of
##'   observations, lower and upper bounds, and alpha parameter.
##'
##' @details The interval score is a proper scoring rule that measures the
##'   accuracy of interval predictions. It takes into account both the coverage
##'   and the width of the prediction interval. A lower score indicates a better
##'   prediction.
##'
##' @param y A numeric vector of observations.
##' @param l A numeric vector of lower bounds for the prediction intervals.
##' @param u A numeric vector of upper bounds for the prediction intervals.
##' @param alpha A numeric value specifying the significance level (e.g., 0.05
##'   for a 95% interval).
##'
##' @return A numeric vector of interval scores.
##'
##' @export
int_score <- function(y, l, u, alpha) {
  stopifnot(alpha > 0 & alpha < 1)
  check_between(y, l, u)
  alpha <- 2 / alpha
  ind_l <- as.numeric(y < l)
  ind_u <- as.numeric(y > u)
  (u - l) + alpha * ((l - y) * ind_l + (y - u) * ind_u)
}

##' @title Predictions based on DRM.
##'
##' @description Considering a new dataset (across the same sites), computes
##'   predictions based on the DRM passed as \code{drm}.
##'
##' @param object An \code{adrm} object containing the output from the
##'   [fit_drm()] function.
##' @param new_data a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions.
##' @param past_data a \code{data.frame} with the dataset from the last year
##'   used in model fitting. Only needed when \code{f_test} is not missing or
##'   when estimating survival.
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (columns) and timepoint (rows).
##' @param type type of predictions to be computed. Admitted values are
##' \itemize{
##' \item `"predictive"` (default): posterior predictive distribution;
##' \item `"expected"`: theoretical mean of the posterior predictive distribution;
##' \item `"latent"`: latent density (i.e., disregarding the observation error);
##' }
##' @param seed a seed used for the predictions. predictions are obtained
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
##' @details The current version of the code reorders the data where predictions
##'   are needed is ordered by "site" and "site" and, in addition, the sites
##'   MUST match ones used to obtain the parameters' estimates from the the
##'   \code{drm} object.
##'
##' @author lcgodoy
##'
##' @name preddrm
##' @seealso [fit_drm()]
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
  ord <- order(new_data[[object$cols$site_col]],
               new_data[[object$cols$time_col]])
  new_data <- new_data[ord, ]
  ## number time points for predictions
  ntime_for <- length(unique(new_data[[object$cols$time_col]]))
  time_for <- new_data[[object$cols$time_col]] -
    min(new_data[[object$cols$time_col]]) + 1
  if (!missing(f_test)) {
    stopifnot(NROW(f_test) == object$data[["n_ages"]])
  }
  ## throw an error if predicing on new sites
  new_sites_factor <- factor(new_data[[object$cols$site_col]], 
                             levels = object$cols$site_levels)
  if (any(is.na(new_sites_factor))) {
    stop("new_data contains sites that were not present in the training data.")
  }
  new_sites_stan <- as.integer(new_sites_factor)
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
    list(n_sites = object$data$n_sites,
         n_ages = object$data$n_ages,
         n_time = ntime_for,
         n_time_train = object$data$n_time,
         time = time_for,
         site = new_sites_stan,
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
         n_edges_adj = object$data$n_edges_adj,
         ages_movement = object$data$ages_movement,
         selectivity_at_age = object$data$selectivity_at_age,
         K_r = object$data$K_r,
         X_r = x_rt)
  pred_data$N <- pred_data$n_sites * pred_data$n_time
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
  gq <- pred_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = pred_data,
                        seed = seed,
                        parallel_chains = cores,
                        ...)
  spt <- data.frame(v1 = new_sites_factor,
                    v2 = pred_data$time + max(object$data$time) + object$data$time_init - 1)
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
}

##' @title Predictions based on SDM.
##'
##' @description Considering a new dataset (across the same sites), computes
##'   predictions based on the SDM passed as \code{sdm}.
##'
##' @param object An \code{sdm} object containing the output of a [fit_sdm] call.
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
##' \item `"latent"`: latent density (i.e., disregarding the observation error);
##' }
##' @param cores number of threads used to compute the predictions. If four
##'   chains were used in the \code{sdm}, then four (or less) threads are
##'   recommended.
##' @param ... additional parameters to be passed to
##'   \code{$generated_quantities}
##'
##' @details The current version of the code assumes the data where predictions
##'   are needed is ordered by "site" and "site" and, in addition, its sites
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
##' @seealso [fit_sdm()]
##' @export
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
  ord <- order(new_data[[object$cols$site_col]],
               new_data[[object$cols$time_col]])
  new_data <- new_data[ord, ]
  ## time points for predictions
  ntime_for <- length(unique(new_data[[object$cols$time_col]]))
  time_for <- new_data[[object$cols$time_col]] -
    min(new_data[[object$cols$time_col]]) + 1
  ## throw an error if predicing on new sites
  new_sites_factor <- factor(new_data[[object$cols$site_col]], 
                             levels = object$cols$site_levels)
  if (any(is.na(new_sites_factor))) {
    stop("new_data contains sites that were not present in the training data.")
  }
  new_sites_stan <- as.integer(new_sites_factor)
  ##--- pars from model fitted ----
  pars <- get_fitted_pars(object$data, "sdm")
  fitted_params <-
    object$stanfit$draws(variables = pars)
  pred_data <-
    list(n_sites = object$data$n_sites,
         n_time = ntime_for,
         n_time_train = object$data$n_time,
         time = time_for,
         site = new_sites_stan,
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
  pred_data$N <- pred_data$n_sites * pred_data$n_time
  pred_data$K_z <- NCOL(pred_data$Z)
  pred_data$K_x <- NCOL(pred_data$X)
  ## load compiled predict_sdm
  pred_comp <-
    instantiate::stan_package_model(name = "predict_sdm",
                                    package = "drmr")
  ## computing predictions
  gq <- pred_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = pred_data,
                        seed = seed,
                        parallel_chains = cores,
                        ...)
  spt <- data.frame(v1 = new_sites_factor,
                    v2 = pred_data$time + max(object$data$time) + object$data$time_init - 1)
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
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
