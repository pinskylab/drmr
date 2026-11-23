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
##'   predictions based on a fitted DRM.
##'
##' @param object The output from the [fit_drm()] function.
##' @param newdata a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions.
##' @param past_data a \code{data.frame} with the dataset from the last year
##'   used in model fitting. Only needed when \code{f_test} is not missing or
##'   when estimating survival.
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (rows) and time point (columns).
##' @param type type of predictions to be computed. Admitted values are
##'   \itemize{ \item `"predictive"` (default): posterior predictive
##'   distribution; \item `"expected"`: theoretical mean of the posterior
##'   predictive distribution; \item `"latent"`: latent density (i.e.,
##'   disregarding the observation error); }
##' @param seed a seed used for the predictions. Predictions are obtained
##'   through Monte Carlo samples from the posterior predictive
##'   distribution. Therefore, a \code{seed} is needed to ensure the results'
##'   reproducibility.
##' @param cores number of threads used for the computations. If four chains
##'   were used to fit the model, then four (or less) threads are recommended.
##' @param ... additional parameters to be passed to \code{$generate_quantities}
##'
##'
##' @details The function sorts \code{newdata} by site, in the site order of
##'   the training data, and then by time. Each site used to fit
##'   \code{object} must appear exactly once in each year of \code{newdata}.
##'   The years in \code{newdata} must be consecutive and must start in the
##'   year after the last training year.
##'
##'   The function builds the design matrices again from \code{newdata}.
##'   Terms that depend on the data, such as \code{scale()} or \code{poly()},
##'   take their values from \code{newdata}, not from the training data. To
##'   standardize a covariate, use the mean and standard deviation from the
##'   training data. A factor must have the same levels as in the training
##'   data.
##'
##' @author Anonymous
##'
##' @name preddrm
##' @seealso [fit_drm()]
##' 
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for predictions.
##'
##' @export
predict.adrm <- function(object,
                         newdata,
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
  ## predict.stan reads the rows as one block per site code, so sort by the
  ## site order of the training data (not by the raw site values)
  ord <- order(factor(newdata[[object$cols$site_col]],
                      levels = object$cols$site_levels),
               newdata[[object$cols$time_col]])
  newdata <- newdata[ord, ]
  ## throw an error if predicing on new sites
  new_sites_factor <- factor(newdata[[object$cols$site_col]],
                             levels = object$cols$site_levels)
  if (any(is.na(new_sites_factor))) {
    stop("newdata contains sites that were not present in the training data.")
  }
  ## the forecast steps one year at a time from the last training year
  last_train <- object$data$time_init + max(object$data$time) - 1
  years <- sort(unique(newdata[[object$cols$time_col]]))
  if (anyNA(newdata[[object$cols$time_col]]) || length(years) == 0 ||
        years[1] != last_train + 1 || any(diff(years) != 1)) {
    stop("The years in newdata must have no missing values, be ",
         "consecutive, and start in ", last_train + 1,
         ", the year after the last training year.")
  }
  ## number time points for predictions
  ntime_for <- length(unique(newdata[[object$cols$time_col]]))
  time_for <- newdata[[object$cols$time_col]] -
    min(newdata[[object$cols$time_col]]) + 1
  if (!missing(f_test)) {
    stopifnot(NROW(f_test) == object$data[["n_ages"]])
  }
  ## every training site must appear exactly once in each year
  if (any(table(new_sites_factor, time_for) != 1)) {
    stop("newdata must have exactly one row for each site and each year.")
  }
  new_sites_stan <- as.integer(new_sites_factor)
  x_tt <- stats::model.matrix(object[["formulas"]][["formula_zero"]],
                              data = newdata)
  x_rt <- stats::model.matrix(object[["formulas"]][["formula_rec"]],
                              data = newdata)
  if (missing(f_test)) {
    f_test <- matrix(0, ncol = ntime_for, nrow = object$data$n_ages)
  }
  ##--- pars from model fitted ----
  ## CmdStan matches the fitted parameters to the parameters block of
  ## predict.stan by position, so pass only the ones it declares (it does
  ## not declare sigma_i, which summary() uses)
  pars <- setdiff(get_fitted_pars(object$data, "drm"), "sigma_i")
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
         n_edges_adj = get_n_edges_adj(object$data),
         areas = get_areas(object$data),
         ages_movement = object$data$ages_movement,
         selectivity_at_age = object$data$selectivity_at_age,
         K_r = object$data$K_r,
         X_r = x_rt)
  pred_data$N <- pred_data$n_sites * pred_data$n_time
  if (length(object$data$K_m) > 0) {
    stopifnot(!missing(past_data))
    ## rows must follow the sites' ordering (site-level random effects)
    past_data <-
      past_data[order(factor(past_data[[object$cols$site_col]],
                             levels = object$cols$site_levels)), ]
    x_mpast <-
      stats::model.matrix(object[["formulas"]][["formula_surv"]],
                          data = past_data)
    x_m <- stats::model.matrix(object[["formulas"]][["formula_surv"]],
                               data = newdata)
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
  check_gq(gq)
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
##'   predictions based on a fitted SDM.
##'
##' @param object The output of a [fit_sdm] call.
##' @inheritParams preddrm
##'
##' @details The function sorts \code{newdata} by site and time. Each site
##'   used to fit \code{object} must appear exactly once in each year of
##'   \code{newdata}. All years in \code{newdata} must come after the last
##'   training year. The years can have gaps: the AR term, if any, also steps
##'   through the missing years.
##'
##'   The function builds the design matrices again from \code{newdata}.
##'   Terms that depend on the data, such as \code{scale()} or \code{poly()},
##'   take their values from \code{newdata}, not from the training data. To
##'   standardize a covariate, use the mean and standard deviation from the
##'   training data. A factor must have the same levels as in the training
##'   data.
##'
##' @author Anonymous
##'
##' @return An object of class \code{CmdStanGQ} (from the \code{instantiate}
##'   package) containing samples for the posterior predictive distribution for
##'   predictions.
##'
##' @rdname predsdm
##' @seealso [fit_sdm()]
##' @export
predict.sdm <- function(object,
                        newdata,
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
  ord <- order(newdata[[object$cols$site_col]],
               newdata[[object$cols$time_col]])
  newdata <- newdata[ord, ]
  ## time points for predictions, counted from the last training year. The
  ## AR term also steps through the years missing from newdata.
  last_train <- object$data$time_init + max(object$data$time) - 1
  time_for <- newdata[[object$cols$time_col]] - last_train
  if (any(time_for < 1)) {
    stop("newdata must contain only years after the last training year (",
         last_train, ").")
  }
  ntime_for <- max(time_for)
  ## throw an error if predicing on new sites
  new_sites_factor <- factor(newdata[[object$cols$site_col]],
                             levels = object$cols$site_levels)
  if (any(is.na(new_sites_factor))) {
    stop("newdata contains sites that were not present in the training data.")
  }
  ## every training site must appear exactly once in each year
  if (any(table(new_sites_factor, time_for) != 1)) {
    stop("newdata must have exactly one row for each site and each year.")
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
         Z = stats::model.matrix.lm(object[["formulas"]][["formula_zero"]],
                                    data = newdata,
                                    na.action = stats::na.fail),
         X = stats::model.matrix.lm(object[["formulas"]][["formula_dens"]],
                                    data = newdata,
                                    na.action = stats::na.fail))
  pred_data$N <- NROW(newdata)
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
  check_gq(gq)
  spt <- data.frame(v1 = new_sites_factor,
                    v2 = newdata[[object$cols$time_col]])
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
}
