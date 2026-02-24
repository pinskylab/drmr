##' @title Fit the dynamic range model.
##' @export
##' @family models
##' @description Fit the DRM Stan model (this function is not stable yet and
##'   has not been fully tested).
##' @param .data A \code{data frame} containing the data for the model.
##' @param y_col A \code{character} specifying the name of the column in `.data`
##'   that contains the response variable.
##' @param time_col A \code{character} specifying the name of the column in
##'   `.data` that contains the time variable.
##' @param site_col A \code{character} specifying the name of the column in
##'   `.data` that contains the site variable.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{ \item
##'   \code{"gamma"} (default): gamma parametrized in terms of its mean; \item
##'   \code{"lognormal"}: log-normal parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its median
##'   (usual parametrization); \item \code{"lognormal_legacy"}: log-normal with
##'   its usual parametrization; }
##' @param formula_zero A \code{formula} specifying the model for the zero
##'   inflation component. Defaults to `~ 1` (intercept only).
##' @param formula_rec A \code{formula} specifying the model for the recruitment
##'   component. Defaults to `~ 1` (intercept only).
##' @param formula_surv A \code{formula} specifying the model for the survival
##'   component. If `NULL` (the default), no survival component is included.
##' @param seed An \code{integer} specifying the random number seed.
##' @param init A scalar specifying the initialization method. The default
##'   ("cmdstan_default") lets \code{cmdstan} initialize parameters.
##'   Other options include: a scalar greater than zero, say \code{x}, which
##'   initializes all parameters uniformly between \code{-x} and \code{x};
##'   \code{0}, which initializes all parameters at \code{0}; or "prior",
##'   which initializes parameters by sampling from their priors.
##' @param algorithm a \code{character} specifying the algorithm used for
##'   inference. Default is \code{nuts} (the default MCMC in Stan). The
##'   remaining options are different flavors of variational bayes algorithms:
##'   "vb" (for ADVI), "pathfinder" (for Pathfinder), "laplace" (normal
##'   approximation centered at the mode of the posterior) or "optimize" for
##'   (penalized) MLEs.
##' @param algo_args a \code{list} with arguments for the sampling
##'   algorithms. For instance, \code{tol_rel_obj} for variational inference.
##' @param ... Passed on to the [make_data()] function used to build the input
##'   \code{list} for our \code{cmdstanr} model.
##' @return An object of class \code{adrm} which is a \code{list} containing the
##'   MCMC draws, the model data, the linear predictors formulas, and the
##'   (response, time, site) column names.
##'   Specifically: \itemize{ \item \code{stanfit}: The MCMC draws from the
##'   fitted model.  \item \code{data}: The data used to fit the model (as a
##'   list).  \item \code{formulas}: The formulas used to create design
##'   matrices.  \item \code{cols}: Important column names.  }
##' @seealso [make_data()]
##' @examples
##' if (instantiate::stan_cmdstan_exists()) {
##'   data(sum_fl)
##'   fit_drm(.data = sum_fl,
##'           y_col = "y",
##'           time_col = "year",
##'           site_col = "patch",
##'           seed = 2025)$stanfit$summary()
##' }
##' @author lcgodoy
fit_drm <- function(.data,
                    y_col,
                    time_col,
                    site_col,
                    family = "gamma",
                    formula_zero = ~ 1,
                    formula_rec = ~ 1,
                    formula_surv = NULL,
                    seed,
                    init = "cmdstan_default",
                    algorithm = "nuts",
                    algo_args = list(),
                    ...) {
  cl <- match.call()
  stopifnot(length(init) == 1)
  stopifnot(length(algorithm) == 1)
  if (is.character(init))
    stopifnot(init %in% c("cmdstan_default", "prior"))
  x_t <- stats::model.matrix.lm(formula_zero, data = .data,
                                na.action = stats::na.fail)
  x_r <- stats::model.matrix.lm(formula_rec, data = .data,
                                na.action = stats::na.fail)
  if (is.null(formula_surv)) {
    model_dat <- make_data(y = .data[[y_col]],
                           time = .data[[time_col]],
                           site = .data[[site_col]],
                           family = family,
                           x_t = x_t,
                           x_r = x_r,
                           ...)
  } else {
    x_m <- stats::model.matrix.lm(formula_surv, data = .data,
                                  na.action = stats::na.fail)
    model_dat <- make_data(y = .data[[y_col]],
                           time = .data[[time_col]],
                           site = .data[[site_col]],
                           family = family,
                           x_t = x_t,
                           x_r = x_r,
                           x_m = x_m,
                           ...)
  }
  model <- instantiate::stan_package_model(
                            name = "drm",
                            package = "drmr"
                        )
  .algo_args <- default_algo(algorithm, algo_args)
  if (init == "cmdstan_default") {
    drm_init <- NULL
  } else if (init == "prior") {
    chains <- ifelse("chains" %in% names(algo_args), algo_args[["chains"]],
              ifelse("num_paths" %in% names(algo_args), algo_args[["num_paths"]],
                     1))
    drm_init <-
      prior_inits(model_dat, chains, "drm")
  } else {
    drm_init <- init
  }

  my_args <- c(list(
      data = model_dat,
      seed = seed,
      init = drm_init
  ), .algo_args)
  if (algorithm == "nuts") {
    draws <- do.call(model$sample, my_args)
  } else if (algorithm == "vb") {
    draws <- do.call(model$variational, my_args)
  } else if (algorithm == "optimize") {    
    draws <- do.call(model$optimize, my_args)    
  } else if (algorithm == "pathfinder") {
    draws <- do.call(model$pathfinder, my_args)
  } else if (algorithm == "laplace") {
    draws <- do.call(model$laplace, my_args)
  } else {
    stop("Algorithm must one of the following: 'nuts', 'vb', 'optimize', 'laplace', or 'pathfinder'")
  }
  output <-
    list("stanfit"  = draws,
         "data"     = model_dat,
         "formulas" = list("formula_zero" = formula_zero,
                           "formula_rec" = formula_rec,
                           "formula_surv" = formula_surv),
         "cols" = list("y_col" = y_col,
                       "time_col" = time_col,
                       "site_col" = site_col),
         "seed" = seed,
         "call" = cl)
  return(new_adrm(output))
}

##' @title Fit a GLM based SDM
##' @export
##' @family models
##' @description Fit the SDM Stan model (this function is not stable yet and
##'   has not been fully tested).
##' @param .data A \code{data frame} containing the data for the model.
##' @param y_col A \code{character} specifying the name of the column in `.data`
##'   that contains the response variable.
##' @param time_col A \code{character} specifying the name of the column in
##'   `.data` that contains the time variable.
##' @param site_col A \code{character} specifying the name of the column in
##'   `.data` that contains the site variable.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{
##'   \item \code{"gamma"} (default): gamma parametrized in terms of its mean;
##'   \item \code{"lognormal"}: log-normal parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its median
##'   (usual parametrization);
##'   \item \code{"lognormal_legacy"}: log-normal with its usual parametrization;
##'    }
##' @param formula_zero A \code{formula} specifying the model for the zero
##'   inflation component. Defaults to `~ 1` (intercept only).
##' @param formula_dens A \code{formula} specifying the model for the non-zero
##'   density component. Defaults to `~ 1` (intercept only).
##' @param seed An \code{integer} specifying the random number seed.
##' @param init A scalar specifying the initialization method. The default
##'   ("cmdstan_default") lets \code{cmdstan} initialize parameters.
##'   Other options include: a scalar greater than zero, say \code{x}, which
##'   initializes all parameters uniformly between \code{-x} and \code{x};
##'   \code{0}, which initializes all parameters at \code{0}; or "prior",
##'   which initializes parameters by sampling from their priors.
##' @param algorithm a \code{character} specifying the algorithm used for
##'   inference. Default is \code{nuts} (the default MCMC in Stan). The
##'   remaining options are different flavors of variational bayes algorithms:
##'   "vb" (for ADVI), "pathfinder" (for Pathfinder), "laplace" (normal
##'   approximation centered at the mode of the posterior) or "optimize" for
##'   (penalized) MLEs.
##' @param algo_args a \code{list} with arguments for the sampling
##'   algorithms. For instance, \code{tol_rel_obj} for variational inference.
##' @param ... Passed on to the [make_data()] function used to build the input
##'   \code{list} for our \code{cmdstanr} model.
##' @return An object of class \code{sdm} which is a \code{list} containing the
##'   MCMC draws, the model data, the linear predictors formulas, and the
##'   (response, time, site) column names.
##'    \itemize{
##'     \item \code{stanfit}: The MCMC draws from the fitted model.
##'     \item \code{data}: The data used to fit the model (as a list).
##'     \item \code{formulas}: The formulas used to create design matrices.
##'     \item \code{cols}: Important column names.
##'   }
##' @seealso [make_data_sdm()]
##' @examples
##' if (instantiate::stan_cmdstan_exists()) {
##'   data(sum_fl)
##'   fit_sdm(.data = sum_fl,
##'           y_col = "y",
##'           time_col = "year",
##'           site_col = "patch",
##'           seed = 2025)$stanfit$summary()
##' }
##' @author lcgodoy
fit_sdm <- function(.data,
                    y_col,
                    time_col,
                    site_col,
                    family = "gamma",
                    formula_zero = ~ 1,
                    formula_dens = ~ 1,
                    seed,
                    init = "cmdstan_default",
                    algorithm = "nuts",
                    algo_args = list(),
                    ...) {
  cl <- match.call()
  stopifnot(length(init) == 1)
  if (is.character(init))
    stopifnot(init %in% c("cmdstan_default", "prior"))
  x_t <- stats::model.matrix.lm(formula_zero, data = .data,
                                na.action = stats::na.fail)
  x_r <- stats::model.matrix.lm(formula_dens, data = .data,
                                na.action = stats::na.fail)
  model_dat <- make_data_sdm(y = .data[[y_col]],
                             time = .data[[time_col]],
                             site = .data[[site_col]],
                             family = family,
                             z = x_t,
                             x = x_r,
                             ...)
  model <- instantiate::stan_package_model(
                            name = "sdm",
                            package = "drmr"                        )
  if (init == "cmdstan_default") {
    sdm_init <- NULL
  } else if (init == "prior") {
    chains <- ifelse("chains" %in% names(algo_args), algo_args[["chains"]],
              ifelse("num_paths" %in% names(algo_args), algo_args[["num_paths"]],
                     1))
    sdm_init <-
      prior_inits(model_dat, chains, "sdm")
  } else {
    sdm_init <- init
  }
  .algo_args <- default_algo(algorithm, algo_args)
  my_args <- c(list(
      data = model_dat,
      seed = seed,
      init = sdm_init
  ), .algo_args)
  if (algorithm == "nuts") {
    draws <- do.call(model$sample, my_args)
  } else if (algorithm == "vb") {
    draws <- do.call(model$variational, my_args)
  } else if (algorithm == "optimize") {    
    draws <- do.call(model$optimize, my_args)    
  } else if (algorithm == "pathfinder") {
    draws <- do.call(model$pathfinder, my_args)
  } else if (algorithm == "laplace") {
    draws <- do.call(model$laplace, my_args)
  } else {
    stop("Algorithm must one of the following: 'nuts', 'vb', 'optimize', 'laplace' or 'pathfinder'")
  }
  output <-
    list("stanfit"  = draws,
         "data"     = model_dat,
         "formulas" = list("formula_zero" = formula_zero,
                           "formula_dens" = formula_dens),
         "cols" = list("y_col" = y_col,
                       "time_col" = time_col,
                       "site_col" = site_col),
         "seed" = seed,
         "call" = cl)
  return(new_sdm(output))
}
