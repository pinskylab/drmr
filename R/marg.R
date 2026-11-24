##' @title Relationships between covariates and demographic processes
##' @description Evaluates and summarizes the relationships between
##'   explanatory variables and recruitment, survival, or absence probability
##'   from a fitted DRM model, or from the prior draws of [prior_pred_drm()].
##'
##' @details The \code{effects_drm} function computes the predicted
##'   relationships across a sequence of values for a focal variable (or
##'   variables), holding all other non-focal variables in the model matrix at
##'   zero.
##'
##'   For \code{process = "pabs"}, the function uses the link function of the
##'   model: the logit link, or the complementary log-log link when the model
##'   has \code{cloglog = 1}. When the model has \code{rho_mu = 1} (the
##'   default), the probability of absence also depends on the expected
##'   density \eqn{\mu} through the term \eqn{\xi \log(\mu)}. The function
##'   leaves out this term, so the curve is the probability of absence when
##'   \eqn{\mu = 1}.
##'
##'   For \code{process = "surv"}, the model must have \code{est_surv = 1}.
##'
##'   For the output of [prior_pred_drm()], the function uses the prior draws
##'   of the coefficients, so the intervals are prior intervals. It uses all
##'   the draws, including those that \code{nonfinite} flags: their parameter
##'   values are valid draws from the priors.
##'
##'   When \code{summary = TRUE}, the function calculates an equal-tailed
##'   credible interval and the median using \code{\link[posterior]{quantile2}},
##'   which is highly optimized for posterior draws.
##'
##'   The function builds the design matrix again from \code{newdata}. Terms
##'   that depend on the data, such as \code{scale()} or \code{poly()}, take
##'   their values from \code{newdata}, not from the training data. To
##'   standardize a covariate, use the mean and standard deviation from the
##'   training data. A factor must have the same levels as in the training
##'   data.
##'
##' @param object The output of [fit_drm()] or of [prior_pred_drm()].
##' @param process A character string indicating the process to evaluate:
##'   \code{"rec"} (recruitment), \code{"surv"} (survival), or \code{"pabs"}
##'   (probability of absence).
##' @param variable A character vector with the name(s) of the focal variable(s)
##'   to examine.
##' @param newdata An optional \code{data.frame} containing the values for the
##'   focal variable(s). If \code{NULL}, a grid is generated automatically based
##'   on the observed range in the model matrix.
##' @param n_pts An integer specifying the number of points to generate for the
##'   sequence of each focal variable when \code{newdata} is
##'   \code{NULL}. Default is 100.
##' @param summary Logical. If \code{TRUE} (the default), returns the quantiles
##'   of the posterior predictions. If \code{FALSE}, returns the raw posterior
##'   draws.
##' @param prob A numeric scalar in \eqn{(0, 1)} specifying the probability mass
##'   of the equal-tailed credible interval. Defaults to \code{0.9}, which
##'   produces a 90% credible interval (i.e., the 5th and 95th percentiles)
##'   together with the median.
##' @param ... Additional arguments passed to methods.
##' 
##' @return A \code{data.frame} with the posterior summaries (or draws) for the
##'   specified process. If \code{summary = TRUE}, it also receives the class
##'   \code{eff_drm} to enable automated plotting.
##' 
##' @name effects_drm
##' @export
##' @author Anonymous
effects_drm <- function(object, ...) {
  UseMethod("effects_drm")
}

##' @rdname effects_drm
##' @export
effects_drm.adrm <- function(object,
                      process = c("rec", "surv", "pabs"),
                      variable, 
                      newdata = NULL,
                      n_pts = 100,
                      summary = TRUE,
                      prob = 0.9, ...) {
  process <- match.arg(process)
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(length(prob) == 1)
  config <- effects_config(process, object$data)
  effects_from_draws(
    formula = object$formulas[[config$form]],
    model_matrix = object$data[[config$X_mat]],
    est_samples = draws(object, variables = config$param, format = "matrix"),
    config = config,
    variable = variable,
    newdata = newdata,
    n_pts = n_pts,
    summary = summary,
    prob = prob
  )
}

##' @rdname effects_drm
##' @export
effects_drm.prior_pred_drm <- function(object,
                                       process = c("rec", "surv", "pabs"),
                                       variable,
                                       newdata = NULL,
                                       n_pts = 100,
                                       summary = TRUE,
                                       prob = 0.9, ...) {
  process <- match.arg(process)
  stopifnot(length(prob) == 1)
  config <- effects_config(process, object$model_data)
  ## the names of the draws are the ones in drm.stan, as in a fit
  est_samples <- posterior::as_draws_matrix(
    posterior::subset_draws(object$pars, variable = config$param)
  )
  effects_from_draws(
    formula = object$formulas[[config$form]],
    model_matrix = object$model_data[[config$X_mat]],
    est_samples = est_samples,
    config = config,
    variable = variable,
    newdata = newdata,
    n_pts = n_pts,
    summary = summary,
    prob = prob
  )
}

##' @title Settings of each process for effects_drm()
##' @description Checks that the model has the process, and returns the
##'   settings of the process. The inverse link for \code{"pabs"} follows the
##'   \code{cloglog} toggle of the model.
##' @param process \code{"rec"}, \code{"surv"}, or \code{"pabs"}.
##' @param model_data the \code{list} that [make_data()] created for the
##'   model.
##' @return a \code{list} with the name of the formula (\code{form}), the
##'   name of the coefficients (\code{param}), the inverse link function
##'   (\code{ilink}), the name of the output column (\code{out}), and the
##'   name of the design matrix in the model data (\code{X_mat}).
##' @keywords internal
##' @noRd
effects_config <- function(process, model_data) {
  if (process == "surv" && !isTRUE(model_data$est_surv == 1)) {
    stop("'process = \"surv\"' needs a model with ",
         "'.toggles = list(est_surv = 1)'.")
  }
  ilink_pabs <- if (isTRUE(model_data$cloglog == 1)) {
    function(x) -expm1(-exp(x))
  } else {
    stats::plogis
  }
  switch(process,
         rec  = list(form = "formula_rec",
                     param = "beta_r",
                     ilink = exp,
                     out = "recruitment",
                     X_mat = "X_r"),
         pabs = list(form = "formula_zero",
                     param = "beta_t",
                     ilink = ilink_pabs,
                     out = "prob_abs",
                     X_mat = "X_t"),
         surv = list(form = "formula_surv",
                     param = "beta_s",
                     ilink = stats::plogis,
                     out = "survival",
                     X_mat = "X_m"))
}

##' @title Relationships between covariates and a process, from draws
##' @description The engine behind the [effects_drm()] methods.
##' @param formula the formula of the process.
##' @param model_matrix the design matrix of the process in the model data.
##'   The function takes the range of each focal variable from it when
##'   \code{newdata} is \code{NULL}.
##' @param est_samples a \code{matrix} of draws of the coefficients, with one
##'   row per draw and one column per column of the design matrix.
##' @param config the output of \code{effects_config()}.
##' @inheritParams effects_drm
##' @return See [effects_drm()].
##' @keywords internal
##' @noRd
effects_from_draws <- function(formula,
                               model_matrix,
                               est_samples,
                               config,
                               variable,
                               newdata,
                               n_pts,
                               summary,
                               prob) {
  if (is.null(newdata)) {
    grid_list <- list()
    for (v in variable) {
      v_min <- min(model_matrix[, v], na.rm = TRUE)
      v_max <- max(model_matrix[, v], na.rm = TRUE)
      grid_list[[v]] <- seq(v_min, v_max, length.out = n_pts)
    }
    newdata <- expand.grid(grid_list)
  }
  all_vars <- all.vars(stats::delete.response(stats::terms(formula)))
  other_vars <- setdiff(all_vars, variable)
  for (v in other_vars) {
    if (!v %in% names(newdata)) newdata[[v]] <- 0
  }
  new_x <- stats::model.matrix(formula, newdata)
  lin_pred <- tcrossprod(est_samples, new_x)
  est_ilink <- config$ilink(lin_pred)
  if (summary) {
    alpha <- (1 - prob) / 2
    probs <- c(alpha, 0.5, 1 - alpha)
    quants <- apply(est_ilink, 2, posterior::quantile2, probs = probs)
    quants_df <- as.data.frame(t(quants))
    
    output <- cbind(newdata, quants_df)
    attr(output, "process") <- config$out
    attr(output, "variable") <- variable
    attr(output, "quant_cols") <- colnames(quants_df)
    attr(output, "prob") <- prob
    class(output) <- c("eff_drm", "data.frame")
    return(output)
  }
  n_sim <- NROW(est_samples)
  idx <- rep(seq_len(NROW(newdata)), each = n_sim)
  output <- newdata[idx, , drop = FALSE]
  rownames(output) <- NULL
  output[[config$out]] <- c(est_ilink)
  return(output)
}

##' @title Internal ggplot2 Backend for eff_drm
##' @description Renders the plot of an \code{eff_drm} object using ggplot2.
##' @keywords internal
##' @noRd
.ploteffects_drm_gg <- function(x,
                         rug_data,
                         focal_var,
                         process_name,
                         col_low,
                         col_est,
                         col_upp, ...) {
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[focal_var]],
                                       y = .data[[col_est]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[col_low]],
                                      ymax = .data[[col_upp]]),
                         alpha = 0.4, fill = "black", color = "transparent") +
    ggplot2::geom_line(ggplot2::aes(y = .data[[col_low]]), linetype = 2) +
    ggplot2::geom_line(ggplot2::aes(y = .data[[col_upp]]), linetype = 2) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(y = process_name, x = focal_var)
  
  if (!is.null(rug_data) && focal_var %in% names(rug_data)) {
    p <- p + ggplot2::geom_rug(data = rug_data, 
                               ggplot2::aes(x = .data[[focal_var]]), 
                               inherit.aes = FALSE, alpha = 0.5)
  }
  
  return(p)
}

##' @title Internal Base R Backend for eff_drm
##' @description Renders the plot of an \code{eff_drm} object using base R
##'   graphics.
##' @keywords internal
##' @noRd
.ploteffects_drm_base <- function(x,
                           rug_data,
                           focal_var,
                           process_name,
                           col_low,
                           col_est,
                           col_upp,
                           ...) {
  x_vals <- x[[focal_var]]
  y_est <- x[[col_est]]
  y_low <- x[[col_low]]
  y_upp <- x[[col_upp]]
  plot(x_vals, y_est, type = "n", 
       ylim = range(c(y_low, y_upp), na.rm = TRUE),
       xlab = focal_var, ylab = process_name, ...)
  graphics::polygon(x = c(x_vals, rev(x_vals)), 
                    y = c(y_low, rev(y_upp)), 
                    col = grDevices::adjustcolor("black", alpha.f = 0.4), 
          border = NA)
  graphics::lines(x_vals, y_low, lty = 2)
  graphics::lines(x_vals, y_upp, lty = 2)
  graphics::lines(x_vals, y_est, lwd = 2)
  if (!is.null(rug_data) && focal_var %in% names(rug_data)) {
    graphics::rug(rug_data[[focal_var]],
                  col = grDevices::adjustcolor("black", alpha.f = 0.5))
  }
  invisible(NULL)
}

##' @title Plot the relationships computed by effects_drm()
##' @description Automatically plots the relationship computed by
##'   \code{effects_drm()}. If \code{ggplot2} is installed, it returns a
##'   ggplot object; otherwise, it falls back to base R graphics.
##' 
##' @param x An object of class \code{eff_drm}, usually the output of
##'   \code{effects_drm()}.
##' @param rug_data An optional \code{data.frame} containing the original data
##'   to add a rug plot to the x-axis.
##' @param ... Additional arguments passed to the underlying plotting functions.
##' 
##' @details This default plotting method is restricted to outputs where exactly
##'   one focal variable was evaluated (\code{length(variable) ==
##'   1}). Visualizing more complicated cases (e.g., 2D interactions) should be
##'   handled manually by the user.
##' 
##' @export
plot.eff_drm <- function(x, rug_data = NULL, ...) {
  process_name <- attr(x, "process")
  focal_var <- attr(x, "variable")
  quant_cols <- attr(x, "quant_cols")
  if (length(focal_var) > 1) {
    stop("More complicated cases should be handled by the user.")
  }
  if (length(quant_cols) < 3) {
    stop("Plotting requires at least 3 probabilities (lower, median, upper).")
  }
  col_low <- quant_cols[1]
  col_est <- quant_cols[2]
  col_upp <- quant_cols[3]
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .ploteffects_drm_gg(x, rug_data, focal_var, process_name,
                        col_low, col_est, col_upp, ...)
  } else {
    .ploteffects_drm_base(x, rug_data, focal_var, process_name,
                          col_low, col_est, col_upp, ...)
  }
}
