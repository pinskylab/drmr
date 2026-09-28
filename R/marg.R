##' @title Effects_Drminal Relationships with Covariates
##' @description Evaluates and summarizes the effects_drminal relationships between
##'   explanatory variables and recruitment, survival, or absence probability
##'   from a fitted DRM model.
##' 
##' @details The \code{effects_drm} function computes the predicted relationships
##'   across a sequence of values for a focal variable (or variables), holding
##'   all other non-focal variables in the model matrix at zero.
##'   
##'   When \code{summary = TRUE}, the function calculates an equal-tailed
##'   credible interval and the median using \code{\link[posterior]{quantile2}},
##'   which is highly optimized for posterior draws.
##' 
##' @param object An object of class \code{adrm}, typically the output of
##'   \code{fit_drm()}.
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
##'   produces a 90\% credible interval (i.e., the 5th and 95th percentiles)
##'   together with the median.
##' @param ... Additional arguments passed to methods.
##' 
##' @return A \code{data.frame} with the posterior summaries (or draws) for the
##'   specified process. If \code{summary = TRUE}, it also receives the class
##'   \code{eff_drm} to enable automated plotting.
##' 
##' @name effects_drm
##' @export
##' @author lcgodoy
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
                      prob = 0.95, ...) {
  process <- match.arg(process)
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  stopifnot(length(prob) == 1)
  config <- switch(process,
                   rec  = list(form = "formula_rec",
                               param = "beta_r",
                               ilink = exp,
                               out = "recruitment",
                               X_mat = "X_r"),
                   pabs = list(form = "formula_zero",
                               param = "beta_t",
                               ilink = stats::plogis,
                               out = "prob_abs",
                               X_mat = "X_t"),
                   surv = list(form = "formula_surv",
                               param = "beta_s",
                               ## ilink = exp,
                               ilink = stats::plogis,
                               out = "survival",
                               X_mat = "X_m"))
  if (process == "surv") {
    stopifnot(!is.null(object$data$K_m))
  }
  my_formula <- object$formulas[[config$form]]
  if (is.null(newdata)) {
    model_matrix <- object$data[[config$X_mat]]
    grid_list <- list()
    for (v in variable) {
      v_min <- min(model_matrix[, v], na.rm = TRUE)
      v_max <- max(model_matrix[, v], na.rm = TRUE)
      grid_list[[v]] <- seq(v_min, v_max, length.out = n_pts)
    }
    newdata <- expand.grid(grid_list)
  }
  all_vars <- all.vars(stats::delete.response(stats::terms(my_formula)))
  other_vars <- setdiff(all_vars, variable)
  for (v in other_vars) {
    if (!v %in% names(newdata)) newdata[[v]] <- 0
  }
  new_x <- stats::model.matrix(my_formula, newdata)
  est_samples <- draws(object, variables = config$param, format = "matrix")
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
##' @description Renders the effects_drminal plot using ggplot2.
##' @keywords internal
##' @noRd
.ploteffects_drm_gg <- function(x,
                         rug_data,
                         focal_var,
                         process_name,
                         col_low,
                         col_est,
                         col_upp, ...) {
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[focal_var]], y = .data[[col_est]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[col_low]], ymax = .data[[col_upp]]), 
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
##' @description Renders the effects_drminal plot using base R graphics.
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
    graphics::rug(rug_data[[focal_var]], col = grDevices::adjustcolor("black", alpha.f = 0.5))
  }
  invisible(NULL)
}

##' @title Plot Effects_Drminal Relationships for ADRM Objects
##' @description Automatically plots the effects_drminal relationship computed by
##'   \code{effects_drm()}. If \code{ggplot2} is installed, it returns a ggplot object;
##'   otherwise, it falls back to base R graphics.
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
    .ploteffects_drm_gg(x, rug_data, focal_var, process_name, col_low, col_est, col_upp, ...)
  } else {
    .ploteffects_drm_base(x, rug_data, focal_var, process_name, col_low, col_est, col_upp, ...)
  }
}
##' @title Stock-Recruitment Curves for DD Models
##' @description Computes expected recruitment across a grid of Spawning Stock Biomass (density)
##'   values, optionally conditioned on environmental regimes.
##' 
##' @param object An object of class \code{adrm}, typically the output of
##'   \code{fit_drm()}.
##' @param density_grid A numeric vector of density values to evaluate. If \code{NULL}, a default 
##'   sequence is generated from 0 to the maximum observed response.
##' @param newdata An optional \code{data.frame} of environmental covariates. Each row 
##'   represents a distinct environmental regime.
##' @param n_pts Integer. Number of points for the default \code{density_grid} when it is
##'   not provided. Defaults to 100.
##' @param summary Logical. If \code{TRUE}, returns posterior quantiles. If \code{FALSE},
##'   returns the raw posterior draws.
##' @param prob Numeric. Credible interval probability mass. Defaults to 0.95.
##' @param ... Additional arguments.
##' 
##' @return A \code{data.frame} of class \code{dd_curve} (if \code{summary = TRUE}) 
##'   containing the Stock-Recruitment evaluations.
##' 
##' @name dd_curves
##' @export
##' @author lcgodoy
dd_curves <- function(object, ...) {
  UseMethod("dd_curves")
}

##' @rdname dd_curves
##' @export
dd_curves.adrm <- function(object, 
                           density_grid = NULL, 
                           newdata = NULL, 
                           n_pts = 100, 
                           summary = TRUE, 
                           prob = 0.95, ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  
  rec_dd_val <- object$data$rec_dd
  if (rec_dd_val == 2) {
    stop("Density dependence must be enabled ('ricker' or 'bh') to plot SR curves.")
  }
  
  my_formula <- object$formulas$formula_rec
  all_vars <- all.vars(stats::delete.response(stats::terms(my_formula)))
  
  if (is.null(newdata)) {
    if (length(all_vars) == 0) {
      newdata <- data.frame(.dummy = 1)
    } else {
      newdata <- as.data.frame(matrix(0, nrow = 1, ncol = length(all_vars)))
      colnames(newdata) <- all_vars
    }
  } else {
    if (nrow(newdata) > 5) {
      warning("newdata has more than 5 rows. The resulting plot may be cluttered.")
    }
    for (v in all_vars) {
      if (!v %in% names(newdata)) newdata[[v]] <- 0
    }
  }
  
  if (is.null(density_grid)) {
    max_y <- max(object$data$y, na.rm = TRUE)
    if (max_y <= 0) max_y <- 100 # fallback
    density_grid <- seq(0, max_y, length.out = n_pts)
  }
  
  new_x <- stats::model.matrix(my_formula, newdata)
  beta_r_draws <- draws(object, variables = "beta_r", format = "matrix")
  kappa_draws <- draws(object, variables = "kappa", format = "matrix")
  
  # log-productivity (alpha)
  lin_pred <- tcrossprod(beta_r_draws, new_x)
  alpha_draws <- exp(lin_pred)
  
  n_draws <- nrow(alpha_draws)
  n_regimes <- nrow(newdata)
  n_ssb <- length(density_grid)
  
  if (summary) {
    alpha_prob <- (1 - prob) / 2
    probs <- c(alpha_prob, 0.5, 1 - alpha_prob)
    
    out_list <- vector("list", n_regimes * n_ssb)
    counter <- 1
    
    for (i in seq_len(n_regimes)) {
      alpha_i <- alpha_draws[, i]
      for (s in density_grid) {
        if (rec_dd_val == 0) {
          # Ricker
          r_draws <- alpha_i * s * exp(-kappa_draws[, 1] * s)
        } else {
          # Beverton-Holt
          r_draws <- (alpha_i * s) / (kappa_draws[, 1] + s)
        }
        
        quants <- posterior::quantile2(r_draws, probs = probs)
        row_data <- cbind(newdata[i, , drop = FALSE], density = s)
        row_data <- cbind(row_data, as.data.frame(t(quants)))
        row_data$.regime <- i
        
        out_list[[counter]] <- row_data
        counter <- counter + 1
      }
    }
    
    output <- do.call(rbind, out_list)
    attr(output, "quant_cols") <- names(posterior::quantile2(1:10, probs = probs))
    attr(output, "model_type") <- ifelse(rec_dd_val == 0, "Ricker", "Beverton-Holt")
    attr(output, "prob") <- prob
    # Remove dummy column if it exists
    if (".dummy" %in% names(output)) output$.dummy <- NULL
    
    class(output) <- c("dd_curve", "data.frame")
    return(output)
  } else {
    # If summary is FALSE, return raw draws
    out_list <- vector("list", n_regimes * n_ssb)
    counter <- 1
    for (i in seq_len(n_regimes)) {
      alpha_i <- alpha_draws[, i]
      for (s in density_grid) {
        if (rec_dd_val == 0) {
          r_draws <- alpha_i * s * exp(-kappa_draws[, 1] * s)
        } else {
          r_draws <- (alpha_i * s) / (kappa_draws[, 1] + s)
        }
        
        row_data <- newdata[rep(i, n_draws), , drop = FALSE]
        row_data$density <- s
        row_data$recruitment <- r_draws
        row_data$.draw <- seq_len(n_draws)
        row_data$.regime <- i
        
        out_list[[counter]] <- row_data
        counter <- counter + 1
      }
    }
    output <- do.call(rbind, out_list)
    if (".dummy" %in% names(output)) output$.dummy <- NULL
    return(output)
  }
}

##' @title Internal ggplot2 Backend for dd_curves
##' @keywords internal
##' @noRd
.plotdd_curves_gg <- function(x, col_low, col_est, col_upp, ...) {
  # determine grouping variables (everything except density, the quantiles, and .regime)
  skip_cols <- c("density", col_low, col_est, col_upp, ".regime")
  group_vars <- setdiff(names(x), skip_cols)
  
  if (length(group_vars) > 0) {
    # Create a single interaction string for coloring/grouping
    x$.group <- apply(x[, group_vars, drop = FALSE], 1, function(row) paste(row, collapse = "-"))
    p <- ggplot2::ggplot(x, ggplot2::aes(x = density, y = .data[[col_est]], color = .group, fill = .group)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[col_low]], ymax = .data[[col_upp]]), alpha = 0.2, color = "transparent") +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(y = "Expected Recruitment", color = "Regime", fill = "Regime")
  } else {
    p <- ggplot2::ggplot(x, ggplot2::aes(x = density, y = .data[[col_est]])) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[col_low]], ymax = .data[[col_upp]]), alpha = 0.2, fill = "black", color = "transparent") +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(y = "Expected Recruitment")
  }
  
  model_type <- attr(x, "model_type")
  return(p)
}

##' @title Internal Base R Backend for dd_curves
##' @keywords internal
##' @noRd
.plotdd_curves_base <- function(x, col_low, col_est, col_upp, ...) {
  regimes <- unique(x$.regime)
  n_regimes <- length(regimes)
  
  y_max <- max(x[[col_upp]], na.rm = TRUE)
  x_max <- max(x$density, na.rm = TRUE)
  
  plot(1, 1, type = "n", xlim = c(0, x_max), ylim = c(0, y_max),
       xlab = "Density", ylab = "Expected Recruitment", 
       ...)
  
  cols <- grDevices::rainbow(n_regimes)
  
  for (i in seq_along(regimes)) {
    sub_x <- x[x$.regime == regimes[i], ]
    x_vals <- sub_x$density
    y_est <- sub_x[[col_est]]
    y_low <- sub_x[[col_low]]
    y_upp <- sub_x[[col_upp]]
    
    col_fill <- grDevices::adjustcolor(cols[i], alpha.f = 0.2)
    graphics::polygon(x = c(x_vals, rev(x_vals)), 
                      y = c(y_low, rev(y_upp)), 
                      col = col_fill, border = NA)
    graphics::lines(x_vals, y_est, col = cols[i], lwd = 2)
  }
  
  if (n_regimes > 1) {
    skip_cols <- c("density", col_low, col_est, col_upp, ".regime")
    group_vars <- setdiff(names(x), skip_cols)
    if (length(group_vars) > 0) {
      leg_labels <- sapply(regimes, function(r) {
        sub_x <- x[x$.regime == r, ][1, group_vars, drop = FALSE]
        paste(sub_x, collapse = "-")
      })
      graphics::legend("topleft", legend = leg_labels, col = cols, lwd = 2, bty = "n")
    }
  }
  invisible(NULL)
}

##' @title Plot Stock-Recruitment Curves for ADRM Objects
##' @description Automatically plots the stock-recruitment curves computed by
##'   \code{dd_curves()}. Uses ggplot2 if available, otherwise base R graphics.
##' 
##' @param x An object of class \code{dd_curve}, usually the output of
##'   \code{dd_curves()}.
##' @param ... Additional arguments passed to the underlying plotting functions.
##' 
##' @export
plot.dd_curve <- function(x, ...) {
  quant_cols <- attr(x, "quant_cols")
  if (length(quant_cols) < 3) {
    stop("Plotting requires at least 3 probabilities (lower, median, upper).")
  }
  col_low <- quant_cols[1]
  col_est <- quant_cols[2]
  col_upp <- quant_cols[3]
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    .plotdd_curves_gg(x, col_low, col_est, col_upp, ...)
  } else {
    .plotdd_curves_base(x, col_low, col_est, col_upp, ...)
  }
}
