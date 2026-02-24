##' @title Marginal Relationships with Covariates
##' @description Evaluates and summarizes the marginal relationships between
##'   explanatory variables and recruitment, survival, or absence probability
##'   from a fitted DRM model.
##' 
##' @details The \code{marg} function computes the predicted relationships
##'   across a sequence of values for a focal variable (or variables), holding
##'   all other non-focal variables in the model matrix at zero.
##'   
##'   When \code{summary = TRUE}, the function calculates credible intervals
##'   using \code{\link[posterior]{quantile2}}, which is highly optimized for
##'   posterior draws.
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
##' @param probs A numeric vector of probabilities to be passed to
##'   \code{posterior::quantile2} when \code{summary = TRUE}. Defaults to
##'   \code{c(0.05, 0.5, 0.95)} (representing the lower bound, median, and upper
##'   bound).
##' @param ... Additional arguments passed to methods.
##' 
##' @return A \code{data.frame} with the posterior summaries (or draws) for the
##'   specified process. If \code{summary = TRUE}, it also receives the class
##'   \code{marg_adrm} to enable automated plotting.
##' 
##' @name marg
##' @export
##' @author lcgodoy
marg <- function(object, ...) {
  UseMethod("marg")
}

##' @rdname marg
##' @export
marg.adrm <- function(object,
                      process = c("rec", "surv", "pabs"),
                      variable, 
                      newdata = NULL,
                      n_pts = 100,
                      summary = TRUE, 
                      probs = c(0.05, 0.5, 0.95), ...) {
  process <- match.arg(process)
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
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
                               ilink = exp,
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
    quants <- apply(est_ilink, 2, posterior::quantile2, probs = probs)
    quants_df <- as.data.frame(t(quants))
    
    output <- cbind(newdata, quants_df)
    attr(output, "process") <- config$out
    attr(output, "variable") <- variable
    attr(output, "quant_cols") <- colnames(quants_df)
    class(output) <- c("marg_adrm", "data.frame")
    return(output)
  }
  n_sim <- NROW(est_samples)
  idx <- rep(seq_len(NROW(newdata)), each = n_sim)
  output <- newdata[idx, , drop = FALSE]
  rownames(output) <- NULL
  output[[config$out]] <- c(est_ilink)
  return(output)
}

##' @title Internal ggplot2 Backend for marg_adrm
##' @description Renders the marginal plot using ggplot2.
##' @keywords internal
##' @noRd
.plotmarg_gg <- function(x,
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

##' @title Internal Base R Backend for marg_adrm
##' @description Renders the marginal plot using base R graphics.
##' @keywords internal
##' @noRd
.plotmarg_base <- function(x,
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

##' @title Plot Marginal Relationships for ADRM Objects
##' @description Automatically plots the marginal relationship computed by
##'   \code{marg()}. If \code{ggplot2} is installed, it returns a ggplot object;
##'   otherwise, it falls back to base R graphics.
##' 
##' @param x An object of class \code{marg_adrm}, usually the output of
##'   \code{marg()}.
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
plot.marg_adrm <- function(x, rug_data = NULL, ...) {
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
    .plotmarg_gg(x, rug_data, focal_var, process_name, col_low, col_est, col_upp, ...)
  } else {
    .plotmarg_base(x, rug_data, focal_var, process_name, col_low, col_est, col_upp, ...)
  }
}
