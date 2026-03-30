##' Update and Re-fit a DRM Model
##'
##' @param object An object of class \code{adrm}.
##' @param ... Arguments to be updated in the new call (e.g., .data, formula_zero, iter_sampling).
##'
##' @return An updated \code{adrm} object.
##' @export
update.adrm <- function(object, ...) {
  call <- object$call
  extras <- list(...)  
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, envir = parent.frame())
}

##' Update and Re-fit a SDM Model
##'
##' @param object An object of class \code{sdm}.
##' @param ... Arguments to be updated in the new call.
##'
##' @return An updated \code{sdm} object.
##' @export
update.sdm <- function(object, ...) {
  call <- object$call
  extras <- list(...)  
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, envir = parent.frame())
}

##' Print method for \code{adrm} and \code{sdm} objects
##'
##' @param x An object of class \code{drmrmodels}.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.drmrmodels <- function(x, ...) {
  is_adrm <- inherits(x, "adrm")
  
  if (is_adrm) {
    cat("Age-structured Dynamic Range Model (ADRM)\n")
  } else {
    cat("Species Distribution Model (SDM)\n")
  }
  cat("\nInference method:\n")
  method <- tryCatch(x$stanfit$metadata()$method, error = function(e) "Unknown") |>
    toupper()
  method <- ifelse(method == "SAMPLE", "NUTS", method)
  cat(paste0("  ", toupper(method), "\n"))
  cat("\nFormulas:\n")
  for (n in names(x$formulas)) {
    if (!is.null(x$formulas[[n]])) {
      if (is_adrm) {
        clean_name <- ifelse(grepl("_zero", n), "Zeros",
                      ifelse(grepl("_rec", n), "Recruitment", "Survival"))
      } else {
        clean_name <- ifelse(grepl("_zero", n), "Zeros", "Non-zero response")
      }
      cat(sprintf("  %s: ", clean_name))
      print(x$formulas[[n]], showEnv = FALSE)
    }
  }
  cat("\nData:\n")
  cat(sprintf("  Timepoints: %d\n", max(x$data$time)))
  cat(sprintf("  Sites: %d\n", length(unique((x$data$site)))))
  if (!is.null(x$data$family)) {
     cat(sprintf("  Non-zero density family: %s\n", x$data$family))
  }
  invisible(x)
}

##' @title Internal helper to generate parameter descriptions
##' @param v A character string representing the parameter name
##' @param covars A list of covariate names extracted from the model object
##' @param model A list of covariate names extracted from the model object
.get_param_description <- function(v, covars,
                                   model = c("drm", "sdm")) {
  model <- match.arg(model)
  base_name <- gsub("\\[.*", "", v)
  if (grepl("\\[", v)) {
    idx_str <- gsub(".*\\[(.*)\\].*", "\\1", v)
    indices <- as.numeric(unlist(strsplit(idx_str, ",")))
    idx <- indices[length(indices)] 
  } else {
    idx <- NA
  }  
  switch(base_name,
         "beta_t"  = sprintf("zero-infl: %s", covars$zero[idx]),
         "beta_r"  = sprintf("%s: %s",
                               ifelse(model == "drm", "rec", "dens"),
                               ifelse(model == "drm",
                                      covars$rec[idx],
                                      covars$dens[idx])),
         "beta_s"  = sprintf("surv: %s", covars$surv[idx]),
         "alpha"   = "Temporal autocorrelation",
         "zeta"    = "Prob. of remaining at the current patch",
         "xi"      = "Relationship between zeros and density",
         "sigma_s" = "Spatial random effect std. dev.",
         "sigma_t" = "Temporal random effect std. dev.",
         "sigma_i" = "IID random effect std. dev.",
         "phi"     = "Dispersion parameter",
         v
  )
}

##' Summary method for \code{adrm} and \code{sdm} objects
##'
##' @param object An object of class \code{drmrmodels}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments (ignored).
##'
##' @export
summary.drmrmodels <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  is_adrm <- inherits(object, "adrm")
  model_type <- ifelse(is_adrm, "drm", "sdm")
  params_out <- get_fitted_pars(object$data, model_type)
  if(object$data$sp_re > 0)
    params_out <- c(params_out, "sigma_s")
  params_out <- params_out[!grepl("^(lambda|z_)", params_out)]
  method <- tryCatch(object$stanfit$metadata()$method, error = function(e) "nuts")  
  if (method == "optimize") {
    res <- object$stanfit$mle(variables = params_out)
    out <- list(estimates = res,
                method = "optimize",
                call = object$call)
  } else {
    summ <- object$stanfit$summary(variables = params_out, 
                                   "mean", "sd",
                                   ~posterior::quantile2(., probs = probs), 
                                   ...)
    out <- list(estimates = summ, method = method, call = object$call, probs = probs)
  }
  covars <- object$covariates
  descriptions <- sapply(out$estimates$variable,
                         .get_param_description,
                         covars = covars,
                         model = model_type)
  out$estimates$description <- unname(descriptions)
  col_order <- c("variable", "description",
                 setdiff(names(out$estimates),
                         c("variable", "description")))
  out$estimates <- out$estimates[, col_order]
  subclass <- ifelse(is_adrm, "summary.adrm", "summary.sdm")
  class(out) <- c(subclass, "summary.drmrmodels", class(out))
  return(out)
}

##' @title Draws method for \code{adrm} and \code{sdm} objects.
##' @description Mirrors the behavior of \code{cmdstanr}'s \code{$draws} method
##'   to retrieve the samples from the posterior distribution of the models'
##'   parameters.
##'
##' @param x An object of class \code{adrm} or \code{sdm}.
##' @param ... currently ignored.
##'
##' @author lcgodoy
##' @return an object of class \code{"drmrdiag"}
##' @name mcmc_diag
##' @export
mcmc_diag <- function(x, ...) UseMethod("mcmc_diag", x)

##' @rdname mcmc_diag
##' @export
mcmc_diag.drmrmodels <- function(x, ...) {
  if (!inherits(x$stanfit, "CmdStanMCMC")) {
    stop("The inference algorithm is not MCMC! Diagnostics are only available for MCMC fits.")
  }
  model_type <- ifelse(inherits(x, "adrm"), "drm", "sdm")
  params_out <- get_fitted_pars(x$data, model_type)
  if (x$data$sp_re > 0) {
    params_out <- c(params_out, "sigma_s")
  }
  params_out <- params_out[!grepl("^(lambda|z_)", params_out)]
  measures <- posterior::default_convergence_measures()
  summ <- x$stanfit$summary(variables = params_out, measures, ...)  
  class(summ) <- c("drmrdiag", class(summ))
  return(summ)
}

##' @export
summary.drmrdiag <- function(object, ...) {
  pars <- object$variable 
  rhats <- object$rhat
  ess_bulk <- object$ess_bulk
  ess_tail <- object$ess_tail
  issues <- list()
  rhat_flag <- !is.na(rhats) & rhats > 1.01
  if (any(rhat_flag)) {
    issues$rhat <- pars[rhat_flag]
  }
  ## rule of thumb
  if (!is.null(ess_bulk)) {
    ess_b_flag <- !is.na(ess_bulk) & ess_bulk < 400
    if (any(ess_b_flag)) {
      issues$ess_bulk <- pars[ess_b_flag]
    }
  }
  if (!is.null(ess_tail)) {
    ess_t_flag <- !is.na(ess_tail) & ess_tail < 400
    if (any(ess_t_flag)) {
      issues$ess_tail <- pars[ess_t_flag]
    }
  }
  if (length(issues) > 0) {
    cat("--- MCMC Diagnostic Issues ---\n")
    if (!is.null(issues$rhat)) {
      cat("\nParameters with R-hat > 1.01 (possible non-convergence):\n")
      print(issues$rhat)
    }
    if (!is.null(issues$ess_bulk)) {
      cat("\nParameters with low Bulk ESS (< 400):\n")
      print(issues$ess_bulk)
    }
    if (!is.null(issues$ess_tail)) {
      cat("\nParameters with low Tail ESS (< 400):\n")
      print(issues$ess_tail)
    }
  } else {
    cat("No issues detected by the MCMC diagnostics.\n")
  }
  invisible(issues)
}

##' @title Print method for summary.drmrmodels
##'
##' @param x An object of class \code{summary.drmrmodels}.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.summary.drmrmodels <- function(x, ...) {
  if (inherits(x, "summary.adrm")) {
    cat("Summary of ADRM Fit\n")
  } else {
    cat("Summary of SDM Fit\n")
  }
  cat("\nEstimates:\n")
  print(x$estimates)
  if (x$method == "optimize") {
    cat("\n(Point estimates from maximum likelihood optimization)\n")
  }
  invisible(x)
}

##' @title Draws method for \code{adrm} and \code{sdm} objects.
##' @description Mirrors the behavior of \code{cmdstanr}'s \code{$draws} method
##'   to retrieve the samples from the posterior distribution of the models'
##'   parameters.
##'
##' @param x An object of class \code{adrm} or \code{sdm}.
##' @param variables a string vector with the name of the variables for which we
##'   want to obtain posteriors samples (defaults to the same variables as the
##'   \code{summary} methods).
##' @param inc_warmup a boolean indicating whether the warmup samples should be
##'   retrieved as well. Defaults to \code{FALSE}.
##' @param format A string. See cmdstanr [$draws
##'   documentation](https://mc-stan.org/cmdstanr/reference/fit-method-draws.html).
##' @param ... currently ignored.
##'
##' @author lcgodoy
##' @return an object of class \code{"draws"} containing the posterior samples
##'   from specified parameters.
##' @name draws
##' @export
draws <- function(x, ...) UseMethod("draws", x)


##' @rdname draws
##' @export
draws.drmrmodels <- function(x,
                             variables = NULL,
                             inc_warmup = FALSE,
                             format = getOption("cmdstanr_draws_format", "draws_array"),
                             ...) {
  if (length(list(...)) > 0) {
    warning("Additional arguments passed to ... are ignored.")
  }
  stopifnot(inherits(x$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                  "CmdStanPathfinder", "CmdStanVB")))
  is_adrm <- inherits(x, "adrm")
  if (is.null(variables)) {
    model_type <- ifelse(is_adrm, "drm", "sdm")
    variables <- get_fitted_pars(x$data, model_type)
    if(x$data$sp_re > 0)
      variables <- c(variables, "sigma_s")
    if (is_adrm) {
      variables <- variables[!grepl("^(lambda|z_)", variables)]
    } else {
      variables <- variables[!grepl("^z_", variables)]
    }
  }
  method <- tryCatch(x$stanfit$metadata()$method, error = function(e) "nuts")  
  if (method == "optimize") {
    stop("Draws are not available for `method = 'optimize'`")
  } else {
    out <- x$stanfit$draws(variables = variables,
                           inc_warmup = inc_warmup,
                           format = format)
  }
  return(out)
}

##' Wrapper for the summary method for \code{CmdStanGQ} objects
##'
##' @param object An object of class \code{CmdStanGQ}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments (ignored).
##'
##' @export
summary.CmdStanGQ <- function(object,
                              probs = c(.05, .5, .95),
                              ...) {
  object$summary(variables = NULL,
                 "mean", "sd", "rhat", "ess_bulk",
                 ~posterior::quantile2(., probs = probs), 
                 ...)
}

##' Summary method for \code{pred_drmr} objects
##'
##' @param object An object of class \code{pred_drmr}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments (ignored).
##'
##' @export
summary.pred_drmr <- function(object,
                              probs = c(.05, .5, .95),
                              ...) {
  output <- cbind(object$spt,
                  summary(object$gq, probs = probs, ...)[, -1])
  return(output)
}

##' Cleaning the variable output from stan for interpretability.
##'
##' @param x A string
##'
##' @export
clean_edens <- function(x) {
  stopifnot(grepl("lambda", x[1]))
  clean_str <- gsub("lambda\\[|\\]", "", x)
  split_list <- strsplit(clean_str, ",")
  return(do.call(rbind, lapply(split_list, as.numeric)))
}

##' Summary method for \code{aesd} objects
##'
##' @param object An object of class \code{aesd}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments (ignored).
##'
##' @export
summary.aesd <- function(object,
                         probs = c(.05, .5, .95),
                         ...) {
  output <- summary(object$gq, probs = probs, ...)
  edens_labels <- as.data.frame(clean_edens(output$variable))
  edens_labels[, 2] <- edens_labels[, 2] + min(object$spt[, 2]) - 1
  colnames(edens_labels) <- c("age", colnames(object$spt)[2], "site_id")
  edens_labels <- merge(edens_labels, object$spt,
                        by = c("site_id", colnames(object$spt)[2]))
  output <- cbind(edens_labels, output[, -1])
  return(output)
}

##' @title Plot Diagnostics for DRM Models
##' @description Provides trace and density plots for the posterior draws of a
##'   DRM (or SDM) model.
##' 
##' @param x An object of class \code{drmrmodels}. Typically, the output of a
##'   \code{fit_drm} or \code{fit_sdm} call.
##' @param variables An optional character vector of parameter names to plot. If
##'   \code{NULL}, plots all.
##' @param type A character vector specifying which plots to draw:
##'   \code{"trace"}, \code{"density"}, or both.
##' @param ask Logical; if \code{TRUE}, the user is prompted before displaying
##'   the next plot. If \code{NULL} (the default), it is calculated dynamically
##'   based on the current graphical layout.
##' @param ... Additional graphical parameters passed to \code{plot}.
##'
##' @importFrom grDevices dev.interactive devAskNewPage palette
##' @importFrom graphics legend lines par
##' @importFrom stats density fitted
##' 
##' @export
plot.drmrmodels <- function(x, 
                            variables = NULL, 
                            type = c("trace", "density"), 
                            ask = NULL,
                            ...) {
  
  if (!inherits(x, "drmrmodels")) {
    stop("Use only with 'drmrmodels' objects.")
  }
  type <- match.arg(type, several.ok = TRUE)
  post <- draws(x) 
  if (length(dim(post)) != 3) {
    stop("Expected draws(x) to return a 3D array: [iterations, chains, parameters].")
  }
  
  n_iter  <- dim(post)[1]
  n_chain <- dim(post)[2]
  p_names <- dimnames(post)[[3]]
  
  if (!is.null(variables)) {
    regex_pattern <- paste0("^(", paste(variables, collapse = "|"), ")(\\[.*\\])?$")
    p_names <- grep(regex_pattern, p_names, value = TRUE)
  }
  if (length(p_names) == 0) stop("No valid parameters to plot.")
  
  show_trace <- "trace" %in% type
  show_dens  <- "density" %in% type
  
  total_plots <- length(p_names) * sum(c(show_trace, show_dens))
  
  if (is.null(ask)) {
    ask <- prod(par("mfcol")) < total_plots && dev.interactive()
  }
  
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  chain_cols <- palette()[1:n_chain]
  for (p in p_names) {
    p_data <- post[, , p, drop = FALSE] 
    if (show_trace) {
      plot(1:n_iter, type = "n", ylim = range(p_data, na.rm = TRUE), 
           xlab = "Iteration", ylab = "Parameter Value",
           main = paste("Trace of", p), ...)
      
      for (c in 1:n_chain) {
        lines(1:n_iter, p_data[, c, ], col = chain_cols[c])
      }
    }
    if (show_dens) {
      dens_list <- lapply(1:n_chain, function(c) density(p_data[, c, ]))
      x_lim <- range(sapply(dens_list, function(d) d$x))
      y_lim <- c(0, max(sapply(dens_list, function(d) d$y)))
      
      plot(1, type = "n", xlim = x_lim, ylim = y_lim,
           xlab = "Parameter Value", ylab = "Density",
           main = paste("Density of", p), ...)
      
      for (c in 1:n_chain) {
        lines(dens_list[[c]], col = chain_cols[c], lwd = 2)
      }
    }
  }
  
  invisible()
}

##' @title Posterior Predictive Checks for DRM Models
##'
##' @description Generates posterior predictive check plots comparing observed
##'   data to simulated data from the posterior predictive distribution.
##' 
##' @param x An object of class \code{drmrmodels}.
##' @param type Character string indicating the type of plot: \code{"density"}
##'   (default) or \code{"ecdf"}.
##' @param npost Integer specifying the number of posterior draws to plot.
##'   Default is 50. Usually, we do not use the total number of draws here
##'   because the plot tends to get too heavy.
##' @param transform Character string specifying the scale for plotting. Either
##'   \code{"none"} (default) or \code{"log1p"}.
##' @param ... Additional graphical parameters passed to \code{plot}.
##' 
##' @name ppc
##' @export
ppc <- function(x, ...) UseMethod("ppc", x)

##' @rdname ppc
##' @export
ppc.drmrmodels <- function(x, 
                           type = c("density", "ecdf"),
                           npost = 50, 
                           transform = c("none", "log1p"),
                           ...) {
  
  if (!inherits(x, "drmrmodels")) {
    stop("Use only with 'drmrmodels' objects.")
  }
  
  type <- match.arg(type)
  transform <- match.arg(transform)
  
  ids <- sort(c(x$data$id_nz, x$data$id_z))
  y <- x$data$y[ids]
  
  ft_vals <- fitted(x)
  y_pp_draws <- ft_vals$gq$draws(variables = sprintf("y_pp[%d]", ids),
                                 format = "draws_matrix")
  
  n_total <- NROW(y_pp_draws)
  if (npost > n_total) npost <- n_total
  y_pp_draws <- y_pp_draws[sample(seq_len(n_total), size = npost), , drop = FALSE]
  
  if (transform == "log1p") {
    y_trans <- log1p(y)
    y_pp_trans <- log1p(y_pp_draws)
    xlab_text <- sprintf("%s (log1p scale)", x$cols$y_col)
  } else {
    y_trans <- y
    y_pp_trans <- y_pp_draws
    xlab_text <- sprintf("%s", x$cols$y_col)
  }
  
  if (type == "density") {
    y_dens <- stats::density(y_trans)
    ft_dens <- apply(y_pp_trans, 1, stats::density)
    
    xlim <- range(c(y_dens$x, sapply(ft_dens, function(d) d$x)))
    ylim <- c(0, max(c(y_dens$y, sapply(ft_dens, function(d) d$y))))
    ylb <- "Density"
    
  } else {
    y_dens <- stats::ecdf(y_trans)
    ft_dens <- apply(y_pp_trans, 1, stats::ecdf)
    xlim <- range(c(y_trans, y_pp_trans))
    ylim <- c(0, 1)
    ylb <- "ECDF"
  }
  
  plot(0, type = "n", xlim = xlim, ylim = ylim,
       ylab = ylb, xlab = xlab_text,
       main = "Posterior Predictive Check", ...)
  
  draw_col <- grDevices::adjustcolor("steelblue", alpha.f = 0.3)
  
  if (type == "density") {
    lapply(ft_dens, function(d) lines(d, col = draw_col))
    lines(y_dens, lwd = 2, col = "black")
    
  } else {
    lapply(ft_dens, function(d) lines(d, col = draw_col, do.points = FALSE, verticals = TRUE))
    lines(y_dens, lwd = 2, col = "black", do.points = FALSE, verticals = TRUE)
  }
  
  legend("topright", 
         legend = c("Observed", "Post. pred. sample"), 
         col = c("black", "steelblue"), 
         lwd = c(2, 1), 
         bty = "n", 
         cex = 0.9)
  
  invisible()
}
