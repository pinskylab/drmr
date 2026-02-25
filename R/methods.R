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
##' @return An updated \code{adrm} object or the unevaluated call.
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
  cat(sprintf("  Patches: %d\n", max(x$data$patch)))
  if (!is.null(x$data$family)) {
     cat(sprintf("  Non-zero density family: %s\n", x$data$family))
  }
  invisible(x)
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
  params_out <- params_out[!grepl("^(lambda|z_)", params_out)]
  method <- tryCatch(object$stanfit$metadata()$method, error = function(e) "nuts")  
  if (method == "optimize") {
    res <- object$stanfit$mle(variables = params_out)
    out <- list(estimates = res,
                method = "optimize",
                call = object$call)
  } else {
    summ <- object$stanfit$summary(variables = params_out, 
                                   "mean", "sd", "rhat", "ess_bulk",
                                   ~posterior::quantile2(., probs = probs), 
                                   ...)
    out <- list(estimates = summ, method = method, call = object$call, probs = probs)
  }
  subclass <- ifelse(is_adrm, "summary.adrm", "summary.sdm")
  class(out) <- c(subclass, "summary.drmrmodels", class(out))
  return(out)
}

##' Print method for summary.spatial_model
##'
##' @param x An object of class \code{summary.spatial_model}.
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
    
    # Maintain your original variable filtering logic
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
