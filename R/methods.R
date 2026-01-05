##' Update and Re-fit a DRM Model
##'
##' @param object An object of class \code{adrm}.
##' @param ... Arguments to be updated in the new call (e.g., .data, formula_zero, iter_sampling).
##'
##' @return An updated \code{adrm} object or the unevaluated call.
##' @export
update.adrm <- function(object, ..., evaluate = TRUE) {
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

##' Print method for adrm objects
##'
##' @param x An object of class \code{adrm}.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.adrm <- function(x, ...) {
  cat("Age-structured Dynamic Range Model (ADRM)\n")
  cat("\nInference method:\n")
  method <- tryCatch(x$stanfit$metadata()$method, error = function(e) "Unknown") |>
    toupper()
  method <- ifelse(method == "SAMPLE", "NUTS", method)
  cat(paste0("  ", toupper(method), "\n"))

  cat("\nFormulas:\n")
  for (n in names(x$formulas)) {
    if (!is.null(x$formulas[[n]])) {
      # Formatting: "formula_zero" -> "Zero"
      clean_name <- ifelse(grepl("_zero", n), "Zeros",
                    ifelse(grepl("_rec", n), "Recruitment",
                           "Survival"))
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

##' Print method for sdm objects
##'
##' @param x An object of class \code{sdm}.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.sdm <- function(x, ...) {
  cat("Species Distribution Model (SDM)\n")
  cat("\nInference method:\n")
  method <- tryCatch(x$stanfit$metadata()$method, error = function(e) "Unknown") |>
    toupper()
  method <- ifelse(method == "SAMPLE", "NUTS", method)
  cat(paste0("  ", toupper(method), "\n"))

  cat("\nFormulas:\n")
  for (n in names(x$formulas)) {
    if (!is.null(x$formulas[[n]])) {
      # Formatting: "formula_zero" -> "Zero"
      clean_name <- ifelse(grepl("_zero", n), "Zeros",
                           "Non-zero response")
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

##' Summary method for adrm objects
##'
##' @param object An object of class \code{adrm}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments passed to the cmdstanr summary method.
##'
##' @export
summary.adrm <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  params_out <- get_fitted_pars(object$data, "drm")
  params_out <- params_out[!grepl("^(lambda|z_)", params_out)]
  method <- tryCatch(object$stanfit$metadata()$method, error = function(e) "nuts")  
  if (method == "optimize") {
    res <- object$stanfit$mle(variables = params_out)
    out <- list(
      estimates = res,
      method = "optimize",
      call = object$call
    )
  } else {
    summ <- object$stanfit$summary(variables = params_out, 
                                   "mean", "sd", "rhat", "ess_bulk",
                                   ~quantile(., probs = probs), 
                                   ...)
    out <- list(
      estimates = summ,
      method = method,
      call = object$call,
      probs = probs
    )
  }
  class(out) <- c("summary.adrm", class(out))
  return(out)
}

##' Summary method for sdm objects
##'
##' @param object An object of class \code{sdm}.
##' @param probs A numeric vector of quantiles to calculate. 
##'   Defaults to c(0.05, 0.5, 0.95).
##' @param ... Additional arguments passed to the cmdstanr summary method.
##'
##' @export
summary.sdm <- function(object, probs = c(0.05, 0.5, 0.95), ...) {
  params_out <- get_fitted_pars(object$data, "sdm")
  params_out <- params_out[!grepl("^(lambda|z_)", params_out)]
  method <- tryCatch(object$stanfit$metadata()$method, error = function(e) "nuts")  
  if (method == "optimize") {
    res <- object$stanfit$mle(variables = params_out)
    out <- list(
      estimates = res,
      method = "optimize",
      call = object$call
    )
  } else {
    summ <- object$stanfit$summary(variables = params_out, 
                                   "mean", "sd", "rhat", "ess_bulk",
                                   ~quantile(., probs = probs), 
                                   ...)
    out <- list(
      estimates = summ,
      method = method,
      call = object$call,
      probs = probs
    )
  }
  class(out) <- c("summary.adrm", class(out))
  return(out)
}

##' Print method for summary.adrm
##'
##' @param x An object of class \code{summary.adrm}.
##' @param digits Number of digits to print.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.summary.adrm <- function(x, digits = 3, ...) {
  cat("Summary of ADRM Fit\n")
  cat("\nEstimates:\n")
  print(x$estimates, digits = digits)
  if (x$method == "optimize") {
    cat("\n(Point estimates from maximum likelihood optimization)\n")
  }
  invisible(x)
}

##' Print method for summary.adrm
##'
##' @param x An object of class \code{summary.sdm}.
##' @param digits Number of digits to print.
##' @param ... Additional arguments (ignored).
##'
##' @export
print.summary.sdm <- function(x, digits = 3, ...) {
  cat("Summary of SDM Fit\n")
  cat("\nEstimates:\n")
  print(x$estimates, digits = digits)
  if (x$method == "optimize") {
    cat("\n(Point estimates from maximum likelihood optimization)\n")
  }
  invisible(x)
}
