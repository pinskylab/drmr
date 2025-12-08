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
