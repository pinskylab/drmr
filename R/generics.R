##' Validate an adrm object
##'
##' Checks if an object is a valid `adrm` object by ensuring it is a \code{list}
##' and contains all required fields.
##'
##' @param x An object to validate.
##'
##' @return Nothing if validation is successful. Throws an error on failure.
##' @keywords internal
validate_drm_m <- function(x) {
  stopifnot(inherits(x, "list"))
  my_fields <- c("stanfit", "data", "formulas", "cols", "seed")
  stopifnot(all(my_fields %in% names(x)))
}

##' Create an adrm object
##'
##' Constructs an object of class `adrm` from a list, after validating its
##' structure. Note that, this function is mostly for internal usage. `adrm`
##' stands for age-structured dynamic range model.
##'
##' @param x A `list` that has the required fields for an 'adrm' object.
##'   Defaults to an empty `list`.
##'
##' @return An object of class `adrm`.
##' @export
new_adrm <- function(x = list()) {
  validate_drm_m(x)
  new_class <- c("adrm", "drmrmodels", class(x))
  structure(x, class = new_class)
}

##' Create an sdm object
##'
##' Constructs an object of class `sdm` from a list, after validating its
##' structure. Note that, this function is mostly for internal usage. `sdm`
##' stands for species distribution model.
##'
##' @param x A `list` that has the required fields for an 'sdm' object.
##'   Defaults to an empty `list`.
##'
##' @return An object of class `sdm`.
##' @export
new_sdm <- function(x = list()) {
  validate_drm_m(x)
  new_class <- c("sdm", "drmrmodels", class(x))
  structure(x, class = new_class)
}

##' Validate a pred_drmr list object
##'
##' Checks if an object is a valid `pred_drmr` object by ensuring it is a
##' \code{list} and contains all required fields.
##'
##' @param x An object to validate.
##'
##' @return Nothing if validation is successful. Throws an error on failure.
##' @keywords internal
validate_pred_drmr <- function(x) {
  stopifnot(inherits(x, "list"))
  # You can adjust these fields based on what you decide to include in your predict output list
  my_fields <- c("gq", "spt") 
  stopifnot(all(my_fields %in% names(x)))
}

##' Create a pred_drmr object
##'
##' Constructs an object of class `pred_drmr` from a list, after validating its
##' structure. Note that this function is mostly for internal usage.
##'
##' @param x A `list` that has the required fields for a 'pred_drmr' object.
##'   Defaults to an empty `list`.
##'
##' @return An object of class `pred_drmr`.
##' @export
new_pred_drmr <- function(x = list()) {
  validate_pred_drmr(x)
  new_class <- c("pred_drmr", class(x))
  structure(x, class = new_class)
}
