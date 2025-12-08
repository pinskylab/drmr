##' Validate an adrm object
##'
##' Checks if an object is a valid `adrm` object by ensuring it is a \code{list}
##' and contains all required fields.
##'
##' @param x An object to validate.
##'
##' @return Nothing if validation is successful. Throws an error on failure.
##' @keywords internal
validate_adrm <- function(x) {
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
  validate_adrm(x)
  new_class <- c("adrm", class(x))
  structure(x, class = new_class)
}

##' Validate an sdm object
##'
##' Checks if an object is a valid `sdm` object by ensuring it is a \code{list}
##' and contains all required fields.
##'
##' @param x An object to validate.
##'
##' @return Nothing if validation is successful. Throws an error on failure.
##' @keywords internal
validate_sdm <- function(x) {
  stopifnot(inherits(x, "list"))
  my_fields <- c("stanfit", "data", "formulas", "cols", "seed")
  stopifnot(all(my_fields %in% names(x)))
}

##' Create an sdm object
##'
##' Constructs an object of class `sdm` from a list, after validating its
##' structure. Note that, this function is mostly for internal usage. `sdm`
##' stands for age-structured dynamic range model.
##'
##' @param x A `list` that has the required fields for an 'sdm' object.
##'   Defaults to an empty `list`.
##'
##' @return An object of class `sdm`.
##' @export
new_sdm <- function(x = list()) {
  validate_sdm(x)
  new_class <- c("sdm", class(x))
  structure(x, class = new_class)
}
