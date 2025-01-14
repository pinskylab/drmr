##' @title Regression coefficient for non-centered variable
##' @description Consider a linear predictor having linear and square terms
##'   associated with a variable \eqn{x}. Assume this variable was centered
##'   before being included in the linear predictor. This functions recovers the
##'   regression coefficient associated to the linear term as if the variable
##'   was not centered.
##' @param beta1 A \code{numeric} regression coefficient associated to the
##'   linear term.
##' @param beta2 A \code{numeric} regression coefficient associated to the
##'   quadratic term.
##' @param offset a \code{numeric} representing the "center" of \eqn{x}.
##' @return a \code{numeric} representing the regression coefficient of the
##'   linear term for the model where \eqn{x} is not centered.
##' @author Lucas Godoy
fix_linbeta <- function(beta1, beta2, offset) {
  beta1 - 2 * offset * beta2
}

##' @title Value of a covariate that maximizes the response variable in a
##'   quadratic model.
##' @description Consider a linear predictor having linear and square terms
##'   associated with a variable \eqn{x}. Assume this variable was centered
##'   before being included in the linear predictor. This functions returns the
##'   value of \eqn{x} (on its original scale) such that the linear predictor is
##'   maximized (or minimized).
##' @param beta1 A \code{numeric} regression coefficient associated to the
##'   linear term.
##' @param beta2 A \code{numeric} regression coefficient associated to the
##'   quadratic term.
##' @param offset a \code{numeric} representing the "center" of \eqn{x}.
##' @return a \code{numeric} representing the uncentered \eqn{x} that maximizes
##'   (or minimizes) the linear predictor.
##' @author Lucas Godoy
##' @export
max_quad_x <- function(beta1, beta2, offset = 0) {
  stopifnot(NROW(beta1) == NROW(beta2))
  stopifnot(NROW(offset) == 1)
  out <- fix_linbeta(beta1, beta2, offset)
  - .5  * out / beta2
}

##' Generates an adjacency matrix for "movement"
##' @title Generates an adjacency matrix
##' @param x an \code{sf} object representing the patches.
##' @return an adjacency \code{matrix}
##' @author lcgodoy
##' @export
gen_adj <- function(x) {
  if (!inherits(x, "sfc")) {
    message("Assuming patches are 1d and ordered in such a way that the first and last patches have only one neighbor.")
    aux <-
      rbind(c(0, 0), c(1, 0),
            c(1, 1), c(0, 0)) |>
      list() |>
      sf::st_polygon() |>
      sf::st_sfc()
    out <- sf::st_make_grid(aux,
                            n = c(1, length(x)),
                            what = "polygons")
  } else
    out <- x
  spdep::poly2nb(out) |>
    spdep::nb2mat(style = "B") |>
    as.matrix()
}

##' Safely modifies a named \code{list}.
##'
##' This function returns an error if any of the names found in the
##' \code{replacements} object are not in the \code{list} to be modified.
##' @title Modifying a named list
##' @param original a named \code{list} with "original" parameters.
##' @param replacements a named \code{list} containing elements to be modified
##'   in the \code{original} \code{list}
##' @return an updated \code{original list}.
##' @author lcgodoy
safe_modify <- function(original, replacements) {
  if (!missing(replacements)) {
    stopifnot(all(names(replacements) %in% names(original)))
    out <- utils::modifyList(original, replacements)
  } else {
    out <- original
  }
  return(out)
}
