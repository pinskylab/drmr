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
