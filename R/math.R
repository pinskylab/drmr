##' @title Regression coefficient for non-centered variable
##' @description Consider a linear predictor having linear and square terms
##'   associated with a variable \eqn{x}. Assume this variable was centered
##'   before being included in the linear predictor. This function recovers the
##'   regression coefficient associated with the linear term as if the variable
##'   was not centered.
##' @param beta1 A \code{numeric} regression coefficient associated with the
##'   linear term.
##' @param beta2 A \code{numeric} regression coefficient associated with the
##'   quadratic term.
##' @param offset a \code{numeric} representing the "center" of \eqn{x}.
##' @return A \code{numeric} value representing the regression coefficient of the
##'   linear term for the model where \eqn{x} is not centered.
##' @author Lucas Godoy
##' @export
fix_linbeta <- function(beta1, beta2, offset) {
  beta1 - 2 * offset * beta2
}

##' @title Value of a covariate that maximizes the response variable in a
##'   quadratic model.
##' @description Consider a linear predictor having linear and square terms
##'   associated with a variable \eqn{x}. Assume this variable was centered
##'   before being included in the linear predictor. This function returns the
##'   value of \eqn{x} (on its original scale) such that the linear predictor is
##'   maximized (or minimized).
##' @param beta1 A \code{numeric} regression coefficient associated with the
##'   linear term.
##' @param beta2 A \code{numeric} regression coefficient associated with the
##'   quadratic term.
##' @param offset a \code{numeric} representing the "center" of \eqn{x}.
##' @return A \code{numeric} value representing the uncentered \eqn{x} that maximizes
##'   (or minimizes) the linear predictor.
##' @author Lucas Godoy
##' @export
max_quad_x <- function(beta1, beta2, offset = 0) {
  stopifnot(NROW(beta1) == NROW(beta2))
  stopifnot(NROW(offset) == 1)
  out <- fix_linbeta(beta1, beta2, offset)
  - .5  * out / beta2
}

##' @title Generalized Inverse
##'
##' @details this function is taken from the package \code{MASS}
##' 
##' @param X A matrix we wish to invert.
##' @param tol A relative tolerance to detect zero singular values.
##' @return The generalized inverse of \code{X}.
ginv <- function (X, tol = sqrt(.Machine$double.eps)) {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}
