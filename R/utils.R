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
##' @export
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

##' @title Check if x is between lb and ub
##'
##' @param x A numeric vector.
##' @param lb A numeric vector of lower bounds.
##' @param ub A numeric vector of upper bounds.
##'
##' @return No return value. This function stops execution if any of the
##'  following conditions are met:
##'  * lb is greater than or equal to ub.
##'  * The lengths of lb and ub are not equal.
##'  * The lengths of lb, ub, and x are not equal.
##'
check_between <- function(x, lb, ub) {
  stopifnot(all(lb <= ub))
  if (length(lb) > 1 | length(ub) > 1) {
    stopifnot(length(lb) == length(ub) &
              length(lb) == length(x))
  }
}

##' @inherit between
between_scalar <- function(x, lb, ub) {
  x >= lb & x <= ub
}

##' @title Check if elements in x are between corresponding elements in lb and
##'   ub
##'
##' @inheritParams check_between
##'
##' @return A logical vector of the same length as x, indicating whether each
##'   element of x is between the corresponding elements of lb and ub.
##'
##' @examples
##' between(1:5, 1, 5)
##' between(1:5, 2, 4)
##'
##' @export
between <- function(x, lb, ub) {
  check_between(x, lb, ub)
  if (length(lb) == 1) {
    sapply(x, between_scalar, lb, ub)
  } else {
    mapply(between_scalar, x, lb, ub)
  }
}

##' @title Calculate the interval score
##'
##' @description This function calculates the interval score for a given set of
##'   observations, lower and upper bounds, and alpha parameter.
##'
##' @details The interval score is a proper scoring rule that measures the
##'   accuracy of interval predictions. It takes into account both the coverage
##'   and the width of the prediction interval. A lower score indicates a better
##'   prediction.
##'
##' @param y A numeric vector of observations.
##' @param l A numeric vector of lower bounds for the prediction intervals.
##' @param u A numeric vector of upper bounds for the prediction intervals.
##' @param alpha A numeric value specifying the penalty parameter for interval
##'   width.
##'
##' @return A numeric vector of interval scores.
##'
##' @export
int_score <- function(y, l, u, alpha) {
  stopifnot(alpha > 0 & alpha < 1)
  check_between(y, l, u)
  alpha <- 2 / alpha
  ind_l <- as.numeric(y < l)
  ind_u <- as.numeric(y > u)
  (u - l) + alpha * ((l - y) * ind_l + (y - u) * ind_u)
}

##' @title Relationships with covariates
##' @details the \code{marg_rec} function works to obtain relationships
##'   regarding recruitment, while \code{marg_surv} and \code{marg_abs} evaluate
##'   relationships with survival and absence probability, respectively.
##' @param drm the output of a [fit_drm()] call.
##' @param newdata a \code{data.frame} with the values of the environmental
##'   variables at which we wish to estimate the recruitment, survival, or
##'   absence probabilities.
##' @return a \code{data.frame} with the samples the recruitment (or survival,
##'   or absence probabilities) posterior distribution.
##' @name marg
##' @export
##' @author lcgodoy
marg_rec <- function(drm, newdata) {
  stopifnot(inherits(drm$stanfit, "CmdStanFit"))
  my_formula <- drm$formulas$formula_rec
  new_x <- stats::model.matrix(my_formula, newdata)
  est_samples <-
    drm$stanfit$draws(variables = c("beta_r"),
                      format = "draws_matrix")
  n_sim <- NROW(est_samples)
  output <- vector(mode = "list", length = NROW(newdata))
  est_samples <-
    tcrossprod(est_samples, new_x)
  est_samples <- exp(c(est_samples))
  rownames(newdata) <- NULL
  for (i in seq_along(output)) {
    output[[i]] <-
      cbind.data.frame(newdata[i, , drop = FALSE], iter = seq_len(n_sim)) |>
      suppressWarnings()
  }
  output <- dplyr::bind_rows(output)
  output <- dplyr::mutate(output, recruitment = est_samples)
  return(output)
}

##' @rdname marg
##' @export
marg_pabs <- function(drm, newdata) {
  stopifnot(inherits(drm$stanfit, "CmdStanFit"))
  my_formula <- drm$formulas$formula_zero
  new_x <- stats::model.matrix(my_formula, newdata)
  est_samples <-
    drm$stanfit$draws(variables = c("beta_t"),
                      format = "draws_matrix")
  n_sim <- NROW(est_samples)
  output <- vector(mode = "list", length = NROW(newdata))
  est_samples <-
    tcrossprod(est_samples, new_x)
  est_samples <- stats::plogis(c(est_samples))
  rownames(newdata) <- NULL
  for (i in seq_along(output)) {
    output[[i]] <-
      cbind.data.frame(newdata[i, , drop = FALSE], iter = seq_len(n_sim)) |>
      suppressWarnings()
  }
  output <- dplyr::bind_rows(output)
  output <- dplyr::mutate(output, prob_abs = est_samples)
  return(output)
}

##' @rdname marg
##' @export
marg_surv <- function(drm, newdata) {
  stopifnot(inherits(drm$stanfit, "CmdStanFit"))
  stopifnot(!is.null(drm$data$K_m))
  my_formula <- drm$formulas$formula_surv
  new_x <- stats::model.matrix(my_formula, newdata)
  est_samples <-
    drm$stanfit$draws(variables = c("beta_s"),
                      format = "draws_matrix")
  n_sim <- NROW(est_samples)
  output <- vector(mode = "list", length = NROW(newdata))
  est_samples <-
    tcrossprod(est_samples, new_x)
  est_samples <- exp(c(est_samples))
  rownames(newdata) <- NULL
  for (i in seq_along(output)) {
    output[[i]] <-
      cbind.data.frame(newdata[i, , drop = FALSE], iter = seq_len(n_sim)) |>
      suppressWarnings()
  }
  output <- dplyr::bind_rows(output)
  output <- dplyr::mutate(output, survival = est_samples)
  return(output)
}
