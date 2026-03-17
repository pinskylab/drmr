##' Safely modifies a named \code{list}.
##'
##' This function returns an error if any of the names found in the
##' \code{replacements} object are not in the \code{list} to be modified.
##' @title Modifying a named list
##' @param original a named \code{list} with "original" parameters.
##' @param replacements a named \code{list} containing elements to be modified
##'   in the \code{original} \code{list}.
##' @return The updated \code{original} list.
##' @author lcgodoy
safe_modify <- function(original, replacements) {
  if (!missing(replacements)) {
    stopifnot(all(names(replacements) %in% names(original)))
    if ("ar_re" %in% names(replacements)) {
      replacements[["ar_re"]] <-
        fix_re(replacements[["ar_re"]])
    }
    if ("iid_re" %in% names(replacements)) {
      replacements[["iid_re"]] <-
        fix_re(replacements[["iid_re"]])
    }
    if ("sp_re" %in% names(replacements)) {
      replacements[["sp_re"]] <-
        fix_re(replacements[["sp_re"]])
    }
    out <- utils::modifyList(original, replacements)
  } else {
    out <- original
  }
  return(out)
}

##' @title Random effects verbose to code
##' @param x a \code{character}
##' @return An \code{integer}.
##' @author lcgodoy
fix_re <- function(x) {
  if (!x %in% c(0:3)) {
    switch(x,
           none = 0,
           rec  = 1,
           surv = 2,
           dens = 3)
  } else x
}

##' @title Estimate phi
##' @inheritParams make_data
##' @return A \code{numeric} scalar representing an estimate for phi.
##' @author lcgodoy
get_phi_hat <- function(y, family) {
  if (family == "gamma") {
    ng0 <- sum(y > 0)
    xbar <- mean(y[y > 0])
    xbar2 <- xbar * xbar
    s2 <- stats::var(y[y > 0]) * (ng0 - 1) / ng0
    return(xbar2 * s2)
  } else if (family == "lognormal") {
    lxbar <- mean(log(y[y > 0]))
    ng0 <- sum(y > 0)
    ls2 <- stats::var(log(y[y > 0])) * (ng0 - 1) / ng0
    muhat <- exp(lxbar + 0.5 * ls2)
    return(muhat * muhat * (exp(ls2) - 1))
  } else {
    return(1)
  }
}

##' @title Zeros and nonzeros
##' @inheritParams make_data
##' @return A \code{list} containing necessary input for stan to identify zeros
##'   and nonzeros in the response variable.
##' @author lcgodoy
get_zeros <- function(y) {
  stopifnot(NCOL(y) == 1)
  stopifnot(all(y[!is.na(y)] >= 0))
  id_z  <- which(y == 0)
  id_nz <- which(y > 0)
  N_z  <- length(id_z)
  N_nz <- length(id_nz)
  return(list(id_z = id_z,
              id_nz = id_nz,
              N_z = N_z,
              N_nz = N_nz))
}
