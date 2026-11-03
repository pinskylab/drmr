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
    bad <- setdiff(names(replacements), names(original))
    if (length(bad) > 0) {
      stop("Unknown name(s): ", paste(bad, collapse = ", "),
           ". Valid names: ", paste(names(original), collapse = ", "), ".",
           call. = FALSE)
    }
    for (re in c("ar_re", "iid_re", "sp_re")) {
      if (re %in% names(replacements)) {
        replacements[[re]] <- fix_re(replacements[[re]], re)
      }
    }
    if ("rec_dd" %in% names(replacements)) {
      replacements[["rec_dd"]] <-
        fix_dd(replacements[["rec_dd"]])
    }
    out <- utils::modifyList(original, replacements)
  } else {
    out <- original
  }
  return(out)
}

##' @title Toggle value to code
##' @description Converts a toggle given as a name, such as \code{"rec"}, to
##'   its code. A valid code passes unchanged. Any other value, including
##'   \code{TRUE} and \code{FALSE}, stops with an error, so that a misspelled
##'   toggle is not ignored.
##' @param x a \code{character} or \code{numeric} scalar.
##' @param codes a named \code{numeric} vector with the valid codes.
##' @param name the name of the toggle, for the error message.
##' @return A \code{numeric} scalar.
##' @author lcgodoy
toggle_code <- function(x, codes, name) {
  if (is.numeric(x) && length(x) == 1 && !is.na(x) && x %in% codes) {
    return(as.numeric(x))
  }
  if (is.character(x) && length(x) == 1 && x %in% names(codes)) {
    return(codes[[x]])
  }
  stop("'", name, "' must be one of ",
       paste0("\"", names(codes), "\"", collapse = ", "), " (or ",
       paste(codes, collapse = ", "), ").", call. = FALSE)
}

##' @title Random effects verbose to code
##' @param x a \code{character} or \code{numeric} scalar.
##' @param name the name of the toggle, for the error message.
##' @return A \code{numeric} scalar.
##' @author lcgodoy
fix_re <- function(x, name = "random effect toggle") {
  toggle_code(x, c(none = 0, rec = 1, surv = 2, dens = 3), name)
}

##' @title Density dependence verbose to code
##' @param x a \code{character} or \code{numeric} scalar.
##' @return A \code{numeric} scalar.
##' @author lcgodoy
fix_dd <- function(x) {
  toggle_code(x, c(ricker = 0, bh = 1, none = 2), "rec_dd")
}

##' @title Number of initial log-densities
##' @description With density dependence, the first time point needs the
##'   densities of ages 1 to \code{n_ages}. Without it, recruitment gives the
##'   density of age 1, so the first time point needs ages 2 to
##'   \code{n_ages} only.
##' @param n_ages number of age classes.
##' @param rec_dd 0 for Ricker, 1 for Beverton-Holt, and 2 for no density
##'   dependence.
##' @return A \code{numeric} scalar.
##' @author lcgodoy
n_init <- function(n_ages, rec_dd) {
  if (rec_dd < 2) n_ages else n_ages - 1
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

##' @title Site areas for the Stan models
##' @param dat the \code{data} element of an \code{adrm} object.
##' @return The \code{areas} element of \code{dat}. A fit from an earlier
##'   version of drmr has no \code{areas}, so the function returns 1 for each
##'   site.
get_areas <- function(dat) {
  if (!is.null(dat$areas)) {
    return(dat$areas)
  }
  return(array(rep(1, dat$n_sites), dim = dat$n_sites))
}

##' @title Number of non-zero entries of the movement matrix
##' @param dat a \code{list} with the elements \code{movement} and
##'   \code{adj_mat}, such as the output of [make_data()] or the \code{data}
##'   element of an \code{adrm} object.
##' @return The number of non-zero entries of the movement matrix that the
##'   Stan models build: one for each non-zero entry of \code{adj_mat} (two
##'   for each pair of neighbors in a symmetric matrix), plus one on the
##'   diagonal for each site with no neighbors. It is 0 without movement.
##'   Fits from earlier versions of drmr stored a smaller count, so the
##'   functions that reuse a fit call this function instead.
get_n_edges_adj <- function(dat) {
  if (!isTRUE(dat$movement == 1)) {
    return(0)
  }
  return(sum(dat$adj_mat != 0) + sum(rowSums(dat$adj_mat) == 0))
}

##' @title Check a generated quantities run
##' @param gq a \code{CmdStanGQ} object.
##' @return \code{gq}, invisibly. Stops if Stan failed in any chain.
check_gq <- function(gq) {
  if (any(gq$return_codes() != 0)) {
    stop("Stan could not compute the generated quantities. ",
         "See the Stan messages above.")
  }
  invisible(gq)
}
