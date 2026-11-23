##' Safely modifies a named \code{list}.
##'
##' This function returns an error if any of the names found in the
##' \code{replacements} object are not in the \code{list} to be modified.
##' @title Modifying a named list
##' @param original a named \code{list} with "original" parameters.
##' @param replacements a named \code{list} containing elements to be modified
##'   in the \code{original} \code{list}.
##' @param re_codes a named \code{numeric} vector with the valid codes of the
##'   random effect toggles (\code{ar_re}, \code{iid_re}, and \code{sp_re}).
##' @return The updated \code{original} list.
##' @author Anonymous
safe_modify <- function(original, replacements,
                        re_codes = c(none = 0, rec = 1, surv = 2, dens = 3)) {
  if (!missing(replacements)) {
    bad <- setdiff(names(replacements), names(original))
    if (length(bad) > 0) {
      stop("Unknown name(s): ", paste(bad, collapse = ", "),
           ". Valid names: ", paste(names(original), collapse = ", "), ".",
           call. = FALSE)
    }
    for (re in c("ar_re", "iid_re", "sp_re")) {
      if (re %in% names(replacements)) {
        replacements[[re]] <- fix_re(replacements[[re]], re, re_codes)
      }
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
##' @author Anonymous
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
##' @inheritParams toggle_code
##' @return A \code{numeric} scalar.
##' @author Anonymous
fix_re <- function(x, name = "random effect toggle",
                   codes = c(none = 0, rec = 1, surv = 2, dens = 3)) {
  toggle_code(x, codes, name)
}

##' @title Match prior hyperparameters to the coefficients
##' @description Gives one hyperparameter for each column of a design
##'   matrix. A single value applies to all the columns.
##' @param x a \code{numeric} vector of hyperparameters.
##' @param k the number of columns of the design matrix.
##' @param default the value to use when \code{x} is empty.
##' @param name the name of the prior, for the error message.
##' @return A \code{numeric} vector with \code{k} elements. When \code{k}
##'   is 1, it is an \code{array} of size 1, as Stan needs.
##' @author Anonymous
fix_beta_prior <- function(x, k, default, name) {
  if (length(x) == 0)
    x <- default
  if (length(x) == 1)
    x <- rep(x, k)
  if (length(x) != k) {
    stop("'", name, "' must have 1 or ", k, " values, one for each ",
         "column of the design matrix.", call. = FALSE)
  }
  if (k == 1)
    x <- array(x, dim = 1)
  return(x)
}

##' @title Zeros and nonzeros
##' @inheritParams make_data
##' @return A \code{list} containing necessary input for stan to identify zeros
##'   and nonzeros in the response variable.
##' @author Anonymous
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
