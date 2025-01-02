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

##' This function creates the \code{list} used as the input for the \code{stan}
##' model.
##'
##' @title Make data for stan models
##' @param y a \code{numeric vector} of species' densities.
##' @param time an \code{vector} indicating the time point associated to each
##'   element of \code{y}.
##' @param site an \code{vector} indicating the sites associated to each element
##'   of \code{y}.
##' @param f_mort an optional \code{matrix} informing the instantaneous fishing
##'   mortality rates at each age (columns) and timepoint (rows).
##' @param m a \code{numeric} value corresponding to the instantaneous natural
##'   mortality rate.
##' @param x_t a design \code{matrix} of variables associated to the probability
##'   of absence at each site/time.
##' @param x_m a design \code{matrix} of variables associated to survival.
##' @param x_r a design \code{matrix} of variables associated to recruitment.
##' @param n_ages an \code{integer} indicating the number of ages for the
##'   underlying population dynamic model.
##' @param age_selectivity an \code{numeric vector} with \code{n_ages} elements,
##'   where each element indicates the selectivity of a respective age. All the
##'   elements of this vector must lie between 0 and 1.
##' @param age_at_maturity a \code{integer} indicating the age at which a
##'   species attains maturity. This is used for movement. That is, every
##'   individual with age greater or equal to \code{age_at_maturity} may move
##'   from one patch to another. Individuals below this age threshold remain
##'   "static".
##' @param adj_mat an adjacency \code{matrix} of dimensions \code{sites}
##'   \eqn{\times} \code{sites}. Its elements are 1 if two sites are neighbors
##'   and zero otherwise.
##' @param .toggles a \code{list} of toggles for model components. The
##'   components are:
##'   \itemize{
##'    \item \code{cloglog}: 1 to use the complementary log-log and 0 for the logit
##'      link function for the absence probabilities.
##'    \item \code{movement}: 1 to allow for (adjacent) moviment; 0 for static.
##'    \item \code{est_mort}: 1 to estimate mortality and 0 otherwise.
##'    \item \code{p_error}: 1 to incorporate an AR(1) process for recruitment.
##'    \item \code{qr_t}: 1 to use QR parametrization for the absence probability
##'     regression coefficients and 0 otherwise.
##'    \item \code{qr_r}: 1 to use QR parametrization for the recruitment
##'     regression coefficients and 0 otherwise.
##'    \item \code{qr_m}: 1 to use QR parametrization for the survival
##'     regression coefficients and 0 otherwise.
##'   }
##' @param .priors a \code{list} of priors hyperparameters.
##' @param reorder a \code{boolean} telling whether the data needs to be
##'   reordered. The default is TRUE and means the data points will be ordered
##'   by site and time, respectively.
##' @return a \code{list} to be used as the input for a \code{stan} model
##' @author lcgodoy
##' @export
make_data <- function(y,
                      time,
                      site,
                      f_mort,
                      m = 0.25,
                      ## x,
                      x_t,
                      x_m,
                      x_r,
                      n_ages,
                      age_selectivity,
                      age_at_maturity,
                      adj_mat = matrix(0, ncol = 1, nrow = 1),
                      ## age_zero = FALSE
                      .toggles,
                      .priors,
                      reorder = TRUE) {
  ## getting the default toggles and using user options
  toggles <- default_toggles() |>
    safe_modify(.toggles)
  ## getting the default "extra quantities" and using user options
  ## extra_qt <- default_qt() |>
  ##   safe_modify(.extra_qt)
  ## getting the default priors and using user options
  priors <- default_priors() |>
    safe_modify(.priors)
  ## additional quantities that can be inferred from data
  n_patches <- length(unique(site))
  ## patches <- seq_len(n_patches)
  n_time <- length(unique(time))
  ## ages <- seq_len(n_ages)
  ## if (age_zero)
  ##   ages <- ages - 1
  if (reorder) {
    my_ord <- order(site, time)
    y <- y[my_ord]
    site <- site[my_ord]
    time <- time[my_ord]
  }
  if (missing(f_mort))
    f_mort <- matrix(0, nrow = n_ages, ncol = n_time)
  if (missing(age_selectivity))
    selectivity_at_age <- rep(1, n_ages)
  if (toggles$movement) {
    stopifnot(ncol(adj_mat) == nrow(adj_mat) &&
              nrow(adj_mat) == n_patches)
    if (missing(age_at_maturity)) {
      age_at_maturity <- array(1L, dim = 1)
    } else {
      age_at_maturity <- array(age_at_maturity, dim = 1)
    }
  } else if (missing(age_at_maturity)) {
    age_at_maturity <- integer(0)
  }
  ## if (toggles$x_dep_movement) {
  ##   stopifnot(ncol(adj_m) == nrow(adj_m) &&
  ##             nrow(adj_m) == n_patches)
  ## } else {
  ##   x <- matrix(0, nrow = 1, ncol = 1)
  ## }
  if (!toggles$est_mort) {
    m <- array(m, dim = 1)
    x_m <- matrix(0, nrow = 1, ncol = 1)
    K_m <- integer(0)
  } else {
    stopifnot(ncol(x_m) >= 1)
    m <- numeric(0)
    if (reorder)
      x_m <- x_m[my_ord, , drop = FALSE]
    K_m <- array(NCOL(x_m), dim = 1)
    if (length(priors$pr_coef_m_mu) == 0)
      priors$pr_coef_m_mu <- rep(0, K_m)
    if (length(priors$pr_coef_m_sd) == 0)
      priors$pr_coef_m_sd <- rep(1, K_m)
  }
  if (missing(x_r)) {
    x_r <- matrix(1, nrow = n_time * n_patches)
    K_r <- 1
  } else {
    K_r <- ncol(x_r)
    if (reorder)
      x_r <- x_r[my_ord, , drop = FALSE]
    if (length(priors$pr_coef_r_mu) < K_r)
      priors$pr_coef_r_mu <- rep(0, K_r)
    if (length(priors$pr_coef_r_sd) < K_r)
      priors$pr_coef_r_sd <- rep(1, K_r)
  }
  if (missing(x_t)) {
    x_t <- matrix(1, nrow = n_time * n_patches)
    K_t <- 1
  } else {
    K_t <- ncol(x_t)
    if (reorder)
      x_t <- x_t[my_ord, , drop = FALSE]
    if (length(priors$pr_coef_t_mu) < K_t)
      priors$pr_coef_t_mu <- rep(0, K_t)
    if (length(priors$pr_coef_t_sd) < K_t)
      priors$pr_coef_t_sd <- rep(1, K_t)
  }
  if (toggles$qr_t) {
    if (K_t == 1) {
      message("turning QR parametrization for theta off!")
      toggles$qr_t <- 0
    }
  }
  if (toggles$qr_r) {
    if (K_r == 1) {
      message("turning QR parametrization for logrec off!")
      toggles$qr_r <- 0
    }
  }
  if (toggles$qr_r) {
    if (!toggles$est_mort || K_m == 1) {
      message("turning QR parametrization for mort off!")
      toggles$qr_m <- 0
    }
  }
  ## if (toggles$sr_rel) {
  ##   x_r <- matrix(0, nrow = 1, ncol = 1)
  ##   K_r <- integer(0)
  ## } else {
  ##   K_r <- array(NCOL(x_r), dim = 1)
  ## }
  output <- list(N = n_time * n_patches,
                 n_ages = n_ages, n_patches = n_patches,
                 n_time = n_time,
                 ## ages = ages,
                 ## patches = patches,
                 y = y,
                 ## area = area,
                 f = f_mort,
                 m = m,
                 X_t = x_t,
                 K_t = K_t,
                 X_m = x_m,
                 K_m = K_m,
                 X_r = x_r,
                 K_r = K_r,
                 adj_mat = adj_mat,
                 age_at_maturity = age_at_maturity,
                 selectivity_at_age = selectivity_at_age) |>
    c(toggles,
      priors)
  return(output)
}
