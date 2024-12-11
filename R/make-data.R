## TODO: better error messages

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

##' Returns default quantities for the DRM model.
##'
##' Based on Alexa's Summer Flounder project.
##' @title Default quantities
##' @return a `list` with (explain values)
##' @author Lucas Godoy
default_qt <- function() {
  list(k = 0.14,
       l_inf = 83.6,
       t0 = -.2,
       sel_100 = 3,
       min_age = 0,
       max_age = 15,
       length_50_sel_guess = 20,
       age_sel = 0,
       h = 0.8)
}

##' Returns default priors' hyperparameters for the DRM model.
##'
##' Based on Alexa's Summer Flounder project.
##' @title Default priors' hyperparameters
##' @return a `list` with (explain values)
##' @author Lucas Godoy
default_priors <- function() {
  list(pr_sigma_obs_mu = 0,
       pr_sigma_obs_sd = 1,
       pr_phi_a = 2,
       pr_phi_b = 1,
       pr_sd_r_u = .1,
       pr_sd_r_alpha = .05,
       pr_rho_u = .9,
       pr_rho_alpha = .1,
       pr_coef_t_mu = 0,
       pr_coef_t_sd = 1,
       pr_coef_m_mu = numeric(0),
       pr_coef_m_sd = numeric(0),
       pr_coef_r_mu = 0,
       pr_coef_r_sd = 1)
}

##' Returns default toggles for the DRM model.
##'
##' Based on Alexa's Summer Flounder project.
##' @title Default toggles
##' @return a `list` with (explain values)
##' @author Lucas Godoy
default_toggles <- function() {
  list(cloglog = 0,
       movement = 0,
       est_mort = 0,
       ## sr_rel = 1,
       ## x_dep_movement = 0,
       ## exp_yn = 0,
       p_error = 1,
       qr_t = 0,
       qr_r = 0,
       qr_m = 0,
       likelihood = 0)
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
##' @author Lucas Godoy
safe_modify <- function(original, replacements) {
  if (!missing(replacements)) {
    stopifnot(all(names(replacements) %in% names(original)))
    out <- modifyList(original, replacements)
  } else {
    out <- original
  }
  return(out)
}

## extra quantities as a list?
make_data <- function(y, ## area,
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
                      M = 30,
                      ## age_zero = FALSE
                      ## wt_at_age,
                      ## dcap = 1 / 3,
                      .toggles,
                      ## .extra_qt,
                      .priors) {
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
  n_patches <- NCOL(y)
  ## patches <- seq_len(n_patches)
  n_time <- NROW(y)
  ## ages <- seq_len(n_ages)
  ## if (age_zero)
  ##   ages <- ages - 1
  if (missing(f_mort))
    f_mort <- matrix(0, nrow = n_time, ncol = n_ages)
  if (missing(age_selectivity))
    selectivity_at_age = rep(1, n_ages)
  if (toggles$movement) {
    stopifnot(ncol(adj_mat) == nrow(adj_mat) &&
              nrow(adj_mat) == n_patches)
    if (missing(age_at_maturity)) {
      age_at_maturity <- array(1L, dim = 1)
    } else {
      age_at_maturity <- array(age_at_maturity, dim = 1)
    }
  } else {
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
                 ## x = x,
                 X_t = x_t,
                 K_t = K_t,
                 X_m = x_m,
                 K_m = K_m,
                 X_r = x_r,
                 K_r = K_r,
                 M = M,
                 ## wt_at_age = wt_at_age,
                 ## dcap = dcap,
                 adj_mat = adj_mat,
                 age_at_maturity = age_at_maturity,
                 selectivity_at_age = selectivity_at_age) |>
    c(## extra_qt,
        toggles,
        priors)
  return(output)
}

## extra quantities as a list?
make_data_vec <- function(y,
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
                          M = 30,
                          ## age_zero = FALSE
                          ## wt_at_age,
                          ## dcap = 1 / 3,
                          .toggles,
                          ## .extra_qt,
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
    f_mort <- matrix(0, nrow = n_time, ncol = n_ages)
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
                 ## x = x,
                 X_t = x_t,
                 K_t = K_t,
                 X_m = x_m,
                 K_m = K_m,
                 X_r = x_r,
                 K_r = K_r,
                 M = M,
                 ## wt_at_age = wt_at_age,
                 ## dcap = dcap,
                 adj_mat = adj_mat,
                 age_at_maturity = age_at_maturity,
                 selectivity_at_age = selectivity_at_age) |>
    c(toggles,
      ## extra_qt,
      priors)
  return(output)
}
