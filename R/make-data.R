##' This function creates the \code{list} used as the input for the \code{stan}
##' model.
##'
##' @title Make data for DRM stan models
##' @param y a \code{numeric vector} of species' densities.
##' @param time an \code{vector} indicating the time point associated to each
##'   element of \code{y}.
##' @param site an \code{vector} indicating the sites associated to each element
##'   of \code{y}.
##' @param init_data an optional vector (of lengh n_ages - 1) to initialize the
##'   population dynamics.
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
##' @param ages_movement a \code{integer} indicating the age at which
##'   individuals from the focal species start moving. In this case, individuals
##'   below this age threshold remain "static". Alternatively, we can input a
##'   \code{vector} of length \code{n_ages}. This vector will have a \code{0}
##'   for age-groups that cannot move and \code{1} for those age-groups that
##'   move. For example, the following vector \code{c(0, 0, 1, 1, 0)} indicates
##'   that ages 1, 2, and 5 are static, while ages 3 and 4 are allowed to move.
##' @param ages_movement An \code{integer} or a \code{numeric vector} specifying
##'   the ages at which individuals of the focal species are assumed to move. If
##'   \code{ages_movement} is an integer, individuals younger than this age are
##'   considered static (non-moving). If \code{ages_movement} is a numeric
##'   vector of length \code{n_ages}, it indicates movement capability for each
##'   age group. A value of \code{0} indicates the corresponding age group is
##'   static, while \code{1} indicates movement is allowed. For example,
##'   \code{c(0, 0, 1, 1, 0)} specifies that age groups 1, 2, and 5 are static,
##'   while 3 and 4 are mobile.
##' @param adj_mat an adjacency \code{matrix} of dimensions \code{sites}
##'   \eqn{\times} \code{sites}. Its elements are 1 if two sites are neighbors
##'   and zero otherwise.
##' @param .toggles a \code{list} of toggles for model components. The
##'   components are: \itemize{ \item \code{cloglog}: 1 to use the complementary
##'   log-log and 0 for the logit link function for the absence probabilities.
##'   \item \code{movement}: 1 to allow for (adjacent) moviment; 0 for static.
##'   \item \code{est_mort}: 1 to estimate mortality and 0 otherwise.  \item
##'   \code{est_init}: 1 to estimate initial values for lambda and 0 otherwise.
##'   \item \code{time_ar}: 1 to incorporate an AR(1) process for recruitment.
##'   \item \code{qr_t}: 1 to use QR parametrization for the absence probability
##'   regression coefficients and 0 otherwise.  \item \code{qr_r}: 1 to use QR
##'   parametrization for the recruitment regression coefficients and 0
##'   otherwise.  \item \code{qr_m}: 1 to use QR parametrization for the
##'   survival regression coefficients and 0 otherwise.  }
##' @param .priors a \code{list} of priors hyperparameters.
##' @param reorder a \code{boolean} telling whether the data needs to be
##'   reordered. The default is TRUE and means the data points will be ordered
##'   by site and time, respectively.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{ \item
##'   \code{"lognormal1"} (default): log-normal with the usual parametrization;
##'   \item \code{"lognormal2"}: log-normal parametrized in terms of its mean;
##'   \item \code{"gamma"}: gamma parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its mean.}
##' @return a \code{list} to be used as the input for a \code{stan} model
##' @author lcgodoy
##' @export
make_data <- function(y,
                      time,
                      site,
                      init_data = numeric(0),
                      f_mort,
                      m = 0.25,
                      x_t,
                      x_m,
                      x_r,
                      n_ages = 2,
                      age_selectivity,
                      ages_movement,
                      adj_mat = matrix(0, ncol = 1, nrow = 1),
                      .toggles,
                      .priors,
                      family = "lognormal1",
                      reorder = TRUE) {
  ## getting the default toggles and using user options
  stopifnot(length(family) == 1)
  stopifnot(family %in% c("lognormal1", "lognormal2",
                          "gamma", "loglogistic",
                          "truncnorm"))
  likelihood <- switch(family,
                       lognormal1  = 0,
                       lognormal2  = 1,
                       gamma       = 2,
                       loglogistic = 3,
                       truncnorm   = 4)
  toggles <- default_toggles() |>
    safe_modify(.toggles) |>
    c(list(likelihood = likelihood))
  ## getting the default priors and using user options
  priors <- default_priors() |>
    safe_modify(.priors)
  ## additional quantities that can be inferred from data
  n_patches <- length(unique(site))
  n_time <- length(unique(time))
  if (reorder) {
    my_ord <- order(site, time)
    y <- y[my_ord]
    site <- site[my_ord]
    time <- time[my_ord]
  }
  if (missing(f_mort))
    f_mort <- matrix(0, nrow = n_ages, ncol = n_time)
  if (missing(age_selectivity)) {
    selectivity_at_age <- rep(1, n_ages)
  } else {
    stopifnot(length(age_selectivity) == n_ages)
    selectivity_at_age <- age_selectivity
  }
  if (toggles$movement) {
    stopifnot(ncol(adj_mat) == nrow(adj_mat) &&
              nrow(adj_mat) == n_patches)
    if (missing(ages_movement)) {
      ages_movement <- rep(1L, n_ages)
    } else if (length(ages_movement) == 1) {
      aux_mov <- seq_len(n_ages)
      ages_movement <- ifelse(aux_mov >= ages_movement, 1, 0)
    } else {
      stopifnot(length(ages_movement) == n_ages)
    }
  } else if (missing(ages_movement)) {
    ages_movement <- integer(0)
  }
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
  if (!toggles$est_init) {
    stopifnot(length(init_data) == n_ages - 1)
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
  if (toggles$qr_m) {
    if (!toggles$est_mort || K_m == 1) {
      message("turning QR parametrization for mort off!")
      toggles$qr_m <- 0
    }
  }
  output <- list(N = n_time * n_patches,
                 n_ages = n_ages,
                 n_patches = n_patches,
                 n_time = n_time,
                 init_data = init_data,
                 y = y,
                 f = f_mort,
                 m = m,
                 X_t = x_t,
                 K_t = K_t,
                 X_m = x_m,
                 K_m = K_m,
                 X_r = x_r,
                 K_r = K_r,
                 adj_mat = adj_mat,
                 ages_movement = ages_movement,
                 selectivity_at_age = selectivity_at_age) |>
    c(toggles,
      priors)
  return(output)
}

##' This function creates the \code{list} used as the input for the \code{stan}
##' model.
##'
##' @title Make data for SDM stan models
##' @param y a \code{numeric vector} of species' densities.
##' @param time an \code{vector} indicating the time point associated to each
##'   element of \code{y}.
##' @param site an \code{vector} indicating the sites associated to each element
##'   of \code{y}.
##' @param z a design \code{matrix} of variables associated to the probability
##'   of absence at each site/time.
##' @param x a design \code{matrix} of variables associated to the non-zero
##'   densities.
##' @param .toggles a \code{list} of toggles for model components. The
##'   components are: \itemize{ \item \code{cloglog}: 1 to use the complementary
##'   log-log and 0 for the logit link function for the absence probabilities.
##'   \item \code{movement}: 1 to allow for (adjacent) moviment; 0 for static.
##'   \item \code{est_mort}: 1 to estimate mortality and 0 otherwise.  \item
##'   \code{time_ar}: 1 to incorporate an AR(1) process for recruitment.  \item
##'   \code{qr_t}: 1 to use QR parametrization for the absence probability
##'   regression coefficients and 0 otherwise.  \item \code{qr_r}: 1 to use QR
##'   parametrization for the recruitment regression coefficients and 0
##'   otherwise.  \item \code{qr_m}: 1 to use QR parametrization for the
##'   survival regression coefficients and 0 otherwise.  }
##' @param .priors a \code{list} of priors hyperparameters.
##' @param reorder a \code{boolean} telling whether the data needs to be
##'   reordered. The default is TRUE and means the data points will be ordered
##'   by site and time, respectively.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{ \item
##'   \code{"lognormal1"} (default): log-normal with the usual parametrization;
##'   \item \code{"lognormal2"}: log-normal parametrized in terms of its mean;
##'   \item \code{"gamma"}: gamma parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its mean.}
##' @return a \code{list} to be used as the input for a \code{stan} model
##' @author lcgodoy
##' @export
make_data_sdm <- function(y,
                          time,
                          site,
                          z,
                          x,
                          ## age_zero = FALSE
                          .toggles,
                          .priors,
                          family = "lognormal1",
                          reorder = TRUE) {
  ## getting the default toggles and using user options
  stopifnot(length(family) == 1)
  stopifnot(family %in% c("lognormal1", "lognormal2",
                          "gamma", "loglogistic"))
  likelihood <- switch(family,
                       lognormal1  = 0,
                       lognormal2  = 1,
                       gamma       = 2,
                       loglogistic = 3)
  toggles <- list(time_ar = 0,
                  cloglog = 0,
                  qr_z = 0,
                  qr_x = 0) |>
    safe_modify(.toggles) |>
    c(list(likelihood = likelihood))
  ## getting the default priors and using user options
  priors <- default_priors() |>
    safe_modify(.priors)
  ## additional quantities that can be inferred from data
  n_patches <- length(unique(site))
  n_time <- length(unique(time))
  if (reorder) {
    my_ord <- order(site, time)
    y <- y[my_ord]
    site <- site[my_ord]
    time <- time[my_ord]
    time <- time - min(time) + 1
  }
  if (missing(x)) {
    x <- matrix(1, nrow = n_time * n_patches)
    K_x <- 1
  } else {
    K_x <- ncol(x)
    if (reorder)
      x <- x[my_ord, , drop = FALSE]
    if (length(priors$pr_coef_r_mu) < K_x)
      priors$pr_coef_r_mu <- rep(0, K_x)
    if (length(priors$pr_coef_r_sd) < K_x)
      priors$pr_coef_r_sd <- rep(1, K_x)
  }
  if (missing(z)) {
    z <- matrix(1, nrow = n_time * n_patches)
    K_z <- 1
  } else {
    K_z <- ncol(z)
    if (reorder)
      z <- z[my_ord, , drop = FALSE]
    if (length(priors$pr_coef_t_mu) < K_z)
      priors$pr_coef_t_mu <- rep(0, K_z)
    if (length(priors$pr_coef_t_sd) < K_z)
      priors$pr_coef_t_sd <- rep(1, K_z)
  }
  if (toggles$qr_z) {
    if (K_z == 1) {
      message("turning QR parametrization for theta off!")
      toggles$qr_z <- 0
    }
  }
  if (toggles$qr_x) {
    if (K_x == 1) {
      message("turning QR parametrization for logrec off!")
      toggles$qr_x <- 0
    }
  }
  output <- list(N = n_time * n_patches,
                 n_patches = n_patches,
                 n_time = n_time,
                 time = time,
                 y = y,
                 X = x,
                 K_x = K_x,
                 Z = z,
                 K_z = K_z) |>
    c(toggles,
      priors)
  return(output)
}
