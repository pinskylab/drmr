##' This function creates the \code{list} used as the input for the \code{stan}
##' model.
##'
##' @title Make data for DRM stan models
##' @details The function sorts the data by site and then by time.
##' @param y a \code{numeric vector} of species' densities.
##' @param time a \code{vector} indicating the time point associated to each
##'   element of \code{y}.
##' @param site a \code{vector} indicating the sites associated to each element
##'   of \code{y}.
##' @param init_data an optional vector with the log-densities at the first
##'   time point. Without density dependence, give ages 2 to \code{n_ages}
##'   (length \code{n_ages - 1}). With density dependence, give ages 1 to
##'   \code{n_ages} (length \code{n_ages}).
##' @param f_mort an optional \code{matrix} informing the instantaneous fishing
##'   mortality rates at each age (rows) and time point (columns).
##' @param m a \code{numeric} value corresponding to the instantaneous natural
##'   mortality rate. The default value for this is \code{-log(.7)}, as it
##'   implies a survival rate of 0.70 between age classes.
##' @param x_t a design \code{matrix} of variables associated to the probability
##'   of absence at each site/time.
##' @param x_m a design \code{matrix} of variables associated to survival.
##' @param x_r a design \code{matrix} of variables associated to recruitment.
##' @param n_ages an \code{integer} indicating the number of ages for the
##'   underlying population dynamic model.
##' @param age_selectivity a \code{numeric vector} with \code{n_ages} elements,
##'   where each element indicates the selectivity of a respective age. All the
##'   elements of this vector must lie between 0 and 1.
##' @param ages_movement An \code{integer} or a \code{numeric vector} specifying
##'   the ages at which individuals of the focal species are assumed to move. If
##'   \code{ages_movement} is an integer, individuals younger than this age are
##'   considered static (non-moving). If \code{ages_movement} is a numeric
##'   vector of length \code{n_ages}, it indicates movement capability for each
##'   age group. A value of \code{0} indicates the corresponding age group is
##'   static, while \code{1} indicates movement is allowed. For example,
##'   \code{c(0, 0, 1, 1, 0)} specifies that age groups 1, 2, and 5 are static,
##'   while 3 and 4 are mobile. If you omit \code{ages_movement}, age 1
##'   (recruits) is static and all other age groups move. This default is the
##'   same as \code{ages_movement = 2}.
##' @param adj_mat an adjacency \code{matrix} of dimensions \code{sites}
##'   \eqn{\times} \code{sites}. Its elements are 1 if two sites are neighbors
##'   and zero otherwise. Its rows and columns follow the order of
##'   \code{levels(factor(site))}. Movement uses the values: entry
##'   \eqn{(i, j)} is the relative rate at which individuals at site \eqn{i}
##'   move to site \eqn{j}, and each row is scaled to sum to 1. The ICAR
##'   random effects (\code{sp_re}) use only which entries are positive, in
##'   either direction. If the sites form more than one connected group, each
##'   group gets its own sum-to-zero constraint and scaling factor, and each
##'   site with no neighbors gets an independent effect (Freni-Sterrantino et
##'   al., 2018). See [icar_data()].
##' @param areas an optional \code{numeric vector} with the area of each site.
##'   All values must be positive. \code{NULL} (the default) means an area of
##'   1 for each site, with a warning when the model uses the areas. Two
##'   parts of the model use \code{areas}:
##'   \itemize{
##'   \item Movement converts density to abundance (density times area),
##'   moves the abundance, and converts it back to density. So movement does
##'   not change the total abundance. Individuals at a site with no neighbors
##'   stay at that site.
##'   \item Density dependence with \code{acc_dd = 1} uses the mean SSB
##'   density of the region, weighted by area.
##'   }
##'   Only the ratios between areas matter, so any unit works. We recommend a
##'   named vector: the function matches the names to the site names and
##'   ignores extra names. An unnamed vector must follow the order of
##'   \code{levels(factor(site))}, which is also the order of the rows of
##'   \code{adj_mat}.
##' @param amat a \code{numeric vector} of length \code{n_ages} specifying the
##'   maturity-at-age.
##' @param weight a \code{numeric vector} of length \code{n_ages} specifying the
##'   weight-at-age.
##' @param .toggles a \code{list} of toggles for model components. The
##'   components are: \itemize{ \item \code{rho_mu}: 1 to explicitly relate rho
##'   to mu and 0 otherwise.  \item \code{cloglog}: 1 to use the complementary
##'   log-log and 0 for the logit link function for the absence probabilities.
##'   \item \code{movement}: 1 to allow for (adjacent) movement; 0 for static.
##'   \item \code{est_surv}: 1 to estimate mortality and 0 otherwise.  \item
##'   \code{est_init}: 1 to estimate initial values for lambda and 0 otherwise.
##'   \item \code{minit}: 1 to use mortality to estimate initial age classes and
##'   0 otherwise. With density-dependent recruitment, 1 starts from the
##'   equilibrium age structure at the first time point, and 2 starts from
##'   the equilibrium age structure times a scale \eqn{d} that the model
##'   estimates, with \eqn{\log d \sim N(}\code{pr_linit_scale_mu},
##'   \code{pr_linit_scale_sd}\eqn{)}. The value 2 needs density-dependent
##'   recruitment. \item
##'   \code{rec_dd}: "ricker" for Ricker, "bh" for
##'   Beverton-Holt, and "none" (default) for no density dependence. \item
##'   \code{acc_dd}: 0 for recruitment to depend on the SSB of the patch
##'   (default), and 1 for recruitment to depend on the "global" (i.e., the
##'   whole study region) SSB. The global SSB is the mean SSB density of the
##'   patches, weighted by \code{areas}. \item \code{ar_re}: a
##'   \code{character}. It
##'   assumes one of the following values: "none" no AR, "rec" AR(1) for
##'   recruitment, "surv" AR(1) for survival (only works when \code{est_surv} is
##'   on), "dens" AR(1) for density. \item \code{iid_re}: a \code{character}. It
##'   assumes one of the following values: "none" no IID random effect, "rec"
##'   IID random effect for recruitment, "surv" IID random effect for survival
##'   (only works when \code{est_surv} is on), "dens" IID random effect for
##'   density.  \item \code{sp_re}: a \code{character}. It assumes one of the
##'   following values: "none" no ICAR random effect, "rec" ICAR random effect
##'   for recruitment, "surv" ICAR random effect for survival (only works when
##'   \code{est_surv} is on), "dens" ICAR random effect for density.}
##' @param .priors a \code{list} of priors hyperparameters.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{ \item
##'   \code{"gamma"} (default): gamma parametrized in terms of its mean; \item
##'   \code{"lognormal"}: log-normal parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its median
##'   (usual parametrization), with shape \code{phi} > 1 so its mean is finite;
##'   \item \code{"lognormal_legacy"}: log-normal with
##'   its usual parametrization; }
##' @param phi_hat a \code{boolean} indicating whether the prior on \code{phi}
##'   should be determined through the data (using [get_phi_hat()]).
##' @return A \code{list} containing the data and settings to be used as input
##'   for the Stan model.
##' @author lcgodoy
##' @export
make_data <- function(y,
                      time,
                      site,
                      init_data = numeric(0),
                      f_mort,
                      m = - log(.7),
                      x_t,
                      x_m,
                      x_r,
                      n_ages = 2,
                      age_selectivity,
                      ages_movement,
                      adj_mat = matrix(0, ncol = 1, nrow = 1),
                      areas = NULL,
                      amat = NULL,
                      weight = NULL,
                      .toggles,
                      .priors,
                      family = "gamma",
                      phi_hat = FALSE) {
  ## getting the default toggles and using user options
  stopifnot(all(!is.na(c(time))))
  stopifnot(all(!is.na(c(site))))
  stopifnot(length(family) == 1)
  stopifnot(family %in% c("lognormal_legacy", "lognormal",
                          "gamma", "loglogistic", "truncnorm"))
  likelihood <- switch(family,
                       lognormal_legacy  = 0,
                       lognormal         = 1,
                       gamma             = 2,
                       loglogistic       = 3,
                       truncnorm         = 4)
  toggles <- default_toggles() |>
    safe_modify(.toggles) |>
    c(list(likelihood = likelihood))
  if (any(c(toggles$ar_re, toggles$iid_re, toggles$sp_re) == 2) &&
        toggles$est_surv == 0) {
    stop("Random effects on survival (\"surv\") need ",
         "'.toggles = list(est_surv = 1)'.", call. = FALSE)
  }
  if (!toggles$minit %in% 0:2) {
    stop("'minit' must be 0, 1, or 2.", call. = FALSE)
  }
  if (toggles$minit == 2 && toggles$rec_dd == 2) {
    stop("'minit = 2' needs density-dependent recruitment. ",
         "Set 'rec_dd' to \"ricker\" or \"bh\".", call. = FALSE)
  }
  ## getting the default priors and using user options
  if (phi_hat) {
    phi_hat <- get_phi_hat(y, family)
    if (missing(.priors)) {
      .priors <- list(pr_phi_a = 2, pr_phi_b = 2 / phi_hat)
    } else if (!"pr_phi_b" %in% names(.priors)) {
      if (!"pr_phi_a" %in% names(.priors)) {
        .priors <- c(.priors, list(pr_phi_a = 2, pr_phi_b = 2 / phi_hat))
      } else {
        .priors <- c(.priors, list(pr_phi_b = .priors$pr_phi_a / phi_hat))
      }
    }
  }
  priors <- default_priors() |>
    safe_modify(.priors)
  ## additional quantities that can be inferred from data
  n_sites <- length(unique(site))
  n_time <- length(unique(time))
  site_factor <- factor(site)
  site_levels <- levels(site_factor)
  site_stan <- as.integer(site_factor)
  time_init <- min(time)
  my_ord <- order(site, time)
  y <- y[my_ord]
  site <- site[my_ord]
  site_stan <- site_stan[my_ord]
  time <- time[my_ord] - time_init + 1
  if (missing(f_mort))
    f_mort <- matrix(0, nrow = n_ages, ncol = n_time)
  if (missing(age_selectivity)) {
    selectivity_at_age <- rep(1, n_ages)
  } else {
    stopifnot(length(age_selectivity) == n_ages)
    selectivity_at_age <- age_selectivity
  }
  if (is.null(amat))
    amat <- rep(1, n_ages)
  if (is.null(weight))
    weight <- rep(1, n_ages)
  if (toggles$movement || toggles$sp_re > 0) {
    check_adj_mat(adj_mat, n_sites)
  }
  if (toggles$movement) {
    if (missing(ages_movement)) {
      ## recruits (age 1) are static by default
      ages_movement <- c(0L, rep(1L, n_ages - 1))
    } else if (length(ages_movement) == 1) {
      aux_mov <- seq_len(n_ages)
      ages_movement <- ifelse(aux_mov >= ages_movement, 1, 0)
    } else {
      stopifnot(length(ages_movement) == n_ages)
    }
  } else if (missing(ages_movement)) {
    ages_movement <- integer(0)
  }
  if (is.null(areas)) {
    areas <- rep(1, n_sites)
    uses_areas <- toggles$movement == 1 ||
      (toggles$rec_dd < 2 && toggles$acc_dd == 1)
    if (uses_areas && n_sites > 1) {
      warning("'areas' is NULL, so all sites have an area of 1. ",
              "Movement and regional density dependence use the areas.",
              call. = FALSE)
    }
  } else {
    if (!is.numeric(areas) || length(dim(areas)) > 1) {
      stop("'areas' must be a numeric vector.")
    }
    if (!is.null(names(areas))) {
      if (anyDuplicated(names(areas)) > 0 ||
          !all(site_levels %in% names(areas))) {
        stop("The names of 'areas' must be unique and include all the ",
             "site names.")
      }
      areas <- areas[site_levels]
    } else if (length(areas) != n_sites) {
      stop("An unnamed 'areas' must have one value per site.")
    }
    if (any(!is.finite(areas)) || any(areas <= 0)) {
      stop("'areas' must contain finite, positive values only.")
    }
    areas <- as.numeric(areas)
  }
  if (n_sites == 1) {
    areas <- array(areas, dim = 1)
  }
  if (!toggles$est_surv) {
    m <- array(m, dim = 1)
    x_m <- matrix(0, nrow = 1, ncol = 1)
    K_m <- integer(0)
  } else {
    stopifnot(ncol(x_m) >= 1)
    m <- numeric(0)
    x_m <- x_m[my_ord, , drop = FALSE]
    K_m <- array(NCOL(x_m), dim = 1)
    if (length(priors$pr_beta_s_mu) < K_m)
      priors$pr_beta_s_mu <- rep(0, K_m)
    if (length(priors$pr_beta_s_sd) < K_m )
      priors$pr_beta_s_sd <- rep(1, K_m)
    if (K_m == 1) {
      priors$pr_beta_s_mu <- array(priors$pr_beta_s_mu,
                                   dim = 1)
      priors$pr_beta_s_sd <- array(priors$pr_beta_s_sd,
                                   dim = 1)
    }
  }
  if (missing(x_r)) {
    x_r <- matrix(1, nrow = n_time * n_sites)
    K_r <- 1
  } else {
    K_r <- ncol(x_r)
    x_r <- x_r[my_ord, , drop = FALSE]
    if (length(priors$pr_beta_r_mu) < K_r)
      priors$pr_beta_r_mu <- rep(0, K_r)
    if (length(priors$pr_beta_r_sd) < K_r)
      priors$pr_beta_r_sd <- rep(1, K_r)
    if (K_r == 1) {
      priors$pr_beta_r_mu <- array(priors$pr_beta_r_mu,
                                   dim = 1)
      priors$pr_beta_r_sd <- array(priors$pr_beta_r_sd,
                                   dim = 1)
    }
  }
  if (missing(x_t)) {
    x_t <- matrix(1, nrow = n_time * n_sites)
    K_t <- 1
  } else {
    K_t <- ncol(x_t)
    x_t <- x_t[my_ord, , drop = FALSE]
    if (length(priors$pr_beta_t_mu) < K_t)
      priors$pr_beta_t_mu <- rep(0, K_t)
    if (length(priors$pr_beta_t_sd) < K_t)
      priors$pr_beta_t_sd <- rep(1, K_t)
    if (K_t == 1) {
      priors$pr_beta_t_mu <- array(priors$pr_beta_t_mu,
                                   dim = 1)
      priors$pr_beta_t_sd <- array(priors$pr_beta_t_sd,
                                   dim = 1)
    }
  }
  if (!toggles$est_init) {
    len_init <- n_init(n_ages, toggles$rec_dd)
    if (length(init_data) > 0) {
      stopifnot(length(init_data) == len_init)
    } else {
      init_data <- log(seq(from = 1, to = .01,
                           length.out = len_init))
      if (length(init_data) == 1)
        init_data <- array(init_data, dim = 1)
    }
  }
  if (toggles$sp_re > 0) {
    check_adj_mat(adj_mat, n_sites)
    icar <- icar_data(adj_mat)
  } else {
    icar <- icar_data_none()
  }
  if (!toggles$movement) {
    adj_mat <- matrix(0, ncol = 1, nrow = 1)
  }
  n_edges_adj <- get_n_edges_adj(list(movement = toggles$movement,
                                      adj_mat = adj_mat))
  zeros <- get_zeros(y)
  if (any(is.na(y)))
    y <- ifelse(is.na(y), -99, y)
  output <- list(N = n_time * n_sites,
                 n_ages = n_ages,
                 n_sites = n_sites,
                 n_time = n_time,
                 time = time,
                 time_init = time_init,
                 site = site_stan,
                 site_levels = site_levels,
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
                 n_edges_adj = n_edges_adj,
                 areas = areas,
                 ages_movement = ages_movement,
                 selectivity_at_age = selectivity_at_age,
                 amat = amat,
                 weight = weight) |>
    c(icar,
      zeros,
      toggles,
      priors)
  return(output)
}

##' This function creates the \code{list} used as the input for the \code{stan}
##' model.
##'
##' @title Make data for SDM stan models
##' @details The function sorts the data by site and then by time.
##' @param y a \code{numeric vector} of species' densities.
##' @param time a \code{vector} indicating the time point associated to each
##'   element of \code{y}.
##' @param site a \code{vector} indicating the sites associated to each element
##'   of \code{y}.
##' @param z a design \code{matrix} of variables associated to the probability
##'   of absence at each site/time.
##' @param x a design \code{matrix} of variables associated to the non-zero
##'   densities.
##' @param adj_mat an adjacency \code{matrix} of dimensions \code{sites}
##'   \eqn{\times} \code{sites}. Its elements are 1 if two sites are neighbors
##'   and zero otherwise. Its rows and columns follow the order of
##'   \code{levels(factor(site))}. The ICAR random effects (\code{sp_re}) use
##'   only which entries are positive, in either direction. If the sites form
##'   more than one connected group, each group gets its own sum-to-zero
##'   constraint and scaling factor, and each site with no neighbors gets an
##'   independent effect (Freni-Sterrantino et al., 2018). See [icar_data()].
##' @param .toggles a \code{list} of toggles for model components. The
##'   components are: \itemize{\item \code{rho_mu}: 1 to explicitly relate
##'   rho to mu and 0 otherwise. \item \code{cloglog}: 1 to use the
##'   complementary log-log and 0 for the logit link function for the absence
##'   probabilities. \item \code{ar_re}: 1 to incorporate an AR(1) process to
##'   density and 0 otherwise. \item \code{iid_re}: 1 to incorporate a site
##'   specific IID random effect to density and 0 otherwise. \item \code{sp_re}:
##'   1 to incorporate a site specific ICAR random effect to density and 0
##'   otherwise.}
##' @param .priors a \code{list} of priors hyperparameters.
##' @param family a \code{character} specifying the family of the probability
##'   distribution assumed for density. The options are: \itemize{ \item
##'   \code{"gamma"} (default): gamma parametrized in terms of its mean; \item
##'   \code{"lognormal"}: log-normal parametrized in terms of its mean; \item
##'   \code{"loglogistic"}: log-logistic parametrized in terms of its median
##'   (usual parametrization), with shape \code{phi} > 1 so its mean is finite;
##'   \item \code{"lognormal_legacy"}: log-normal with
##'   its usual parametrization; }
##' @param phi_hat a \code{boolean} indicating whether the prior on \code{phi}
##'   should be determined through the data (using [get_phi_hat()]).
##' @return A \code{list} containing the data and settings to be used as input
##'   for the Stan model.
##' @author lcgodoy
##' @export
make_data_sdm <- function(y,
                          time,
                          site,
                          z,
                          x,
                          adj_mat = matrix(0, ncol = 1, nrow = 1),
                          .toggles,
                          .priors,
                          family = "gamma",
                          phi_hat = FALSE) {
  ## getting the default toggles and using user options
  stopifnot(all(!is.na(c(time))))
  stopifnot(all(!is.na(c(site))))
  stopifnot(length(family) == 1)
  stopifnot(family %in% c("lognormal_legacy", "lognormal",
                          "gamma", "loglogistic", "truncnorm"))
  likelihood <- switch(family,
                       lognormal_legacy = 0,
                       lognormal        = 1,
                       gamma            = 2,
                       loglogistic      = 3,
                       truncnorm        = 4)
  toggles <- list(rho_mu = 1,
                  ar_re = 0,
                  sp_re = 0,
                  iid_re = 0,
                  cloglog = 0) |>
    safe_modify(.toggles) |>
    c(list(likelihood = likelihood))
  for (re in c("ar_re", "iid_re", "sp_re")) {
    if (toggles[[re]] > 1) {
      stop("In an SDM, '", re, "' must be 0 or 1 (\"none\" or \"rec\").",
           call. = FALSE)
    }
  }
  ## getting the default priors and using user options
  if (phi_hat) {
    phi_hat <- get_phi_hat(y, family)
    if (missing(.priors)) {
      .priors <- list(pr_phi_a = 2, pr_phi_b = 2 / phi_hat)
    } else if (!"pr_phi_b" %in% names(.priors)) {
      if (!"pr_phi_a" %in% names(.priors)) {
        .priors <- c(.priors, list(pr_phi_a = 2, pr_phi_b = 2 / phi_hat))
      } else {
        .priors <- c(.priors, list(pr_phi_b = .priors$pr_phi_a / phi_hat))
      }
    }
  }
  priors <- default_priors() |>
    safe_modify(.priors)
  ## additional quantities that can be inferred from data
  n_sites <- length(unique(site))
  n_time <- length(unique(time))
  site_factor <- factor(site)
  site_levels <- levels(site_factor)
  site_stan <- as.integer(site_factor)
  time_init <- min(time)
  my_ord <- order(site, time)
  y <- y[my_ord]
  site <- site[my_ord]
  site_stan <- site_stan[my_ord]
  time <- time[my_ord]
  time <- time - time_init + 1
  if (missing(x)) {
    x <- matrix(1, nrow = n_time * n_sites)
    K_x <- 1
  } else {
    K_x <- ncol(x)
    x <- x[my_ord, , drop = FALSE]
    if (length(priors$pr_beta_r_mu) < K_x)
      priors$pr_beta_r_mu <- rep(0, K_x)
    if (length(priors$pr_beta_r_sd) < K_x)
      priors$pr_beta_r_sd <- rep(1, K_x)
    if (K_x == 1) {
      priors$pr_beta_r_mu <- array(priors$pr_beta_r_mu,
                                   dim = 1)
      priors$pr_beta_r_sd <- array(priors$pr_beta_r_sd,
                                   dim = 1)
    }
  }
  if (missing(z)) {
    z <- matrix(1, nrow = n_time * n_sites)
    K_z <- 1
  } else {
    K_z <- ncol(z)
    z <- z[my_ord, , drop = FALSE]
    if (length(priors$pr_beta_t_mu) < K_z)
      priors$pr_beta_t_mu <- rep(0, K_z)
    if (length(priors$pr_beta_t_sd) < K_z)
      priors$pr_beta_t_sd <- rep(1, K_z)
    if (K_z == 1) {
      priors$pr_beta_t_mu <- array(priors$pr_beta_t_mu,
                                   dim = 1)
      priors$pr_beta_t_sd <- array(priors$pr_beta_t_sd,
                                   dim = 1)
    }
  }
  if (toggles$sp_re > 0) {
    check_adj_mat(adj_mat, n_sites)
    icar <- icar_data(adj_mat)
  } else {
    icar <- icar_data_none()
  }
  zeros <- get_zeros(y)
  if (any(is.na(y)))
    y <- ifelse(is.na(y), -99, y)
  output <- list(N = n_time * n_sites,
                 n_sites = n_sites,
                 n_time = n_time,
                 time = time,
                 time_init = time_init,
                 site = site_stan,
                 site_levels = site_levels,
                 y = y,
                 X = x,
                 K_x = K_x,
                 Z = z,
                 K_z = K_z) |>
    c(icar,
      zeros,
      toggles,
      priors)
  return(output)
}
