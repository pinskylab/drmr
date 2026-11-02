##' @title Simulate a zero-mean AR(1) process
##' @description Follows \code{drm.stan}: \eqn{z_1 = \sigma w_1} and
##'   \eqn{z_t = \alpha z_{t - 1} + \sigma w_t}, where the \eqn{w_t} are
##'   independent standard normal draws. The first value does not come from
##'   the stationary distribution.
##' @param n_time number of time points.
##' @param alpha autocorrelation parameter, between -1 and 1.
##' @param sigma conditional standard deviation.
##' @return a \code{numeric} vector of length \code{n_time}.
##' @author lcgodoy
sim_ar <- function(n_time, alpha, sigma) {
  w <- stats::rnorm(n_time)
  out <- numeric(n_time)
  out[1] <- sigma * w[1]
  for (i in seq_len(n_time)[-1]) {
    out[i] <- alpha * out[i - 1] + sigma * w[i]
  }
  return(out)
}

##' @title Simulate an IID random effect as in drm.stan
##' @description Draws a sum-to-zero vector with marginal standard deviation
##'   \code{sigma}. This is the prior on \code{z_i} in \code{drm.stan}.
##' @param n_sites number of sites (at least 2).
##' @param sigma marginal standard deviation.
##' @return a \code{numeric} vector of length \code{n_sites} that sums to
##'   zero.
##' @author lcgodoy
sim_iid <- function(n_sites, sigma) {
  stopifnot(n_sites > 1)
  e <- stats::rnorm(n_sites, sd = sigma * sqrt(n_sites / (n_sites - 1)))
  return(e - mean(e))
}

##' @title Simulate an ICAR random effect as in drm.stan
##' @description Draws \code{w} as in the prior on \code{w_s} in
##'   \code{drm.stan}. Within each connected group of sites, \code{w} has the
##'   Gaussian density proportional to
##'   \eqn{\exp\{-0.5 \sum_k (w_{n_{1k}} - w_{n_{2k}})^2\}}, where
##'   \eqn{(n_{1k}, n_{2k})} is column \eqn{k} of \code{neighbors}, and sums
##'   to zero. A site with no neighbors gets an independent standard normal
##'   value. The function returns \code{sigma * w / sqrt(scaling)}, which is
##'   \code{z_s} in \code{drm.stan}.
##' @param neighbors a \code{matrix} with 2 rows. Each column is a pair of
##'   neighbors, and each pair appears once (see [icar_data()]).
##' @param n_sites number of sites.
##' @param sigma standard deviation parameter (\code{sigma_s}).
##' @param scaling one scaling factor per site, from [icar_data()].
##' @return a \code{numeric} vector of length \code{n_sites}.
##' @author lcgodoy
sim_icar <- function(neighbors, n_sites, sigma, scaling) {
  prec <- matrix(0, nrow = n_sites, ncol = n_sites)
  for (k in seq_len(NCOL(neighbors))) {
    i <- neighbors[1, k]
    j <- neighbors[2, k]
    prec[i, i] <- prec[i, i] + 1
    prec[j, j] <- prec[j, j] + 1
    prec[i, j] <- prec[i, j] - 1
    prec[j, i] <- prec[j, i] - 1
  }
  ## the null space of prec holds one constant vector per connected group,
  ## so drawing in its complement gives one sum-to-zero constraint per group
  w <- numeric(n_sites)
  eig <- eigen(prec, symmetric = TRUE)
  pos <- eig$values > 1e-8 * max(eig$values, 1)
  if (any(pos)) {
    w <- as.numeric(eig$vectors[, pos, drop = FALSE] %*%
                      (stats::rnorm(sum(pos)) / sqrt(eig$values[pos])))
  }
  single <- setdiff(seq_len(n_sites), c(neighbors))
  w[single] <- stats::rnorm(length(single))
  return(sigma * w / sqrt(scaling))
}

##' @title Density-dependent recruitment as in drm.stan
##' @param psi log-recruitment linear predictor at one time point, one value
##'   per site.
##' @param ssb spawning stock biomass at the previous time point, one value
##'   per site.
##' @param kappa density-dependence parameter.
##' @param rec_dd 0 for Ricker and 1 for Beverton-Holt.
##' @return a \code{numeric} vector of expected recruits. Recruits are zero
##'   where \code{ssb} is not above \code{1e-10}.
##' @author lcgodoy
dd_rec <- function(psi, ssb, kappa, rec_dd) {
  out <- numeric(length(ssb))
  ## a NaN SSB gives zero recruits, as in Stan
  ok <- !is.na(ssb) & ssb > 1e-10
  s <- ssb[ok]
  log_dd <- if (rec_dd == 0) -kappa * s else -log(kappa + s)
  out[ok] <- exp(psi[ok] + log(s) + log_dd)
  return(out)
}

##' @title Simulate population dynamics
##'
##' @description Computes the expected density by age, time, and site. The
##'   function follows the Stan functions \code{simplest},
##'   \code{simplest_movement}, \code{pop_rec_dd}, and
##'   \code{pop_rec_dd_movement} in \code{age_struct.stanfunctions}.
##'
##' @details These details follow the Stan code on purpose:
##'   \itemize{
##'   \item At the first time point, recruitment is
##'   \code{exp(log_rec[1, ])}, also with density dependence.
##'   \item \code{minit} has no effect with density dependence.
##'   \item Movement is \code{zeta * x + (1 - zeta) * solve(A) \%*\% t(W)
##'   \%*\% A \%*\% x}, where \code{W} is the row-standardized adjacency
##'   matrix and \code{A} is \code{diag(areas)}. So the \code{(1 - zeta)}
##'   share of the abundance at each site moves to its neighbors in
##'   proportion to the values in its row of \code{adj_mat} (in equal parts
##'   for a 0/1 matrix), and the total abundance (density times area) stays
##'   the same.
##'   For a site with no neighbors, \code{W} has a 1 on the diagonal, so its
##'   individuals stay there.
##'   \item With \code{acc_dd = 1}, the SSB is the mean SSB density of the
##'   sites, weighted by \code{areas}.
##'   }
##'
##' @param log_rec a \code{n_time} by \code{n_sites} \code{matrix} with the
##'   log-recruitment linear predictor (\eqn{\psi_{t, i}}).
##' @param log_surv a \code{n_time} by \code{n_sites} \code{matrix} with the
##'   log-survival (\eqn{\log s_{t, i}}).
##' @param f a \code{n_ages} by \code{n_time} \code{matrix} with the
##'   instantaneous fishing mortality rates.
##' @param init a \code{numeric} vector of length \code{n_ages - 1} with the
##'   log-densities of ages 2 to \code{n_ages} at the first time point.
##' @param minit 1 to start from the age structure that survival implies at
##'   the first time point, and 0 to use \code{init}.
##' @param rec_dd 0 for Ricker, 1 for Beverton-Holt, and 2 for no density
##'   dependence.
##' @param kappa density-dependence parameter.
##' @param amat maturity-at-age vector.
##' @param weight weight-at-age vector.
##' @param acc_dd 1 to use the area-weighted mean SSB density of the whole
##'   region, and 0 to use the SSB density of each site.
##' @param zeta probability of staying at the same site. \code{NULL} means no
##'   movement.
##' @param adj_mat a \code{n_sites} by \code{n_sites} adjacency
##'   \code{matrix}.
##' @param ages_movement a 0/1 vector of length \code{n_ages}. A value of 1
##'   means the age group moves.
##' @param areas a positive \code{numeric} vector with the area of each site.
##'   \code{NULL} means 1 for each site.
##' @return an \code{array} of expected densities with dimensions
##'   \code{n_ages}, \code{n_time}, and \code{n_sites}.
##' @author lcgodoy
pop_dyn <- function(log_rec, log_surv, f, init,
                    minit = 0,
                    rec_dd = 2,
                    kappa = NULL,
                    amat = NULL,
                    weight = NULL,
                    acc_dd = 0,
                    zeta = NULL,
                    adj_mat = NULL,
                    ages_movement = NULL,
                    areas = NULL) {
  n_time <- NROW(log_rec)
  n_sites <- NCOL(log_rec)
  n_ages <- NROW(f)
  use_minit <- minit == 1 && rec_dd == 2
  stopifnot(NROW(log_surv) == n_time,
            NCOL(log_surv) == n_sites,
            NCOL(f) == n_time,
            use_minit || length(init) == n_ages - 1)
  if (is.null(areas)) {
    areas <- rep(1, n_sites)
  }
  stopifnot(length(areas) == n_sites,
            all(is.finite(areas) & areas > 0))
  movement <- !is.null(zeta)
  if (movement) {
    stopifnot(NROW(adj_mat) == n_sites,
              NCOL(adj_mat) == n_sites,
              length(ages_movement) == n_ages)
    rs <- rowSums(adj_mat)
    adj_std <- adj_mat
    adj_std[rs > 0, ] <- adj_mat[rs > 0, ] / rs[rs > 0]
    ## a site with no neighbors keeps its individuals
    for (k in which(rs <= 0)) {
      adj_std[k, ] <- 0
      adj_std[k, k] <- 1
    }
    ## element (i, j) is W[j, i] * areas[j] / areas[i]
    adj_std <- t(adj_std) * outer(1 / areas, areas)
  }
  lambda <- array(0, dim = c(n_ages, n_time, n_sites))
  lambda[1, 1, ] <- exp(log_rec[1, ])
  for (a in seq_len(n_ages)[-1]) {
    lambda[a, 1, ] <- if (use_minit) {
      lambda[a - 1, 1, ] * exp(log_surv[1, ] - f[a - 1, 1])
    } else {
      exp(init[a - 1])
    }
  }
  mat_weight <- amat * weight
  for (i in seq_len(n_time)[-1]) {
    if (rec_dd == 2) {
      lambda[1, i, ] <- exp(log_rec[i, ])
    } else {
      ssb <- colSums(matrix(lambda[, i - 1, ], nrow = n_ages) * mat_weight)
      if (acc_dd == 1) {
        ssb <- rep(sum(ssb * areas) / sum(areas), n_sites)
      }
      lambda[1, i, ] <- dd_rec(log_rec[i, ], ssb, kappa, rec_dd)
    }
    for (a in seq_len(n_ages)[-1]) {
      x <- lambda[a - 1, i - 1, ] * exp(log_surv[i - 1, ] - f[a - 1, i - 1])
      if (movement && ages_movement[a] == 1) {
        x <- zeta * x + (1 - zeta) * as.numeric(adj_std %*% x)
      }
      lambda[a, i, ] <- x
    }
  }
  return(lambda)
}

##' @title Simulate observed densities as in drm.stan
##' @description Follows \code{drmsdm_rng} in \code{lpdfs.stanfunctions}.
##'   Each value is zero with probability \code{rho}. Otherwise, it comes
##'   from the family that \code{likelihood} selects.
##' @param mu a \code{numeric} vector. It is the mean for codes 0, 1, and 2,
##'   the median for code 3, and the location before truncation for code 4.
##' @param rho a \code{numeric} vector with the probabilities of a zero.
##' @param likelihood 0 (log-normal, legacy), 1 (log-normal), 2 (gamma),
##'   3 (log-logistic), or 4 (truncated normal).
##' @param phi \code{sigma_obs} for code 0, and \code{phi} otherwise.
##' @return a \code{numeric} vector with the same length as \code{mu}.
##' @author lcgodoy
sim_obs <- function(mu, rho, likelihood, phi) {
  y <- numeric(length(mu))
  pos <- which(stats::rbinom(length(mu), size = 1, prob = rho) == 0)
  m <- mu[pos]
  k <- length(pos)
  y[pos] <- switch(
    as.character(likelihood),
    "0" = stats::rlnorm(k, meanlog = log(m) - phi^2 / 2, sdlog = phi),
    "1" = stats::rlnorm(k, meanlog = log(m) - phi, sdlog = sqrt(2 * phi)),
    "2" = stats::rgamma(k, shape = phi, rate = phi / m),
    "3" = {
      u <- stats::runif(k)
      m * (u / (1 - u))^(1 / phi)
    },
    "4" = m - phi * stats::qnorm(stats::runif(k, 0, stats::pnorm(m / phi))),
    stop("Unknown likelihood code: ", likelihood, ".")
  )
  return(y)
}

##' @title Check a data list before a simulation
##' @param dat a \code{list} created by [make_data()].
##' @return \code{dat}, invisibly. The function stops if it cannot simulate
##'   from \code{dat}.
##' @author lcgodoy
check_sim_dat <- function(dat) {
  panel <- identical(as.numeric(dat$time),
                     as.numeric(rep(seq_len(dat$n_time), dat$n_sites))) &&
    identical(as.numeric(dat$site),
              as.numeric(rep(seq_len(dat$n_sites), each = dat$n_time)))
  if (!panel) {
    stop("The data must have exactly one row for each site and time ",
         "point, and the time points must be consecutive integers.")
  }
  if (NROW(dat$f) != dat$n_ages || NCOL(dat$f) != dat$n_time) {
    stop("'f_mort' must have one row per age and one column per time ",
         "point.")
  }
  if (length(dat$amat) != dat$n_ages || length(dat$weight) != dat$n_ages) {
    stop("'amat' and 'weight' must have one value per age.")
  }
  re_codes <- c(dat$ar_re, dat$iid_re, dat$sp_re)
  if (any(re_codes == 2) && dat$est_surv == 0) {
    stop("Random effects on survival need ",
         "'.toggles = list(est_surv = 1)'.")
  }
  return(invisible(dat))
}

##' @title Check the parameters of a simulation
##' @param pars a named \code{list} of parameters. See [sim_drm()].
##' @param dat a \code{list} created by [make_data()].
##' @return \code{pars}, invisibly. The function stops if \code{pars} does
##'   not fit the model that \code{dat} defines.
##' @author lcgodoy
check_sim_pars <- function(pars, dat) {
  need <- c("beta_r", "beta_t",
            if (dat$likelihood == 0) "sigma_obs" else "phi",
            if (dat$est_surv == 1) "beta_s",
            if (dat$rho_mu == 1) "xi",
            if (dat$rec_dd < 2) "kappa",
            if (dat$movement == 1) "zeta",
            if (dat$est_init == 1) "log_init",
            if (dat$ar_re > 0 && is.null(pars$z_t)) c("alpha", "sigma_t"),
            if (dat$iid_re > 0 && is.null(pars$z_i)) "sigma_i",
            if (dat$sp_re > 0 && is.null(pars$z_s)) "sigma_s")
  miss <- setdiff(need, names(pars))
  if (length(miss) > 0) {
    stop("'pars' must contain: ", paste(miss, collapse = ", "), ".")
  }
  scalars <- c("xi", "sigma_obs", "phi", "kappa", "zeta", "alpha",
               "sigma_t", "sigma_i", "sigma_s")
  for (p in intersect(scalars, names(pars))) {
    if (length(pars[[p]]) != 1) {
      stop("'", p, "' must be a single value.")
    }
  }
  stopifnot(
    length(pars$beta_r) == dat$K_r,
    length(pars$beta_t) == dat$K_t,
    dat$est_surv == 0 || length(pars$beta_s) == dat$K_m[1],
    dat$rho_mu == 0 || (length(pars$xi) == 1 && pars$xi < 0),
    dat$likelihood > 0 || pars$sigma_obs > 0,
    dat$likelihood == 0 || pars$phi > as.numeric(dat$likelihood == 3),
    dat$rec_dd == 2 || pars$kappa > 0,
    dat$movement == 0 || (pars$zeta >= 0 && pars$zeta <= 1),
    dat$est_init == 0 || length(pars$log_init) == dat$n_ages - 1,
    is.null(pars$z_t) || length(pars$z_t) == dat$n_time,
    dat$ar_re == 0 || !is.null(pars$z_t) ||
      (abs(pars$alpha) < 1 && pars$sigma_t > 0),
    dat$iid_re == 0 || !is.null(pars$z_i) || pars$sigma_i > 0,
    dat$sp_re == 0 || !is.null(pars$z_s) || pars$sigma_s > 0
  )
  toggle_of <- c(z_t = "ar_re", z_i = "iid_re", z_s = "sp_re")
  for (z in intersect(names(toggle_of), names(pars))) {
    if (dat[[toggle_of[[z]]]] == 0) {
      stop("'", z, "' has no effect unless '", toggle_of[[z]],
           "' is on.")
    }
    if (!is.numeric(pars[[z]]) || any(!is.finite(pars[[z]]))) {
      stop("'", z, "' must contain finite numbers only.")
    }
  }
  sums_to_zero <- function(v) abs(sum(v)) <= 1e-8 * sum(abs(v))
  v <- pars$z_i
  if (!is.null(v) && (length(v) != dat$n_sites || !sums_to_zero(v))) {
    stop("'z_i' must have one value per site and sum to zero.")
  }
  v <- pars$z_s
  if (!is.null(v)) {
    ok <- length(v) == dat$n_sites
    if (ok && is.null(dat$comp_size_s)) {
      ## data from an earlier version: one constraint over all sites
      ok <- sums_to_zero(v)
    }
    pos <- 1
    for (k in dat$comp_size_s) {
      ok <- ok && sums_to_zero(v[dat$comp_order_s[pos:(pos + k - 1)]])
      pos <- pos + k
    }
    if (!ok) {
      stop("'z_s' must have one value per site and sum to zero within ",
           "each connected group of sites. Sites with no neighbors are ",
           "not constrained.")
    }
  }
  return(invisible(pars))
}

##' @title Draw the random effects of a simulation
##' @description Draws \code{z_t}, \code{z_i}, and \code{z_s} as in
##'   \code{drm.stan}, unless \code{pars} already contains them.
##' @inheritParams check_sim_pars
##' @return \code{pars} with the random effects added.
##' @author lcgodoy
draw_re <- function(pars, dat) {
  if (dat$ar_re > 0 && is.null(pars$z_t)) {
    pars$z_t <- sim_ar(dat$n_time, pars$alpha, pars$sigma_t)
  }
  if (dat$iid_re > 0 && is.null(pars$z_i)) {
    pars$z_i <- sim_iid(dat$n_sites, pars$sigma_i)
  }
  if (dat$sp_re > 0 && is.null(pars$z_s)) {
    pars$z_s <- sim_icar(dat$neighbors, dat$n_sites,
                         pars$sigma_s, dat$scaling)
  }
  return(pars)
}

##' @title Random effects for one target process
##' @inheritParams check_sim_pars
##' @param target 1 for recruitment, 2 for survival, and 3 for density.
##' @return a \code{numeric} vector with one value for each row of
##'   \code{dat}.
##' @author lcgodoy
re_at <- function(pars, dat, target) {
  out <- numeric(dat$N)
  if (dat$ar_re == target) {
    out <- out + pars$z_t[dat$time]
  }
  if (dat$iid_re == target) {
    out <- out + pars$z_i[dat$site]
  }
  if (dat$sp_re == target) {
    out <- out + pars$z_s[dat$site]
  }
  return(out)
}

##' @title Draw parameters from the priors in drm.stan
##' @description Draws one set of parameters from the priors in the model
##'   block of \code{drm.stan}. The names and scales are the ones that
##'   [sim_drm()] expects.
##' @param dat a \code{list} created by [make_data()].
##' @return a named \code{list} of parameters.
##' @seealso [pp_sim()], [sim_drm()]
##' @author lcgodoy
prior_pars <- function(dat) {
  rn <- function(n, mu, sd) stats::rnorm(n, mean = mu, sd = sd)
  ## normal truncated below at 0; sampling the upper tail stays finite when
  ## most of the mass is below 0
  rn_pos <- function(mu, sd) {
    p_up <- stats::pnorm(0, mean = mu, sd = sd, lower.tail = FALSE)
    stats::qnorm(stats::runif(1, min = 0, max = p_up),
                 mean = mu, sd = sd, lower.tail = FALSE)
  }
  out <- list(beta_r = rn(dat$K_r, dat$pr_beta_r_mu, dat$pr_beta_r_sd),
              beta_t = rn(dat$K_t, dat$pr_beta_t_mu, dat$pr_beta_t_sd))
  if (dat$est_surv == 1) {
    out$beta_s <- rn(dat$K_m[1], dat$pr_beta_s_mu, dat$pr_beta_s_sd)
  }
  if (dat$rho_mu == 1) {
    out$xi <- -exp(rn(1, dat$pr_lmxi_mu, dat$pr_lmxi_sd))
  }
  if (dat$likelihood == 0) {
    out$sigma_obs <- rn_pos(dat$pr_sigma_obs_mu, dat$pr_sigma_obs_sd)
  } else {
    ## the log-logistic requires phi > 1, as in drm.stan
    lower <- as.numeric(dat$likelihood == 3)
    p_up <- stats::pgamma(lower, shape = dat$pr_phi_a, rate = dat$pr_phi_b,
                          lower.tail = FALSE)
    out$phi <- stats::qgamma(stats::runif(1, min = 0, max = p_up),
                             shape = dat$pr_phi_a, rate = dat$pr_phi_b,
                             lower.tail = FALSE)
  }
  if (dat$rec_dd < 2) {
    out$kappa <- rn_pos(dat$pr_kappa_mu, dat$pr_kappa_sd)
  }
  if (dat$movement == 1) {
    out$zeta <- stats::rbeta(1, dat$pr_zeta_a, dat$pr_zeta_b)
  }
  if (dat$est_init == 1) {
    out$log_init <- stats::rnorm(dat$n_ages - 1)
  }
  if (dat$ar_re > 0) {
    out$alpha <- tanh(rn(1, dat$pr_alpha_mu, dat$pr_alpha_sd))
    out$sigma_t <- exp(rn(1, dat$pr_lsigma_t_mu, dat$pr_lsigma_t_sd))
  }
  if (dat$iid_re > 0) {
    out$sigma_i <- exp(rn(1, dat$pr_lsigma_i_mu, dat$pr_lsigma_i_sd))
  }
  if (dat$sp_re > 0) {
    out$sigma_s <- exp(rn(1, dat$pr_lsigma_s_mu, dat$pr_lsigma_s_sd))
  }
  return(out)
}

## The rlang import below is for the ggplot2 code in R/marg.R.

##' @title Simulate from a DRM given a data list and parameters
##'
##' @description The engine behind [sim_drm()] and [pp_sim()]. It follows the
##'   transformed parameters block of \code{drm.stan} and the random number
##'   generator \code{drmsdm_rng}.
##'
##' @param dat a \code{list} created by [make_data()]. The function ignores
##'   the response \code{y}.
##' @param pars a named \code{list} of parameters. See [sim_drm()].
##' @return A \code{list} with the elements below. Vectors follow the row
##'   order of \code{dat}: by site, then by time. \itemize{
##'   \item \code{y}: simulated densities.
##'   \item \code{mu}: expected density, \eqn{\mu_{t, i}}.
##'   \item \code{rho}: probability of a zero, \eqn{\rho_{t, i}}.
##'   \item \code{log_rec}: log-recruitment linear predictor,
##'   \eqn{\psi_{t, i}}, with its random effects.
##'   \item \code{surv}: survival before fishing mortality,
##'   \eqn{s_{t, i}}.
##'   \item \code{lambda}: expected density by age, time, and site.
##'   \item \code{pars}: \code{pars} plus the random effects drawn.
##'   }
##' @importFrom rlang .data
##' @author lcgodoy
model_sim <- function(dat, pars) {
  check_sim_dat(dat)
  check_sim_pars(pars, dat)
  pars <- draw_re(pars, dat)
  log_rec <- as.numeric(dat$X_r %*% pars$beta_r) + re_at(pars, dat, 1)
  if (dat$est_surv == 1) {
    lp_surv <- as.numeric(dat$X_m %*% pars$beta_s) + re_at(pars, dat, 2)
    log_surv <- stats::plogis(lp_surv, log.p = TRUE)
  } else {
    log_surv <- rep(-dat$m[1], dat$N)
  }
  init <- if (dat$est_init == 1) pars$log_init else dat$init_data
  lambda <- pop_dyn(log_rec = matrix(log_rec, nrow = dat$n_time),
                    log_surv = matrix(log_surv, nrow = dat$n_time),
                    f = dat$f,
                    init = init,
                    minit = dat$minit,
                    rec_dd = dat$rec_dd,
                    kappa = pars$kappa,
                    amat = dat$amat,
                    weight = dat$weight,
                    acc_dd = dat$acc_dd,
                    zeta = if (dat$movement == 1) pars$zeta,
                    adj_mat = dat$adj_mat,
                    ages_movement = dat$ages_movement,
                    areas = get_areas(dat))
  mu <- apply(lambda, c(2, 3),
              function(x) sum(x * dat$selectivity_at_age))
  mu <- as.numeric(mu)
  mu[is.nan(mu)] <- 0
  mu <- mu * exp(re_at(pars, dat, 3))
  eta <- as.numeric(dat$X_t %*% pars$beta_t)
  if (dat$rho_mu == 1) {
    eta <- eta + pars$xi * log(mu)
  }
  rho <- if (dat$cloglog == 1) -expm1(-exp(eta)) else stats::plogis(eta)
  phi <- if (dat$likelihood == 0) pars$sigma_obs else pars$phi
  return(list(y = sim_obs(mu, rho, dat$likelihood, phi),
              mu = mu,
              rho = rho,
              log_rec = log_rec,
              surv = exp(log_surv),
              lambda = lambda,
              pars = pars))
}

##' @title Simulate from the prior predictive distribution of a DRM
##' @description Draws one set of parameters from the priors in
##'   \code{drm.stan}. Then it simulates data with [model_sim()].
##' @param dat a \code{list} created by [make_data()], or the \code{data}
##'   element of a [fit_drm()] output.
##' @return The output of [model_sim()], plus \code{hyper_pars}: the prior
##'   hyperparameters.
##' @seealso [prior_pars()], [model_sim()]
##' @author lcgodoy
pp_sim <- function(dat) {
  out <- model_sim(dat, prior_pars(dat))
  out$hyper_pars <- dat[grepl("^pr_", names(dat))]
  return(out)
}

##' @title Simulate data from a dynamic range model
##' @export
##' @description Simulates densities from the DRM that [fit_drm()] fits, for
##'   known parameter values. Use it for parameter-recovery studies and other
##'   simulation studies. The data, column, family, formula, and \code{...}
##'   arguments are the same as in [fit_drm()]. To fit the simulated data,
##'   call [fit_drm()] on \code{sim$data} with the same arguments, plus
##'   \code{y_col} and \code{seed}, and without \code{pars}.
##'
##' @details The function builds the model data with [make_data()]. Then it
##'   follows \code{drm.stan} step by step. Call [set.seed()] first to make
##'   the results reproducible.
##'
##'   Give \code{pars} on the natural scale: standard deviations, not their
##'   logarithms, and \code{alpha}, not its inverse hyperbolic tangent. The
##'   exception is \code{log_init}, which holds the log-densities of ages 2
##'   to \code{n_ages} at the first time point. The elements that you need
##'   depend on the toggles:
##'   \itemize{
##'   \item \code{beta_r} and \code{beta_t}: always.
##'   \item \code{phi}, or \code{sigma_obs} for \code{"lognormal_legacy"}:
##'   always.
##'   \item \code{beta_s}: when \code{est_surv = 1}.
##'   \item \code{xi} (negative): when \code{rho_mu = 1}.
##'   \item \code{kappa}: when \code{rec_dd} is \code{"ricker"} or
##'   \code{"bh"}.
##'   \item \code{zeta}: when \code{movement = 1}.
##'   \item \code{log_init}: when \code{est_init = 1}. Otherwise, the
##'   simulation uses \code{init_data} from [make_data()].
##'   \item \code{alpha} and \code{sigma_t}: when \code{ar_re} is not
##'   \code{"none"}.
##'   \item \code{sigma_i}: when \code{iid_re} is not \code{"none"}.
##'   \item \code{sigma_s}: when \code{sp_re} is not \code{"none"}.
##'   }
##'   To fix a random effect instead of drawing it, give \code{z_t} (one value
##'   per time point), \code{z_i}, or \code{z_s} (one value per site; the
##'   values of \code{z_i} must sum to zero, and those of \code{z_s} must sum
##'   to zero within each connected group of sites; sites with no neighbors
##'   are not constrained).
##'
##'   The data must have one row for each site and time point, and the time
##'   points must be consecutive integers. To study missing surveys, set some
##'   simulated values to \code{NA} before you fit the model.
##'
##' @inheritParams fit_drm
##' @param pars a named \code{list} with the true parameter values. See
##'   Details.
##' @param y_col name of the column for the simulated densities. The function
##'   replaces an existing column with this name.
##' @param ... passed on to [make_data()]. For example: \code{n_ages},
##'   \code{f_mort} (one row per age and one column per time point),
##'   \code{m}, \code{init_data}, \code{age_selectivity},
##'   \code{ages_movement}, \code{adj_mat}, \code{areas}, \code{amat},
##'   \code{weight}, and \code{.toggles}.
##' @return A \code{list} with four elements:
##'   \itemize{
##'   \item \code{data}: \code{.data} sorted by site and time, with the
##'   simulated densities in column \code{y_col}.
##'   \item \code{latent}: a \code{data.frame} with the true latent values
##'   for each row of \code{data}: \code{psi} (log-recruitment linear
##'   predictor), \code{surv} (survival before fishing mortality), \code{mu}
##'   (expected density), and \code{rho} (probability of a zero).
##'   \item \code{lambda}: an \code{array} with the true expected density by
##'   age, time, and site. \code{lambda[1, , ]} is the recruitment.
##'   \item \code{pars}: \code{pars} plus the random effects drawn.
##'   }
##' @seealso [fit_drm()], [ages_edens()]
##' @examples
##' set.seed(2026)
##' grid <- expand.grid(site = 1:5, time = 1:30)
##' grid$temp <- rnorm(nrow(grid))
##' sim <- sim_drm(grid,
##'                time_col = "time",
##'                site_col = "site",
##'                formula_rec = ~ 1 + temp,
##'                formula_surv = ~ 1,
##'                n_ages = 4,
##'                .toggles = list(est_surv = 1, rho_mu = 0),
##'                pars = list(beta_r = c(1, 0.5), beta_s = 1,
##'                            beta_t = -2, phi = 5))
##' head(sim$data)
##' @author lcgodoy
sim_drm <- function(.data,
                    time_col,
                    site_col,
                    pars,
                    family = "gamma",
                    formula_zero = ~ 1,
                    formula_rec = ~ 1,
                    formula_surv = NULL,
                    y_col = "y",
                    ...) {
  est_surv <- list(...)[[".toggles"]][["est_surv"]]
  if (isTRUE(est_surv == 1) && is.null(formula_surv)) {
    stop("'.toggles = list(est_surv = 1)' needs 'formula_surv'.")
  }
  mm <- function(f) {
    stats::model.matrix.lm(f, data = .data, na.action = stats::na.fail)
  }
  args <- list(y = rep(1, NROW(.data)),
               time = .data[[time_col]],
               site = .data[[site_col]],
               family = family,
               x_t = mm(formula_zero),
               x_r = mm(formula_rec))
  if (!is.null(formula_surv)) {
    args$x_m <- mm(formula_surv)
  }
  dat <- do.call(make_data, c(args, list(...)))
  if (!is.null(formula_surv) && dat$est_surv == 0) {
    stop("'formula_surv' needs '.toggles = list(est_surv = 1)'. ",
         "Otherwise, make_data() fixes survival at exp(-m).")
  }
  sim <- model_sim(dat, pars)
  ord <- order(.data[[site_col]], .data[[time_col]])
  out <- .data[ord, , drop = FALSE]
  rownames(out) <- NULL
  out[[y_col]] <- sim$y
  latent <- out[c(site_col, time_col)]
  latent$psi <- sim$log_rec
  latent$surv <- sim$surv
  latent$mu <- sim$mu
  latent$rho <- sim$rho
  dimnames(sim$lambda) <-
    list(age = seq_len(dat$n_ages),
         time = sort(unique(.data[[time_col]])),
         site = dat$site_levels)
  return(list(data = out,
              latent = latent,
              lambda = sim$lambda,
              pars = sim$pars))
}

##' @title Turn an array of density per age, time, and patch/site into a
##'   \code{data.frame}
##' @param lbd a 3-dimensional \code{array}
##' @return return a \code{data.frame}.
##' @author lcgodoy
lambda2df <- function(lbd) {
  dim_lambdas <- dim(lbd)
  dimnames(lbd) <-
    list("age"  = seq_len(dim_lambdas[1]),
         "time" = seq_len(dim_lambdas[2]),
         "site" = seq_len(dim_lambdas[3]))
  out <-
    array2DF(lbd,
             responseName = "density")
  out <- out |>
    transform(age = as.integer(out[["age"]]),
              time = as.integer(out[["time"]]),
              site = as.integer(out[["site"]]))
  return(out)
}

## ##' @title Generate samples from the prior predictive distribution of model
## ##'   parameters
## ##' @description Generates samples from the prior distributions of the model
## ##'   parameters. This is primarily used for prior predictive checks or to
## ##'   generate initial values for MCMC.
## ##' @param n a \code{integer} representing the number of samples form the
## ##'   posterior predictive distribution.
## ##' @inheritParams model_sim
## ##' @param ... parameters passed to [pop_dyn()].
## ##' @return A \code{list} containing samples drawn from the prior distributions
## ##'   of the model parameters.
## ##' @seealso [prior_inits()] [pps_aux()]
## ##' @author lcgodoy
## prior_predictive_sim <- function(n, dat, model,
##                                  selectivity,
##                                  ...) {
##   out <-
##     replicate(n, pps_aux(dat, model,
##                          selectivity, ...),
##               simplify = FALSE)
##   return(out)
## }
