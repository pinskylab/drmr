##--- generating data and computing log-lik for truncated t and normal ----

dtt <- function(x, ncp = 0, df = 3, range = c(0, Inf)) {
  ll <- pt(q = min(range), ncp = ncp, df = df)
  hh <- pt(q = max(range), ncp = ncp, df = df)
  dt(x, df = df, ncp = ncp) / (hh - ll)
}

rtt <- function(n, ncp = 0, df = 3, range = c(0, Inf)) {
  ll <- pt(q = min(range), ncp = ncp, df = df)
  hh <- pt(q = max(range), ncp = ncp, df = df)
  u  <- runif(n, min = ll, max = hh)
  qt(p = u, df = df, ncp = ncp)
}

dtn <- function(x, mean = 0, sd = 1, range = c(0, Inf)) {
  ll <- pnorm(q = min(range), mean = mean, sd = sd)
  hh <- pnorm(q = max(range), mean = mean, sd = sd)
  dnorm(x, mean = mean, sd = sd) / (hh - ll)
}

rtn <- function(n, mean = 0, sd = 1, range = c(0, Inf)) {
  ll <- pnorm(q = min(range), mean = mean, sd = sd)
  hh <- pnorm(q = max(range), mean = mean, sd = sd)
  u  <- runif(n, min = ll, max = hh)
  qnorm(p = u, mean = mean, sd = sd)
}

##--- initial values for DRM ----

## given the list used to fit the model and the number of chains, generate
## initial values from the prior distribution
inits_drm <- function(dat, chains) {
  names_check <-
    c("N",
      "n_ages", "n_patches", "n_time",
      "y",
      "f",
      "m",
      "X_t",
      "K_t",
      "X_m",
      "K_m",
      "X_r",
      "K_r",
      "adj_mat",
      "age_at_maturity",
      "selectivity_at_age",
      "cloglog",
      "movement",
      "est_mort",
      "p_error",
      "qr_t",
      "qr_r",
      "qr_m",
      "likelihood",
      "pr_sigma_obs_mu",
      "pr_sigma_obs_sd",
      "pr_phi_a",
      "pr_phi_b",
      "pr_sd_r_u",
      "pr_sd_r_alpha",
      "pr_rho_u",
      "pr_rho_alpha",
      "pr_coef_t_mu",
      "pr_coef_t_sd",
      "pr_coef_m_mu",
      "pr_coef_m_sd",
      "pr_coef_r_mu",
      "pr_coef_r_sd")
  stopifnot(all(names(dat) %in% names_check))
  out <- vector(mode = "list", length = chains)
  for (k in seq_len(chains)) {
    out[[k]] <-
      list(coef_r0 = rnorm(dat$K_r,
                           mean = dat$pr_coef_r_mu,
                           sd = dat$pr_coef_r_sd),
           coef_t0 = rnorm(dat$K_t,
                           mean = dat$pr_coef_t_mu,
                           sd = dat$pr_coef_t_sd))
    if (dat$likelihood == 0) {
      out[[k]] <-
        list(sigma_obs = array(rtn(1,
                                   mean = dat$pr_sigma_obs_mu,
                                   sd = dat$pr_sigma_obs_sd,
                                   range = c(0, Inf)),
                               dim = 1))
    } else {
      out[[k]] <-
        list(phi = array(rgamma(1,
                                shape = dat$pr_phi_a,
                                scale = dat$pr_phi_b),
                         dim = 1))
    }
    if (dat$movement) {
      out[[k]] <- c(out[[k]],
                    list(logit_not_mov_prob = array(rnorm(1),
                                                    dim = 1)))
    }
    if (dat$est_m) {
      out[[k]] <- c(out[[k]],
                    list(coef_m0 = rnorm(dat$K_m,
                                         mean = dat$pr_coef_m_mu,
                                         sd = dat$pr_coef_m_sd)))
    }
    if (dat$p_error) {
      out[[k]] <-
        c(out[[k]],
          list(log_sigma_r = array(rnorm(1),
                                   dim = 1),
               rho   = array(runif(1, min = -1, max = 1),
                             dim = 1),
               raw     = rnorm(dat$n_time)))
    }
  }
  return(out)
}
