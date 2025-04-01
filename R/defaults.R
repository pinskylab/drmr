##' Returns default priors' hyperparameters for the DRM model.
##'
##' @title Default priors' hyperparameters
##' @return a `list` with (explain values)
##' @author lcgodoy
##' @export
default_priors <- function() {
  list(pr_sigma_obs_mu = 0,
       pr_sigma_obs_sd = 1,
       pr_phi_a = 2,
       pr_phi_b = 1,
       pr_logsd_r_mu = -2,
       pr_logsd_r_sd = .25,
       pr_alpha_a = .5,
       pr_alpha_b = .5,
       pr_zeta_a = .5,
       pr_zeta_b = .5,
       pr_beta_t_mu = 0,
       pr_beta_t_sd = 1,
       pr_beta_m_mu = numeric(0),
       pr_beta_m_sd = numeric(0),
       pr_beta_r_mu = 0,
       pr_beta_r_sd = 1)
}

##' Returns default toggles for the DRM model.
##'
##' @title Default toggles
##' @return a `list` with (explain values)
##' @author lcgodoy
##' @export
default_toggles <- function() {
  list(cloglog = 0,
       movement = 0,
       est_mort = 0,
       est_init = 1,
       time_ar = 0)
}
