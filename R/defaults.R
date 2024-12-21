##' Returns default priors' hyperparameters for the DRM model.
##'
##' Based on Alexa's Summer Flounder project.
##' @title Default priors' hyperparameters
##' @return a `list` with (explain values)
##' @author lcgodoy
##' @export
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
##' @author lcgodoy
##' @export
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
