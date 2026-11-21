##' @title Names of the fitted parameters that a function needs
##'
##' @description Lists the parameters to take from the draws of a fitted DRM
##'   (or SDM). \code{summary()}, \code{mcmc_diag()}, and \code{draws()}
##'   (which the plot method uses) start from the list that
##'   \code{get_fitted_pars()} gives. \code{predict()} uses the same list,
##'   but for a DRM it drops \code{sigma_i}, which \code{predict.stan} does
##'   not need. \code{fitted_pars_lambda()} gives the inputs that
##'   \code{ages_edens()} and \code{get_adults()} need, and
##'   \code{fitted_pars_ll()} gives the inputs that \code{log_lik()} and
##'   \code{fitted()} need.
##'
##' @details \code{fitted_pars_lambda()} also needs the element \code{proj} in
##'   \code{dat}: 1 if the densities are also computed for new data, and 0
##'   otherwise.
##'
##' @inheritParams prior_sample
##' @return a \code{character} vector with the names of the parameters.
##' @author lcgodoy
get_fitted_pars <- function(dat, model = "drm") {
  stopifnot(inherits(dat, "list"))
  stopifnot(length(model) == 1)
  if (model == "drm") {
    output <- fitted_pars_drm(dat)
  } else
    output <- fitted_pars_sdm(dat)
  return(output)
}

##' @rdname get_fitted_pars
fitted_pars_drm <- function(dat) {
  output <- c("beta_t", "beta_r",
              "lambda")
  if (dat$rho_mu == 1)
    output <- c(output, "xi")
  if (dat$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (dat$movement == 1)
    output <- c(output, "zeta")
  if (dat$rec_dd < 2)
    output <- c(output, "kappa")
  if (dat$est_surv == 1)
    output <- c(output, "beta_s")
  if (dat$ar_re > 0)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (dat$iid_re > 0)
    output <- c(output, "z_i", "sigma_i")
  if (dat$sp_re > 0)
    output <- c(output, "z_s")
  return(output)
}

##' @rdname get_fitted_pars
fitted_pars_sdm <- function(dat) {
  output <- c("beta_t", "beta_r")
  if (dat$rho_mu == 1)
    output <- c(output, "xi")
  if (dat$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (dat$ar_re == 1)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (dat$iid_re == 1)
    output <- c(output, "z_i")
  if (dat$sp_re == 1)
    output <- c(output, "z_s")
  return(output)
}

##' @rdname get_fitted_pars
fitted_pars_lambda <- function(dat) {
  output <- c("log_rec")
  if (dat$proj == 1)
    output <- c(output, "beta_r")
  if (dat$ar_re > 0)
    output <- c(output, "z_t")
  if (dat$proj == 1 && dat$ar_re > 0)
    output <- c(output, "alpha", "sigma_t")
  if (dat$iid_re > 0)
    output <- c(output, "z_i")
  if (dat$sp_re > 0)
    output <- c(output, "z_s")
  if (dat$est_surv == 1)
    output <- c(output, "beta_s")
  if (dat$movement == 1)
    output <- c(output, "zeta")
  if (dat$rec_dd < 2)
    output <- c(output, "kappa")
  if (dat$est_init == 1)
    output <- c(output, "log_init")
  if (dat$minit == 2 && dat$rec_dd < 2)
    output <- c(output, "log_init_scale")
  return(output)
}

##' @rdname get_fitted_pars
fitted_pars_ll <- function(dat) {
  output <- c("mu", "beta_t")
  if (dat$rho_mu == 1)
    output <- c(output, "lxi")
  if (dat$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  return(output)
}
