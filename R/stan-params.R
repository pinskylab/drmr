##' @title Retrieve parameters needed for predictions
##'
##' @description This function identifies the parameters necessary for carrying
##'   out predictions based on the data used to fit a DRM (or SDM).
##'
##' @param data_list a \code{list} used as input for model fitting. Typically,
##'   the output from the [make_data] function.
##' @param model a \code{character} indicating which model predictions are
##'   sought for. This input admits two possible entries: "drm" (default) or
##'   "sdm".
##' @return a \code{character} vector of labels indicating the parameters
##'   necessary for the predictions.
##' @author lcgodoy
get_fitted_pars <- function(data_list, model = "drm") {
  stopifnot(inherits(data_list, "list"))
  stopifnot(length(model) == 1)
  if (model == "drm") {
    output <- fitted_pars_drm(data_list)
  } else
    output <- fitted_pars_sdm(data_list)
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_drm <- function(data_list) {
  output <- c("beta_t", "beta_r",
              "lambda")
  if (data_list$rho_mu == 1)
    output <- c(output, "xi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$movement == 1)
    output <- c(output, "zeta")
  if (data_list$est_surv == 1)
    output <- c(output, "beta_s")
  if (data_list$ar_re > 0)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (data_list$iid_re > 0)
    output <- c(output, "z_i", "sigma_i")
  if (data_list$sp_re > 0)
    output <- c(output, "z_s")
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_sdm <- function(data_list) {
  output <- c("beta_t", "beta_r")
  if (data_list$rho_mu == 1)
    output <- c(output, "xi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  if (data_list$ar_re == 1)
    output <- c(output, "z_t", "alpha", "sigma_t")
  if (data_list$iid_re == 1)
    output <- c(output, "z_i")
  if (data_list$sp_re == 1)
    output <- c(output, "z_s")
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_lambda <- function(data_list) {
  output <- c("log_rec")
  if (data_list$est_surv == 1)
    output <- c(output, "beta_s")
  if (data_list$movement == 1)
    output <- c(output, "zeta")
  if (data_list$est_init == 1)
    output <- c(output, "log_init")
  return(output)
}

##' @inherit get_fitted_pars
fitted_pars_ll <- function(data_list) {
  output <- c("mu", "beta_t")
  if (data_list$rho_mu == 1)
    output <- c(output, "lxi")
  if (data_list$likelihood == 0) {
    output <- c(output, "sigma_obs")
  } else {
    output <- c(output, "phi")
  }
  return(output)
}
