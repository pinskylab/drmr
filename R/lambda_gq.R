##' @title Age-specific expected densities based on DRM.
##'
##' @param drm A \code{list} object containing the output from the [fit_drm]
##'   function.
##' @param new_data a \code{data.frame} with the dataset at which we wish to
##'   obtain predictions.
##' @param f_test a \code{matrix} informing the instantaneous fishing mortality
##'   rates at each age (columns) and timepoint (rows).
##' @param cores number of threads used for the forecast. If four chains were
##'   used in the \code{drm}, then four (or less) threads are recommended.
##'
##' @author lcgodoy
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for forecasting.
##'
##' @name age_dens
##' 
##' @export
ages_edens <- function(drm,
                       new_data = NULL,
                       f_test = NULL,
                       cores = 1) {
  stopifnot(inherits(drm$stanfit,
                     c("CmdStanFit", "CmdStanLaplace",
                       "CmdStanPathfinder", "CmdStanVB")))  
  ##--- list for lambda_gq ----
  gq_data <- drm$data
  gq_data$proj <- as.integer(!is.null(new_data))
  
  if (gq_data$proj == 1) {
    ## Ensure new_data is ordered by site and then time
    ord <- order(new_data[[drm$cols$site_col]],
                 new_data[[drm$cols$time_col]])
    new_data <- new_data[ord, ]

    ## number time points for forecasting
    ntime_for <- length(unique(new_data[[drm$cols$time_col]]))

    time_for <- new_data[[drm$cols$time_col]] -
      min(new_data[[drm$cols$time_col]]) + 1
    
    ## throw an error if predicing on new sites
    new_sites_factor <- factor(new_data[[drm$cols$site_col]], 
                               levels = drm$cols$site_levels)
    if (any(is.na(new_sites_factor))) {
      stop("new_data contains sites that were not present in the training data.")
    }
    new_sites_stan <- as.integer(new_sites_factor)
    
    if (is.null(f_test)) {
      f_test <- matrix(0, ncol = ntime_for, nrow = drm$data$n_ages)
    }
    
    gq_data$n_proj <- array(ntime_for, dim = 1)
    gq_data$time_proj <- time_for
    gq_data$site_proj <- new_sites_stan
    gq_data$f_proj <- f_test
    gq_data$X_r <- stats::model.matrix(drm$formulas$formula_rec,
                                       data = new_data)
    if (drm$data$est_surv == 1) {
      gq_data$X_mproj <- stats::model.matrix(drm$formulas$formula_surv,
                                             data = new_data)
    } else {
      gq_data$X_mproj <- matrix(0, nrow = 1, ncol = 1)
    }
  } else {
    gq_data$n_proj <- integer(0)
    gq_data$time_proj <- integer(0)
    gq_data$site_proj <- integer(0)
    gq_data$f_proj <- matrix(0, nrow = 0, ncol = 0)
    gq_data$X_r <- matrix(0, nrow = 0, ncol = gq_data$K_r)
    gq_data$X_mproj <- matrix(0, nrow = 1, ncol = 1)
  }

  ##--- pars from model fitted ----
  pars <- fitted_pars_lambda(gq_data)
  fitted_params <-
    drm$stanfit$draws(variables = pars)
  ## load compiled lambda
  lambda_comp <-
    instantiate::stan_package_model(name = "lambda_gq",
                                    package = "drmr")
  output <- lambda_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = gq_data,
                        parallel_chains = cores)
  spt <- data.frame(v1 = drm$cols$site_levels[drm$data$site],
                    v2 = drm$data$time + drm$data$time_init - 1,
                    site_id = drm$data$site,
                    time = drm$data$time)
  colnames(spt)[1:2] <- rev(unname(unlist(drm$cols[2:3])))
  spt <- unique(spt)
  
  spt_proj <- NULL
  if (gq_data$proj == 1) {
    spt_proj <- data.frame(v1 = new_sites_factor,
                           v2 = time_for + max(drm$data$time) + drm$data$time_init - 1,
                           site_id = new_sites_stan,
                           time = time_for)
    colnames(spt_proj)[1:2] <- rev(unname(unlist(drm$cols[2:3])))
    spt_proj <- unique(spt_proj)
  }
  
  output <- list("gq" = output,
                 "spt" = spt,
                 "spt_proj" = spt_proj)
  return(new_aesd(output))
}

##' @rdname age_dens
##' @param ... params to be passed to \code{age_edens}
##' @export
lambda_drm <- function(...) {
  lifecycle::deprecate_warn("0.3.1", "lambda_drm()", "ages_edens()")
  ages_edens(...)
}

##' Extract and Summarize Age-Specific Densities
##'
##' Parses the output of \code{lambda_drm} to return a data frame of
##' age-specific densities, indexed by age, time, and patch.
##'
##' @param lambda_obj A \code{CmdStanGQ} object returned by \code{lambda_drm}.
##' @param ages An optional vector of integers. If provided, the output is filtered
##'   to include only these specific ages.
##' @param probs A numeric vector of quantiles to calculate in the summary.
##'   Defaults to c(0.05, 0.5, 0.95).
##'
##' @return A \code{data.frame} containing the summary statistics and parsed indices
##'   (age, time, patch).
##' @export
summarise_adens <- function(lambda_obj, ages = NULL,
                            probs = c(0.05, 0.5, 0.95)) {
  summ <- lambda_obj$summary(NULL, "mean", "sd",
                             ~quantile(., probs = probs))
  summ <- as.data.frame(summ)
  summ <- summ[grep("^lambda\\[", summ$variable), ]
  vars <- summ$variable
  age  <- as.numeric(sub(".*\\[([0-9]+),.*", "\\1", vars))
  time <- as.integer(sub(".*,([0-9]+),.*", "\\1", vars))
  patch <- as.integer(sub(".*,([0-9]+)].*", "\\1", vars))
  out <- cbind(
    data.frame(age = age, time = time, patch = patch),
    summ
  )
  if (!is.null(ages)) {
    out <- out[out$age %in% ages, ]
  }
  rownames(out) <- NULL 
  return(out)
}
