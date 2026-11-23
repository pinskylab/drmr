##' @title Age-specific expected densities based on DRM.
##'
##' @inheritParams preddrm
##'
##' @details The function builds the design matrices again from
##'   \code{newdata}. Terms that depend on the data, such as \code{scale()} or
##'   \code{poly()}, take their values from \code{newdata}, not from the
##'   training data. To standardize a covariate, use the mean and standard
##'   deviation from the training data. A factor must have the same levels as
##'   in the training data.
##'
##' @author Anonymous
##'
##' @return an object of class \code{"CmdStanGQ"} containing samples for the
##'   posterior predictive distribution for forecasting.
##'
##' @name age_dens
##' 
##' @export
ages_edens <- function(object,
                       newdata = NULL,
                       f_test = NULL,
                       cores = 1) {
  stopifnot(inherits(object$stanfit,
                     c("CmdStanFit", "CmdStanLaplace",
                       "CmdStanPathfinder", "CmdStanVB")))  
  ##--- list for lambda_gq ----
  gq_data <- object$data
  gq_data$areas <- get_areas(gq_data)
  gq_data$n_edges_adj <- get_n_edges_adj(gq_data)
  gq_data$proj <- as.integer(!is.null(newdata))
  
  if (gq_data$proj == 1) {
    ## lambda_gq.stan reads the rows as one block per site code, so sort by
    ## the site order of the training data (not by the raw site values),
    ## and then by time
    ord <- order(factor(newdata[[object$cols$site_col]],
                        levels = object$cols$site_levels),
                 newdata[[object$cols$time_col]])
    newdata <- newdata[ord, ]

    ## number time points for forecasting
    ntime_for <- length(unique(newdata[[object$cols$time_col]]))

    time_for <- newdata[[object$cols$time_col]] -
      min(newdata[[object$cols$time_col]]) + 1
    
    ## throw an error if predicing on new sites
    new_sites_factor <- factor(newdata[[object$cols$site_col]], 
                               levels = object$cols$site_levels)
    if (any(is.na(new_sites_factor))) {
      stop("newdata contains sites that were not present in the training data.")
    }
    new_sites_stan <- as.integer(new_sites_factor)
    
    if (is.null(f_test)) {
      f_test <- matrix(0, ncol = ntime_for, nrow = object$data$n_ages)
    }
    
    gq_data$n_proj <- array(ntime_for, dim = 1)
    gq_data$time_proj <- time_for
    gq_data$site_proj <- new_sites_stan
    gq_data$f_proj <- f_test
    gq_data$X_r <- stats::model.matrix(object$formulas$formula_rec,
                                       data = newdata)
    if (object$data$est_surv == 1) {
      gq_data$X_mproj <- stats::model.matrix(object$formulas$formula_surv,
                                             data = newdata)
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
    object$stanfit$draws(variables = pars)
  ## load compiled lambda
  lambda_comp <-
    instantiate::stan_package_model(name = "lambda_gq",
                                    package = "drmr")
  output <- lambda_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = gq_data,
                        parallel_chains = cores)
  check_gq(output)
  spt <- data.frame(v1 = object$cols$site_levels[object$data$site],
                    v2 = object$data$time + object$data$time_init - 1,
                    site_id = object$data$site,
                    time = object$data$time)
  colnames(spt)[1:2] <- rev(unname(unlist(object$cols[2:3])))
  spt <- unique(spt)
  
  spt_proj <- NULL
  if (gq_data$proj == 1) {
    spt_proj <- data.frame(v1 = new_sites_factor,
                           v2 = time_for + max(object$data$time) + object$data$time_init - 1,
                           site_id = new_sites_stan,
                           time = time_for) |>
      unique()
    colnames(spt_proj)[1:2] <- rev(unname(unlist(object$cols[2:3])))
  }
  
  output <- list("gq" = output,
                 "spt" = spt,
                 "spt_proj" = spt_proj)
  return(new_aesd(output))
}

##' Extract and Summarize Age-Specific Densities
##'
##' Parses the \code{gq} element of the [ages_edens()] output to return a data
##' frame of age-specific densities, indexed by age, time, and patch.
##'
##' @inheritParams summary.CmdStanGQ object
##' @inheritParams summary.drmrmodels probs
##' @param ages An optional vector of integers. If provided, the output is filtered
##'   to include only these specific ages.
##'
##' @return A \code{data.frame} containing the summary statistics and parsed indices
##'   (age, time, patch).
##' @export
summarise_adens <- function(object, ages = NULL,
                            probs = c(0.05, 0.5, 0.95)) {
  summ <- object$summary(NULL, "mean", "sd",
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
