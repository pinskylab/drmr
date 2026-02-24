##' Returns default priors' hyperparameters for the DRM model.
##'
##' @title Default priors' hyperparameters
##' @return A \code{list} containing the default hyperparameter values for the priors.
##' @author lcgodoy
##' @export
default_priors <- function() {
  list(pr_sigma_obs_mu = 0,
       pr_sigma_obs_sd = 1,
       pr_phi_a = 2,
       pr_phi_b = 1,
       pr_lsigma_t_mu = -2,
       pr_lsigma_t_sd = .25,
       pr_lsigma_i_mu = -2,
       pr_lsigma_i_sd = .25,
       pr_lsigma_s_mu = -2,
       pr_lsigma_s_sd = .25,
       pr_alpha_a = .5,
       pr_alpha_b = .5,
       pr_zeta_a = .5,
       pr_zeta_b = .5,
       pr_lmxi_mu = 0,
       pr_lmxi_sd = 2.5,
       pr_beta_t_mu = 0,
       pr_beta_t_sd = 1,
       pr_beta_s_mu = numeric(0),
       pr_beta_s_sd = numeric(0),
       pr_beta_r_mu = 0,
       pr_beta_r_sd = 1)
}

##' Returns default toggles for the DRM model.
##'
##' @title Default toggles
##' @return A \code{list} containing the default values for model toggles.
##' @author lcgodoy
##' @export
default_toggles <- function() {
  list(rho_mu = 1,
       cloglog = 0,
       movement = 0,
       est_surv = 0,
       est_init = 0,
       minit = 0,
       ar_re = 0,
       iid_re = 0,
       sp_re = 0)
}

##' Returns default for NUTS
##'
##' @title Default NUTS arguments
##' @return A \code{list} containing the default arguments for the NUTS algorithm.
##' @author lcgodoy
##' @export
default_nuts <- function() {
  list(chains = 4,
       parallel_chains = getOption("mc.cores", 1),
       threads_per_chain = NULL,
       opencl_ids = NULL,
       iter_warmup = NULL,
       iter_sampling = NULL,
       save_warmup = FALSE,
       thin = NULL,
       max_treedepth = NULL,
       adapt_engaged = TRUE,
       adapt_delta = NULL,
       step_size = NULL,
       metric = NULL,
       metric_file = NULL,
       inv_metric = NULL,
       init_buffer = NULL,
       term_buffer = NULL,
       window = NULL,
       fixed_param = FALSE,
       show_messages = TRUE,
       show_exceptions = TRUE,
       diagnostics = c("divergences", "treedepth", "ebfmi"),
       save_metric = NULL,
       save_cmdstan_config = NULL)
}

##' Returns default for Variational Bayes
##'
##' @title Default VB arguments
##' @return A \code{list} containing the default arguments for the Variational Bayes algorithm.
##' @author lcgodoy
##' @export
default_vb <- function() {
  list(threads = NULL,
       opencl_ids = NULL,
       algorithm = NULL,
       iter = NULL,
       grad_samples = NULL,
       elbo_samples = NULL,
       eta = NULL,
       adapt_engaged = NULL,
       adapt_iter = NULL,
       tol_rel_obj = NULL,
       eval_elbo = NULL,
       output_samples = NULL,
       draws = NULL,
       show_messages = TRUE,
       show_exceptions = TRUE,
       save_cmdstan_config = NULL)
}

##' Returns default for Pathfinder
##'
##' @title Default Pathfinder arguments
##' @return A \code{list} containing the default arguments for the Pathfinder algorithm.
##' @author lcgodoy
##' @export
default_pf <- function() {
  list(opencl_ids = NULL,
       num_threads = NULL,
       init_alpha = NULL,
       tol_obj = NULL,
       tol_rel_obj = NULL,
       tol_grad = NULL,
       tol_rel_grad = NULL,
       tol_param = NULL,
       history_size = NULL,
       single_path_draws = NULL,
       draws = NULL,
       num_paths = 4,
       max_lbfgs_iters = NULL,
       num_elbo_draws = NULL,
       save_single_paths = NULL,
       psis_resample = NULL,
       calculate_lp = NULL,
       show_messages = TRUE,
       show_exceptions = TRUE,
       save_cmdstan_config = NULL)
}

##' Returns default for Laplace
##'
##' @title Default Laplace arguments
##' @return A \code{list} containing the default arguments for the Laplace algorithm.
##' @author lcgodoy
##' @export
default_laplace <- function() {
  list(threads = NULL,
       opencl_ids = NULL,
       mode = NULL,
       opt_args = NULL,
       jacobian = TRUE,
       draws = NULL,
       show_messages = TRUE,
       show_exceptions = TRUE,
       save_cmdstan_config = NULL)
}

##' Returns default for Optimization
##'
##' @title Default Optimization arguments
##' @return A \code{list} containing the default arguments for the Optimization algorithm.
##' @author lcgodoy
##' @export
default_opt <- function() {
  list(threads = NULL,
       opencl_ids = NULL,
       algorithm = NULL,
       jacobian = FALSE,
       init_alpha = NULL,
       iter = NULL,
       tol_obj = NULL,
       tol_rel_obj = NULL,
       tol_grad = NULL,
       tol_rel_grad = NULL,
       tol_param = NULL,
       history_size = NULL,
       show_messages = TRUE,
       show_exceptions = TRUE,
       save_cmdstan_config = NULL)
}

##' Returns default for a given algorithm
##'
##' @title Default arguments for inference algorithm
##' @inheritParams fit_drm
##' @return A \code{list} containing the default arguments for the specified algorithm.
##' @author lcgodoy
##' @export
default_algo <- function(algorithm = "nuts",
                         algo_args = list()) {
    switch(algorithm,
           nuts       = default_nuts(),
           vb         = default_vb(),
           pathfinder = default_pf(),
           laplace    = default_laplace(),
           optimize   = default_opt()) |>
      safe_modify(algo_args)
}
