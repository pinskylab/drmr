## Valid arguments for prior_pred_drm(): 5 sites and 10 time points. The
## covariate is random, so call set.seed() first.
pp_args <- function(n_draws = 5) {
  grid <- expand.grid(site = 1:5, time = 1:10)
  grid$temp <- stats::rnorm(nrow(grid))
  list(.data = grid,
       time_col = "time",
       site_col = "site",
       n_draws = n_draws,
       formula_rec = ~ 1 + temp,
       formula_surv = ~ 1,
       n_ages = 3,
       .toggles = list(est_surv = 1, rho_mu = 0))
}

## pp_args() with all the toggles on, and with priors that differ from the
## defaults.
all_on_args <- function(n_draws = 5) {
  args <- pp_args(n_draws)
  args$.toggles <- list(est_surv = 1, rho_mu = 1, movement = 1,
                        est_init = 1, ar_re = "rec", iid_re = "dens",
                        sp_re = "surv")
  args$adj_mat <- line_adj(5)
  args$areas <- rep(1, 5)
  args$.priors <- list(pr_beta_r_mu = c(1, -1), pr_beta_r_sd = c(0.5, 2),
                       pr_beta_t_mu = 0.5, pr_beta_t_sd = 0.3,
                       pr_beta_s_mu = -1, pr_beta_s_sd = 0.7,
                       pr_lmxi_mu = -1, pr_lmxi_sd = 0.5,
                       pr_phi_a = 3, pr_phi_b = 2,
                       pr_zeta_a = 2, pr_zeta_b = 5,
                       pr_alpha_mu = 0.2, pr_alpha_sd = 0.5,
                       pr_sigma_t_loc = 0.5, pr_sigma_t_scale = 2,
                       pr_sigma_i_loc = 0, pr_sigma_i_scale = 0.3,
                       pr_sigma_s_loc = -1, pr_sigma_s_scale = 1)
  args
}

## CDF of a Cauchy distribution truncated to positive values.
ptcauchy <- function(q, location, scale) {
  p0 <- stats::pcauchy(0, location, scale)
  (stats::pcauchy(q, location, scale) - p0) / (1 - p0)
}

## The names of the parameters of the random effects in the draws.
re_vars <- list(
  ar_re = c("alpha[1]", "sigma_t[1]", sprintf("z_t[%d]", 1:10)),
  iid_re = c("sigma_i[1]", sprintf("z_i[1,%d]", 1:5)),
  sp_re = c("sigma_s[1]", sprintf("z_s[%d]", 1:5))
)

test_that("prior_pred_drm() returns the documented output", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, pp_args(n_draws = 7))
  expect_s3_class(pp, "prior_pred_drm")
  expect_named(pp, c("y_rep", "pars", "data", "nonfinite", "priors",
                     "formulas", "model_data"))
  expect_named(pp$formulas, c("formula_zero", "formula_rec", "formula_surv"))
  expect_equal(pp$formulas$formula_rec, ~ 1 + temp, ignore_attr = TRUE)
  expect_equal(pp$model_data$N, 50)
  expect_equal(pp$model_data$est_surv, 1)
  expect_true(is.matrix(pp$y_rep))
  expect_equal(dim(pp$y_rep), c(7, 50))
  expect_s3_class(pp$pars, "draws_df")
  expect_equal(posterior::ndraws(pp$pars), 7)
  expect_setequal(posterior::variables(pp$pars),
                  c("beta_r[1]", "beta_r[2]", "beta_t[1]", "beta_s[1,1]",
                    "phi[1]"))
  expect_equal(pp$nonfinite, rep(FALSE, 7))
  expect_true(all(grepl("^pr_", names(pp$priors))))
  expect_equal(pp$priors$pr_phi_a, 2)
  expect_equal(nrow(pp$data), 50)
})

test_that("print() summarizes a prior_pred_drm object", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, pp_args(n_draws = 7))
  expect_output(print(pp), "Prior predictive draws from a DRM")
  expect_output(print(pp), "Draws: 7 (0 with Inf, NaN, or NA values)",
                fixed = TRUE)
  expect_output(print(pp), "Recruitment: ~1 + temp", fixed = TRUE)
  expect_output(print(pp), "Survival: ~1", fixed = TRUE)
  expect_output(print(pp), "Time points: 10", fixed = TRUE)
  expect_output(print(pp), "Sites: 5", fixed = TRUE)
  expect_output(print(pp), "Non-zero density family: gamma", fixed = TRUE)
  expect_output(print(pp), "Parameters: beta_r, beta_t, beta_s, phi",
                fixed = TRUE)
  utils::capture.output(res <- withVisible(print(pp)))
  expect_false(res$visible)
  expect_identical(res$value, pp)
})

test_that("prior_pred_drm() works with one draw", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, pp_args(n_draws = 1))
  expect_equal(dim(pp$y_rep), c(1, 50))
  expect_equal(posterior::ndraws(pp$pars), 1)
})

test_that("the columns of y_rep follow the rows of data", {
  ## sites 1 and 2 have rho close to 1, so their densities are 0. The other
  ## sites have rho close to 0, so their densities are positive.
  set.seed(2026)
  args <- pp_args(n_draws = 20)
  args$.data$x <- ifelse(args$.data$site <= 2, 20, -20)
  args$.data <- args$.data[sample(nrow(args$.data)), ]
  args$formula_zero <- ~ 1 + x
  args$.priors <- list(pr_beta_t_mu = c(0, 1), pr_beta_t_sd = c(1e-3, 1e-3),
                       pr_phi_a = 100, pr_phi_b = 10)
  pp <- do.call(prior_pred_drm, args)
  expect_equal(order(pp$data$site, pp$data$time), seq_len(50))
  zero <- pp$data$x > 0
  expect_true(all(pp$y_rep[, zero] == 0))
  expect_true(all(pp$y_rep[, !zero] > 0))
})

test_that("prior_pred_drm() gives the same results with the same seed", {
  set.seed(2026)
  args <- all_on_args(n_draws = 5)
  set.seed(1)
  a <- suppressWarnings(do.call(prior_pred_drm, args))
  set.seed(1)
  b <- suppressWarnings(do.call(prior_pred_drm, args))
  expect_identical(a, b)
})

test_that("prior_pred_drm() matches prior_pars() followed by sim_drm()", {
  set.seed(2026)
  args <- all_on_args(n_draws = 1)
  dat_args <- args[setdiff(names(args), "n_draws")]
  dat <- do.call(build_sim_data,
                 c(dat_args, list(family = "gamma", formula_zero = ~ 1)))
  set.seed(10)
  pp <- do.call(prior_pred_drm, args)
  set.seed(10)
  sim <- do.call(sim_drm, c(dat_args, list(pars = prior_pars(dat))))
  expect_equal(pp$y_rep[1, ], sim$data$y)
  sim_data <- sim$data
  sim_data$y <- NULL
  expect_equal(pp$data, sim_data)
  pars <- stan_par_names(sim$pars)
  draw <- unlist(as.data.frame(pp$pars)[1, names(pars)])
  expect_equal(draw, pars)
})

test_that("prior_pred_drm() works with each random effect", {
  for (re in names(re_vars)) {
    for (target in c("rec", "surv", "dens")) {
      set.seed(2026)
      args <- pp_args(n_draws = 3)
      args$.toggles[[re]] <- target
      if (re == "sp_re") {
        args$adj_mat <- line_adj(5)
      }
      pp <- do.call(prior_pred_drm, args)
      expect_true(all(re_vars[[re]] %in% posterior::variables(pp$pars)),
                  label = paste(re, "=", target))
    }
  }
})

test_that("prior_pred_drm() allows phi below 1 for the log-logistic", {
  set.seed(2026)
  args <- pp_args(n_draws = 50)
  args$family <- "loglogistic"
  pp <- do.call(prior_pred_drm, args)
  expect_true(any(posterior::extract_variable(pp$pars, "phi[1]") < 1))
})

test_that("prior_pred_drm() flags the draws with non-finite values", {
  ## exp(800) overflows, so every draw has infinite densities
  set.seed(2026)
  args <- pp_args(n_draws = 10)
  args$.priors <- list(pr_beta_r_mu = c(800, 0), pr_beta_r_sd = c(1, 1e-3))
  expect_warning(pp <- do.call(prior_pred_drm, args),
                 "10 of 10 draws have Inf, NaN, or NA values.",
                 fixed = TRUE)
  expect_true(all(pp$nonfinite))
  expect_equal(dim(pp$y_rep), c(10, 50))
  expect_equal(posterior::ndraws(pp$pars), 10)
})

test_that("prior_pred_drm() gives one warning when some draws overflow", {
  ## log-recruitment spread around log(.Machine$double.xmax), about 709.8,
  ## and phi near 10: about a third of the draws overflow
  set.seed(2026)
  args <- pp_args(n_draws = 20)
  args$.priors <- list(pr_beta_r_mu = c(700, 0), pr_beta_r_sd = c(20, 1e-3),
                       pr_phi_a = 100, pr_phi_b = 10)
  msgs <- character(0)
  pp <- withCallingHandlers(
    do.call(prior_pred_drm, args),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(pp$nonfinite))
  expect_false(all(pp$nonfinite))
  expect_identical(msgs,
                   sprintf(paste("%d of 20 draws have Inf, NaN, or NA",
                                 "values. See 'nonfinite' in the output."),
                           sum(pp$nonfinite)))
  expect_true(all(is.finite(pp$y_rep[!pp$nonfinite, ])))
})

test_that("the prior draws follow the priors in drm.stan", {
  ## Kolmogorov-Smirnov tests with a fixed seed. The priors differ from the
  ## defaults, so the tests also check that .priors reaches each draw.
  set.seed(2026)
  pp <- suppressWarnings(do.call(prior_pred_drm, all_on_args(1000)))
  x <- function(v) posterior::extract_variable(pp$pars, v)
  ks_p <- function(v, ...) suppressWarnings(stats::ks.test(v, ...)$p.value)
  p_values <- c(
    "beta_r[1]" = ks_p(x("beta_r[1]"), "pnorm", 1, 0.5),
    "beta_r[2]" = ks_p(x("beta_r[2]"), "pnorm", -1, 2),
    "beta_t[1]" = ks_p(x("beta_t[1]"), "pnorm", 0.5, 0.3),
    "beta_s[1,1]" = ks_p(x("beta_s[1,1]"), "pnorm", -1, 0.7),
    "log(-xi[1])" = ks_p(log(-x("xi[1]")), "pnorm", -1, 0.5),
    "phi[1]" = ks_p(x("phi[1]"), "pgamma", shape = 3, rate = 2),
    "zeta[1]" = ks_p(x("zeta[1]"), "pbeta", 2, 5),
    "log_init[1]" = ks_p(x("log_init[1]"), "pnorm"),
    "log_init[2]" = ks_p(x("log_init[2]"), "pnorm"),
    "atanh(alpha[1])" = ks_p(atanh(x("alpha[1]")), "pnorm", 0.2, 0.5),
    "sigma_t[1]" = ks_p(x("sigma_t[1]"), ptcauchy, 0.5, 2),
    "sigma_i[1]" = ks_p(x("sigma_i[1]"), ptcauchy, 0, 0.3),
    "sigma_s[1]" = ks_p(x("sigma_s[1]"), ptcauchy, -1, 1),
    ## z_t[1] = sigma_t w_1, and z_i[1, k] has marginal SD sigma_i
    "z_t[1] / sigma_t[1]" = ks_p(x("z_t[1]") / x("sigma_t[1]"), "pnorm"),
    "z_i[1,1] / sigma_i[1]" = ks_p(x("z_i[1,1]") / x("sigma_i[1]"), "pnorm")
  )
  for (v in names(p_values)) {
    expect_gt(p_values[[v]], 1e-3, label = v)
  }
  ## z_i sums to zero, and so does z_s (the 5 sites form one group)
  z_i <- as.matrix(as.data.frame(pp$pars)[sprintf("z_i[1,%d]", 1:5)])
  z_s <- as.matrix(as.data.frame(pp$pars)[sprintf("z_s[%d]", 1:5)])
  expect_true(all(abs(rowSums(z_i)) <= 1e-8 * rowSums(abs(z_i))))
  expect_true(all(abs(rowSums(z_s)) <= 1e-8 * rowSums(abs(z_s))))
})

test_that("prior_pred_drm() stops when n_draws is not valid", {
  set.seed(2026)
  bad <- list(0, -1, 2.5, "a", NA, NA_real_, c(1, 2), numeric(0))
  for (n in bad) {
    args <- pp_args()
    args["n_draws"] <- list(n)
    expect_error(do.call(prior_pred_drm, args),
                 "'n_draws' must be a single whole number of 1 or more.",
                 fixed = TRUE)
  }
})
