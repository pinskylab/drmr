## Valid arguments for prior_pred_drm() with the covariate temp in each
## formula: 5 sites, 10 time points, and 20 draws. The covariate is random,
## so call set.seed() first.
eff_args <- function(toggles = list(est_surv = 1, rho_mu = 0)) {
  grid <- expand.grid(site = 1:5, time = 1:10)
  grid$temp <- stats::rnorm(nrow(grid))
  list(.data = grid,
       time_col = "time",
       site_col = "site",
       n_draws = 20,
       formula_zero = ~ 1 + temp,
       formula_rec = ~ 1 + temp,
       formula_surv = if (isTRUE(toggles$est_surv == 1)) ~ 1 + temp,
       n_ages = 3,
       .toggles = toggles)
}

## A mock of a fitted DRM whose draws are the prior draws in pp. draws()
## reads the draws through the stanfit element, as for a real fit.
mock_fit <- function(pp) {
  stanfit <- list(
    draws = function(variables = NULL, inc_warmup = FALSE, format = NULL) {
      posterior::as_draws_matrix(
        posterior::subset_draws(pp$pars, variable = variables)
      )
    },
    metadata = function() list(method = "sample")
  )
  class(stanfit) <- "CmdStanFit"
  out <- list(stanfit = stanfit, data = pp$model_data,
              formulas = pp$formulas)
  class(out) <- c("adrm", "drmrmodels")
  out
}

## The draws of the coefficients named par, as a matrix with one row per
## draw.
coef_draws <- function(pp, par) {
  unclass(posterior::as_draws_matrix(
    posterior::subset_draws(pp$pars, variable = par)
  ))
}

newdata <- data.frame(temp = c(-1, 0, 1))

test_that("effects_drm() on prior draws follows each process", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, eff_args())
  x <- stats::model.matrix(~ 1 + temp, newdata)
  expected <- list(
    rec = exp(tcrossprod(coef_draws(pp, "beta_r"), x)),
    surv = stats::plogis(tcrossprod(coef_draws(pp, "beta_s"), x)),
    pabs = stats::plogis(tcrossprod(coef_draws(pp, "beta_t"), x))
  )
  out <- c(rec = "recruitment", surv = "survival", pabs = "prob_abs")
  for (p in names(expected)) {
    eff <- effects_drm(pp, process = p, variable = "temp",
                       newdata = newdata, summary = FALSE)
    expect_equal(eff[[out[[p]]]], c(expected[[p]]), label = p)
    expect_equal(eff$temp, rep(newdata$temp, each = 20), label = p)
  }
})

test_that("effects_drm() summarizes the prior draws", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, eff_args())
  eff <- effects_drm(pp, process = "rec", variable = "temp", n_pts = 7)
  expect_s3_class(eff, "eff_drm")
  expect_equal(nrow(eff), 7)
  expect_equal(attr(eff, "quant_cols"), c("q5", "q50", "q95"))
  expect_equal(attr(eff, "process"), "recruitment")
  expect_equal(range(eff$temp), range(pp$model_data$X_r[, "temp"]))
})

test_that("effects_drm() gives the same results for a fit and prior draws", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, eff_args())
  fit <- mock_fit(pp)
  for (p in c("rec", "surv", "pabs")) {
    for (s in c(TRUE, FALSE)) {
      expect_equal(effects_drm(fit, process = p, variable = "temp",
                               summary = s),
                   effects_drm(pp, process = p, variable = "temp",
                               summary = s),
                   label = paste(p, "with summary =", s))
    }
  }
})

test_that("effects_drm() uses the cloglog link for pabs when the model does", {
  set.seed(2026)
  args <- eff_args(list(est_surv = 1, rho_mu = 0, cloglog = 1))
  pp <- do.call(prior_pred_drm, args)
  x <- stats::model.matrix(~ 1 + temp, newdata)
  expected <- -expm1(-exp(tcrossprod(coef_draws(pp, "beta_t"), x)))
  eff <- effects_drm(pp, process = "pabs", variable = "temp",
                     newdata = newdata, summary = FALSE)
  expect_equal(eff$prob_abs, c(expected))
  eff_fit <- effects_drm(mock_fit(pp), process = "pabs", variable = "temp",
                         newdata = newdata, summary = FALSE)
  expect_equal(eff_fit$prob_abs, c(expected))
})

test_that("effects_drm() stops for survival when est_surv is off", {
  set.seed(2026)
  pp <- do.call(prior_pred_drm, eff_args(list(est_surv = 0, rho_mu = 0)))
  msg <- paste("'process = \"surv\"' needs a model with",
               "'.toggles = list(est_surv = 1)'.")
  expect_error(effects_drm(pp, process = "surv", variable = "temp"),
               msg, fixed = TRUE)
  expect_error(effects_drm(mock_fit(pp), process = "surv", variable = "temp"),
               msg, fixed = TRUE)
})

test_that("plot() works on effects_drm() output from prior draws", {
  skip_if_not_installed("ggplot2")
  set.seed(2026)
  pp <- do.call(prior_pred_drm, eff_args())
  eff <- effects_drm(pp, process = "rec", variable = "temp")
  expect_true(inherits(plot(eff), "ggplot"))
})
