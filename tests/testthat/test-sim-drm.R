## Valid arguments for sim_drm(), from the example in its help page. The
## covariate is random, so call set.seed() first.
sim_args <- function() {
  grid <- expand.grid(site = 1:5, time = 1:30)
  grid$temp <- stats::rnorm(nrow(grid))
  list(.data = grid,
       time_col = "time",
       site_col = "site",
       formula_rec = ~ 1 + temp,
       formula_surv = ~ 1,
       n_ages = 4,
       .toggles = list(est_surv = 1, rho_mu = 0),
       pars = list(beta_r = c(1, 0.5), beta_s = 1,
                   beta_t = -2, phi = 5))
}

## sim_args() with survival fixed at exp(-m).
fixed_surv_args <- function() {
  args <- sim_args()
  args$.toggles$est_surv <- 0
  args$formula_surv <- NULL
  args$pars$beta_s <- NULL
  args
}

## sim_args() with movement between the 5 sites.
movement_args <- function() {
  args <- sim_args()
  args$.toggles$movement <- 1
  args$adj_mat <- line_adj(5)
  args$areas <- rep(1, 5)
  args$pars$zeta <- 0.5
  args
}

## sim_args() with other toggles and the parameters that they need.
toggle_args <- function(toggles, pars) {
  args <- sim_args()
  args$.toggles <- utils::modifyList(args$.toggles, toggles)
  args$pars <- utils::modifyList(args$pars, pars)
  args
}

test_that("sim_drm() works with valid inputs", {
  set.seed(2026)
  args <- sim_args()
  sim <- do.call(sim_drm, args)
  expect_equal(dim(sim$data), c(150, 4))
  expect_equal(dim(sim$lambda), c(4, 30, 5))
  expect_type(do.call(sim_drm, fixed_surv_args()), "list")
  expect_type(do.call(sim_drm, movement_args()), "list")
  ## the toggles that the tests below use
  ok <- list(
    list(list(rec_dd = "ricker"), list(kappa = 1)),
    list(list(rec_dd = "bh"), list(kappa = 1)),
    list(list(rho_mu = 1), list(xi = -0.5)),
    list(list(est_init = 1), list(log_init = c(0, 0, 0))),
    list(list(ar_re = "rec"), list(alpha = 0.5, sigma_t = 0.3)),
    list(list(ar_re = "rec"), list(z_t = stats::rnorm(30))),
    list(list(iid_re = "rec"), list(sigma_i = 0.5)),
    list(list(iid_re = "rec"), list(z_i = c(-2, -1, 0, 1, 2)))
  )
  for (x in ok) {
    expect_type(do.call(sim_drm, toggle_args(x[[1]], x[[2]])), "list")
  }
  ## one site, the log-logistic family, and phi below 1 with the gamma
  one <- args
  one$.data <- args$.data[args$.data$site == 1, ]
  expect_type(do.call(sim_drm, one), "list")
  args$family <- "loglogistic"
  expect_type(do.call(sim_drm, args), "list")
  args <- sim_args()
  args$pars$phi <- 0.5
  expect_type(do.call(sim_drm, args), "list")
})

## .data
test_that("sim_drm() stops when .data is not a data.frame", {
  set.seed(2026)
  args <- sim_args()
  args$.data <- as.matrix(sim_args()$.data)
  expect_error(do.call(sim_drm, args))
  args$.data <- as.list(sim_args()$.data)
  expect_error(do.call(sim_drm, args))
  args[".data"] <- list(NULL)
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when a site-time row is missing or duplicated", {
  set.seed(2026)
  args <- sim_args()
  grid <- args$.data
  msg <- "exactly one row for each site and time point"
  args$.data <- grid[-7, ]
  expect_error(do.call(sim_drm, args), msg)
  args$.data <- rbind(grid, grid[7, ])
  expect_error(do.call(sim_drm, args), msg)
})

test_that("sim_drm() stops when the time points are not consecutive", {
  set.seed(2026)
  args <- sim_args()
  grid <- args$.data
  args$.data$time[grid$time > 15] <- grid$time[grid$time > 15] + 1
  expect_error(do.call(sim_drm, args), "consecutive integers")
  args$.data$time <- 2 * grid$time
  expect_error(do.call(sim_drm, args), "consecutive integers")
  args$.data$time <- grid$time / 2
  expect_error(do.call(sim_drm, args), "consecutive integers")
})

test_that("sim_drm() stops when .data has no rows", {
  set.seed(2026)
  args <- sim_args()
  args$.data <- args$.data[0, ]
  expect_error(suppressWarnings(do.call(sim_drm, args)))
})

## time_col and site_col
test_that("sim_drm() stops when time_col is not a column of .data", {
  set.seed(2026)
  args <- sim_args()
  args$time_col <- "tm"
  expect_error(suppressWarnings(do.call(sim_drm, args)))
})

test_that("sim_drm() stops when time_col has two names", {
  set.seed(2026)
  args <- sim_args()
  args$time_col <- c("time", "site")
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when site_col is not a column of .data", {
  set.seed(2026)
  args <- sim_args()
  args$site_col <- "patch"
  expect_error(do.call(sim_drm, args))
})

## formulas
test_that("sim_drm() stops when a covariate has NA", {
  set.seed(2026)
  args <- sim_args()
  args$.data$temp[3] <- NA
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when formula_rec uses a missing column", {
  set.seed(2026)
  args <- sim_args()
  args$formula_rec <- ~ 1 + sst
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when est_surv = 1 and formula_surv is NULL", {
  set.seed(2026)
  args <- sim_args()
  args$formula_surv <- NULL
  expect_error(do.call(sim_drm, args), "needs 'formula_surv'",
               fixed = TRUE)
})

test_that("sim_drm() stops when formula_surv is given and est_surv = 0", {
  set.seed(2026)
  args <- sim_args()
  args$.toggles$est_surv <- 0
  expect_error(do.call(sim_drm, args),
               "'formula_surv' needs '.toggles = list(est_surv = 1)'",
               fixed = TRUE)
})

## pars
test_that("sim_drm() stops when pars lacks a required parameter", {
  set.seed(2026)
  for (p in c("beta_r", "phi", "beta_s")) {
    args <- sim_args()
    args$pars[[p]] <- NULL
    expect_error(do.call(sim_drm, args),
                 paste0("'pars' must contain: ", p, "."), fixed = TRUE)
  }
})

test_that("sim_drm() stops when pars lacks a parameter for a toggle", {
  set.seed(2026)
  args <- sim_args()
  args$.toggles$rec_dd <- "bh"
  expect_error(do.call(sim_drm, args), "'pars' must contain: kappa.",
               fixed = TRUE)
  args <- movement_args()
  args$pars$zeta <- NULL
  expect_error(do.call(sim_drm, args), "'pars' must contain: zeta.",
               fixed = TRUE)
})

test_that("sim_drm() stops when pars is not a list", {
  set.seed(2026)
  args <- sim_args()
  args["pars"] <- list(NULL)
  expect_error(do.call(sim_drm, args), "'pars' must contain:",
               fixed = TRUE)
  args$pars <- unlist(sim_args()$pars)
  expect_error(do.call(sim_drm, args), "'pars' must contain: beta_r.",
               fixed = TRUE)
  args$pars <- c(beta_r = 1, beta_r = 0.5, beta_s = 1, beta_t = -2,
                 phi = 5)
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when a beta vector has the wrong length", {
  set.seed(2026)
  args <- sim_args()
  args$pars$beta_r <- 1
  expect_error(do.call(sim_drm, args))
  args <- sim_args()
  args$pars$beta_t <- c(-2, 1)
  expect_error(do.call(sim_drm, args))
  args <- sim_args()
  args$pars$beta_s <- c(1, 2)
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when a parameter is not numeric", {
  set.seed(2026)
  args <- sim_args()
  args$pars$phi <- "5"
  expect_error(do.call(sim_drm, args))
  args <- sim_args()
  args$pars$beta_r <- c("1", "0.5")
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when beta_r has NA", {
  skip("No input check yet")
  set.seed(2026)
  args <- sim_args()
  args$pars$beta_r <- c(NA, 0.5)
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when beta_r names do not match the design", {
  skip("No input check yet")
  set.seed(2026)
  args <- sim_args()
  args$pars$beta_r <- c(temp = 0.5, "(Intercept)" = 1)
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when pars has an unknown name", {
  skip("No input check yet")
  set.seed(2026)
  args <- sim_args()
  args$pars$foo <- 1
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when pars has a parameter the model ignores", {
  skip("No input check yet")
  set.seed(2026)
  args <- sim_args()
  args$pars$sigma_t <- 0.3
  expect_error(do.call(sim_drm, args))
  args <- fixed_surv_args()
  args$pars$beta_s <- 3
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when phi is not positive", {
  set.seed(2026)
  args <- sim_args()
  args$pars$phi <- 0
  expect_error(do.call(sim_drm, args))
  args$pars$phi <- NA
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when phi is below 1 for the log-logistic", {
  set.seed(2026)
  args <- sim_args()
  args$family <- "loglogistic"
  args$pars$phi <- 0.5
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when phi has two values", {
  set.seed(2026)
  args <- sim_args()
  args$pars$phi <- c(5, 5)
  expect_error(do.call(sim_drm, args), "'phi' must be a single value.",
               fixed = TRUE)
})

test_that("sim_drm() stops when xi is positive", {
  set.seed(2026)
  args <- toggle_args(list(rho_mu = 1), list(xi = 0.5))
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when kappa is negative", {
  set.seed(2026)
  args <- toggle_args(list(rec_dd = "ricker"), list(kappa = -1))
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when zeta is above 1", {
  set.seed(2026)
  args <- movement_args()
  args$pars$zeta <- 1.5
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when log_init has the wrong length", {
  set.seed(2026)
  args <- toggle_args(list(est_init = 1), list(log_init = c(0, 0)))
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when z_t has the wrong length", {
  set.seed(2026)
  args <- toggle_args(list(ar_re = "rec"), list(z_t = stats::rnorm(10)))
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when alpha or sigma_t is out of range", {
  set.seed(2026)
  args <- toggle_args(list(ar_re = "rec"), list(alpha = 0.5, sigma_t = -1))
  expect_error(do.call(sim_drm, args))
  args <- toggle_args(list(ar_re = "rec"), list(alpha = 1.2, sigma_t = 0.3))
  expect_error(do.call(sim_drm, args))
})

test_that("sim_drm() stops when z_i does not sum to zero", {
  set.seed(2026)
  args <- toggle_args(list(iid_re = "rec"), list(z_i = 1:5))
  expect_error(do.call(sim_drm, args),
               "'z_i' must have one value per site and sum to zero.",
               fixed = TRUE)
})

test_that("sim_drm() stops when iid_re is on with a single site", {
  set.seed(2026)
  args <- toggle_args(list(iid_re = "rec"), list(sigma_i = 0.5))
  args$.data <- args$.data[args$.data$site == 1, ]
  expect_error(do.call(sim_drm, args))
})

## arguments passed on to make_data(). test-make-data.R tests the checks
## that make_data() runs on them.
test_that("sim_drm() stops when f_mort has the wrong dimensions", {
  set.seed(2026)
  args <- sim_args()
  msg <- "'f_mort' must have one row per age and one column per time point."
  args$f_mort <- matrix(0, 4, 29)
  expect_error(do.call(sim_drm, args), msg, fixed = TRUE)
  args$f_mort <- matrix(0, 30, 4)
  expect_error(do.call(sim_drm, args), msg, fixed = TRUE)
})

test_that("sim_drm() stops when amat or weight has the wrong length", {
  set.seed(2026)
  msg <- "'amat' and 'weight' must have one value per age."
  args <- sim_args()
  args$amat <- c(0, 1, 1)
  expect_error(do.call(sim_drm, args), msg, fixed = TRUE)
  args <- sim_args()
  args$weight <- 1:5
  expect_error(do.call(sim_drm, args), msg, fixed = TRUE)
})

test_that("sim_drm() stops when ... has an unknown argument", {
  set.seed(2026)
  args <- sim_args()
  args$foo <- 1
  expect_error(do.call(sim_drm, args))
})

## y_col
test_that("sim_drm() stops when y_col is not a column name", {
  skip("No input check yet")
  set.seed(2026)
  args <- sim_args()
  args$y_col <- 1
  expect_error(do.call(sim_drm, args))
})
