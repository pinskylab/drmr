## Valid arguments for make_data_sdm() with ICAR random effects.
sdm_sp_args <- function() {
  c(sdm_args(), list(adj_mat = line_adj(), .toggles = list(sp_re = 1)))
}

## Valid inputs
test_that("make_data_sdm() works with valid inputs", {
  expect_no_error(out <- do.call(make_data_sdm, sdm_args()))
  expect_equal(out$N, 350)
  expect_equal(out$n_sites, 10)
  expect_equal(out$n_time, 35)
  expect_equal(dim(out$X), c(350, 2))
  expect_equal(dim(out$Z), c(350, 2))
  expect_equal(out$N_z + out$N_nz, 350)
  expect_no_error(out <- do.call(make_data_sdm, sdm_sp_args()))
  expect_equal(out$sp_re, 1)
})

test_that("make_data_sdm() accepts NA in y", {
  args <- sdm_args()
  args$y[1] <- NA
  expect_no_error(out <- do.call(make_data_sdm, args))
  expect_equal(out$y[1], -99)
  expect_equal(out$N_z + out$N_nz, 349)
})

test_that("make_data_sdm() gives the same output for rows in any order", {
  args <- sdm_args()
  idx <- rev(seq_along(args$y))
  args$y <- args$y[idx]
  args$time <- args$time[idx]
  args$site <- args$site[idx]
  args$x_t <- args$x_t[idx, ]
  args$x <- args$x[idx, ]
  expect_equal(do.call(make_data_sdm, args),
               do.call(make_data_sdm, sdm_args()))
})

## y
test_that("make_data_sdm() stops when y is not numeric", {
  skip("No input check yet")
  args <- sdm_args()
  args$y <- as.character(sdm_args()$y)
  expect_error(do.call(make_data_sdm, args))
  args$y <- sdm_args()$y > 0
  expect_error(do.call(make_data_sdm, args))
  args$y <- as.list(sdm_args()$y)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y is a factor", {
  args <- sdm_args()
  args$y <- factor(args$y)
  expect_error(suppressWarnings(do.call(make_data_sdm, args)))
})

test_that("make_data_sdm() stops when y has a negative value", {
  args <- sdm_args()
  args$y[1] <- -1
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y has an infinite value", {
  skip("No input check yet")
  args <- sdm_args()
  args$y[1] <- Inf
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y is a matrix", {
  skip("No input check yet")
  args <- sdm_args()
  args$y <- cbind(args$y, args$y)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y has the wrong length", {
  skip("No input check yet")
  args <- sdm_args()
  args$y <- sdm_args()$y[-1]
  expect_error(do.call(make_data_sdm, args))
  args$y <- c(sdm_args()$y, 1)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when all values of y are NA", {
  skip("No input check yet")
  args <- sdm_args()
  args$y <- rep(NA_real_, length(args$y))
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y is NULL", {
  args <- sdm_args()
  args["y"] <- list(NULL)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when y is missing", {
  args <- sdm_args()
  args$y <- NULL
  expect_error(do.call(make_data_sdm, args))
})

## time
test_that("make_data_sdm() stops when time is not numeric", {
  args <- sdm_args()
  args$time <- as.character(sdm_args()$time)
  expect_error(do.call(make_data_sdm, args))
  args$time <- factor(sdm_args()$time)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time is a Date", {
  skip("No input check yet")
  args <- sdm_args()
  args$time <- as.Date(paste0(args$time, "-01-01"))
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time has non-integer values", {
  skip("No input check yet")
  args <- sdm_args()
  args$time <- args$time + 0.5 * (args$site %% 2)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time points are not consecutive", {
  skip("No input check yet")
  args <- sdm_args()
  args$time <- 2 * args$time
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time has an NA", {
  args <- sdm_args()
  args$time[1] <- NA
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time has the wrong length", {
  args <- sdm_args()
  args$time <- args$time[-1]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when a site-time pair is duplicated", {
  skip("No input check yet")
  args <- sdm_args()
  args$time[2] <- args$time[1]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when a site-time pair has no row", {
  skip("No input check yet")
  args <- sdm_args()
  for (v in c("y", "time", "site")) args[[v]] <- args[[v]][-5]
  args$x_t <- args$x_t[-5, ]
  args$x <- args$x[-5, ]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when time is missing", {
  args <- sdm_args()
  args$time <- NULL
  expect_error(do.call(make_data_sdm, args))
})

## site
test_that("make_data_sdm() stops when site has an NA", {
  args <- sdm_args()
  args$site[1] <- NA
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when site has the wrong length", {
  args <- sdm_args()
  args$site <- args$site[-1]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when site is missing", {
  args <- sdm_args()
  args$site <- NULL
  expect_error(do.call(make_data_sdm, args))
})

## x
test_that("make_data_sdm() stops when x has too few rows", {
  args <- sdm_args()
  args$x <- args$x[-1, ]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x has too many rows", {
  skip("No input check yet")
  args <- sdm_args()
  args$x <- rbind(args$x, args$x[1, ])
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x is a data.frame", {
  skip("No input check yet")
  args <- sdm_args()
  args$x <- as.data.frame(args$x)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x is a vector", {
  args <- sdm_args()
  args$x <- args$x[, 2]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x has an NA", {
  skip("No input check yet")
  args <- sdm_args()
  args$x[1, 2] <- NA
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x has no columns", {
  skip("No input check yet")
  args <- sdm_args()
  args$x <- args$x[, 0, drop = FALSE]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x is a character matrix", {
  skip("No input check yet")
  args <- sdm_args()
  args$x <- matrix(as.character(args$x), ncol = ncol(args$x))
  expect_error(do.call(make_data_sdm, args))
})

## x_t
test_that("make_data_sdm() stops when x_t has too few rows", {
  args <- sdm_args()
  args$x_t <- args$x_t[-1, ]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x_t is a data.frame", {
  skip("No input check yet")
  args <- sdm_args()
  args$x_t <- as.data.frame(args$x_t)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x_t is a vector", {
  args <- sdm_args()
  args$x_t <- args$x_t[, 2]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x_t has an NA", {
  skip("No input check yet")
  args <- sdm_args()
  args$x_t[1, 2] <- NA
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x_t has no columns", {
  skip("No input check yet")
  args <- sdm_args()
  args$x_t <- args$x_t[, 0, drop = FALSE]
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when x_t is a character matrix", {
  skip("No input check yet")
  args <- sdm_args()
  args$x_t <- matrix(as.character(args$x_t), ncol = ncol(args$x_t))
  expect_error(do.call(make_data_sdm, args))
})

## adj_mat
test_that("make_data_sdm() stops when sp_re is on and adj_mat is not given", {
  args <- sdm_sp_args()
  args$adj_mat <- NULL
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must have one row and one column per site")
})

test_that("make_data_sdm() stops when adj_mat has the wrong dimensions", {
  args <- sdm_sp_args()
  args$adj_mat <- line_adj(9)
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must have one row and one column per site")
  args$adj_mat <- line_adj(10)[, -1]
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must have one row and one column per site")
  args$adj_mat <- c(line_adj(10))
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must have one row and one column per site")
})

test_that("make_data_sdm() stops when adj_mat has negative or NA values", {
  args <- sdm_sp_args()
  args$adj_mat[2, 1] <- -1
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must contain finite, non-negative values only")
  args$adj_mat[2, 1] <- NA
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must contain finite, non-negative values only")
})

test_that("make_data_sdm() stops when adj_mat is a character matrix", {
  args <- sdm_sp_args()
  args$adj_mat <- matrix(as.character(args$adj_mat), nrow = 10)
  expect_error(do.call(make_data_sdm, args),
               "'adj_mat' must contain finite, non-negative values only")
})

## .toggles
test_that("make_data_sdm() stops when a toggle name is unknown", {
  args <- sdm_args()
  args$.toggles <- list(foo = 1)
  expect_error(do.call(make_data_sdm, args), "Unknown name(s): foo",
               fixed = TRUE)
  args$.toggles <- list(movement = 1)
  expect_error(do.call(make_data_sdm, args), "Unknown name(s): movement",
               fixed = TRUE)
})

test_that("make_data_sdm() stops when ar_re is not 0 or 1", {
  args <- sdm_args()
  args$.toggles <- list(ar_re = "surv")
  expect_error(do.call(make_data_sdm, args),
               "'ar_re' must be one of \"none\", \"rec\" (or 0, 1).",
               fixed = TRUE)
  args$.toggles <- list(ar_re = 2)
  expect_error(do.call(make_data_sdm, args),
               "'ar_re' must be one of \"none\", \"rec\" (or 0, 1).",
               fixed = TRUE)
})

test_that("make_data_sdm() stops when sp_re is an unknown name", {
  args <- sdm_args()
  args$.toggles <- list(sp_re = "foo")
  expect_error(do.call(make_data_sdm, args),
               "'sp_re' must be one of \"none\", \"rec\" (or 0, 1).",
               fixed = TRUE)
})

test_that("make_data_sdm() stops when iid_re is logical", {
  args <- sdm_args()
  args$.toggles <- list(iid_re = TRUE)
  expect_error(do.call(make_data_sdm, args), "'iid_re' must be one of")
})

test_that("make_data_sdm() stops when rho_mu is not 0 or 1", {
  skip("No input check yet")
  args <- sdm_args()
  args$.toggles <- list(rho_mu = 2)
  expect_error(do.call(make_data_sdm, args))
  args$.toggles <- list(rho_mu = NA)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when cloglog is not 0 or 1", {
  skip("No input check yet")
  args <- sdm_args()
  args$.toggles <- list(cloglog = "yes")
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when .toggles is not a list", {
  args <- sdm_args()
  args$.toggles <- c(ar_re = 1)
  expect_error(do.call(make_data_sdm, args))
  args$.toggles <- "ar_re"
  expect_error(do.call(make_data_sdm, args))
  args[".toggles"] <- list(NULL)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when .toggles is an unnamed list", {
  skip("No input check yet")
  args <- sdm_args()
  args$.toggles <- list(1)
  expect_error(do.call(make_data_sdm, args))
})

## .priors
test_that("make_data_sdm() stops when a prior name is unknown", {
  args <- sdm_args()
  args$.priors <- list(foo = 1)
  expect_error(do.call(make_data_sdm, args), "Unknown name(s): foo",
               fixed = TRUE)
})

test_that("make_data_sdm() stops when a prior sd is negative", {
  skip("No input check yet")
  args <- sdm_args()
  args$.priors <- list(pr_alpha_sd = -1)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when a prior is not numeric", {
  skip("No input check yet")
  args <- sdm_args()
  args$.priors <- list(pr_alpha_sd = "a")
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when a prior is NA", {
  skip("No input check yet")
  args <- sdm_args()
  args$.priors <- list(pr_phi_a = NA)
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when .priors is not a list", {
  args <- sdm_args()
  args$.priors <- c(pr_alpha_sd = 1)
  expect_error(do.call(make_data_sdm, args))
  args$.priors <- 1
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when pr_beta_r_mu is longer than ncol(x)", {
  args <- sdm_args()
  args$.priors <- list(pr_beta_r_mu = c(0, 0, 0))
  expect_error(do.call(make_data_sdm, args),
               "'pr_beta_r_mu' must have 1 or 2 values", fixed = TRUE)
})

test_that("make_data_sdm() uses one pr_beta_t_sd value for all coefficients", {
  args <- sdm_args()
  args$.priors <- list(pr_beta_t_sd = 5)
  expect_equal(do.call(make_data_sdm, args)$pr_beta_t_sd, c(5, 5))
})

test_that("make_data_sdm() stops when a prior is used only by the DRM", {
  skip("No input check yet")
  args <- sdm_args()
  args$.priors <- list(pr_beta_s_mu = 1)
  expect_error(do.call(make_data_sdm, args))
})

## family
test_that("make_data_sdm() stops when family is not supported", {
  args <- sdm_args()
  args$family <- "poisson"
  expect_error(do.call(make_data_sdm, args))
  args$family <- "lognormal-legacy"
  expect_error(do.call(make_data_sdm, args))
  args$family <- NA
  expect_error(do.call(make_data_sdm, args))
})

test_that("make_data_sdm() stops when family is not of length 1", {
  args <- sdm_args()
  args$family <- c("gamma", "lognormal")
  expect_error(do.call(make_data_sdm, args))
  args["family"] <- list(NULL)
  expect_error(do.call(make_data_sdm, args))
})
