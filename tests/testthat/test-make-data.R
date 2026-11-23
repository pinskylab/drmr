## Valid arguments with the movement, survival, or ICAR toggles on.
mov_args <- function() {
  args <- drm_args()
  args$adj_mat <- line_adj()
  args$areas <- patch_areas()
  args$.toggles <- list(movement = 1)
  args
}

surv_args <- function() {
  args <- drm_args()
  args$x_m <- stats::model.matrix(~ c_stemp, data = sum_fl)
  args$.toggles <- list(est_surv = 1)
  args
}

icar_args <- function() {
  args <- drm_args()
  args$adj_mat <- line_adj()
  args$.toggles <- list(sp_re = "dens")
  args
}

test_that("make_data() works with valid inputs", {
  out <- do.call(make_data, drm_args())
  expect_type(out, "list")
  expect_equal(out$N, 350)
  expect_no_error(do.call(make_data, mov_args()))
  expect_no_error(do.call(make_data, surv_args()))
  expect_no_error(do.call(make_data, icar_args()))
  args <- drm_args()
  args$y[1] <- NA
  expect_no_error(do.call(make_data, args))
  idx <- rev(seq_len(nrow(sum_fl)))
  args <- drm_args()
  args[c("y", "time", "site")] <- lapply(args[c("y", "time", "site")],
                                         function(v) v[idx])
  args$x_t <- args$x_t[idx, ]
  args$x_r <- args$x_r[idx, ]
  expect_identical(do.call(make_data, args), out)
  args <- drm_args()
  args$site <- factor(args$site)
  expect_identical(do.call(make_data, args), out)
  args <- drm_args()
  args$site <- as.character(args$site)
  expect_no_error(do.call(make_data, args))
  args <- drm_args()
  args$site <- factor(args$site, levels = c(1:5, 99, 6:10))
  expect_equal(do.call(make_data, args)$n_sites, 10)
  args <- drm_args()
  args$n_ages <- 1
  expect_no_error(do.call(make_data, args))
})

## y
test_that("make_data() stops when y is missing", {
  args <- drm_args()
  args$y <- NULL
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y is NULL", {
  args <- drm_args()
  args["y"] <- list(NULL)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y is a factor", {
  args <- drm_args()
  args$y <- factor(args$y)
  expect_error(suppressWarnings(do.call(make_data, args)))
})

test_that("make_data() stops when y has a negative value", {
  args <- drm_args()
  args$y[1] <- -1
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y is not numeric", {
  skip("No input check yet")
  y <- drm_args()$y
  for (value in list(as.character(y), y > 0, as.list(y))) {
    args <- drm_args()
    args$y <- value
    expect_error(do.call(make_data, args))
  }
})

test_that("make_data() stops when y is a matrix with two columns", {
  skip("No input check yet")
  args <- drm_args()
  args$y <- cbind(args$y, args$y)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y has an infinite value", {
  skip("No input check yet")
  args <- drm_args()
  args$y[1] <- Inf
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y has the wrong length", {
  skip("No input check yet")
  args <- drm_args()
  args$y <- args$y[-1]
  expect_error(do.call(make_data, args))
  args <- drm_args()
  args$y <- c(args$y, 5)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when y is all NA", {
  skip("No input check yet")
  args <- drm_args()
  args$y <- rep(NA_real_, length(args$y))
  expect_error(do.call(make_data, args))
})

## time
test_that("make_data() stops when time is missing", {
  args <- drm_args()
  args$time <- NULL
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time is character", {
  args <- drm_args()
  args$time <- as.character(args$time)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time is a factor", {
  args <- drm_args()
  args$time <- factor(args$time)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time has an NA", {
  args <- drm_args()
  args$time[1] <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time has the wrong length", {
  args <- drm_args()
  args$time <- args$time[-1]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time is a Date", {
  skip("No input check yet")
  args <- drm_args()
  args$time <- as.Date(paste0(args$time, "-01-01"))
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time is not a whole number", {
  skip("No input check yet")
  args <- drm_args()
  args$time <- args$time * 0.5
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when time skips some years", {
  skip("No input check yet")
  args <- drm_args()
  args$time <- ifelse(args$time > 1999, args$time + 5, args$time)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when a site has the same time twice", {
  skip("No input check yet")
  args <- drm_args()
  args$time[2] <- args$time[1]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when a site-time pair is missing", {
  skip("No input check yet")
  args <- drm_args()
  args[c("y", "time", "site")] <- lapply(args[c("y", "time", "site")],
                                         function(v) v[-1])
  args$x_t <- args$x_t[-1, ]
  args$x_r <- args$x_r[-1, ]
  expect_error(do.call(make_data, args))
})

## site
test_that("make_data() stops when site is missing", {
  args <- drm_args()
  args$site <- NULL
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when site has an NA", {
  args <- drm_args()
  args$site[1] <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when site has the wrong length", {
  args <- drm_args()
  args$site <- args$site[-1]
  expect_error(do.call(make_data, args))
})

## x_t
test_that("make_data() stops when x_t has too few rows", {
  args <- drm_args()
  args$x_t <- args$x_t[-1, ]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_t is a vector", {
  args <- drm_args()
  args$x_t <- sum_fl$c_stemp
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_t has too many rows", {
  skip("No input check yet")
  args <- drm_args()
  args$x_t <- rbind(args$x_t, args$x_t)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_t has no columns", {
  skip("No input check yet")
  args <- drm_args()
  args$x_t <- args$x_t[, 0, drop = FALSE]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_t has an NA", {
  skip("No input check yet")
  args <- drm_args()
  args$x_t[1, 2] <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_t is character", {
  skip("No input check yet")
  args <- drm_args()
  args$x_t <- matrix(as.character(args$x_t), nrow = nrow(args$x_t))
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when the rows of x_t do not follow y", {
  skip("No input check yet")
  args <- drm_args()
  args$x_t <- args$x_t[rev(seq_len(nrow(args$x_t))), ]
  expect_error(do.call(make_data, args))
})

## x_r
test_that("make_data() stops when x_r has too few rows", {
  args <- drm_args()
  args$x_r <- args$x_r[-1, ]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_r is a vector", {
  args <- drm_args()
  args$x_r <- sum_fl$c_btemp
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_r has too many rows", {
  skip("No input check yet")
  args <- drm_args()
  args$x_r <- rbind(args$x_r, args$x_r)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_r has no columns", {
  skip("No input check yet")
  args <- drm_args()
  args$x_r <- args$x_r[, 0, drop = FALSE]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_r has an NA", {
  skip("No input check yet")
  args <- drm_args()
  args$x_r[1, 2] <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_r is character", {
  skip("No input check yet")
  args <- drm_args()
  args$x_r <- matrix(as.character(args$x_r), nrow = nrow(args$x_r))
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when the rows of x_r do not follow y", {
  skip("No input check yet")
  args <- drm_args()
  args$x_r <- args$x_r[rev(seq_len(nrow(args$x_r))), ]
  expect_error(do.call(make_data, args))
})

## x_m (with est_surv = 1)
test_that("make_data() stops when x_m is missing and est_surv is 1", {
  args <- surv_args()
  args$x_m <- NULL
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m has too few rows", {
  args <- surv_args()
  args$x_m <- args$x_m[-1, ]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m is a vector", {
  args <- surv_args()
  args$x_m <- sum_fl$c_stemp
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m has no columns", {
  args <- surv_args()
  args$x_m <- args$x_m[, 0, drop = FALSE]
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m has too many rows", {
  skip("No input check yet")
  args <- surv_args()
  args$x_m <- rbind(args$x_m, args$x_m)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m has an NA", {
  skip("No input check yet")
  args <- surv_args()
  args$x_m[1, 2] <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when x_m is character", {
  skip("No input check yet")
  args <- surv_args()
  args$x_m <- matrix(as.character(args$x_m), nrow = nrow(args$x_m))
  expect_error(do.call(make_data, args))
})

## f_mort
test_that("make_data() stops when f_mort has the wrong dimensions", {
  skip("No input check yet")
  args <- drm_args()
  args$f_mort <- matrix(0, nrow = 3, ncol = 35)
  expect_error(do.call(make_data, args))
  args$f_mort <- matrix(0, nrow = 2, ncol = 34)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when f_mort is a vector", {
  skip("No input check yet")
  args <- drm_args()
  args$f_mort <- rep(0, 70)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when f_mort is negative", {
  skip("No input check yet")
  args <- drm_args()
  args$f_mort <- matrix(-0.1, nrow = 2, ncol = 35)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when f_mort has NA values", {
  skip("No input check yet")
  args <- drm_args()
  args$f_mort <- matrix(NA_real_, nrow = 2, ncol = 35)
  expect_error(do.call(make_data, args))
})

## m
test_that("make_data() stops when m is negative", {
  skip("No input check yet")
  args <- drm_args()
  args$m <- -0.5
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when m has more than one value", {
  skip("No input check yet")
  args <- drm_args()
  args$m <- c(0.3, 0.4)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when m is character", {
  skip("No input check yet")
  args <- drm_args()
  args$m <- "0.3"
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when m is NA", {
  skip("No input check yet")
  args <- drm_args()
  args$m <- NA_real_
  expect_error(do.call(make_data, args))
})

## n_ages
test_that("make_data() stops when n_ages is not positive", {
  args <- drm_args()
  args$n_ages <- 0
  expect_error(do.call(make_data, args))
  args$n_ages <- -1
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when n_ages is character", {
  args <- drm_args()
  args$n_ages <- "3"
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when n_ages has more than one value", {
  args <- drm_args()
  args$n_ages <- c(2, 3)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when n_ages is not a whole number", {
  skip("No input check yet")
  args <- drm_args()
  args$n_ages <- 2.5
  expect_error(do.call(make_data, args))
})

## age_selectivity
test_that("make_data() stops when age_selectivity has the wrong length", {
  args <- drm_args()
  args$age_selectivity <- 1
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when age_selectivity is outside [0, 1]", {
  skip("No input check yet")
  args <- drm_args()
  args$age_selectivity <- c(0.5, 1.5)
  expect_error(do.call(make_data, args))
  args$age_selectivity <- c(-0.5, 1)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when age_selectivity is character", {
  skip("No input check yet")
  args <- drm_args()
  args$age_selectivity <- c("0.5", "1")
  expect_error(do.call(make_data, args))
})

## ages_movement (with movement = 1)
test_that("make_data() stops when ages_movement has the wrong length", {
  args <- mov_args()
  args$ages_movement <- c(0, 1, 1)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when ages_movement is above n_ages", {
  skip("No input check yet")
  args <- mov_args()
  args$ages_movement <- 5
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when ages_movement is not 0 or 1", {
  skip("No input check yet")
  args <- mov_args()
  args$ages_movement <- c(0, 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when ages_movement is character", {
  skip("No input check yet")
  args <- mov_args()
  args$ages_movement <- "2"
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when ages_movement is given without movement", {
  skip("No input check yet")
  args <- drm_args()
  args$ages_movement <- 2
  expect_error(do.call(make_data, args))
})

## amat and weight
test_that("make_data() stops when amat has the wrong length", {
  skip("No input check yet")
  args <- drm_args()
  args$amat <- 1
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when amat is outside [0, 1]", {
  skip("No input check yet")
  args <- drm_args()
  args$amat <- c(-1, 1)
  expect_error(do.call(make_data, args))
  args$amat <- c(0.5, 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when weight has the wrong length", {
  skip("No input check yet")
  args <- drm_args()
  args$weight <- c(1, 1, 1)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when weight is negative", {
  skip("No input check yet")
  args <- drm_args()
  args$weight <- c(-1, 1)
  expect_error(do.call(make_data, args))
})

## init_data
test_that("make_data() stops when init_data has the wrong length", {
  args <- drm_args()
  args$init_data <- c(0, 0)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when init_data is character", {
  skip("No input check yet")
  args <- drm_args()
  args$init_data <- "0"
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when init_data is NA", {
  skip("No input check yet")
  args <- drm_args()
  args$init_data <- NA_real_
  expect_error(do.call(make_data, args))
})

## adj_mat (with movement = 1 and with sp_re = "dens")
test_that("make_data() stops when adj_mat has the wrong size", {
  msg <- "'adj_mat' must have one row and one column per site"
  for (base in list(mov_args(), icar_args())) {
    args <- base
    args$adj_mat <- line_adj(9)
    expect_error(do.call(make_data, args), msg)
    args$adj_mat <- line_adj()[, -1]
    expect_error(do.call(make_data, args), msg)
  }
})

test_that("make_data() stops when adj_mat has bad values", {
  msg <- "'adj_mat' must contain finite, non-negative values only"
  for (base in list(mov_args(), icar_args())) {
    args <- base
    args$adj_mat[1, 2] <- -1
    expect_error(do.call(make_data, args), msg)
    args$adj_mat[1, 2] <- NA
    expect_error(do.call(make_data, args), msg)
  }
})

test_that("make_data() stops when adj_mat is character", {
  msg <- "'adj_mat' must contain finite, non-negative values only"
  for (base in list(mov_args(), icar_args())) {
    args <- base
    args$adj_mat <- matrix(as.character(args$adj_mat), nrow = 10)
    expect_error(do.call(make_data, args), msg)
  }
})

## areas (with movement = 1)
test_that("make_data() stops when areas is not positive", {
  args <- mov_args()
  args$areas[1] <- -1
  expect_error(do.call(make_data, args),
               "'areas' must contain finite, positive values only")
  args$areas[1] <- 0
  expect_error(do.call(make_data, args),
               "'areas' must contain finite, positive values only")
})

test_that("make_data() stops when areas has an NA", {
  args <- mov_args()
  args$areas[1] <- NA
  expect_error(do.call(make_data, args),
               "'areas' must contain finite, positive values only")
})

test_that("make_data() stops when an unnamed areas has the wrong length", {
  args <- mov_args()
  args$areas <- unname(args$areas)[-1]
  expect_error(do.call(make_data, args),
               "An unnamed 'areas' must have one value per site")
})

test_that("make_data() stops when the names of areas miss a site", {
  args <- mov_args()
  names(args$areas) <- paste0("p", names(args$areas))
  expect_error(do.call(make_data, args),
               "The names of 'areas' must be unique and include all")
  args <- mov_args()
  args$areas <- args$areas[-1]
  expect_error(do.call(make_data, args),
               "The names of 'areas' must be unique and include all")
})

test_that("make_data() stops when areas is a matrix", {
  args <- mov_args()
  args$areas <- matrix(patch_areas(), nrow = 2)
  expect_error(do.call(make_data, args), "'areas' must be a numeric vector")
  args$areas <- matrix(patch_areas(), ncol = 1)
  expect_error(do.call(make_data, args), "'areas' must be a numeric vector")
})

test_that("make_data() stops when areas is not numeric", {
  args <- mov_args()
  args$areas <- as.character(patch_areas())
  expect_error(do.call(make_data, args), "'areas' must be a numeric vector")
  args$areas <- as.list(patch_areas())
  expect_error(do.call(make_data, args), "'areas' must be a numeric vector")
})

## .toggles
test_that("make_data() stops when .toggles has an unknown name", {
  args <- drm_args()
  args$.toggles <- list(foo = 1)
  expect_error(do.call(make_data, args), "Unknown name(s): foo",
               fixed = TRUE)
})

test_that("make_data() stops when .toggles is not a list", {
  args <- drm_args()
  args$.toggles <- "movement"
  expect_error(do.call(make_data, args))
  args$.toggles <- c(rho_mu = "0")
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when .toggles is NULL", {
  args <- drm_args()
  args[".toggles"] <- list(NULL)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when minit is not 0, 1, or 2", {
  args <- drm_args()
  args$.toggles <- list(minit = 3)
  expect_error(do.call(make_data, args), "'minit' must be 0, 1, or 2.",
               fixed = TRUE)
})

test_that("make_data() stops when minit is 2 without density dependence", {
  args <- drm_args()
  args$.toggles <- list(minit = 2)
  expect_error(do.call(make_data, args),
               "'minit = 2' needs density-dependent recruitment")
})

test_that("make_data() stops when rec_dd is not a known value", {
  args <- drm_args()
  args$.toggles <- list(rec_dd = "foo")
  expect_error(do.call(make_data, args), "'rec_dd' must be one of")
})

test_that("make_data() stops when ar_re is not a known value", {
  args <- drm_args()
  args$.toggles <- list(ar_re = "foo")
  expect_error(do.call(make_data, args), "'ar_re' must be one of")
})

test_that("make_data() stops when a survival effect has est_surv = 0", {
  args <- drm_args()
  msg <- "Random effects on survival (\"surv\") need"
  args$.toggles <- list(ar_re = "surv")
  expect_error(do.call(make_data, args), msg, fixed = TRUE)
  args$.toggles <- list(iid_re = "surv")
  expect_error(do.call(make_data, args), msg, fixed = TRUE)
  args$.toggles <- list(sp_re = "surv")
  expect_error(do.call(make_data, args), msg, fixed = TRUE)
})

test_that("make_data() stops when movement is character", {
  args <- mov_args()
  args$.toggles <- list(movement = "1")
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when movement is not 0 or 1", {
  skip("No input check yet")
  args <- mov_args()
  args$.toggles <- list(movement = 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when movement is TRUE", {
  skip("No input check yet")
  args <- mov_args()
  args$.toggles <- list(movement = TRUE)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when est_surv is not 0 or 1", {
  skip("No input check yet")
  args <- surv_args()
  args$.toggles <- list(est_surv = 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when est_init is not 0 or 1", {
  skip("No input check yet")
  args <- drm_args()
  args$.toggles <- list(est_init = 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when rho_mu is not 0 or 1", {
  skip("No input check yet")
  args <- drm_args()
  args$.toggles <- list(rho_mu = 2)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when rho_mu is NA", {
  skip("No input check yet")
  args <- drm_args()
  args$.toggles <- list(rho_mu = NA)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when cloglog is character", {
  skip("No input check yet")
  args <- drm_args()
  args$.toggles <- list(cloglog = "yes")
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when acc_dd is not 0 or 1", {
  skip("No input check yet")
  args <- drm_args()
  args$.toggles <- list(acc_dd = 5)
  expect_error(do.call(make_data, args))
})

## .priors
test_that("make_data() stops when .priors has an unknown name", {
  args <- drm_args()
  args$.priors <- list(pr_foo = 1)
  expect_error(do.call(make_data, args), "Unknown name(s): pr_foo",
               fixed = TRUE)
})

test_that("make_data() stops when .priors is not a list", {
  args <- drm_args()
  args$.priors <- c(pr_alpha_mu = 0)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when a scale or shape prior is negative", {
  skip("No input check yet")
  args <- drm_args()
  args$.priors <- list(pr_alpha_sd = -1)
  expect_error(do.call(make_data, args))
  args$.priors <- list(pr_phi_a = -1)
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when a prior is character", {
  skip("No input check yet")
  args <- drm_args()
  args$.priors <- list(pr_alpha_mu = "0")
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when a scalar prior has two values", {
  skip("No input check yet")
  args <- drm_args()
  args$.priors <- list(pr_alpha_mu = c(0, 1))
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when pr_beta_t_mu is longer than K_t", {
  args <- drm_args()
  args$.priors <- list(pr_beta_t_mu = c(0, 0, 0))
  expect_error(do.call(make_data, args),
               "'pr_beta_t_mu' must have 1 or 2 values", fixed = TRUE)
})

test_that("make_data() uses one pr_beta_t_mu value for all coefficients", {
  args <- drm_args()
  args$.priors <- list(pr_beta_t_mu = 2)
  expect_equal(do.call(make_data, args)$pr_beta_t_mu, c(2, 2))
})

## family
test_that("make_data() stops when family is not a known family", {
  args <- drm_args()
  args$family <- "poisson"
  expect_error(do.call(make_data, args))
  args$family <- "Gamma"
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when family is NA", {
  args <- drm_args()
  args$family <- NA
  expect_error(do.call(make_data, args))
})

test_that("make_data() stops when family does not have one value", {
  args <- drm_args()
  args$family <- c("gamma", "lognormal")
  expect_error(do.call(make_data, args))
  args["family"] <- list(NULL)
  expect_error(do.call(make_data, args))
})
