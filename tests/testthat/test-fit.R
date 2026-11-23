## Input checks of fit_drm() and fit_sdm(): how they read .data and check
## their own arguments. test-make-data.R and test-make-data-sdm.R test the
## checks that make_data() and make_data_sdm() run, such as the checks on
## y, time, site, family and .toggles.

fits <- c("fit_drm", "fit_sdm")

## Valid arguments for fit_drm() and fit_sdm().
fit_args <- function() {
  list(.data = sum_fl, y_col = "y", time_col = "year", site_col = "patch",
       seed = 2025)
}

## Valid formulas with one covariate each, named by the formula arguments of
## each function (without formula_surv).
cov_formulas <- list(
  fit_drm = list(formula_zero = ~ c_stemp, formula_rec = ~ c_btemp),
  fit_sdm = list(formula_zero = ~ c_stemp, formula_dens = ~ c_btemp)
)

## All the formula arguments of each function.
all_formulas <- list(
  fit_drm = c("formula_zero", "formula_rec", "formula_surv"),
  fit_sdm = c("formula_zero", "formula_dens")
)

## The iid_re toggle that turns on the IID random effect in each model.
iid_on <- list(fit_drm = "rec", fit_sdm = 1)

## sum_fl with the column col set to value.
with_col <- function(col, value) {
  d <- sum_fl
  d[[col]] <- value
  d
}

test_that("fit_drm() and fit_sdm() accept valid inputs", {
  local_no_stan()
  valid_data <- list(
    sum_fl,
    as.data.frame(sum_fl),
    sum_fl[rev(seq_len(nrow(sum_fl))), ],
    with_col("y", replace(sum_fl$y, 1, NA)),
    with_col("patch", factor(sum_fl$patch)),
    with_col("patch", as.character(sum_fl$patch)),
    sf::st_as_sf(sum_fl, coords = c("lon", "lat")),
    sum_fl[sum_fl$patch == 1, ]
  )
  for (fit in fits) {
    for (d in valid_data) {
      args <- fit_args()
      args$.data <- d
      expect_reaches_stan(do.call(fit, args))
    }
    expect_reaches_stan(do.call(fit, c(fit_args(), cov_formulas[[fit]])))
    args <- fit_args()
    args$.toggles <- list(iid_re = iid_on[[fit]])
    expect_reaches_stan(do.call(fit, args))
  }
  args <- fit_args()
  args$formula_surv <- ~ c_btemp
  args$.toggles <- list(est_surv = 1)
  expect_reaches_stan(do.call(fit_drm, args))
  args <- fit_args()
  args$.toggles <- list(est_init = 1)
  expect_reaches_stan(do.call(fit_drm, args))
})

## .data
for (fit in fits) {
  test_that(paste0(fit, "() stops when .data is a matrix"), {
    local_no_stan()
    args <- fit_args()
    args$.data <- as.matrix(sum_fl)
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when .data is a list"), {
    local_no_stan()
    args <- fit_args()
    args$.data <- as.list(sum_fl)
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when .data is NULL"), {
    local_no_stan()
    args <- fit_args()
    args[".data"] <- list(NULL)
    expect_input_error(suppressWarnings(do.call(fit, args)))
  })

  test_that(paste0(fit, "() stops when .data is a list and the formulas ",
                   "have covariates"), {
    skip("No input check yet")
    local_no_stan()
    args <- c(fit_args(), cov_formulas[[fit]])
    args$.data <- as.list(sum_fl)
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when .data has no rows"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$.data <- sum_fl[0, ]
    expect_input_error(suppressWarnings(do.call(fit, args)))
  })
}

## y_col, time_col and site_col
cols <- c(y_col = "y", time_col = "year", site_col = "patch")
for (fit in fits) {
  for (col in names(cols)) {
    test_that(paste0(fit, "() stops when ", col, " is not in .data"), {
      local_no_stan()
      args <- fit_args()
      args[[col]] <- "nope"
      expect_input_error(suppressWarnings(do.call(fit, args)))
    })

    test_that(paste0(fit, "() stops when ", col, " is not one name"), {
      local_no_stan()
      for (value in list(NULL, c(cols[[col]], "dens"), NA_character_)) {
        args <- fit_args()
        args[col] <- list(value)
        expect_input_error(do.call(fit, args))
      }
    })

    test_that(paste0(fit, "() stops when ", col, " is missing"), {
      local_no_stan()
      args <- fit_args()
      args[[col]] <- NULL
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", col, " is a column number"), {
      skip("No input check yet")
      local_no_stan()
      args <- fit_args()
      args[[col]] <- match(cols[[col]], names(sum_fl))
      expect_input_error(do.call(fit, args))
    })
  }

  test_that(paste0(fit, "() stops when site_col is the same as time_col"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$site_col <- "year"
    expect_input_error(do.call(fit, args))
  })
}

## Formulas
for (fit in fits) {
  for (arg in all_formulas[[fit]]) {
    test_that(paste0(fit, "() stops when ", arg, " uses a column not in ",
                     ".data"), {
      local_no_stan()
      args <- fit_args()
      args[[arg]] <- ~ nope
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " is a string"), {
      local_no_stan()
      for (value in c("~ c_stemp", "c_stemp")) {
        args <- fit_args()
        args[[arg]] <- value
        expect_input_error(do.call(fit, args))
      }
    })

    test_that(paste0(fit, "() stops when ", arg, " uses a covariate with ",
                     "NaN"), {
      local_no_stan()
      args <- fit_args()
      args[[arg]] <- ~ ssalin
      expect_input_error(do.call(fit, args))
    })
  }

  for (arg in names(cov_formulas[[fit]])) {
    test_that(paste0(fit, "() stops when ", arg, " is NULL"), {
      local_no_stan()
      args <- fit_args()
      args[arg] <- list(NULL)
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " uses a short variable ",
                     "outside .data"), {
      local_no_stan()
      ext_short <- seq_len(10)
      args <- fit_args()
      args[[arg]] <- ~ ext_short
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " is two-sided"), {
      skip("No input check yet")
      local_no_stan()
      args <- fit_args()
      args[[arg]] <- y ~ c_stemp
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " uses a variable outside ",
                     ".data"), {
      skip("No input check yet")
      local_no_stan()
      ext_cov <- seq_len(nrow(sum_fl))
      args <- fit_args()
      args[[arg]] <- ~ ext_cov
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " has an offset"), {
      skip("No input check yet")
      local_no_stan()
      args <- fit_args()
      args[[arg]] <- ~ offset(log(area_km2))
      expect_input_error(do.call(fit, args))
    })

    test_that(paste0(fit, "() stops when ", arg, " uses the response"), {
      skip("No input check yet")
      local_no_stan()
      args <- fit_args()
      args[[arg]] <- ~ y
      expect_input_error(do.call(fit, args))
    })
  }
}

test_that("fit_drm() stops when est_surv is on and formula_surv is NULL", {
  local_no_stan()
  args <- fit_args()
  args$.toggles <- list(est_surv = 1)
  expect_input_error(do.call(fit_drm, args))
})

test_that("fit_drm() stops when formula_surv is given and est_surv is off", {
  skip("No input check yet")
  local_no_stan()
  args <- fit_args()
  args$formula_surv <- ~ c_btemp
  expect_input_error(do.call(fit_drm, args))
})

## seed
for (fit in fits) {
  test_that(paste0(fit, "() stops when seed is missing"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$seed <- NULL
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when seed is not a number"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$seed <- "abc"
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when seed has length 2"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$seed <- c(1, 2)
    expect_input_error(do.call(fit, args))
  })
}

## init
for (fit in fits) {
  test_that(paste0(fit, "() stops when init is an unknown string"), {
    local_no_stan()
    args <- fit_args()
    args$init <- "foo"
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when init has more than one value"), {
    local_no_stan()
    for (value in list(c(1, 2), rep(list(list(beta_t = 0)), 4))) {
      args <- fit_args()
      args$init <- value
      expect_input_error(do.call(fit, args))
    }
  })

  test_that(paste0(fit, "() stops when init is negative"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$init <- -1
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when init is a list"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$init <- list(list(beta_t = 0))
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when init is a function"), {
    local_no_stan()
    args <- fit_args()
    args$init <- function() list(beta_t = 0)
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when init is NA"), {
    local_no_stan()
    args <- fit_args()
    args$init <- NA
    expect_input_error(do.call(fit, args))
  })
}

## algorithm
test_that("fit_drm() stops when algorithm does not have length 1", {
  local_no_stan()
  for (value in list(c("nuts", "vb"), NULL)) {
    args <- fit_args()
    args["algorithm"] <- list(value)
    expect_input_error(do.call(fit_drm, args))
  }
})

test_that("fit_sdm() stops when algorithm does not have length 1", {
  local_no_stan()
  for (value in list(c("nuts", "vb"), NULL)) {
    args <- fit_args()
    args["algorithm"] <- list(value)
    expect_input_error(do.call(fit_sdm, args))
  }
})

for (fit in fits) {
  test_that(paste0(fit, "() stops when algorithm is unknown"), {
    local_no_stan()
    args <- fit_args()
    args$algorithm <- "foo"
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when algorithm is NA"), {
    local_no_stan()
    args <- fit_args()
    args$algorithm <- NA
    expect_input_error(do.call(fit, args))
  })
}

## algo_args
for (fit in fits) {
  test_that(paste0(fit, "() stops when algo_args has an unknown name"), {
    local_no_stan()
    args <- fit_args()
    args$algo_args <- list(foo = 1)
    expect_input_error(do.call(fit, args), "Unknown name(s): foo",
                       fixed = TRUE)
  })

  test_that(paste0(fit, "() stops when algo_args is not a list"), {
    local_no_stan()
    for (value in list("foo", c(chains = 2))) {
      args <- fit_args()
      args$algo_args <- value
      expect_input_error(do.call(fit, args))
    }
  })

  test_that(paste0(fit, "() stops when algo_args has an unnamed value"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$algo_args <- list(500)
    expect_input_error(do.call(fit, args))
  })
}

## Arguments passed on to make_data() and make_data_sdm()
for (fit in fits) {
  test_that(paste0(fit, "() stops when ... has an unknown argument"), {
    local_no_stan()
    args <- fit_args()
    args$foo <- 1
    expect_input_error(do.call(fit, args))
  })

  test_that(paste0(fit, "() stops when iid_re is on with a single site"), {
    skip("No input check yet")
    local_no_stan()
    args <- fit_args()
    args$.data <- sum_fl[sum_fl$patch == 1, ]
    args$.toggles <- list(iid_re = iid_on[[fit]])
    expect_input_error(do.call(fit, args))
  })
}

test_that("fit_sdm() stops when ... has a DRM-only argument", {
  local_no_stan()
  args <- fit_args()
  args$n_ages <- 3
  expect_input_error(do.call(fit_sdm, args))
})

test_that("fit_drm() stops when an argument in ... is a partial name", {
  skip("No input check yet")
  local_no_stan()
  args <- fit_args()
  args$n_age <- 3
  expect_input_error(do.call(fit_drm, args))
})

test_that("fit_drm() stops when init_data is given and est_init is on", {
  skip("No input check yet")
  local_no_stan()
  args <- fit_args()
  args$.toggles <- list(est_init = 1)
  args$init_data <- -1
  expect_input_error(do.call(fit_drm, args))
})
