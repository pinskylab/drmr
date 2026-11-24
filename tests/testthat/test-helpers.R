## Six square patches on a 3 x 2 grid, as an sf object.
sf_patches <- function() {
  box <- rbind(c(0, 0), c(3, 0), c(3, 2), c(0, 2), c(0, 0))
  box <- sf::st_sfc(sf::st_polygon(list(box)))
  sf::st_sf(id = 1:6, geometry = sf::st_make_grid(box, n = c(3, 2)))
}

test_that("the helpers work with valid inputs", {
  adj <- gen_adj(sf_patches())
  expect_equal(dim(adj), c(6, 6))
  expect_true(all(adj %in% c(0, 1)))
  expect_message(adj <- gen_adj(1:4), "Assuming patches are 1d")
  expect_equal(adj, line_adj(4), ignore_attr = TRUE)
  expect_equal(default_algo("nuts", list(chains = 2))$chains, 2)
  dat <- do.call(make_data, drm_args())
  expect_length(prior_inits(dat, chains = 2), 2)
  expect_type(prior_sample(dat), "list")
  expect_no_error(check_pars(dat, "drm"))
  dsdm <- do.call(make_data_sdm, sdm_args())
  expect_length(prior_inits(dsdm, chains = 2, model = "sdm"), 2)
})

## gen_adj()
test_that("gen_adj() stops when x is a data.frame without geometry", {
  skip("No input check yet")
  x <- sf::st_drop_geometry(sf_patches())
  x$x <- 1:6
  x$y <- 1
  expect_error(suppressMessages(gen_adj(x)))
})

test_that("gen_adj() stops when x is a matrix of coordinates", {
  skip("No input check yet")
  x <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(sf_patches())))
  expect_error(suppressMessages(gen_adj(x)))
})

test_that("gen_adj() stops when the geometries are not polygons", {
  x <- sf_patches()
  pts <- sf::st_sf(id = x$id,
                   geometry = sf::st_centroid(sf::st_geometry(x)))
  expect_error(gen_adj(pts))
  lines <- sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(1, 1))),
                      sf::st_linestring(rbind(c(1, 1), c(2, 0))))
  expect_error(gen_adj(sf::st_sf(id = 1:2, geometry = lines)))
})

test_that("gen_adj() stops when there are no patches", {
  expect_error(gen_adj(sf_patches()[0, ]))
  expect_error(suppressMessages(gen_adj(NULL)))
  expect_error(suppressMessages(gen_adj(integer(0))))
})

## default_algo()
test_that("default_algo() stops when algorithm is unknown", {
  expect_error(default_algo("hmc"))
  expect_error(default_algo("NUTS"))
})

test_that("default_algo() stops when algorithm is not a single value", {
  expect_error(default_algo(NULL))
  expect_error(default_algo(c("nuts", "vb")))
})

test_that("default_algo() stops when algorithm is a number", {
  skip("No input check yet")
  expect_error(default_algo(2))
})

test_that("default_algo() stops when algo_args has an unknown name", {
  expect_error(default_algo("nuts", list(chain = 2)),
               "Unknown name(s): chain.", fixed = TRUE)
})

test_that("default_algo() stops when algo_args is not a list", {
  expect_error(default_algo("nuts", c(chains = 2)))
  expect_error(default_algo("nuts", NULL))
})

test_that("default_algo() stops when algo_args has an unnamed value", {
  skip("No input check yet")
  expect_error(default_algo("nuts", list(2)))
})

test_that("default_algo() stops when an algo_args value is the wrong type", {
  skip("No input check yet")
  expect_error(default_algo("nuts", list(chains = "a")))
})

## prior_inits(), prior_sample(), and check_pars(): chains
test_that("prior_inits() stops when chains is missing", {
  dat <- do.call(make_data, drm_args())
  expect_error(prior_inits(dat))
})

test_that("prior_inits() stops when chains is negative", {
  dat <- do.call(make_data, drm_args())
  expect_error(prior_inits(dat, chains = -1))
})

test_that("prior_inits() stops when chains is not numeric", {
  dat <- do.call(make_data, drm_args())
  expect_error(suppressWarnings(prior_inits(dat, chains = "a")))
})

test_that("prior_inits() stops when chains is zero", {
  skip("No input check yet")
  dat <- do.call(make_data, drm_args())
  expect_error(prior_inits(dat, chains = 0))
})

test_that("prior_inits() stops when chains is not a whole number", {
  skip("No input check yet")
  dat <- do.call(make_data, drm_args())
  expect_error(prior_inits(dat, chains = 2.7))
})

## prior_sample() and prior_pars(): phi
test_that("the prior draws of phi are not truncated for the log-logistic", {
  ## with the default Gamma(2, 1) prior, P(phi < 1) is about 0.26
  set.seed(2026)
  dat <- do.call(make_data, c(drm_args(), list(family = "loglogistic")))
  expect_true(any(replicate(200, prior_sample(dat)$phi) < 1))
  expect_true(any(replicate(200, prior_pars(dat)$phi) < 1))
})

## dat
test_that("the helpers stop when dat is a data.frame", {
  expect_error(prior_inits(sum_fl, chains = 2))
  expect_error(prior_sample(sum_fl))
  expect_error(check_pars(sum_fl, "drm"))
})

test_that("prior_inits() stops when dat is a matrix", {
  expect_error(prior_inits(as.matrix(sum_fl[1:3]), chains = 2))
})

test_that("prior_inits() stops when dat is empty", {
  expect_error(suppressWarnings(prior_inits(list(), chains = 2)))
  expect_error(suppressWarnings(prior_inits(NULL, chains = 2)))
})

test_that("prior_inits() stops when dat is for the other model", {
  dat <- do.call(make_data, drm_args())
  dsdm <- do.call(make_data_sdm, sdm_args())
  expect_error(prior_inits(dsdm, chains = 2))
  expect_error(prior_inits(dat, chains = 2, model = "sdm"))
})

test_that("prior_inits() stops when dat has an extra element", {
  dat <- do.call(make_data, drm_args())
  dat$extra <- 1
  expect_error(prior_inits(dat, chains = 2))
})

test_that("prior_inits() stops when dat lacks an element", {
  dat <- do.call(make_data, drm_args())
  dat$K_r <- NULL
  expect_error(suppressWarnings(prior_inits(dat, chains = 2)))
})

## model
test_that("prior_inits() stops when model is not drm or sdm", {
  skip("No input check yet")
  dsdm <- do.call(make_data_sdm, sdm_args())
  expect_error(prior_inits(dsdm, chains = 2, model = "foo"))
})

test_that("check_pars() stops when model is in the wrong case", {
  dat <- do.call(make_data, drm_args())
  expect_error(check_pars(dat, "DRM"))
})

test_that("check_pars() stops when model is not a single value", {
  dat <- do.call(make_data, drm_args())
  expect_error(check_pars(dat, c("drm", "sdm")))
  expect_error(check_pars(dat, NULL))
})

test_that("check_pars() stops when model is missing", {
  dat <- do.call(make_data, drm_args())
  expect_error(check_pars(dat))
})
