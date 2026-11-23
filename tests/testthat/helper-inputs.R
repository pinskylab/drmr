## Valid inputs and expectations that all the test files use. Each test
## starts from a valid input and changes one argument.

## Valid arguments for make_data(). sum_fl has one row per patch and year,
## sorted by patch and then by year.
drm_args <- function() {
  list(y = sum_fl$y,
       time = sum_fl$year,
       site = sum_fl$patch,
       x_t = stats::model.matrix(~ c_stemp, data = sum_fl),
       x_r = stats::model.matrix(~ c_btemp, data = sum_fl))
}

## Valid arguments for make_data_sdm().
sdm_args <- function() {
  list(y = sum_fl$y,
       time = sum_fl$year,
       site = sum_fl$patch,
       x_t = stats::model.matrix(~ c_stemp, data = sum_fl),
       x = stats::model.matrix(~ c_btemp, data = sum_fl))
}

## Adjacency matrix for n sites along a line.
line_adj <- function(n = 10) {
  adj <- matrix(0, nrow = n, ncol = n)
  adj[cbind(seq_len(n - 1), seq_len(n)[-1])] <- 1
  adj + t(adj)
}

## Named area of each patch in sum_fl.
patch_areas <- function() {
  areas <- tapply(sum_fl$area_km2, sum_fl$patch, unique)
  stats::setNames(as.numeric(areas), names(areas))
}

## Replaces instantiate::stan_package_model() with a stub that stops with
## an error of class "drmr_test_stan". So fit_drm() and fit_sdm() stop
## after all their input checks, and no test runs Stan.
local_no_stan <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    stan_package_model = function(...) {
      rlang::abort("The inputs passed all the checks.",
                   class = "drmr_test_stan")
    },
    .package = "instantiate",
    .env = env
  )
}

## Expects that a call passes all the input checks and reaches Stan. Use it
## after local_no_stan().
expect_reaches_stan <- function(object) {
  testthat::expect_error(object, class = "drmr_test_stan")
}

## Expects that an input check stops the call, before it reaches Stan.
## regexp is matched against the error message, if given.
expect_input_error <- function(object, regexp = NULL, ...) {
  err <- rlang::inject(
    testthat::expect_error(!!rlang::enquo(object), regexp = regexp, ...)
  )
  testthat::expect(!inherits(err, "drmr_test_stan"),
                   "The call passed all the input checks and reached Stan.")
  invisible(err)
}
