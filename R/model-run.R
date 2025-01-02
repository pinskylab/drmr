##' @title Fit the dynamic range model.
##' @export
##' @family models
##' @description Fit the DRM Stan model and return posterior summary.
##' @return A data frame of posterior summaries.
##' @param data_list a \code{list} created using the [make_data()] function.
##' @param ... Passed on to the \code{sample()} method of CmdStan model objects:
##'   <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
##' @seealso [make_data()]
##' @examples
##' if (instantiate::stan_cmdstan_exists()) {
##'   data(sum_fl)
##'   make_data(y = sum_fl$y,
##'             time = sum_fl$year,
##'             site = sum_fl$patch,
##'             n_ages = 8,
##'             age_at_maturity = integer(0),
##'             m = 0.25) |>
##'   run_drm()
##' }
run_drm <- function(data_list, ...) {
  model <- instantiate::stan_package_model(
    name = "drm",
    package = "drmr"
  )
  fit <- model$sample(data = data_list, ...)
  fit$summary()
}
