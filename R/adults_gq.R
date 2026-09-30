##' @title Density of adults in reproductive age
##'
##' @description Computes posterior samples of the density (or abundance,
##'   depending on the response variable) of **adults in reproductive age** at
##'   each site and time point used to fit an \code{adrm} object. When the
##'   weight at age (\code{weight}) is provided at the time of fitting, the
##'   output is the **biomass** of adults in reproductive age instead, a
##'   quantity known as the Spawning Stock Biomass (SSB) in fisheries science.
##'
##' @details For each site \eqn{s} and time \eqn{t}, the adults in reproductive
##'   age are given by \deqn{\sum_{a = 1}^{A} m_a w_a \lambda_{a, t, s},} where
##'   \eqn{\lambda_{a, t, s}} is the expected density at age \eqn{a} (see
##'   [ages_edens()]), \eqn{m_a} is the maturity at age (\code{amat}), and
##'   \eqn{w_a} is the weight at age (\code{weight}). If \code{amat} was not
##'   provided at the time of fitting, all ages are assumed to be in
##'   reproductive age. Similarly, if \code{weight} was not provided,
##'   \eqn{w_a = 1} for all ages. See [make_data()] for details.
##'
##'   When density dependence is enabled, this is the quantity driving
##'   recruitment: the adults at time \eqn{t} produce the recruits at time
##'   \eqn{t + 1}. If recruitment depends on the adults from the whole region
##'   (\code{acc_dd = 1}), the relevant quantity at each time point is the sum
##'   of the adults across sites. Contrasting these samples with the curves
##'   from [dd_curves()] helps assess the range of densities (or biomass)
##'   informing the density-dependence relationship.
##'
##'   The adults are computed from the latent age-specific densities.
##'   Therefore, they do not account for the selectivity at age, the
##'   observation error, or random effects acting directly on the density
##'   (i.e., \code{"dens"}).
##'
##' @param object An object of class \code{adrm}, typically the output of
##'   [fit_drm()].
##' @param cores number of threads used to compute the samples. If four chains
##'   were used in the \code{fit_drm}, then four (or less) threads are
##'   recommended.
##' @param ... additional parameters to be passed to
##'   \code{$generate_quantities}
##'
##' @return An object of class \code{pred_drmr}, which is a \code{list} with
##'   the elements: \itemize{ \item \code{gq}: a \code{"CmdStanGQ"} object
##'   containing the posterior samples of the adults in reproductive age
##'   (variable \code{adults}).  \item \code{spt}: a \code{data.frame}
##'   identifying the site and time associated with each element of
##'   \code{adults}.  } Its \code{summary} method returns a \code{data.frame}
##'   with posterior summaries for each site and time point.
##'
##' @seealso [dd_curves()], [ages_edens()], [make_data()]
##' @examples
##' \donttest{
##' if (instantiate::stan_cmdstan_exists()) {
##'   data(sum_fl)
##'   fit <- fit_drm(.data = sum_fl,
##'                  y_col = "y",
##'                  time_col = "year",
##'                  site_col = "patch",
##'                  .toggles = list(rec_dd = "ricker"),
##'                  seed = 2025)
##'   summary(get_adults(fit))
##' }
##' }
##' @name get_adults
##' @export
##' @author lcgodoy
get_adults <- function(object, ...) {
  UseMethod("get_adults")
}

##' @rdname get_adults
##' @export
get_adults.adrm <- function(object,
                            cores = 1,
                            ...) {
  stopifnot(inherits(object$stanfit, c("CmdStanFit", "CmdStanLaplace",
                                       "CmdStanPathfinder", "CmdStanVB")))
  ##--- pars from model fitted ----
  pars <- fitted_pars_lambda(c(object$data, list(proj = 0)))
  fitted_params <-
    object$stanfit$draws(variables = pars)
  ## load compiled model
  adults_comp <-
    instantiate::stan_package_model(name = "adults_drm",
                                    package = "drmr")
  ## computing the adults
  gq <- adults_comp$
    generate_quantities(fitted_params = fitted_params,
                        data = object$data,
                        parallel_chains = cores,
                        ...)
  spt <- data.frame(v1 = object$cols$site_levels[object$data$site],
                    v2 = object$data$time + object$data$time_init - 1)
  colnames(spt) <- rev(unname(unlist(object$cols[2:3])))
  output <- list("gq" = gq,
                 "spt" = spt)
  return(new_pred_drmr(output))
}
