% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-run.R
\name{fit_drm}
\alias{fit_drm}
\title{Fit the dynamic range model.}
\usage{
fit_drm(
  .data,
  y_col,
  time_col,
  site_col,
  family = "lognormal1",
  formula_zero = ~1,
  formula_rec = ~1,
  formula_surv = NULL,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4,
  seed,
  init = "cmdstan_default",
  ...
)
}
\arguments{
\item{.data}{A \code{data frame} containing the data for the model.}

\item{y_col}{A \code{character} specifying the name of the column in \code{.data}
that contains the response variable.}

\item{time_col}{A \code{character} specifying the name of the column in
\code{.data} that contains the time variable.}

\item{site_col}{A \code{character} specifying the name of the column in
\code{.data} that contains the site variable.}

\item{family}{a \code{character} specifying the family of the probability
distribution assumed for density. The options are: \itemize{
\item \code{"lognormal1"} (default): log-normal with the usual parametrization;
\item \code{"lognormal2"}: log-normal parametrized in terms of its mean;
\item \code{"gamma"}: gamma parametrized in terms of its mean;
\item \code{"loglogistic"}: log-logistic parametrized in terms of its mean.
}}

\item{formula_zero}{A \code{formula} specifying the model for the zero
inflation component. Defaults to \code{~ 1} (intercept only).}

\item{formula_rec}{A \code{formula} specifying the model for the recruitment
component. Defaults to \code{~ 1} (intercept only).}

\item{formula_surv}{A \code{formula} specifying the model for the survival
component. If \code{NULL} (the default), no survival component is included.}

\item{iter_warmup}{An \code{integer} specifying the number of warmup
iterations for the MCMC sampler. Defaults to 1000.}

\item{iter_sampling}{An \code{integer} specifying the number of sampling
iterations for the MCMC sampler. Defaults to 1000.}

\item{chains}{An \code{integer} specifying the number of MCMC chains.
Defaults to 4.}

\item{parallel_chains}{An \code{integer} specifying the number of chains to
run in parallel. Defaults to 4.}

\item{seed}{An \code{integer} specifying the random number seed.}

\item{init}{A \code{character} specifying the initialization method.  Can be
"cmdstan_default" (the default), "prior" (to initialize the model
parameters using samples from their prior) or "pathfinder".}

\item{...}{Passed on to the \code{\link[=make_data]{make_data()}} function used to build the input
\code{list} for our \code{cmdstanr} model.}
}
\value{
A \code{list} containing the MCMC draws and the model data.
Specifically: \itemize{
\item \code{draws}: The MCMC draws from the fitted model.
\item\code{data}: The data used to fit the model (as a list).
}
}
\description{
Fit the DRM Stan model (this function is not stable yet and
have not been fully tested).
}
\examples{
if (instantiate::stan_cmdstan_exists()) {
  data(sum_fl)
  fit_drm(.data = sum_fl,
          y_col = "y",
          time_col = "year",
          site_col = "patch",
          seed = 2025)$draws$summary()
}
}
\seealso{
\code{\link[=make_data]{make_data()}}

Other models: 
\code{\link{fit_sdm}()}
}
\author{
lcgodoy
}
\concept{models}
