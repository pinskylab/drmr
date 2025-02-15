% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict_drm}
\alias{predict_drm}
\title{Forecasts based on DRM.}
\usage{
predict_drm(
  drm,
  drm_data,
  ntime_for,
  x_tt,
  x_rt,
  x_mt = matrix(1, ncol = 1),
  x_mpast = matrix(1, ncol = 1),
  f_test,
  seed = 1,
  cores = 1
)
}
\arguments{
\item{drm}{A \code{CmdStanFit} object containing samples from the posterior
distribution.}

\item{drm_data}{a \code{list} used as input for model fitting. Typically,
the output from the \link{make_data} function.}

\item{ntime_for}{an \code{integer} denoting the number of timepoints for the
forecast.}

\item{x_tt}{a design \code{matrix} of variables associated to the
probability of absence at each site/time.}

\item{x_rt}{a design \code{matrix} of variables associated to recruitment.}

\item{x_mt}{a design \code{matrix} of variables associated to survival.}

\item{x_mpast}{a design \code{matrix} of variables associated to survival
from the last year of the "training" perior (i.e., the period using to fit
the model and obtain the parameters' estimates).}

\item{f_test}{a \code{matrix} informing the instantaneous fishing mortality
rates at each age (columns) and timepoint (rows).}

\item{seed}{a seed used for the forecasts. Forecasts are obtained through
Monte Carlo samples from the posterior predictive distribution. Therefore,
a \code{seed} is needed to ensure the results' reproducibility.}

\item{cores}{number of threads used for the forecast. If four chains were
used in the \code{drm}, then four (or less) threads are recommended.}
}
\value{
an object of class \code{"CmdStanGQ"} containing samples for the
posterior predictive distribution for forecasting.
}
\description{
Considering a new dataset (across the same patches), computes
forecasts based on the DRM passed as \code{drm}.
}
\details{
It is important that the rows of design matrices \code{x_tt},
\code{x_rt} and \code{x_mt} are associated to the same patch/site and
timepoints. In addition, the current version of the code assumes the data
where forecasts are needed is ordered by "patch" and "site" and, in
addition, its patches MUST be the same as the ones used to obtain the
parameters' estimates from the the \code{drm} object.
}
\author{
lcgodoy
}
