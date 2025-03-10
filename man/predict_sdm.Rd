% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict_sdm}
\alias{predict_sdm}
\title{Forecasts based on SDM.}
\usage{
predict_sdm(sdm, sdm_data, ntime_for, z_t, x_t, time_for, seed = 1, cores = 1)
}
\arguments{
\item{sdm}{A \code{CmdStanFit} object containing samples from the posterior
distribution.}

\item{sdm_data}{a \code{list} used as input for model fitting. Typically,
the output from the \link{make_data} function.}

\item{ntime_for}{an \code{integer} denoting the number of timepoints for the
forecast.}

\item{z_t}{a design \code{matrix} of variables associated to the probability
of absence at each site/time.}

\item{x_t}{a design \code{matrix} of variables associated to the non-zero
densities.}

\item{time_for}{an \code{integer vector} indicating timepoints for
forecasting.}

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
forecasts based on the SDM passed as \code{sdm}.

Consider a linear predictor having linear and square terms
associated with a variable \eqn{x}. Assume this variable was centered
before being included in the linear predictor. This functions returns the
value of \eqn{x} (on its original scale) such that the linear predictor is
maximized (or minimized).
}
\details{
It is important ethat the rows of design matrices \code{z_t} and
\code{x_t} are associated to the same patch/site and timepoints. In
addition, the current version of the code assumes the data where forecasts
are needed is ordered by "patch" and "site" and, in addition, its patches
MUST be the same as the ones used to obtain the parameters' estimates from
the the \code{sdm} object.
}
\author{
lcgodoy
}
