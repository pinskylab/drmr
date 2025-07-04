% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make-data.R
\name{make_data_sdm}
\alias{make_data_sdm}
\title{Make data for SDM stan models}
\usage{
make_data_sdm(
  y,
  time,
  site,
  z,
  x,
  .toggles,
  .priors,
  family = "gamma",
  reorder = TRUE,
  phi_hat = FALSE
)
}
\arguments{
\item{y}{a \code{numeric vector} of species' densities.}

\item{time}{an \code{vector} indicating the time point associated to each
element of \code{y}.}

\item{site}{an \code{vector} indicating the sites associated to each element
of \code{y}.}

\item{z}{a design \code{matrix} of variables associated to the probability
of absence at each site/time.}

\item{x}{a design \code{matrix} of variables associated to the non-zero
densities.}

\item{.toggles}{a \code{list} of toggles for model components. The
components are: \itemize{\item \code{rho_mu}: 1 to use
explicitly relate rho to mu and 0 otherwise. \item \code{cloglog}: 1 to
use the complementary log-log and 0 for the logit link function for the
absence probabilities.  \item \code{movement}: 1 to allow for (adjacent)
moviment; 0 for static.  \item \code{est_surv}: 1 to estimate survival
rates and 0 otherwise.  \item \code{ar_re}: "rec" to incorporate an AR(1)
process density. The only other accepted option is "none"}}

\item{.priors}{a \code{list} of priors hyperparameters.}

\item{family}{a \code{character} specifying the family of the probability
distribution assumed for density. The options are: \itemize{ \item
\code{"gamma"} (default): gamma parametrized in terms of its mean; \item
\code{"lognormal"}: log-normal parametrized in terms of its mean; \item
\code{"loglogistic"}: log-logistic parametrized in terms of its mean.
\item \code{"lognormal_legacy"} (default): log-normal with its usual
parametrization; }}

\item{reorder}{a \code{boolean} telling whether the data needs to be
reordered. The default is TRUE and means the data points will be ordered
by site and time, respectively.}

\item{phi_hat}{a \code{boolean} indicating whether the prior on \code{phi}
should be determined through the data.}
}
\value{
a \code{list} to be used as the input for a \code{stan} model
}
\description{
This function creates the \code{list} used as the input for the \code{stan}
model.
}
\author{
lcgodoy
}
