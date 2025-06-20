% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-gen.R
\name{sim_log_rec}
\alias{sim_log_rec}
\title{Simulate log-recruitment}
\usage{
sim_log_rec(n_patches, n_time, x_rec, pars, ar_time = TRUE)
}
\arguments{
\item{n_patches}{number of patches}

\item{n_time}{number of timepoints}

\item{x_rec}{matrix of environmental factors affecting recruitment}

\item{pars}{a named \code{list} of parameters used to simulate
log-recruitment. It must contain a vector named \code{"beta_r"} with
length equal to the number of columnts in \code{x_rec}. In addition, if
\code{ar_time = TRUE}, the list must also contain a named \code{vector}
called "ar". This named \code{vector} must contain an element called
\code{alpha} (the autocorrelation parameter) and another called \code{tau}
(the conditional SD).}

\item{ar_time}{a \code{boolean} indicating whether an AR(1) term should be
included to the log-recruitment.}
}
\value{
a \code{matrix} with \code{n_patches} columns and \code{n_time} rows
representing the log-recruitment at each patch/site and time.
}
\description{
Simulate log-recruitment
}
\author{
lcgodoy
}
