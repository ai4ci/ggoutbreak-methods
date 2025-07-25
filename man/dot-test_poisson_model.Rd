% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-standalone-test-utils.R
\name{.test_poisson_model}
\alias{.test_poisson_model}
\title{Generate an outbreak case count series defined by growth rates using a poisson model.}
\usage{
.test_poisson_model(
  changes = tibble::tibble(time = c(0, 20, 40, 60, 80), r = c(0.1, 0, -0.1, 0, 0.1)),
  kappa = 1,
  initial = 100,
  max_time = 104,
  seed = Sys.time()
)
}
\arguments{
\item{changes}{a dataframe holding change time points (`time`) and
growth rate per week (`r`) columns}

\item{kappa}{dispersion parameter. 1 is no dispersion, smaller values mean more dispersion.}

\item{initial}{the size of the initial outbreak}

\item{max_time}{the desired length of the time series}

\item{seed}{a random seed}
}
\value{
A dataframe containing the following columns:

* count (positive_integer) - Positive case counts associated with the specified timeframe
* time (ggoutbreak::time_period + group_unique) - A (usually complete) set of singular observations per unit time as a `time_period`

Ungrouped.
}
\description{
Generate an outbreak case count series defined by growth rates using a poisson model.
}
\examples{
if (interactive()) {

ggplot2::ggplot(
  .test_poisson_model(kappa=0.1, seed=100),
  ggplot2::aes(x=time,y=count)
)+
ggplot2::geom_point()

}
}
\keyword{internal}
