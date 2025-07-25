% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-standalone-test-utils.R
\name{.test_multinomial}
\alias{.test_multinomial}
\title{Generate a multinomial outbreak defined by per class growth rates and a
poisson model}
\usage{
.test_multinomial(
  changes = tibble::tibble(time = c(0, 20, 40, 60, 80), variant1 = c(0.1, 0, -0.1, 0,
    0.1), variant2 = c(0.15, 0.05, -0.05, -0.01, 0.05), variant3 = c(0, 0.05, -0.05,
    +0.05, -0.05), ),
  initial = c(100, 1, 100),
  ...
)
}
\arguments{
\item{changes}{a list of time points and growth rates per week per class.}

\item{initial}{the size of the inital outbreak per class}

\item{...}{
  Arguments passed on to \code{\link[=.test_poisson_model]{.test_poisson_model}}
  \describe{
    \item{\code{seed}}{a random seed}
    \item{\code{kappa}}{dispersion parameter. 1 is no dispersion, smaller values mean more dispersion.}
    \item{\code{max_time}}{the desired length of the time series}
  }}
}
\value{
a case count time series including `class`, `count` and `time` columns
}
\description{
Generate a multinomial outbreak defined by per class growth rates and a
poisson model
}
\examples{
if (interactive()) {

ggplot2::ggplot(
  .test_multinomial(),
  ggplot2::aes(x=time,y=count,colour=class)
)+ggplot2::geom_point()

}
}
