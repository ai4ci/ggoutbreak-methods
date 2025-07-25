% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-standalone-ggplot-utils.R
\name{.gg_transform_logit}
\alias{.gg_transform_logit}
\title{logit scale}
\usage{
.gg_transform_logit(n = 5, ...)
}
\arguments{
\item{n}{number of breas}

\item{...}{not used}
}
\value{
A scales object
}
\description{
Perform logit scaling with right axis formatting. To not be
used directly but with ggplot (e.g. ggplot2::scale_y_continuous(trans =
"logit")
}
\examples{
try({
tibble::tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) \%>\%
 ggplot2::ggplot(ggplot2::aes(fold_change , pvalue)) +
 ggplot2::geom_point() +
 ggplot2::scale_y_continuous(transform= .gg_transform_logit())
})
}
\concept{ggplot}
\keyword{internal}
