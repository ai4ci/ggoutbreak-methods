% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-standalone-file-output-utils.R
\name{.convert_pdf_to_pngs}
\alias{.convert_pdf_to_pngs}
\title{Convert a single pdf to multiple pngs one per page}
\usage{
.convert_pdf_to_pngs(pdfFile)
}
\arguments{
\item{pdfFile}{an input pdf file}
}
\value{
a list of png files one per page.
}
\description{
Convert a single pdf to multiple pngs one per page
}
\examples{
try({
 hux = iris \%>\% huxtable::as_hux() \%>\% huxtable::theme_mondrian(font="Arial")
 html = hux \%>\% huxtable::to_html()
 tmp = .print_html_with_chrome(html,maxWidth = std_size$A4$width,maxHeight = std_size$A4$height)
 pngs = .convert_pdf_to_pngs(tmp)
})
}
\concept{output}
\keyword{internal}
