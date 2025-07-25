% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-standalone-file-output-utils.R
\name{.is_html_output}
\alias{.is_html_output}
\title{Check the current input and output type}
\usage{
.is_html_output(fmt = pandoc_to(), excludes = NULL)
}
\arguments{
\item{fmt}{A character vector of output formats to be checked against. If not
provided, \code{is_html_output()} uses \code{pandoc_to()}, and
\code{pandoc_to()} returns the output format name.}

\item{excludes}{A character vector of output formats that should not be
considered as HTML format. Options are: markdown, epub, epub2, html, html4, html5,
revealjs, s5, slideous, slidy, and gfm.}
}
\description{
The function \code{is_latex_output()} returns \code{TRUE} when the output
format is LaTeX; it works for both \file{.Rnw} and R Markdown documents (for
the latter, the two Pandoc formats \code{latex} and \code{beamer} are
considered LaTeX output). The function \code{is_html_output()} only works for
R Markdown documents and will test for several Pandoc HTML based output
formats (by default, these formats are considered as HTML formats:
\code{c('markdown', 'epub', 'epub2', 'html', 'html4', 'html5', 'revealjs', 's5',
'slideous', 'slidy', 'gfm')}).
}
\details{
The function \code{pandoc_to()} returns the Pandoc output format, and
\code{pandoc_from()} returns Pandoc input format. \code{pandoc_to(fmt)}
allows to check the current output format against a set of format names. Both
are to be used with R Markdown documents.

These functions may be useful for conditional output that depends on the
output format. For example, you may write out a LaTeX table in an R Markdown
document when the output format is LaTeX, and an HTML or Markdown table when
the output format is HTML. Use \code{pandoc_to(fmt)} to test a more specific
Pandoc format.

Internally, the Pandoc output format of the current R Markdown document is
stored in \code{knitr::\link[knitr]{opts_knit}$get('rmarkdown.pandoc.to')}, and the
Pandoc input format in
\code{knitr::\link[knitr]{opts_knit}$get('rmarkdown.pandoc.from')}
}
\note{
See available Pandoc formats, in
  \href{https://pandoc.org/MANUAL.html}{Pandoc's Manual}
}
\examples{
# check for output formats type
knitr::is_latex_output()
knitr::is_html_output()
knitr::is_html_output(excludes = c("markdown", "epub"))
# Get current formats
knitr::pandoc_from()
knitr::pandoc_to()
# Test if current output format is 'docx'
knitr::pandoc_to("docx")
}
\concept{output}
\keyword{internal}
