\name{mix.dens}
\alias{mix.dens}
\title{ Estimated densities }
\description{
Plot the estimated density or log-density (univariate case, p=1).
}
\usage{
mix.dens(y, model, log=FALSE, ylab=NULL, xlab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response vector }
  \item{model}{ a variable returned by \code{\link{smsn.mix}} }
  \item{log}{ Logical, plot log-density if TRUE (default = FALSE) }
  \item{ylab}{ Title of the ylab, if NULL default is selected }
  \item{xlab}{ Title of the xlab, if NULL default is selected }
  \item{main}{ Main Title, if NULL default is selected }
  \item{...}{further arguments to \code{\link{plot}}}
}

\examples{ ## see \code{\link{bmi}} and \code{\link{smsn.mix}}
}

\seealso{ \code{\link{smsn.mix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }

