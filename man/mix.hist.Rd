\name{mix.hist}
\alias{mix.hist}
\title{ Estimated densities }
\description{
Plot the histogram along with the estimated density (univariate case, p=1).
}
\usage{
mix.hist(y, model, breaks, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response vector }
  \item{model}{ a variable returned by \code{\link{smsn.mix}} }
  \item{breaks}{ the same option in \cite{\link{histogram}} }
  \item{...}{further arguments to \code{\link{hist}}}
}

\examples{ ## see \code{\link{bmi}} and \code{\link{smsn.mix}}
}

\seealso{ \code{\link{smsn.mix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }

