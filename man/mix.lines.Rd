\name{mix.lines}
\alias{mix.lines}
\title{ Plot lines of smsn densities }
\description{
Add lines of smsn estimated denisty or log-density in mix.dens plots (univariate case, p=1).
}
\usage{
mix.lines(y, model, log=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response vector }
  \item{model}{ a variable returned by \code{\link{smsn.mix}} }
  \item{log}{ Logical, plot log-density if TRUE (default = FALSE) }
  \item{...}{further arguments to \code{\link{lines}}}
}

\examples{ ## see \code{\link{bmi}} and \code{\link{smsn.mix}}
}

\seealso{ \code{\link{smsn.mix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lines }

