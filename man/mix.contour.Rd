\name{mix.contour}
\alias{mix.contour}
\title{ Print the selected groups with contours }
\description{
Plot the contour of the observations with the group selection.
}
\usage{
mix.contour(y, model,
            slice=100, ncontour=10,
            x.min=1, x.max=1,
            y.min=1,y.max=1,
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response matrix (dimension nx2) }
  \item{model}{ a variable returned by \code{\link{smsn.mmix}} }
  \item{slice}{ number of slices in the sequenceo the contour }
  \item{ncontour}{number of contours to be ploted }
  \item{x.min}{ value to be subtracted of the smallest observation in the x-axis }
  \item{x.max}{ value to be added of the biggest observation in the x-axis }
  \item{y.min}{ value to be subtracted of the smallest observation in the y-axis }
  \item{y.max}{ value to be added of the biggest observation in the y-axis }
  \item{...}{further arguments to \code{\link{contour}}}
}

\examples{ ## see \code{\link{smsn.mmix}}
}

\seealso{ \code{\link{smsn.mmix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }

