\name{rmmix}
\alias{rmmix}
\title{ Random multivariate FM-SMSN generator }
\description{
Random generator of multivariate FM-SMSN distributions.}
\usage{
rmmix(n, pii, family, arg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations  }
  \item{pii}{ a vector of weights for the mixture (dimension of the number \code{g} of clusters). Must sum to one!}
  \item{family}{ distribution family to be used in fitting ("Skew.t", "Skew.nc", "Skew.slash", "Skew.normal", "Normal") }
  \item{arg}{ a list of \code{g} lists with each list containing the necessary parameters of the selected \code{family} }
}

\examples{ ## see \code{\link{smsn.mmix}}}

\author{Marcos Prates \email{marcos.prates@uconn.edu}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{ \code{\link{smsn.mmix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

