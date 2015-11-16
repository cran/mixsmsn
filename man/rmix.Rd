\name{rmix}
\alias{rmix}
\title{ Random univariate FM-SMSN generator }
\description{
Random generator of univariate FM-SMSN distributions.}
\usage{
rmix(n, pii, family, arg, cluster=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations  }
  \item{pii}{ a vector of weights for the mixture (dimension of the number \code{g} of clusters). Must sum to one!}
  \item{family}{ distribution family to be used in fitting ("t", "Skew.t", "Skew.cn", "Skew.slash", "Skew.normal", "Normal") }
  \item{arg}{ a list with each entry containing a vector of size equal to the number of clusters of the necessary parameters from a \code{family} }
  \item{cluster}{TRUE or FALSE if the true observations clusters must be returned. }
}

\examples{ ## see \code{\link{smsn.mix}}
}

\author{Marcos Prates \email{marcosop@est.ufmg.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{ \code{\link{smsn.mix}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

