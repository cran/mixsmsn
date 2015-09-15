\name{im.smsn}
\alias{im.smsn}
\title{ Information matrix }

\description{
Calculate the information matrix of returned analysis based on the 
model \code{family} choice (univariate case, p=1).
}
\usage{
im.smsn(y, model)
}
\arguments{
  \item{y}{ the response vector }
  \item{model}{ a variable returned by \code{\link{smsn.mix}} }
}

\value{
  Estimate the Information Matrix of the parameters.
}

\examples{ ## see \code{\link{bmi}}
}
\author{Marcos Prates \email{marcosop@est.ufmg.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{ \code{\link{smsn.mix}}}

\keyword{interface}
\keyword{models}
