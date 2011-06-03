\name{imm.smsn}
\alias{imm.smsn}
\title{ Information matrix }

\description{
Calculate the information matrix of returned analysis based on the 
model \code{family} choice (multivariate case, p>=2).
}
\usage{
imm.smsn(y, model)
}
\arguments{
  \item{y}{ the response vector (p>2)}
  \item{model}{ a variable returned by \code{\link{smsn.mmix}} }
}

\value{
  Estimate the Information Matrix of the parameters.
  Note: In the Information Matrix the scale parameters 
        estimates are relative to the entries of square 
        root matrix of Sigma.

}

\examples{ ## see \code{\link{faithful}}
}
\author{Marcos Prates \email{marcosop@ime.unicamp.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{ \code{\link{smsn.mmix}}}

\keyword{interface}
\keyword{models}
