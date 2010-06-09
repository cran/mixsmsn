\name{smsn.mix.search}
\alias{smsn.mix.search}
\title{ Find the best number of cluster for a determined data set. }

\description{
Search for the best fitting for number of cluster from \code{g.min} to \code{g.max}
for a selected \code{family} and \code{criteria}.}
\usage{
smsn.mix.search(y, nu,
                g.min = 1, g.max = 3,
                family = "Skew.normal", criteria = "bic",
                error = 0.00001, iter.max = 100, calc.im = FALSE)
}
\arguments{
  \item{y}{ the response vector }
  \item{nu}{ the parameter of the scale variable (vector os scalar) of the SMSN family (hyperparameter). It is necessary to all distributions. For the "Skew.cn" must be a vector of length 2 and values in (0,1) }
  \item{g.min}{ the minimum number of cluster to be modeled }
  \item{g.max}{ the maximum number of cluster to be modeled }
  \item{family}{ distribution famility to be used in fitting ("Skew.t", "Skew.nc", "Skew.slash", "Skew.normal", "Normal") }
  \item{criteria}{ the selection criteria method to be used ("aic", "bic", "edc") } 
  \item{error}{ the covergence maximum error }
  \item{iter.max}{ the maximum number of iterations of the EM algorithm }
  \item{calc.im}{ if TRUE, the infomation matrix is calculated and the starndard erros are reported }
}

\value{
  Estimated values of the location, scale, skewness and hyperparameter from the optimum number of clusters.
}

\examples{ ## see \code{\link{bmi}}}

\author{Marcos Prates \email{marcos.prates@uconn.edu}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{\code{\link{smsn.mix}}}

\keyword{interface}
\keyword{models}
