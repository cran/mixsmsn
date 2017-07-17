\name{smsn.search}
\alias{smsn.search}
\title{ Find the best number of cluster for a determined data set. }

\description{
Search for the best fitting for number of cluster from \code{g.min} to \code{g.max}
for a selected \code{family} and \code{criteria} for both univariate and multivariate
distributions.}
\usage{
smsn.search(y, nu,
            g.min = 1, g.max = 3,
            family = "Skew.normal", criteria = "bic",
            error = 0.0001, iter.max = 100, 
            calc.im = FALSE, uni.Gama = FALSE, kmeans.param = NULL, ...)
}
\arguments{
  \item{y}{ the response vector(matrix) }
  \item{nu}{ the parameter of the scale variable (vector or scalar) of the SMSN family (kurtosis parameter). It is necessary to all distributions. For the "Skew.cn" must be a vector of length 2 and values in (0,1) }
  \item{g.min}{ the minimum number of cluster to be modeled }
  \item{g.max}{ the maximum number of cluster to be modeled }
  \item{family}{ distribution famility to be used in fitting ("t", "Skew.t", "Skew.nc", "Skew.slash", "Skew.normal", "Normal") }
  \item{criteria}{ the selection criteria method to be used ("aic", "bic", "edc", "icl") } 
  \item{error}{ the covergence maximum error }
  \item{iter.max}{ the maximum number of iterations of the EM algorithm }
  \item{calc.im}{ if TRUE, the infomation matrix is calculated and the starndard erros are reported }
  \item{uni.Gama}{ if TRUE, the Gamma parameters are restricted to be the same for all clusters (Only valid in the multivariate case, p>1) }
  \item{kmeans.param}{ a list with alternative parameters for the kmeans function when generating initial values, 
                      list(iter.max = 10, n.start = 1, algorithm = "Hartigan-Wong")}
  \item{...}{ other parameters for the hist function }
}

\value{
  Estimated values of the location, scale, skewness and kurtosis parameter from the optimum number of clusters.
}

\examples{ ## see \code{\link{bmi}} and \code{\link{faithful}}
}

\author{Marcos Prates \email{marcosop@est.ufmg.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{\code{\link{smsn.mix}} and \code{\link{smsn.mmix}}}

\keyword{interface}
\keyword{models}
