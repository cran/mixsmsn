\name{smsn.mmix}
\alias{smsn.mmix}
\title{ Fit multivariate FM-SMSN distributions.}

\encoding{UTF-8}

\description{
Return EM algorithm output for multivariate FM-SMSN distributions.
}

\usage{
smsn.mmix(y, nu=1,
          mu = NULL, Sigma = NULL, shape = NULL, pii = NULL,
          g = NULL, get.init = TRUE, criteria = TRUE,
          group = FALSE, family = "Skew.normal", 
          error = 0.0001, iter.max = 100, uni.Gama = FALSE,
          calc.im=FALSE, obs.prob = FALSE, kmeans.param = NULL)
}
\arguments{
  \item{y}{ the response matrix (dimension nxp) }
  \item{nu}{ the parameter of the scale variable (vector or scalar) of the SMSN family (kurtosis parameter). It is necessary to all distributions. For the "Skew.cn" must be a vector of length 2 and values in (0,1) }
  \item{mu}{ a list of \code{g} arguments of vectors of initial values (dimension p) for the location parameters }
  \item{Sigma}{ a list of \code{g} arguments of matrices of initial values (dimension pxp) for the scale parameters }
  \item{shape}{ a list of \code{g} arguments of vectors of initial values (dimension p)for the skewness parameters }
  \item{pii}{ the vector of initial values (dimension g) for the weights for each cluster. Must sum one! }
  \item{g}{ the number of cluster to be considered in fitting}
  \item{get.init}{ if TRUE, the initial values are generated via k-means}
  \item{criteria}{ if TRUE, log-likelihood (logLik), AIC, DIC, EDC and ICL will be calculated }
  \item{group}{ if TRUE, the vector with the classification of the response is returned }
  \item{family}{ distribution famility to be used in fitting ("Skew.t", "t", "Skew.cn", "Skew.slash", "slash", "Skew.normal", "Normal") }
  \item{error}{ the covergence maximum error }
  \item{iter.max}{ the maximum number of iterations of the EM algorithm. Default = 100 }
  \item{uni.Gama}{ if TRUE, the Gamma parameters are restricted to be the same for all clusters}
  \item{calc.im}{ if TRUE, the information matrix is calculated and the starndard erros are reported }
  \item{obs.prob}{ if TRUE, the posterior probability of each observation belonging to one of the g groups is reported}  
  \item{kmeans.param}{ a list with alternative parameters for the kmeans function when generating initial values, 
                      list(iter.max = 10, n.start = 1, algorithm = "Hartigan-Wong")}
}

\value{
  Estimated values of the location, scale, skewness and kurtosis parameter.
  Note: The scale parameters estimated are relative to the entries of 
  the squae root matrix of Sigma.
}

\examples{
mu1 <- c(0,0)
Sigma1 <- matrix(c(3,1,1,3), 2,2)
shape1 <-c(4,4)
nu1 <- 4

mu2 <- c(5,5)
Sigma2 <- matrix(c(2,1,1,2), 2,2)
shape2 <-c(2,2)
nu2 <- 4

pii<-c(0.6,0.4)

arg1 = list(mu=mu1, Sigma=Sigma1, shape=shape1, nu=nu1)
arg2 = list(mu=mu2, Sigma=Sigma2, shape=shape2, nu=nu2)
y <- rmmix(n= 500, p = pii, "Skew.t", list(arg1,arg2))

\dontrun{

## Normal fit giving intial values
mu <- list(mu1,mu2)
Sigma <- list(Sigma1,Sigma2)
shape <- list(shape1,shape2)
pii <- c(0.6,0.4)

Norm.analysis <- smsn.mmix(y, nu=3, mu=mu, Sigma=Sigma, shape=shape, pii = pii,
                           criteria = TRUE, g=2, get.init = FALSE, group = TRUE,
                           family = "Normal")
mix.contour(y,Norm.analysis)

## Normal fit 
Norm.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Normal")
mix.contour(y,Norm.analysis)

## Normal fit with a unique Gamma
Norm.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Normal", uni.Gama = TRUE)
mix.contour(y,Norm.analysis)


## Skew Normal fit
Snorm.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                            group = TRUE, family = "Skew.normal")
mix.contour(y,Snorm.analysis)

## t fit
t.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                         group = TRUE, family = "t")
mix.contour(y,t.analysis)

## Skew t fit
St.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                         group = TRUE, family = "Skew.t")
mix.contour(y,St.analysis)

## Skew Contaminated Normal fit
Scn.analysis <- smsn.mmix(y, nu=c(0.1,0.1), g=2, get.init = TRUE, criteria = TRUE, 
                          group = TRUE, family = "Skew.cn",error=0.01)
mix.contour(y,Scn.analysis)

## Skew Contaminated Normal fit
Sslash.analysis <- smsn.mmix(y, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                             group = TRUE, family = "Skew.slash", error=0.1)
mix.contour(y,Sslash.analysis)

}
}

\references{
  Cabral, C. R. B., Lachos, V. H. and Prates, M. O. (2012). "Multivariate Mixture Modeling Using Skew-Normal Independent Distributions". \emph{Computational Statistics & Data Analysis}, 56, 126-142, doi:10.1016/j.csda.2011.06.026.

  Marcos Oliveira Prates, Celso Romulo Barbosa Cabral, Victor Hugo Lachos (2013)."mixsmsn: Fitting Finite Mixture of Scale Mixture of Skew-Normal Distributions". Journal of Statistical Software, 54(12), 1-20., URL https://doi.org/10.18637/jss.v054.i12.
}


\author{Marcos Prates \email{marcosop@est.ufmg.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{\code{\link{mix.contour}},  \code{\link{rmmix}} and \code{\link{smsn.search}}}

\keyword{interface}
\keyword{models}
