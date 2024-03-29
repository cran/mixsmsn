\name{smsn.mix}
\alias{smsn.mix}
\title{ Fit univariate FM-SMSN distribution }

\encoding{UTF-8}

\description{
Return EM algorithm output for FM-SMSN distributions (univaritate case, p=1).
}

\usage{
smsn.mix(y, 
         nu, mu = NULL, sigma2 = NULL, shape = NULL, pii = NULL,
         g = NULL, get.init = TRUE,
         criteria = TRUE, group = FALSE, family = "Skew.normal",
         error = 0.00001, iter.max = 100, calc.im = TRUE, obs.prob = FALSE,
         kmeans.param = NULL)
}
\arguments{
  \item{y}{ the response vector }
  \item{nu}{ the parameter of the scale variable (vector or scalar) of the SMSN family (kurtosis parameter). It is necessary to all distributions. For the "Skew.cn" must be a vector of length 2 and values in (0,1) }
  \item{mu}{ the vector of initial values (dimension g) for the location parameters }
  \item{sigma2}{ the vector of initial values (dimension g) for the scale parameters }
  \item{shape}{ the vector of initial values (dimension g) for the skewness parameters }
  \item{pii}{ the vector of initial values (dimension g) for the weights for each cluster. Must sum one! }
  \item{g}{ the number of cluster to be considered in fitting}
  \item{get.init}{ if TRUE, the initial values are generated via k-means}
  \item{criteria}{ if TRUE, AIC, DIC, EDC and ICL will be calculated }
  \item{group}{ if TRUE, the vector with the classification of the response is returned }
  \item{family}{ distribution family to be used in fitting ("Skew.t", "t", "Skew.cn", "Skew.slash", "slash", "Skew.normal", "Normal") }
  \item{error}{ the covergence maximum error }
  \item{iter.max}{ the maximum number of iterations of the EM algorithm. Default = 100 }
  \item{calc.im}{ if TRUE, the information matrix is calculated and the standard errors are reported}
  \item{obs.prob}{ if TRUE, the posterior probability of each observation belonging to one of the g groups is reported}
  \item{kmeans.param}{ a list with alternative parameters for the kmeans function when generating initial values, 
                      list(iter.max = 10, n.start = 1, algorithm = "Hartigan-Wong")}
}

\value{
  Estimated values of the location, scale, skewness and kurtosis parameter.
}

\examples{
mu1 <- 5; mu2 <- 20; mu3 <- 35
sigma2.1 <- 9; sigma2.2 <- 16; sigma2.3 <- 9
lambda1 <- 5; lambda2 <- -3; lambda3 <- -6
nu = 5

mu <- c(mu1,mu2,mu3)
sigma2 <- c(sigma2.1,sigma2.2,sigma2.3)
shape <- c(lambda1,lambda2,lambda3)
pii <- c(0.5,0.2,0.3)

arg1 = c(mu1, sigma2.1, lambda1, nu)
arg2 = c(mu2, sigma2.2, lambda2, nu)
arg3 = c(mu3, sigma2.3, lambda3, nu)
y <- rmix(n=1000, p=pii, family="Skew.t", arg=list(arg1,arg2,arg3))

\dontrun{
par(mfrow=c(1,2))
## Normal fit
Norm.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                          group = TRUE, family = "Normal", calc.im=FALSE)
mix.hist(y,Norm.analysis)
mix.print(Norm.analysis)
mix.dens(y,Norm.analysis)

## Skew Normal fit
Snorm.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Skew.normal", calc.im=FALSE)
mix.hist(y,Snorm.analysis)
mix.print(Snorm.analysis)
mix.dens(y,Snorm.analysis)

## t fit
t.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                        group = TRUE, family = "t", calc.im=FALSE)
mix.hist(y,t.analysis)
mix.print(t.analysis)
mix.dens(y,t.analysis)

## Skew t fit
St.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                        group = TRUE, family = "Skew.t", calc.im=FALSE)
mix.hist(y,St.analysis)
mix.print(St.analysis)
mix.dens(y,St.analysis)

## Skew Contaminated Normal fit
Scn.analysis <- smsn.mix(y, nu = c(0.3,0.3), g = 3, get.init = TRUE, criteria = TRUE, 
                         group = TRUE, family = "Skew.cn", calc.im=FALSE)
mix.hist(y,Scn.analysis)
mix.print(Scn.analysis)
mix.dens(y,Scn.analysis)

par(mfrow=c(1,1))
mix.dens(y,Norm.analysis)
mix.lines(y,Snorm.analysis,col="green")
mix.lines(y,t.analysis,col="red")
mix.lines(y,St.analysis,col="blue")
mix.lines(y,Scn.analysis,col="grey")
}
}


\references{
  Rodrigo M. Basso, Victor H. Lachos, Celso R. B. Cabral, Pulak Ghosh (2010). "Robust mixture modeling based on scale mixtures of skew-normal
  distributions". Computational Statistics and Data Analysis, 54, 2926-2941. doi: 10.1016/j.csda.2009.09.031

  Marcos Oliveira Prates, Celso Romulo Barbosa Cabral, Victor Hugo Lachos (2013)."mixsmsn: Fitting Finite Mixture of Scale Mixture of Skew-Normal Distributions". Journal of Statistical Software, 54(12), 1-20., URL https://doi.org/10.18637/jss.v054.i12.
}

\author{Marcos Prates \email{marcosop@est.ufmg.br}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{\code{\link{mix.hist}}, \code{\link{im.smsn}} and \code{\link{smsn.search}}}

\keyword{interface}
\keyword{models}
