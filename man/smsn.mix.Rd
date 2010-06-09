\name{smsn.mix}
\alias{smsn.mix}
\title{ Fit finite mixture FM-SMSN distribution.}

\encoding{UTF-8}

\description{
Return EM algorithm output for FM-SMSN distributions.
}

\usage{
smsn.mix(y, 
         nu, mu = NULL, sigma2 = NULL, shape = NULL, pii = NULL,
         g = NULL, get.init = TRUE,
         criteria = TRUE, group = FALSE, family = "Skew.normal",
         error = 0.00001, iter.max = 100, calc.im = TRUE)
}
\arguments{
  \item{y}{ the response vector }
  \item{nu}{ the parameter of the scale variable (vector os scalar) of the SMSN family (hyperparameter). It is necessary to all distributions. For the "Skew.cn" must be a vector of length 2 and values in (0,1) }
  \item{mu}{ the vector of initial values (dimension g) for the location parameters }
  \item{sigma2}{ the vector of initial values (dimension g) for the scale parameters }
  \item{shape}{ the vector of initial values (dimension g) for the skewness parameters }
  \item{pii}{ the vector of initial values (dimension g) for the weights for each cluster. Must sum one! }
  \item{g}{ the number of cluster to be considered in fitting}
  \item{get.init}{ if TRUE, the initial values are generated via k-means}
  \item{criteria}{ if TRUE, AIC, DIC and EDC will be calculated }
  \item{group}{ if TRUE, the vector with the classification of the response is returned }
  \item{family}{ distribution famility to be used in fitting ("Skew.t", "Skew.nc", "Skew.slash", "Skew.normal", "Normal") }
  \item{error}{ the covergence maximum error }
  \item{iter.max}{ the maximum number of iterations of the EM algorithm. Default = 100 }
  \item{calc.im}{ if TRUE, the infomation matrix is calculated and the standard errors are reported}
}

\value{
  Estimated values of the location, scale, skewness and degree of freedom.
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
par(mfrow=c(2,2))
## Normal fit
Norm.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                          group = TRUE, family = "Normal", calc.im=FALSE)
mix.hist(y,Norm.analysis)
mix.print(Norm.analysis)

## Skew Normal fit
Snorm.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Skew.normal", calc.im=FALSE)
mix.hist(y,Snorm.analysis)
mix.print(Snorm.analysis)

## Skew t fit
St.analysis <- smsn.mix(y, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, 
                        group = TRUE, family = "Skew.t", calc.im=FALSE)
mix.hist(y,St.analysis)
mix.print(St.analysis)

## Skew Contaminated Normal fit
Scn.analysis <- smsn.mix(y, nu = c(0.3,0.3), g = 3, get.init = TRUE, criteria = TRUE, 
                         group = TRUE, family = "Skew.cn", calc.im=FALSE)
mix.hist(y,Scn.analysis)
mix.print(Scn.analysis)
}
}

\references{
  Rodrigo M. Basso, Víctor H. Lachos, Celso R. B. Cabral, Pulak Ghosh (2009). "Robust mixture modeling based on scale mixtures of skew-normal
  distributions". \emph{Computational Statistics and Data Analysis} (in press). doi: 10.1016/j.csda.2009.09.031
}

\author{Marcos Prates \email{marcos.prates@uconn.edu}, Victor Lachos \email{hlachos@ime.unicamp.br} and Celso Cabral \email{celsoromulo@gmail.com}}

\seealso{\code{\link{mix.hist}}, \code{\link{im.smsn}} and \code{\link{smsn.mix.search}}}

\keyword{interface}
\keyword{models}
