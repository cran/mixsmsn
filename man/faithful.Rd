\name{faithful}
\docType{data}
\alias{faithful}
\title{ Old Faithful Geyser Data }
\description{
Waiting time between eruptions and the duration of the eruption for the Old Faithful geyser in Yellowstone National Park, Wyoming, USA.
}
\usage{data(faithful)}
\format{
A data frame with 272 observations on 2 variables (p=2)
}

\encoding{latin1}

\source{
Härdle, W. (1991) "Smoothing Techniques with Implementation in S". \emph{New York: Springer}. 

Azzalini, A. and Bowman, A. W. (1990). "A look at some data on the Old Faithful geyser". \emph{Applied Statistics} 39, 357--365. 
}

\examples{
\dontrun{
data(faithful)

## Maximum likelihood estimaton (MLE) for the multivariate FM-SMSN distribution
## with generated values
## Normal
Norm.analysis <- smsn.mmix(faithful, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Normal")
mix.contour(faithful,Norm.analysis,x.min=1,x.max=1,y.min=15,y.max=10,
            levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))

## Calculate the information matrix (when the calc.im option in smsn.mmix is set FALSE)
Norm.im <-  imm.smsn(faithful, Norm.analysis)

## Skew-Normal
Snorm.analysis <- smsn.mmix(faithful, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Skew.normal")
mix.contour(faithful,Snorm.analysis,x.min=1,x.max=1,y.min=15,y.max=10,
            levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))

## Calculate the information matrix (when the calc.im option in smsn.mmix is set FALSE)
Snorm.im <-  imm.smsn(faithful, Snorm.analysis)

## Skew-t
St.analysis <- smsn.mmix(faithful, nu=3, g=2, get.init = TRUE, criteria = TRUE, 
                           group = TRUE, family = "Skew.t")
mix.contour(faithful,St.analysis,x.min=1,x.max=1,y.min=15,y.max=10,
            levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))

## Calculate the information matrix (when the calc.im option in smsn.mmix is set FALSE)
St.im <-  imm.smsn(faithful, St.analysis)

## Passing initial values to MLE and automaticaly calculate the information matrix
mu1 <- c(5,77)
Sigma1 <- matrix(c(0.18,0.60,0.60,41), 2,2)
shape1 <- c(0.69,0.64)

mu2 <- c(2,52)
Sigma2 <- matrix(c(0.15,1.15,1.15,40), 2,2)
shape2 <- c(4.3,2.7)

pii<-c(0.65,0.35)

mu <- list(mu1,mu2)
Sigma <- list(Sigma1,Sigma2)
shape <- list(shape1,shape2)

Snorm.analysis <- smsn.mmix(faithful, nu=3, mu=mu, Sigma=Sigma, shape=shape, pii=pii,
                            g=2, get.init = FALSE, group = TRUE,
                            family = "Skew.normal", calc.im=TRUE)
                            
mix.contour(faithful,Snorm.analysis,x.min=1,x.max=1,y.min=15,y.max=10,
            levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))

## Search for the best number of clusters from g=1 to g=3
faithful.analysis <- smsn.search(faithful, nu = 3, g.min = 1, g.max=3)
mix.contour(faithful,faithful.analysis$best.model,x.min=1,x.max=1,
            y.min=15,y.max=10,levels = c(0.1, 0.015, 0.005, 0.0009,
            0.00015)) 
}
}

\references{
  Marcos Oliveira Prates, Celso Romulo Barbosa Cabral, Victor Hugo Lachos (2013)."mixsmsn: Fitting Finite Mixture of Scale Mixture of 
         Skew-Normal Distributions". Journal of Statistical Software, 54(12), 1-20., URL http://www.jstatsoft.org/v54/i12/.
}

\keyword{datasets}
