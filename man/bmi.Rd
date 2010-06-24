\name{bmi}
\docType{data}
\alias{bmi}
\title{ Body Mass Index }
\description{
The data set has the measure of the Body Mass Index (bmi) for 2107 people.
}
\usage{data(bmi)}
\format{
  A data frame with 2107 observations of bmi
}

\encoding{latin1}

\source{
Rodrigo M. Basso, Victor H. Lachos, Celso R. B. Cabral, Pulak Ghosh (2009). "Robust mixture modeling based on scale mixtures of skew-normal
distributions". \emph{Computational Statistics and Data Analysis} (in press). doi: 10.1016/j.csda.2009.09.031
}

\examples{
\dontrun{
data(bmi)
y <-bmi$bmi

hist(y,breaks=40)

## Maximum likelihood estimaton (MLE) with generated values
bmi.analysis <- smsn.mix(y, nu = 3, g = 2, get.init = TRUE, criteria = TRUE, 
                         group = TRUE, calc.im=FALSE)
mix.hist(y,bmi.analysis)

## Passing initial values to MLE
mu1 <- 20; mu2 <- 35
sigma2.1 <- 9; sigma2.2 <- 9;
lambda1 <- 0; lambda2 <- 0;
pii<- c(0.5,0.5)

mu <- c(mu1,mu2)
sigma2 <- c(sigma2.1,sigma2.2)
shape <- c(lambda1,lambda2)

bmi.analysis <- smsn.mix(y, nu = 3, mu, sigma2 , shape, pii, get.init = FALSE,
                         criteria = TRUE, group = TRUE, calc.im=FALSE)
mix.hist(y,bmi.analysis)

## Calculate the information matrix (when the calc.im option in smsn.mix is set FALSE)
bmi.im <-  im.smsn(y, bmi.analysis)

## Search for the best number of clusters from g=1 to g=5
bmi.analysis <- smsn.search(y, nu = 3, g.min = 1, g.max=5)
mix.hist(y,bmi.analysis$best.model)
}
}
\keyword{datasets}
