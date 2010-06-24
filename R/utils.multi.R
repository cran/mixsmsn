######################################################
#     algoritmo para graficar os contornos

mix.contour <- function(y, model, slice = 100, ncontour = 10, x.min=1, x.max=1, y.min=1,y.max=1, ...){
   # Essa função só serve para graficar os contornos de misturas finitas de SMSN BIVARIADO!!!
   # dat: o cosliceunto de dat a ser plotado no R^2
   # model: deve ser um objeto resultante da função EMmulti.MIXSNI
   # slice e ncountor são parametros passados para a função countor
   # ?contour para detalhes
   dat <- y
   y <- NULL
   n <- nrow(dat)
   p <- ncol(dat)

   if(p != 2) stop("The mix.contour function is only appropriate for the bivariate analysis.\n")
   if((class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.",sep=" "))
   if(length(model$group) == 0) stop("The groups were not save in the model.\n")
   if ((x.min < 0) || (x.max < 0) || (y.min < 0) || (y.max < 0)) stop("All limits must be non negative.\n")
   g <- length(model$pii)
   if (class(model) == "Normal"){
      mixed.Normal <- function(x, y, pii, mu, Sigma) {
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]*dmvnorm(cbind(x, y), mu[[j]], Sigma[[j]])
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)    #faithfull
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)  #faithfull
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <- mixed.Normal(x[i], y[j], model$pii, model$mu, model$Sigma)
      contour(x, y, f, , nlevels = ncontour, ...) 
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.normal"){
      mixed.SN <- function(x, y, pii, mu, Sigma, lambda) {
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]*2*dmvnorm(cbind(x, y), mu[[j]], Sigma[[j]])*pnorm(t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%t(cbind(x, y) - mu[[j]]))
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.SN(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.t"){
      mixed.ST <- function(x, y, pii, mu, Sigma, lambda, nu) {
#        n <- nrow(dat)
#        p <- ncol(dat)
        dens <- 0
        for (j in 1:g) {
          denst <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Sigma[[j]])^(-1/2)*(1 + mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]])/nu)^(-(p+nu)/2)
          dens <- dens + pii[j] * 2*(denst)*pt(sqrt((p + nu)/(mahalanobis(cbind(x,y), mu[[j]], Sigma[[j]]) + nu))*t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]]), df = nu + p)
        }
      }
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.ST(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape, model$nu)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.cn"){
      mixed.SNC <- function(x, y, pii, mu, Sigma, lambda, nu) {
#        n <- nrow(y)
#        p <- ncol(y)
        dens <- 0
        for (j in 1:g) dens <- dens + pii[j]* 2*(nu[1]*dmvnorm(cbind(x,y), mu[[j]], Sigma[[j]]/nu[2])*pnorm(sqrt(nu[2])*t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]]) ) + (1 - nu[1])*dmvnorm(cbind(x,y), mu[[j]], Sigma[[j]])*pnorm(t(lambda[[j]])%*%solve(matrix.sqrt(Sigma[[j]]))%*%(c(x,y) - mu[[j]])) )
      }
#      xc <- yc <- 0
#      for (j in 1:g) {
#        xc <- xc + model$mu[[j]][1]
#        yc <- yc + model$mu[[j]][2]
#      }
#      xc <- xc / g
#      yc <- yc / g
#      x <- seq(xc - x.min, xc + x.max, length = slice)
#      y <- seq(yc - y.min, yc + y.max, length = slice)
#      f <- matrix(0,slice,slice)
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  mixed.SNC(x[i], y[j], model$pii, model$mu, model$Sigma, model$shape, model$nu)
      contour(x, y, f, nlevels = ncontour, ...)
      points(dat[,1], dat[,2], col = (model$group+1))
   }

   if (class(model) == "Skew.slash"){
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)
      #y <- seq(min(dat[,2])-15,max(dat[,2])+10, length = slice)
      #x <- seq(min(dat[,1])-1,max(dat[,1])+1, length = slice)     # Swiss
      #y <- seq(min(dat[,2])-1,max(dat[,2])+1, length = slice)     # Swiss
      x <- seq(min(dat[,1])-x.min,max(dat[,1])+x.max, length = slice)    
      y <- seq(min(dat[,2])-y.min,max(dat[,2])+y.max, length = slice)  
      f <- matrix(0,slice,slice)
      for (i in 1:slice) for (j in 1:slice) f[i,j] <-  d.mixedmvSS(cbind(x[i],y[j]), model$pii, model$mu, model$Sigma, model$shape, model$nu)
      #contour(x, y, f, nlevels = 5, levels = c(0.1, 0.015, 0.005, 0.0009, 0.00015))
      contour(x, y, f, , nlevels = ncontour, ...)
      #contour(x, y, f, , levels = c(0.3,0.15, 0.075,0.0375,0.01875)) # Swiss
      points(dat[,1], dat[,2], col = (model$group+1))
   }
   title(main=paste("Contour plot for",class(model)), font = 2)


}

