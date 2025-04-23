##################################################################
##########     Plot histogram                #####################

mix.hist <- function(y, model, breaks=40, main=paste("Histogram of" , class(model),"fit"), col.hist="grey", col.dens="red", ...){
  if(!inherits(model, 't') && !inherits(model, 'Skew.t') && !inherits(model, 'Skew.cn') &&
     !inherits(model, 'Skew.slash') && !inherits(model, 'Skew.normal') &&
     !inherits(model, 'Normal')) stop(paste("Class of family",class(model),"not recognized.",sep=" "))
  
  y <- as.matrix(y) 
  if (dim(y)[2] != 1) stop("The mix.hist function is only appropriate for the univariate analysis.\n")
   
  if (inherits(model, 'Skew.t') || inherits(model, 't')){
      #### grafico ajuste
      hist(y, breaks = breaks,probability=T,col=col.hist,main=main,...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedST(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col=col.dens)
     }
  if (inherits(model, 'Skew.cn')){
      hist(y, breaks = breaks,probability=T,col=col.hist,main=main,...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSNC(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col=col.dens)
     }
  if (inherits(model, 'Skew.slash')){
      hist(y, breaks = breaks,probability=T,col=col.hist,main=main,...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSS(xx, model$pii, model$mu, model$sigma2, model$shape, model$nu),col=col.dens)
     }
  if (inherits(model, 'Skew.normal') || inherits(model, 'Normal')){
      hist(y, breaks = breaks,probability=T,col=col.hist,main=main,...)
      xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      lines(xx,d.mixedSN(xx, model$pii, model$mu, model$sigma2, model$shape),col=col.dens)
     }
}
