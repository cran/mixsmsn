##################################################################
##########     Funções para gerar SNI e mistura       ############


  rmix <- function(n, p, family, arg) {
    #Função para gerar misturas de g populações
    #n: numero de amostras geradas
    #p: vetor de proporções das misturas (tamanho g)
    #arg: deve ser um tipo list com cada entrada contendo um vetor de tamanho g de agrumentos a ser passado para rF1
 if((family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) stop(paste("Family",family,"not recognized.",sep=" "))

    if((family == "Normal") || (family == "Skew.normal") ) {
                                                              rF1 <- gen.Skew.normal
                                                              for (i in 1:length(arg)) if(length(arg[[i]]) != 4 && length(arg[[i]]) != 3) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                                           }
                                                           
    if (family == "Skew.t"){
                             rF1 <- gen.Skew.t
                             for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                           }
    if (family == "Skew.cn"){
                              rF1 <- gen.Skew.cn
                              for (i in 1:length(arg)) if(length(arg[[i]]) != 5) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                            }
    if (family == "Skew.slash") {
                                 rF1 <- gen.Skew.slash
                                 for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                }

    x1 <- vector(mode = "numeric", length = n)
    g <- length(p)
    interval <- c(0)
    for (j in 1:g-1) interval <- cbind(interval, interval[j] + p[j])
    interval <- cbind(interval, 1)
    for(i in 1:n) {
      u <- runif(1)
      x1[i] <- do.call("rF1", c(list(1), arg[[findInterval(u, interval)]]))
    }
    return(x1)
  }


  gen.Skew.normal <- function(n, mu, sigma2, shape, nu=NULL){
    #Função para gerar valores aleatórios de uma Skew-Normal
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locação, escala e assimetria, respectivamente
    delta <- shape / sqrt(1 + shape^2)
    y <- mu*rep(1,n) + sqrt(sigma2)*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
    return(y)
  }

  gen.Skew.t <- function(n, mu, sigma2, shape, nu ){
    #Função para gerar Skew-t
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locação, escala e assimetria, respectivamente
    y <- mu + (rgamma(n, nu/2, nu/2))^(-1/2)*gen.Skew.normal(n, 0, sigma2, shape)
  }


  gen.Skew.cn <- function(n, mu, sigma2, shape, nu){
    #Função para gerar Skew Normal Contaminada
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locação, escala e assimetria, respectivamente
    rmix(n, nu[1], gen.Skew.normal, list(c(mu,sigma2/nu[2],shape), c(mu,sigma2,shape)))
  }


  gen.Skew.slash <- function(n, mu, sigma2, shape, nu){
    # Função para gerar Skew Slash
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locação, escala e assimetria, respectivamente
    u1 <- runif(n)
    u2 <- u1^(1/(nu))   # formula 10 do artigo e método da inversão
    ys <- mu + (u2)^(-1/2)*gen.Skew.normal(n, 0, sigma2, shape)
    return(ys)
  }


##########  FIM   Funções para gerar SNI e mistura     ###########
##################################################################

