library("mvtnorm")
library(Rcpp)
library("RcppArmadillo")
library(fda.usc)



########################################################7
rbvunif_escenario_2 <- function(n,rho) {
  x <- runif(n)
  if ((rho > 1.0) || (rho < -1.0)) {
    stop('rbvunif::rho not in [-1,+1]')
  }
  else if (rho==1.0) xy <- cbind(x,x)
  else if (rho==-1.0) xy <- cbind(x,1-x)
  else if (rho==0.0) xy <- cbind(x,runif(n))
  else {
    a <- (sqrt((49+rho)/(1+rho))-5)/2
    u <- rbeta(n,a,1.0)
    y <- runif(n)
    y <- ifelse(y<0.5 ,abs(u-x), 1-abs(1-u-x) )
    xy <- cbind(x,y)
  }
  return(xy)
}




generar_objeto_escenario_2 = function(n=300, p=2, a1=0.1, b1=0.005, rho=0.3, caso1=1, caso2=2) {
  X = matrix(rnorm(n*p), ncol=p)
  #X= matrix(runif(n*p), ncol=p)
  p1= 1/(1+exp(-1+X[,1]+X[,2]))
  etiqueta= c()
  for(i in 1:n){
    etiqueta= c(etiqueta,rbinom(1,1,prob= p1[i]))
  }
  m= glm(as.factor(etiqueta)~X[,1]+X[,2], family = binomial())
  pnew= m$fitted.values
  # w= etiqueta/(p*n)
  wteorico= etiqueta/(pnew*n)
  wteorico = wteorico/sum(wteorico)
  n1=sum(wteorico>10^-5)
  wteorico = wteorico[wteorico>10^-5]

  if(caso1==1){
    edad1= runif(n1,30,50)
  }

  if(caso1==2){
    edad1= runif(n1,50,70)
  }

  if(caso2==1){
    edad2= runif(n1,30,50)
  }

  if(caso2==2){
    edad2= runif(n1,50,70)
  }

  lregilla= seq(0,1,length=100)
  Q0= 70+230*lregilla
  a0=b0=0

  V1=rbvunif_escenario_2(n1, rho)
  V2=rbvunif_escenario_2(n1, rho)

  V1primera= -20+V1[,1]*40
  V1segunda= -20+V1[,2]*40
  V2primera=0.8+V2[,1]*0.4
  V2segunda=0.8+V2[,2]*0.4



  X1= matrix(0, nrow = n1, ncol= 100)
  Y1= matrix(0, nrow = n1, ncol= 100)


  a0=5
  b0=1


  for(i in 1:n1){
    X1[i,]= V1primera[i]*(a0+a1*edad1[i])+V2primera[i]*(b0+b1*edad1[i])*Q0
    Y1[i,]= V1segunda[i]*(a0+a1*edad2[i])+V2segunda[i]*(b0+b1*edad2[i])*Q0

  }


  # X2 = datos[(n1+1):(n1+n2), c(1,2)]
  # Y2 = datos[(n1+n2+1):(n1+n2+n3), c(3,4)]

  X1 = as.data.frame(X1)
  Y1 = as.data.frame(Y1)
  # X2 = as.data.frame(X2)
  # Y2 = as.data.frame(Y2)

  # result <- list(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2)
  result <- list(X1 = X1, Y1 = Y1, w=wteorico)
  return (result)
}




#######################################################3
Rcpp::sourceCpp('computar_missing.cpp')

calcular_pvalor_escenario_2 = function(n=300, p=2, a1=0.1, b1=0.005, rho=0.3, caso1=1, caso2=2) {
  simular = generar_objeto_escenario_2(n=n, p=p, a1=a1, b1=b1, rho=rho, caso1=caso1, caso2=caso2)

  X1= simular$X1
  Y1= simular$Y1
  w= simular$w

  # indicesconjuntos= 1:2
  # n1= dim(objeto$X1)[1]

  # print("X:")
  # print(objeto$X1)
  # print("Y:")
  # print(objeto$Y1)
  # print("w:")
  # print(objeto$w)

  pvalor1 = energypairedwildboostrapmissing(as.matrix(X1)/100, as.matrix(Y1)/100, w, 1)
  pvalor2 = energypairedwildboostrapmissing(as.matrix(X1)/100, as.matrix(Y1)/100, w, 2)

  lista = list(
    objeto=simular,
    pv_rbf=pvalor1, pv_energy=pvalor2)

  return(lista)
}


# pvalores= c()
# for(i in 1:50){
#   #pvalores= c(pvalores,calcular_pvalor()$pv_rbf)
#   pvalor = calcular_pvalor()
#   #print(pvalor$pv_rbf)
# }






