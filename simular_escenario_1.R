library("mvtnorm")
library(Rcpp)
library("RcppArmadillo")
library(fda.usc)


rbvunif_escenario_1 <- function(n,rho) {
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


generar_objeto_escenario_1 = function(n=300, n1=150, n2=150, n3=150, a1=0.1, b1=0.005, rho=0.3, caso1=1, caso2=1) {
  p = 2

  n1=150
  n2=150
  n2=n3

  # print(caso1)
  # print(caso2)
  X = matrix(rnorm(n*p), ncol=p)
  #X= matrix(runif(n*p), ncol=p)

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



  if(caso1==1){
    edad3= runif(n2,30,50)
  }

  if(caso1==2){
    edad3= runif(n2,50,70)
  }

  if(caso2==1){
    edad4= runif(n3,30,50)
  }

  if(caso2==2){
    edad4= runif(n3,50,70)
  }








  lregilla= seq(0,1,length=100)
  Q0= 70+230*lregilla
  a0=b0=0

  V1=rbvunif_escenario_1(n1, rho)
  V2=rbvunif_escenario_1(n1, rho)

  V1primera= -20+V1[,1]*40
  V1segunda= -20+V1[,2]*40
  V2primera=0.8+V2[,1]*0.4
  V2segunda=0.8+V2[,2]*0.4


  V11= runif(n2,-20,20)
  V12= runif(n3,-20,20)

  V21= runif(n2,0.8,1.2)
  V22= runif(n3,0.8,1.2)

  # V11= -20+V11*40
  # V12= -20*V12*40
  #
  # V21= 0.8*V21*0.4
  # V22= 0.8+V22*0.4
  #


  X1= matrix(0, nrow = n1, ncol= 100)
  Y1= matrix(0, nrow = n1, ncol= 100)

  X1segundo= matrix(0, nrow = n2, ncol= 100)
  Y1segundo= matrix(0, nrow = n3, ncol= 100)

  a0=5
  b0=1

  for(i in 1:n1){
    X1[i,]= V1primera[i]*(a0+a1*edad1[i])+V2primera[i]*(b0+b1*edad1[i])*Q0
    Y1[i,]= V1segunda[i]*(a0+a1*edad2[i])+V2segunda[i]*(b0+b1*edad2[i])*Q0

  }

  for(i in 1:n2){
    X1segundo[i,]= V11[i]*(a0+a1*edad3[i])+V21[i]*(b0+b1*edad3[i])*Q0
    Y1segundo[i,]= V12[i]*(a0+a1*edad4[i])+V22[i]*(b0+b1*edad4[i])*Q0

  }

  # X2 = datos[(n1+1):(n1+n2), c(1,2)]
  # Y2 = datos[(n1+n2+1):(n1+n2+n3), c(3,4)]

  X1 = as.data.frame(X1)
  Y1 = as.data.frame(Y1)
   X2 = as.data.frame(X1segundo)
   Y2 = as.data.frame(Y1segundo)

   result <- list(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2)
 # result <- list(X1 = X1, Y1 = Y1, w=wteorico)
  return (result)
}




calcular_pvalor_escenario_1 = function(n, n1, n2, n3, a1, b1, rho, caso1, caso2, permutaciones=10) {
  sourceCpp('computar.cpp')

   objeto = generar_objeto_escenario_1(n, n1, n2, n3, a1, b1, rho, caso1, caso2)

  indices_conjuntos = 0:(n2+n3-1)
  perm = matrix(0, nrow=permutaciones, ncol=length(indices_conjuntos))
  nsize = length(indices_conjuntos)
  for(j in 1:permutaciones) {
    perm[j,] = sample(indices_conjuntos, size=nsize, replace=FALSE)
  }


  pvalor1 = energypairedwildboostrap2(as.matrix(objeto$X1)/100, as.matrix(objeto$Y1)/100,
                                     as.matrix(objeto$X2)/100, as.matrix(objeto$Y2)/100, perm, 1)

 # pvalor2 = energypairedwildboostrap2(as.matrix(objeto$X1)/100, as.matrix(objeto$Y1)/100,
                                     #as.matrix(objeto$X2)/100, as.matrix(objeto$Y2)/100, perm, 2)

  lista = list(objeto=objeto, pv_rbf=pvalor1)

  # print(paste(pvalor1))


  return(lista)
}


