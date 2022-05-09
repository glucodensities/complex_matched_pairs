#include <RcppArmadillo.h>
#include <cmath>
#include <stdio.h>
// [[Rcpp::depends(RcppArmadillo)]]
//using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

//Dino Sejdinovic, 2013
//D. Sejdinovic, A. Gretton and W. Bergsma.  A KERNEL TEST FOR THREE-VARIABLE INTERACTIONS, 2013.

//median heuristic bandwidth selection

float median_heur(mat Z) {

  mat Zmed(100,Z.n_cols);

  uvec indtemp;
  int size1=Z.n_rows;



  if (size1>100){//choose 100 random samples from Z if it contains
    //more than 100 points
    uvec ind=sort_index(randu(size1));
    Zmed = Z.rows(ind.subvec(0,99));
    size1 = 100;
  }else{
    Zmed = Z;
  }

  mat G = sum((Zmed%Zmed),1);
  mat Q = repmat(G,1,size1);
  mat R = repmat(G.t(),size1,1);
  mat dists = Q + R - 2*(Zmed*Zmed.t());
  dists = dists-trimatl(dists);
  vec distsr = reshape(dists,pow(size1,2),1);
  distsr = sqrt(distsr);
  float sig = sqrt(0.5*median(distsr(find(distsr>0))));
  return(sig);

}


class Kernel {
public:
  virtual mat getKernelMatrix(const mat &X, const mat &Y) const = 0;
};


class RBFKernel : public Kernel {

private:
  double s; //s is the kernel parameter

public:
  RBFKernel(double sigma = 1) {
    s = sigma;
  }

  mat getKernelMatrix(const mat &X, const mat &Y) const {
    mat G = sum((X%X),1);
    mat H = sum((Y%Y),1);

    mat Q = repmat(G,1,Y.n_rows);
    mat R = repmat(H.t(),X.n_rows,1);

    H = Q + R - 2*X*Y.t();


    mat H2=exp(-(H/2)/pow(s,2));

    return(H2);
  }
};



class EnergyKernel : public Kernel {

private:
  double alpha; //s is the kernel parameter

public:
  EnergyKernel(double _alpha) {
    alpha = _alpha;
  }

  mat getKernelMatrix(const mat &X, const mat &Y) const {
    mat R(X.n_rows, Y.n_rows, fill::zeros);
    for (uword i=0; i < X.n_rows; i++) {
      for (uword j=0; j <= i; j++) {
        // R(i,j) = arma::sum(X.row(i) * Y.row(j));
        // R(i,j) = pow(norm(X.row(i)),2) + pow(norm(Y.row(j)),2) - 2 * pow(norm(X.row(i) - Y.row(j)),2);
        R(i,j) = 0.5 * (pow(norm(X.row(i)),2) + pow(norm(Y.row(j)),2) - pow(norm(X.row(i) - Y.row(j)),2));
        R(j,i) = R(i,j);
      }
    }
    return(R);
  }
};




// [[Rcpp::export]]
arma::mat wild(double ln, int permutaciones, int n1){

  arma::mat resultado=arma::randn(permutaciones,n1);
  arma::mat auxiliar=arma::randn(permutaciones,1);

  resultado=   resultado*sqrt(1-exp(-2/ln));

  arma::mat resultado2(permutaciones,n1);


  for(int i=0;i<n1;i++){

    if(i==0){

      resultado2.col(i)=auxiliar*exp(-1/ln)+resultado.col(i);

    }else{

      resultado2.col(i)=resultado2.col(i-1)*exp(-1/ln)+resultado.col(i);

    }


  }



  return resultado2;

}

// [[Rcpp::export]]
double energypairedwildboostrap(arma::mat d11, arma::mat d12, arma::mat d13, arma::mat d21, arma::mat d22, arma::mat d23,arma::umat permutaciones){

  int n1= d11.n_rows;
  int n2= d21.n_rows;
  int n3= d22.n_rows;



  double acum=0;
  double acum2=0;
  double acum3=0;
  double acum4=0;





  for(int i=0;i<n1;i++){
    for(int j=0;j<n1;j++){
      acum= acum+d11(i,j)+d12(i,j)-2*d13(i,j);

    }
  }

  acum= acum*(1/double(n1))*(1/double(n1));





  for(int i=0;i<n2;i++){
    for(int j=0;j<n2;j++){
      acum2= acum2+d21(i,j);

    }
  }

  acum2= acum2*(1/double(n2))*(1/double(n2));




  for(int i=0;i<n3;i++){
    for(int j=0;j<n3;j++){
      acum3= acum3+d22(i,j);

    }
  }

  acum3= acum3*(1/double(n3))*(1/double(n3));



  for(int i=0;i<n2;i++){
    for(int j=0;j<n3;j++){
      acum4= acum4+d23(i,j);

    }
  }

  acum4=  2*acum4*(1/double(n2))*(1/double(n3));



  double resfinal= acum;
  double ln= log(n1);
  int npermutaciones=2000;

  arma::mat matrizW=wild(ln,npermutaciones,n1);

  arma::vec resultados(npermutaciones);

  arma::rowvec vectorw(npermutaciones);

  // Rcout << "The value of v : " << vectorw << "\n";



  for(int r=0;r<npermutaciones;r++){

    acum=0;
    // Rcout << "The value of v : " << matrizW.row(r) << "\n";

    vectorw=matrizW.row(r);
    // Rcout << "The value of v : " << vectorw << "\n";

    for(int i=0;i<n1;i++){
      for(int j=0;j<n1;j++){
        acum= acum+  (vectorw(i)*vectorw(j))*(d11(i,j)+d12(i,j)-2*d13(i,j));

      }
    }

    acum= acum*(1/double(n1))*(1/double(n1));
    resultados(r)= acum;
    // Rcout << "The value of v : " << resultados(r) << "\n";

  }

  int contar=0;

  for(int i=0; i<npermutaciones;i++){

    if(resultados(i)>=resfinal){
      contar= contar+1;

    }

  }


  double pvalor= double(contar)/double(npermutaciones);



  return pvalor;

  // return List::create(Named("coefficients") = coef,Named("stderr") = stderrest);
}




// [[Rcpp::export]]
double energypairedwildboostrapmissing(arma::mat X1, arma::mat Y1, arma::vec w, int ker_type){

  int n1= X1.n_rows;
  // Rcpp::Rcout << "The value is " << contador2 << std::endl;

  arma::mat X= arma::join_cols(X1,Y1);

  // Rcpp::Rcout << "The value is funciona " << contador2 << std::endl;

  float sx=median_heur(X);
  // Kernel kernel_X(sx);
  Kernel *kernel_X;
  if (ker_type == 1) {
    kernel_X = new RBFKernel(sx);
  } else {
    kernel_X = new EnergyKernel(1);
  }

  // Rcpp::Rcout << "The value is kernel " << contador2 << std::endl;

  // arma::mat d11=  kernel_X(X1,X1);
  // arma::mat d12=  kernel_X(Y1,Y1);
  // arma::mat d13=  kernel_X(X1,Y1);
  arma::mat d11 =  (*kernel_X).getKernelMatrix(X1,X1);
  arma::mat d12 =  (*kernel_X).getKernelMatrix(Y1,Y1);
  arma::mat d13 =  (*kernel_X).getKernelMatrix(X1,Y1);

  double acum=0;
  w= w/sum(w);
  for(int i=0;i<n1;i++){
    for(int j=0;j<n1;j++){
      acum= acum+ (w(i)*w(j))*(d11(i,j)+d12(i,j)-2*d13(i,j));

    }
  }


  double resfinal= acum;
  //double ln= log(n1);
  //double ln=sqrt(n1);
  double ln=0.01;
  int npermutaciones=2000;

  arma::mat matrizW=wild(ln,npermutaciones,n1);

  arma::vec resultados(npermutaciones);

  arma::rowvec vectorw(npermutaciones);

  // Rcout << "The value of v : " << vectorw << "\n";

  for(int r=0;r<npermutaciones;r++){
    acum=0;
    // Rcout << "The value of v : " << matrizW.row(r) << "\n";

    vectorw=matrizW.row(r);
    // Rcout << "The value of v : " << vectorw << "\n";

    for(int i=0;i<n1;i++){
      for(int j=0;j<n1;j++){
        acum= acum+(w(i)*w(j))*(vectorw(i)*vectorw(j))*(d11(i,j)+d12(i,j)-2*d13(i,j));
      }
    }

    acum= acum;
    // Rcout << "The value of v : " << resultados(r) << "\n";

    // arma::mat d11=  kernel_X(X1,X1);
    // arma::mat d12=  kernel_X(Y1,Y1);
    // arma::mat d13=  kernel_X(X1,Y1);
    // Rcpp::Rcout << "The value is permutaciones " << permutaciones << std::endl;

    resultados(r)= acum;
  }

  int contarn=0;
  for(int i=0; i<npermutaciones;i++){
    if(resultados(i)>=resfinal){
      contarn= contarn+1;
    }
  }

  double pvalor= double(contarn)/double(npermutaciones);

  return pvalor;
  // return List::create(Named("coefficients") = coef,Named("stderr") = stderrest);
}

