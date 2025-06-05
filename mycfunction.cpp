// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
mat FunC(mat mx, mat mY) {
  mat mZ = mx * mY;
  mat mZInv = mZ.i();
  return mZInv;
}