// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double fRosen_cpp(vec vX) {
  double pi = atan(1)*4; //pi
  return(-pow((pi - vX[0]), 2) - 100 * pow((vX[1] - pow(vX[0], 2)), 2));
}

// [[Rcpp::export]]
vec fRosen_grad_cpp(vec vX) {
  vec vOutput(2);
  double pi = atan(1)*4; //pi
  vOutput[0] = 2 * pi - 2 * vX[0] - 400 * pow(vX[0], 3) + 400 * vX[0] * vX[1];
  vOutput[1] = 200 * (pow(vX[0], 2) - vX[1]);
  return(vOutput);
}