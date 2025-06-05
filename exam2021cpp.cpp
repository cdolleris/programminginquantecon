// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
vec transToOrig_cpp(vec vInput) {
  vec vOut = zeros<vec>(vInput.size());
  vOut[0] = exp(vInput[0]);
  vOut[1] = exp(vInput[1]) / (1 + exp(vInput[1]) + exp(vInput[2]));
  vOut[2] = exp(vInput[2]) / (1 + exp(vInput[1]) + exp(vInput[2]));
  vOut[3] = exp(vInput[3]);
  vOut[4] = exp(vInput[4]) / (1 + exp(vInput[4]) + exp(vInput[5]));
  vOut[5] = exp(vInput[5]) / (1 + exp(vInput[4]) + exp(vInput[5]));
  vOut[6] = -1 + 2 * (exp(vInput[6])) / (1 + exp(vInput[6]));
  
  return vOut;
}