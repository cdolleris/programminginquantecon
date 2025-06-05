// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double factorial_cpp(int iN) {
  if (iN < 0) {
    Rcpp::stop("Input to factorial must be a non-negative integer.");
  }
  if (iN == 0) {
    return 1.0;
  }
  double result = 1.0;
  for (int i = 1; i <= iN; ++i) {
    result *= i;
  }
  return result;
}

// [[Rcpp::export]]
double fLnL_cpp(vec vParams, vec vY, mat mX) {
  double dSum = 0.0;
  int iN = mX.n_rows;
  
  for (int i = 0; i < iN; i++) {
    //dSum = dSum + vY[i] * as_scalar(mX.row(i) * vParams) - exp(as_scalar(mX.row(i) * vParams)) - log(factorial_cpp(vY[i]));
    dSum = dSum + vY[i] * as_scalar(mX.row(i) * vParams) - exp(as_scalar(mX.row(i) * vParams)) - lgamma(vY[i] + 1.0);
  }
  
  return(dSum / iN);
}

// [[Rcpp::export]]
vec fScore_2_cpp(vec vParams, vec vY, mat mX){
  vec vP = exp(mX*vParams);
  vec vScore = trans(mX)*(vY - vP) / vY.n_elem;
  
  return vScore;
}

// [[Rcpp::export]]
vec fScore_cpp(vec vParams, vec vY, mat mX){
  
  int n_obs = mX.n_rows; // Number of observations
  int p = mX.n_cols;     // Number of parameters
  
  vec vOut = zeros<vec>(p);
  
  for (int i = 0; i < n_obs; i++) {
    vOut += (vY[i] - exp(as_scalar(mX.row(i) * vParams))) * mX.row(i).t();
  }
  
  return vOut / n_obs;
}

