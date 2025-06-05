// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double log_n_choose_k_cpp(int n, int k) {
  if (k == 0) {
    return(0);
  } else {
    double dSumPart1 = 0;
    double dSumPart2 = 0;
    for (int i = (n - k + 1); i <= n; i++) {
      dSumPart1 = dSumPart1 + log(i);
    }
    for (int i = 1; i <= k; i++) {
      dSumPart2 = dSumPart2 + log(i);
    }
    return(dSumPart1 - dSumPart2);
  }
}

// [[Rcpp::export]]
double fLnL_cpp(double p, int n, vec d) {
  int len = d.n_elem;
  vec result(len);
  
  for (int i = 0; i < len; i++) {
    result[i] = log_n_choose_k_cpp(n, d[i]) + d[i] * log(p) + (n - d[i]) * log(1 - p);
  }
  return(mean(result));
}

// [[Rcpp::export]]
double fScore_cpp(double p, int n, vec d) {
  return(mean(d / p - (n - d) / (1 - p)));
}

// [[Rcpp::export]]
double fSD_cpp(double p, int n, vec d) {
  return(mean(-(d / pow(p, 2)) - (n - d) / pow((1 - p), 2)));
}