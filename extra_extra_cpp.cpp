// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double poisson_log_lik_cpp(int k, double lambda) {
  return (k * log(lambda) - lambda - tgamma(k + 1));
}

//[[Rcpp::export]]
double total_poisson_log_lik_cpp(Rcpp::IntegerVector counts, double lambda) {
  double dSum = 0.0;
  for (int i = 0; i < counts.size(); i++) {
    dSum += poisson_log_lik_cpp(counts[i], lambda);
  }
  return dSum;
}