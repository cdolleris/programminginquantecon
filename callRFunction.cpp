#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
SEXP callRFunction(NumericVector vX, Function f) {
  SEXP res = f(vX);
  return res;
}