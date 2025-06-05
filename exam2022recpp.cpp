// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec CppGfunc(Function g, arma::vec x) {
  
  arma::vec vOut = zeros<vec>(x.size());
  
  for (arma::uword i = 0; i < x.size(); i++) {
    SEXP s;
    s = g(x[i]);
    vOut[i] = as<double>(s) + 2.0;
  }
  
  return(vOut);
}