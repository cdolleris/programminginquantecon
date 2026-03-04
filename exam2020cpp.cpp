// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List bisection_cpp(Function f, double dXLeft = -1.0, double dXRight = 1.0, double dTol = 0.0001, int maxIter = 1000) {
  if (dXLeft >= dXRight) {
    stop("Starting conditions not met");
  }
  double fLeft = as<double>(f(dXLeft));
  double fRight = as<double>(f(dXRight));
  if (fLeft == 0) {
    List lOut;
    lOut["root"] = dXLeft;
    lOut["func_at_root"] = fLeft;
    lOut["num_ite"] = 0;
    return lOut;
  } else if (fRight == 0) {
    List lOut;
    lOut["root"] = dXRight;
    lOut["func_at_root"] = fRight;
    lOut["num_ite"] = 0;
    return lOut;
  } else if (fLeft * fRight > 0) {
    stop("error: f(x.l)*f(x.r) > 0");
  }
  
  int iter = 0;
  while ((dXRight - dXLeft) > dTol && (iter < maxIter)) {
    double dXMid = (dXLeft + dXRight)/2;
    double fMid = as<double>(f(dXMid));
    if (fMid == 0) {
      return(dXMid);
    } else if (fLeft * fMid < 0) {
      dXRight = dXMid;
      fRight = fMid;
    } else {
      dXLeft = dXMid;
      fLeft = fMid;
    }
    iter = iter + 1;
  }
  
  List lOut;
  lOut["root"] = (dXLeft + dXRight)/2;
  lOut["func_at_root"] = as<double>(f((dXLeft + dXRight)/2));
  lOut["num_ite"] = iter;
  return lOut;
}