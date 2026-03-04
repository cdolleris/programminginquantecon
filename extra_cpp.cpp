// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List arma_matrix_ops(mat mA, mat mB) {
  List lOut;
  lOut["mProd"] = mA * mB;
  arma::mat mInvA;
  bool success = inv(mInvA, mA); // Try to invert
  if (!success) {
    Rcpp::Rcout << "Matrix mA is singular." << std::endl;
    // mInvA will be empty or you can fill with NAs if returning to R
  } else {
    lOut["mInvA"] = mA.i();
  }
  lOut["vEigenvalsA"] = eig_gen(mA);
  return lOut;
}

// [[Rcpp::export]]
List simulate_garch_cpp(int iT, double dOmega, double dAlpha, double dBeta, double dSigmaInit) {
  
  vec vY(iT); //observations
  vec vSigma2(iT); // conditional variances
  
  // initialize at the unconditional value
  vSigma2(0) = dSigmaInit;
  
  // sample the first observation
  vY(0) = pow(vSigma2(0), 0.5) * Rf_rnorm(0.0, 1.0);
  
  for (int t = 1; t < iT; t++) {
    vSigma2(t) = dOmega + dAlpha * pow(vY(t - 1), 2.0) + dBeta * vSigma2(t - 1);
    vY(t) = pow(vSigma2(t), 0.5) * Rf_rnorm(0.0, 1.0);
  }
  
  List lOut;
  lOut["vSigma2"] = vSigma2;
  lOut["vR"] = vY;
  
  return lOut;
}

//[[Rcpp::export]]
vec cpp_vector_scalar_mult(vec v, double s) {
  return v * s;
}

/*** R
vX <- 1:10
dS <- 5
res <- cpp_vector_scalar_mult(vX, dS)
print(res)
*/