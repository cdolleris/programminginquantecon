// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// Helper function to get a shock
inline double get_shock(Function f, const vec& vInputs) {
  SEXP s;
  int n_extra_args = vInputs.n_elem;
  
  if (n_extra_args == 0) { // e.g. f takes only n, like a pre-wrapped user function or rnorm(n) with defaults
    s = f(1);
  } else if (n_extra_args == 1) { // e.g. rt(n, df) or rnorm(n, mean)
    s = f(1, vInputs[0]);
  } else if (n_extra_args == 2) { // e.g. rnorm(n, mean, sd)
    s = f(1, vInputs[0], vInputs[1]);
  } else {
    // Fallback or error for more arguments, or could extend this if/else
    Rcpp::stop("vInputs has too many elements for this simplified handler. Max 2 supported.");
  }
  if (TYPEOF(s) != REALSXP || Rf_length(s) != 1) {
    Rcpp::stop("The user-supplied function 'f' did not return a single numeric value.");
  }
  return Rcpp::as<double>(s);
}


// [[Rcpp::export]]
List GarchSim(int iT_in, vec vParams, Function f, vec vInputs) {
  int iT = iT_in;
  if (iT <= 0) { iT = 1; } // Simplified default handling
  
  double dOmega = vParams[0];
  double dAlpha = vParams[1];
  double dBeta  = vParams[2];
  
  vec vY(iT, fill::zeros);
  vec vSigma2(iT, fill::zeros);
  
  if (dOmega <= 0 || dAlpha < 0 || dBeta < 0 || (dAlpha + dBeta) >= 1.0) {
    Rcpp::stop("Invalid GARCH parameters.");
  }
  
  double unconditional_variance = dOmega / (1.0 - dAlpha - dBeta);
  vSigma2(0) = unconditional_variance;
  
  double current_eps = get_shock(f, vInputs); // Call helper
  vY(0) = sqrt(vSigma2(0)) * current_eps;
  
  for (int t = 1; t < iT; t++) {
    vSigma2(t) = dOmega + dAlpha * pow(current_eps, 2) + dBeta * vSigma2(t-1);
    if (vSigma2(t) <= 0) vSigma2(t) = unconditional_variance; 
    
    current_eps = get_shock(f, vInputs); // Call helper
    vY(t) = sqrt(vSigma2(t)) * current_eps;
  }
  
  return List::create(_["vSigma2"] = vSigma2, _["vY"] = vY);
}