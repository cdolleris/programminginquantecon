// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double probit_neg_log_lik_cpp(vec beta, mat X, vec Y) {
  int n = Y.size();
  double dSum = 0.0;
  vec linear_term = X * beta;
  vec dProb = zeros<vec>(X.n_rows);
    
  for (int i = 0; i < n; i++) {
    dProb(i) = R::pnorm(linear_term(i), 0.0, 1.0, 1, 0);
    dSum += Y[i] * log(dProb[i]) + (1 - Y[i]) * log(1 - dProb[i]);
  }
  return(-dSum);
}

// [[Rcpp::export]]
double fGARCH11NegLogLik_cpp(vec vParams, vec vReturns, double sigma2_initial) {
  double dOmega = vParams[0];
  double dAlpha = vParams[1];
  double dBeta = vParams[2];
  int n = vReturns.size();
  
  double dSum = 0;
  vec vSigma2 = zeros<vec>(vReturns.size());
  vSigma2[0] = sigma2_initial;
  double dPi = atan(1)*4;
  
  for (int t = 1; t < n; t++) {
    vSigma2[t] = dOmega + dAlpha * pow(vReturns[t-1], 2) + dBeta * vSigma2[t-1];
    dSum = dSum - 0.5*(log(2*dPi) + log(vSigma2[t]) + pow(vReturns[t], 2) / vSigma2[t]);
  }
  return(-dSum);
}

// [[Rcpp::export]]
double fVAR1NegLogLik(vec vParams, mat mY) {
  vec vC = zeros<vec>(2);
  vC[0] = vParams[0];
  vC[1] = vParams[1];
  mat mA(2,2);
  mA(0,0) = vParams[2];
  mA(1,0) = vParams[3];
  mA(0,1) = vParams[4];
  mA(1,1) = vParams[5];
  mat mSigma(2,2);
  mSigma(0,0) = exp(vParams[6]);
  mSigma(0,1) = 0;
  mSigma(1,0) = vParams[7];
  mSigma(1,1) = exp(vParams[8]);
  mSigma = mSigma * mSigma.t();
  double iK = mY.n_cols;
  double dPi = atan(1)*4;
  
  double dSum = 0.0;
  for (arma::uword t = 1; t < mY.n_rows; t++) {
    vec u_t = mY.row(t).t() - vC - mA * mY.row(t-1).t();
    dSum = dSum - iK / 2.0 * log(2.0 * dPi) - 0.5 * log(det(mSigma)) - 0.5 * as_scalar(trans(u_t) * mSigma.i() * (u_t));
  }
  
  return(-dSum);
}

// [[Rcpp::export]]
double logistic_neg_log_lik_cpp(vec vBeta, mat mX, vec vY) {
  double dSum = 0.0;
  vec linear_term = mX * vBeta;
  int n = vY.size();
  for (int i = 0; i < n; i++) {
    dSum += (vY(i) * linear_term(i)) - log(1 + exp(linear_term(i)));
  }
  return(dSum);
}

