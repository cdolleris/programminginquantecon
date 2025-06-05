// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double f_Y1(double y1, vec vParams) {
  double dPi = atan(1)*4;
  double dMu = vParams[0];
  double dPhi = vParams[1];
  double dSigma2 = vParams[2];
  
  double dOut = 1 / pow(2 * dPi * dSigma2 / (1 - pow(dPhi, 2)), 0.5) * exp(-0.5 * pow((y1 - dMu / (1 - dPhi)), 2) / (dSigma2 / (1 - dPhi)));
  
  return dOut;
}

// [[Rcpp::export]]
double f_cond(double y1, double yt_1, vec vParams) {
  double dPi = atan(1)*4;
  double dMu = vParams[0];
  double dPhi = vParams[1];
  double dSigma2 = vParams[2];
  
  double dOut = 1 / pow(2 * dPi * dSigma2, 0.5) * exp(-0.5 * pow((y1 - dMu - dPhi * yt_1), 2) / dSigma2);
  
  return dOut;
}

// [[Rcpp::export]]
arma::vec simulateAR1_buggy(double mu_true = 0.5, double phi_true = 0.7, double sigma2_true = 1.5, double T_val = 500) {
  
  if (abs(phi_true) >= 1) {
    stop("phi must be between -1 and 1 for stationarity.");
  }
  if (sigma2_true <= 0) {
    stop("sigma2 must be positive.");
  }
  
  arma::vec y = zeros<vec>(T_val);
  
  double mean_y1 = mu_true / (1 - phi_true);
  double var_y1 = sigma2_true / (1 - pow(phi_true, 2));
  double sd_y1 = sqrt(var_y1);
  
  y[0] = Rf_rnorm(mean_y1, sd_y1);
  
  double sd_epsilon = sqrt(sigma2_true);
  double epsilon_t = 0.0;
  for (int t = 1; t < T_val; t++) {
    epsilon_t = Rf_rnorm(0, sd_epsilon);
    y[t] = mu_true + phi_true * y[t-1] + epsilon_t;
  }
  
  return(y);
}