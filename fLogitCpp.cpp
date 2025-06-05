//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//initialize auxiliary functions: prb, LL, score, hessian

//
vec fP(vec vBeta, mat mX){
  return(1 / (1 + exp( -mX * vBeta )));
}
//
double fLnL(vec vBeta, vec vY, mat mX){
  vec vP = fP(vBeta,mX);
  double vLnL = mean( vY % log(vP) + (ones<vec>(vY.n_elem) - vY) % log(ones<vec>(vP.n_elem) - vP));
  return vLnL;
}
//
vec fScore(vec vBeta, vec vY, mat mX){
  vec vP = fP(vBeta,mX);
  vec vScore = trans(mX)*(vY - vP) / vY.n_elem;
  
  return vScore;
}
//
mat fHessian(vec vBeta, vec vY, mat mX){
  int iN = vY.n_elem;
  vec vP = fP(vBeta,mX);
  mat mP = zeros<mat>(iN,iN);
  
  mP.diag(0) = vP % (ones<vec>(iN) - vP);
  mat mHessian = -trans(mX)*mP*mX / iN;
  
  return mHessian;
}

//[[Rcpp::export]]
List fLogitCpp(vec vY, mat mX, vec vBeta0, 
               bool constant = true, 
               double dTol = 1e-9, 
               int imax = 200){
  
  List lOut;
  int iN = vY.n_elem;
  
  if (constant == true){
    vec vOnes = ones<vec>(iN);
    vBeta0 = join_cols(zeros<vec>(1),vBeta0);
    mX = join_rows(vOnes,mX);
    //Rcout << "Note: A constant has been added to data matrix mX" << std::endl;
  }
  
  vec vBeta = vBeta0;
  vec vScore = fScore(vBeta,vY,mX);
  mat mHessian = fHessian(vBeta,vY,mX);
  
  int i = 0;
  while ((max(abs(vScore)) > dTol) & (i < imax)) {
    //Rcout << i << "\n";
    vBeta = vBeta - solve(fHessian(vBeta,vY,mX),fScore(vBeta,vY,mX));
    vScore = fScore(vBeta,vY,mX);
    mHessian = fHessian(vBeta,vY,mX);
    
    i++;
  }
  if (i == imax){
    Rcout << "newton failed to converge" << std::endl;
    return lOut;
  } else {
   // Rcout << "Convergence achieved after " << i << " iterations" << std::endl
   //       << "due to max(abs(score)) <" << dTol << std::endl;
    lOut["coefficients"] = vBeta;
    lOut["logLik0"] = fLnL(zeros<vec>(vBeta.n_elem),vY,mX);
    lOut["logLik"] = fLnL(vBeta,vY,mX);
    lOut["score"] = vScore;
    lOut["hessian"] = mHessian;
    lOut["phat"] = fP(vBeta,mX);
    
    return lOut;
  }
}


