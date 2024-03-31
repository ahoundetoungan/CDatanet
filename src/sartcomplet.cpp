// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// simulate tobit 
// complete information
//[[Rcpp::export]]
int fyTobit(arma::vec& yst,
            arma::vec& y,
            arma::vec& Gy,
            arma::vec& Ztlamda,
            List& G,
            const arma::vec& eps,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::vec& psi,
            const int& n,
            const double& lambda,
            const double& tol,
            const int& maxit) {
  int n1, n2, t = 0;
  arma::vec vzero        = arma::zeros(n);
  
  computeyst: ++t;
  arma::vec y0           = y;
  Ztlamda.subvec(0, n-1) = lambda*Gy + psi;
  yst.subvec(0, n-1)     = Ztlamda + eps;
  y.subvec(0, n-1)       = arma::max(arma::join_rows(vzero, yst), 1);
  for (int m(0); m < ngroup; ++ m) {
    n1                = igroup(m,0);
    n2                = igroup(m,1);
    arma::mat Gm      = G[m];
    Gy.subvec(n1, n2) = Gm*y.subvec(n1, n2);
  }
  double dist         = max(arma::abs((y - y0)/(y0 + 1e-50)));
  if (dist > tol && t < maxit) goto computeyst;
  return t; 
}

// [[Rcpp::export]]
double foptimTobit(const arma::vec& theta,
                   const arma::mat& X,
                   arma::vec& logdetA2,
                   arma::vec& alphatilde,
                   List& G2,
                   List& I2,
                   const int& K,
                   const arma::vec& y,
                   const arma::vec& Gy,
                   const arma::uvec& idpos,
                   const arma::uvec& idzero,
                   const int& npos,
                   const int& ngroup,
                   List& I,
                   List& W,
                   const int& n, 
                   const arma::mat igroup,
                   const bool print = false){
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1.0);
  double sigma         = exp(theta(K + 1));
  arma::vec tmp        = alpha*Gy + xb;
  arma::vec tmpzer     = tmp.elem(idzero);
  NumericVector tzcpp  = wrap(tmpzer);
  tmp                  = y - tmp;
  if(print){
    NumericVector thetacpp = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::Rcout<<"---------------"<< endl;
    Rcpp::Rcout<<"Estimate: "<< endl;
    thetacpp(0)          = alpha;
    thetacpp(K + 1)      = sigma;
    Rcpp::print(thetacpp);
  }
  
  if(alphatilde(0) != theta(0)) {
    logdetA2(0)     = 0;
    for (int i(0); i < ngroup; ++ i) {
      double vali, signi;
      arma::mat G2i = G2[i];
      arma::mat I2i = I2[i];
      log_det(vali, signi, I2i - alpha*G2i);
      logdetA2(0)  += vali;
      logdetA2(0)  += log(signi);
    }
  }
  
  alphatilde(0)     = theta(0);
  
  
  double llh        = sum(Rcpp::pnorm(tzcpp, 0, sigma, false, true)) -
   npos*(0.5*log(2*acos(-1)) + log(sigma)) + logdetA2(0) - 0.5*sum(pow(tmp.elem(idpos)/sigma, 2));
  if(llh < -1e293) {
    llh           = -1e293;
  }
  if(print) Rcpp::Rcout <<"Likelihood: "<< llh << endl;
  return -llh;
}


//[[Rcpp::export]]
arma::vec fgradvecTobit(arma::vec& theta,
                        arma::mat& X,
                        const arma::vec& logdetA2,
                        const arma::vec& alphatilde,
                        List& G2,
                        List& I2,
                        const int& K,
                        const arma::vec& y,
                        const arma::vec& Gy,
                        const arma::uvec& idpos,
                        const arma::uvec& idzero,
                        const int& ngroup,
                        List& I,
                        List& W,
                        const int& n, 
                        const arma::vec& indzero,
                        const arma::vec& indpos,
                        const arma::mat igroup) {
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1.0);
  
  double sigma         = exp(theta(K + 1)) + 1e-8;
  arma::vec tmp        = y - alpha*Gy - xb;
  NumericVector tmpcpp = wrap(tmp);
  
  NumericVector irmcpp = exp(Rcpp::dnorm(tmpcpp/sigma, 0, 1, true) - Rcpp::pnorm(tmpcpp/sigma, 0, 1, true, true));
  arma::vec     irm    = as<arma::vec>(irmcpp);
  
  arma::vec rvec(n);
  
  for (int i(0); i < ngroup; ++ i) {
    arma::mat Wi   = W[i];
    arma::mat Ii   = I[i];
    arma::mat Rmat = arma::inv(Ii - alpha*Wi)*Wi;
    rvec.subvec(igroup(i,0), igroup(i,1)) = arma::diagvec(Rmat);
  }
  
  arma::mat out(n, K + 2);
  
  // derivation with respect alpha
  out.col(0)       = ((indzero%(irm)/sigma - indpos%tmp/pow(sigma, 2))%Gy + rvec)*exp(theta(0))/pow(1 + exp(theta(0)), 2);
  out.cols(1, K)   = arma::repmat(indzero%irm/sigma - indpos%tmp/pow(sigma,2), 1, K)%X;
  out.col(K + 1)   = (-indzero%irm%(alpha*Gy + xb)/pow(sigma, 2) + indpos%(1 - pow(tmp/sigma, 2))/sigma)*sigma;
  
  NumericVector outcpp = wrap(arma::trans(sum(out, 0)));
  outcpp.attr("dim") = R_NilValue;
  
  return outcpp;
}


//[[Rcpp::export]]
List fcovSTC(const arma::vec& theta,
                  const arma::mat& X,
                  List& G2,
                  List& I,
                  List& W,
                  const int& K,
                  const int& n,
                  const arma::vec& y,
                  const arma::vec& Gy,
                  const arma::vec& indzero,
                  const arma::vec& indpos,
                  const arma::mat& igroup,
                  const int& ngroup,
                  const bool& ccov) {
  List out;
  double lambda        = 1.0/(exp(-theta(0)) + 1);
  double sigma         = exp(theta(K + 1));
  arma::vec ZtL        = lambda*Gy + X*theta.subvec(1, K); 
  NumericVector ZtLst  = wrap(ZtL/sigma);
  NumericVector PhiZtLst(Rcpp::pnorm(ZtLst, 0, 1, true, false));
  double avPhiZtLst    = sum(PhiZtLst)/n;
  arma::vec lbeta      = arma::join_cols(arma::ones(1)*lambda, theta.subvec(1, K));
  arma::mat meffects   = as<arma::vec>(PhiZtLst)*lbeta.t();
  
  if(ccov) {
    arma::vec tmp        = y - ZtL;
    NumericVector tmpcpp = wrap(tmp);
    NumericVector irmcpp = exp(Rcpp::dnorm(tmpcpp/sigma, 0, 1, true) - Rcpp::pnorm(tmpcpp/sigma, 0, 1, true, true));
    arma::vec     irm    = as<arma::vec>(irmcpp);
    
    arma::vec rvec(n);
    
    for (int i(0); i < ngroup; ++ i) {
      arma::mat Wi   = W[i];
      arma::mat Ii   = I[i];
      arma::mat Rmat = arma::inv(Ii - lambda*Wi)*Wi;
      rvec.subvec(igroup(i,0), igroup(i,1)) = arma::diagvec(Rmat);
    }
    
    arma::mat qvec(n, K + 2);
    
    qvec.col(0)      = ((-indzero%(irm)/sigma + indpos%tmp/pow(sigma, 2))%Gy - rvec);
    qvec.cols(1, K)  = arma::repmat(-indzero%irm/sigma + indpos%tmp/pow(sigma, 2), 1, K)%X;
    qvec.col(K + 1)  = (indzero%irm%ZtL/pow(sigma, 2) - indpos%(1 - pow(tmp/sigma, 2))/sigma);
    arma::mat covt   = arma::inv(arma::cov(qvec))/n;
    
    // cov marginal effects
    NumericVector phiZtLst   = Rcpp::dnorm4(ZtLst, 0, 1, false);
    arma::mat Z              = arma::join_rows(Gy, X);
    arma::rowvec ZavphiZtLst = arma::mean(Z.each_col()%as<arma::vec>(phiZtLst), 0);
    
    arma::mat tmp1 = arma::eye<arma::mat>(K + 1, K + 1)*avPhiZtLst + lbeta*ZavphiZtLst/sigma;
    arma::vec tmp2 = -lbeta*mean(ZtLst*phiZtLst)/sigma;
    arma::mat tmp3 = arma::join_rows(tmp1, tmp2);
    arma::mat covm = tmp3*covt*tmp3.t();
    
    out = List::create(Named("meff")  = meffects,
                       Named("covt")  = covt,
                       Named("covm") = covm);
  } else {
    out = List::create(Named("meff") = meffects);
  }

  return out;
}
