/* GENERAL NOTATIONS
 * The notations used are consistent with the paper.
 * The index used if the comments are consistent with cpp spirit.
 * for example X(0) is used for the first entry of X and X(1) the second entry.
 * This code is free and can be customized for other uses.
 * y       : is the vector of the outcome values
 * X       : is the matrix of the explanatory variables. Add the intercept 
 *           if it is included in the model. The intercept will not be added automatically. 
 * G       : is the network matrix List. That is G[s] is the subnetwork of the group s. 
 *           I work with the row normalization version.
 *           Gs(i,j) = 1/ni if i knows j, where ni is the number of friend i has
 *           and Gs(i,j) = 0 otherwise. 
 * igroup  : is the matrix of groups indexes. The data should be ordered by the group. 
 *           igroup[s,] is a 2-dimension vector ans gives the firt and the last rows
 *           for the group s
 * ngroup  : is the number of groups.
 * theta   : is the vector of parameters ordered as follow: peer effects, explanatory variables
 *           coefficients, delta2, delta3 ..., deltaRbar.
 * n       : The sample size.
 * tol     : A tolerance value for the iterative method to compute P convergence. 
 * maxit   : The maximum number of iterations of the iterative method to compute P. If this
 *           number is reached, the algorithm stops and the last P is used as the solution. maxit
 *           is important for numerical reasons if tol is too small. For example a tol = 1e-16
 *           may not be reached and the algorithm will stop after maxit iterations.
 * yb      : is the vector of equilibrium outcome expectation.
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::vec fLTBT(const NumericVector& ZtLambda,
                 const double sigma) {
  return ZtLambda*Rcpp::pnorm5(ZtLambda/sigma, 0, 1, true, false) + 
    sigma*Rcpp::dnorm4(ZtLambda/sigma, 0, 1, false);;
}

// fyb: Takes an initial value of yb and finds the equilibrium
// Gyb is G*yb
//[[Rcpp::export]]
int fybtbit(arma::vec& yb,
            arma::vec& Gyb,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::vec& psi,
            const double& lambda,
            const double& sigma,
            const int& n, 
            const double& tol,
            const int& maxit) {
  int n1, n2, t = 0;
  arma::vec ZtLambda(n);
  
  computeL: ++t;
  ZtLambda       = lambda*Gyb + psi;
  arma::vec ybst = fLTBT(wrap(ZtLambda), sigma);
  double dist    = arma::accu(arma::abs(yb - ybst));
  yb             = ybst;
  
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1);
    arma::mat Gm       = G[m];
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
  if (dist > tol && t < maxit) goto computeL;
  return t; 
}

// foptimREM: compute the likelihood given the patameters
//[[Rcpp::export]]
double foptimRE_TBT(arma::vec& yb,
                 arma::vec& Gyb,
                 const arma::vec& theta,
                 const arma::vec& yidpos,
                 const arma::mat& X,
                 List& G,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const int& npos,
                 const arma::uvec& idpos,
                 const arma::uvec& idzero,
                 const int& K,
                 const int& n,
                 const double& tol = 1e-13,
                 const int& maxit  = 1e3) {

  double lambda   = 1.0/(exp(-theta(0)) + 1);
  double sigma    = exp(theta(K + 1));
  arma::vec psi   = X * theta.subvec(1, K);
  
  // compute ybar
  fybtbit(yb, Gyb, G, igroup, ngroup, psi, lambda, sigma, n, tol, maxit);
  
  arma::vec ZtLambda = lambda*Gyb + psi;
  arma::vec ZtL0     = ZtLambda.elem(idzero);
  NumericVector tmp  = wrap(ZtL0);
  double llh         = sum(Rcpp::pnorm5(tmp, 0, sigma, false, true)) -
    npos*(0.5*log(2*acos(-1)) + log(sigma))  - 0.5*sum(pow((yidpos - ZtLambda.elem(idpos))/sigma, 2));
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}

//// NPL solution
//[[Rcpp::export]]
double foptimTBT_NPL(const arma::vec& yidpos,
                     const arma::vec& Gyb,
                     const arma::mat& X,
                     const arma::vec& theta,
                     const int& npos,
                     const arma::uvec& idpos,
                     const arma::uvec& idzero,
                     const int& K) {
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  double sigma       = exp(theta(K + 1));
  arma::vec ZtLambda = lambda*Gyb + X*theta.subvec(1, K);
  arma::vec ZtL0     = ZtLambda.elem(idzero);
  NumericVector tmp  = wrap(ZtL0);
  double llh         = sum(Rcpp::pnorm5(tmp, 0, sigma, false, true)) -
    npos*(0.5*log(2*acos(-1)) + log(sigma))  - 0.5*sum(pow((yidpos - ZtLambda.elem(idpos))/sigma, 2));
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}

//[[Rcpp::export]]
void fLTBT_NPL(arma::vec& yb,
            arma::vec& Gyb,
            List& G,
            const arma::mat& X,
            const arma::vec& theta,
            const arma::mat& igroup,
            const int& ngroup,
            const int& n,
            const int& K) {
  int n1, n2;
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  double sigma       = exp(theta(K + 1));
  arma::vec ZtLambda = lambda*Gyb + X*theta.subvec(1, K);
  
  // new yb
  yb.subvec(0, n-1)  = fLTBT(wrap(ZtLambda), sigma);
  
  // new Gyb
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1);
    arma::mat Gm       = G[m];
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
}

//[[Rcpp::export]]
void fnewybTBT(arma::vec& yb,
               arma::vec& Gyb,
               List& G,
               const arma::mat& igroup,
               const int& ngroup,
               const arma::mat& X,
               const arma::vec& theta,
               const int& K,
               const int& n, 
               const double& tol,
               const int& maxit) {
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  double sigma       = exp(theta(K + 1));
  arma::vec psi      = X*theta.subvec(1, K);

  fybtbit(yb, Gyb, G, igroup, ngroup, psi, lambda, sigma, n, tol, maxit);
}



// variance
//[[Rcpp::export]]
List fcovSTI(const int& n,
             const arma::vec& Gyb,
             const arma::vec& theta,
             const arma::mat& X,
             const int& K,
             List& G,
             const arma::mat& igroup,
             const int& ngroup,
             const bool& ccov) {
  List out;
  double lambda          = 1.0/(exp(-theta(0)) + 1);
  double sigma           = exp(theta(K + 1));
  arma::vec ZtL          = lambda*Gyb + X*theta.subvec(1, K); 
  NumericVector ZtLst    = wrap(ZtL/sigma);
  NumericVector PhiZtLst = Rcpp::pnorm(ZtLst, 0, 1, true, false);
  double avPhiZtLst      = sum(PhiZtLst)/n;
  arma::vec lbeta        = arma::join_cols(arma::ones(1)*lambda, theta.subvec(1, K));
  arma::vec meffects     = avPhiZtLst*lbeta;
  
  if(ccov){
    arma::mat Z                = arma::join_rows(Gyb, X);
    NumericVector phiZtLst     = Rcpp::dnorm4(ZtLst, 0, 1, false); 
    NumericVector Zst1phiZtLst = ZtLst*phiZtLst;
    NumericVector Zst2phiZtLst = ZtLst*Zst1phiZtLst;
    NumericVector Zst3phiZtLst = ZtLst*Zst2phiZtLst;
    NumericVector RphiPhi      = phiZtLst/(1 - PhiZtLst);
    
    // compute GinvSW
    arma::mat GinvSW(n, K + 2);
    {// Compute b and d
      arma::vec d = as<arma::vec>(PhiZtLst); 
      arma::vec b = as<arma::vec>(phiZtLst); 
      
      // G*inv(S)*W
      arma::mat W = arma::join_rows(Z.each_col() % d, b/sigma);
      for (int m(0); m < ngroup; ++ m) {
        int n1              = igroup(m,0);
        int n2              = igroup(m,1);
        int nm              = n2 - n1 + 1;
        arma::mat Gm        = G[m];
        arma::mat Sm        = Gm.each_col() % d.subvec(n1, n2);
        Sm                  = arma::eye<arma::mat>(nm, nm) - lambda*Sm;
        GinvSW.rows(n1, n2) = Gm * arma::solve(Sm, W.rows(n1, n2));
      }
    }
  
    // Compute Sigma0 and Omega0
    arma::mat Sigma(K + 2, K + 2), Omega(K + 2, K + 2);
    {
    NumericVector tmp1 = Zst1phiZtLst - phiZtLst*RphiPhi - PhiZtLst; 
    NumericVector tmp2 = -phiZtLst - Zst2phiZtLst + RphiPhi*Zst1phiZtLst; 
    arma::mat Ztmp1    = arma::trans(Z.each_col()%as<arma::vec>(tmp1))/sigma;
    arma::vec Ztmp2    = Z.t()*as<arma::vec>(tmp2)/sigma;

    Sigma.submat(0, 0, K, K)         = Ztmp1*Z;
    Omega.rows(0, K)                 = Ztmp1*GinvSW;
    
    Sigma.submat(0, K + 1, K, K + 1) = Ztmp2;
    Sigma.submat(K + 1, 0, K + 1, K) = Ztmp2.t();
    Omega.row(K + 1)                 = arma::trans(as<arma::vec>(tmp2))*GinvSW;
    
    Sigma(K + 1, K + 1)              = sum(Zst1phiZtLst + Zst3phiZtLst - Zst2phiZtLst*RphiPhi - 2*PhiZtLst);
    }
    
    Omega          = Omega*lambda/sigma;
    // covt
    arma::mat covt = arma::inv(Sigma + Omega);
    covt           = -covt*Sigma*covt.t();
    // covm
    arma::rowvec ZavphiZtLst = arma::mean(Z.each_col()%as<arma::vec>(phiZtLst), 0);
    arma::mat tmp1 = arma::eye<arma::mat>(K + 1, K + 1)*avPhiZtLst + lbeta*ZavphiZtLst/sigma;
    arma::vec tmp2 = -lbeta*mean(ZtLst*phiZtLst);
    arma::mat tmp3 = arma::join_rows(tmp1, tmp2);
    arma::mat covm = tmp3*covt*tmp3.t();

    out            = List::create(Named("meffects")    = meffects,
                                  Named("covtheta")    = covt,
                                  Named("covmeffects") = covm,
                                  Named("var.comp")    = List::create(Named("Sigma") = Sigma/n, Named("Omega") = Omega/n));
  } else{
    out            = List::create(Named("meffects")    = meffects);
  }
  
  return out;
}
