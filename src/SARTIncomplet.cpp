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
