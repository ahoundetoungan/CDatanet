/* GENERAL NOTATIONS
 * The notations used are consistent with the paper.
 * The index used if the comments are consistent with cpp spirit.
 * for example X(0) is used for the first entry of X and X(1) the second entry.
 * This code is free and can be customized for other uses.
 * y       : is the vector of the outcome values
 * X       : is the matrix of the explanatory variables. Add the intercept 
 *           if it is included in the model. The intercept will not be added automatically. 
 * G       : is the network matrix List. That is G[s] is the subnetwork of the group s. 
 *           I work with the row norlalization version.
 *           Gs(i,j) = 1/ni if i knows j, where ni is the number of friend i has
 *           and Gs(i,j) = 0 otherwise. 
 * igroup  : is the matrix of groups indexes. The data should be ordered by the group. 
 *           igroup[s,] is a 2-dimension vector ans gives the firt and the last rows
 *           for the group s
 * ngroup  : is the number of groups.
 * theta   : is the vector of parameters ordered as follow: peer effects, explanatory variables
 *           coefficients, sigma (not sigma^2).
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

/* LOCAL NOTATIONS
 *  the notation m in front of objects refers to the group m
 *  psi    : is X*beta + cte 
 *  lambda  : is the peer effect
 *  sigma  : is sigma as in the paper
 */

// fL: is the mapping L as in the paper.
//[[Rcpp::export]]
arma::vec fL(const arma::vec& um,
             const double& lambda,
             const double& sigma,
             const arma::vec& psim,
             const arma::mat& Gm,
             const int& nm) {
  double lF0, lFr, lu2, tmpi;
  int r;
  bool cont;
  
  arma::vec tmp = lambda*Gm*um + psim;
  arma::vec logu(nm);
  
  for (int i(0); i < nm; ++ i) {
    tmpi    = tmp(i);
    lF0     = R::pnorm(tmpi, 0, sigma, true, true);
    lu2     = 1;
    r       = 1;
    cont    = true;
    while(cont) {
      lFr   = R::pnorm(tmpi - r, 0, sigma, true, true);
      lu2  += exp(lFr - lF0); 
      ++ r;
      cont  = (lFr > -1000);
    }
    logu(i) = lF0 + log(lu2);
  }
  
  return exp(logu);
}


// fyb: Takes an initial value of yb and finds the equilibrium
// Gyb is G*yb
//[[Rcpp::export]]
arma::vec fyb(arma::vec& yb,
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
  int tm, nm;
  arma::vec Gybm, ybm, psim, t(ngroup); 
  //loop over group
  for (int m(0); m < ngroup; ++ m) {
    tm                                   = 0;
    nm                                   = igroup(m,1) - igroup(m,0) + 1;
    ybm                                  = yb.subvec(igroup(m,0), igroup(m,1));
    arma::mat Gm                         = G[m];
    psim                                 = psi.subvec(igroup(m,0), igroup(m,1));
    
    newstep: arma::vec tmp1 = fL(ybm, lambda, sigma, psim, Gm, nm);
    
    double tmp2                          = arma::accu(arma::abs(tmp1 - ybm));
    
    ++tm;
    ybm                                  = tmp1;
    if (tmp2 > tol && tm < maxit) goto newstep;
    
    yb.rows(igroup(m,0), igroup(m,1))    = ybm;
    Gyb.rows(igroup(m,0), igroup(m,1))   = Gm * ybm;

    t(m)                                 = tm;
  }
  return t;  // the number of iteration if needed
}

arma::vec flogp(const arma::vec& y,
                const arma::vec& yb,
                const arma::vec& Gyb,
                const arma::vec& psi,
                const arma::vec& h1,
                const double& lambda,
                const double& sigma,
                const int& n){
  
  arma::vec tmp   = lambda*Gyb + psi;
  arma::vec out(n);
  for (int i(0); i < n; ++ i) {
    out(i) = R::pnorm(tmp(i) - h1(y(i)), 0, sigma, true, true) + 
      log(1 - exp(R::pnorm(tmp(i) - h1(y(i) + 1), 0, sigma, true, true) -
      R::pnorm(tmp(i) - h1(y(i)), 0, sigma, true, true)));
  }
  return out;
}


// foptimREM: compute the likelihood given the patameters
//[[Rcpp::export]]
double foptimREM(arma::vec& yb,
                 arma::vec& Gyb,
                 const arma::vec& theta,
                 const arma::mat& X,
                 List& G,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const arma::vec& h1,
                 const int& K,
                 const int& n,
                 const arma::vec& y,
                 const double& tol = 1e-13,
                 const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  Rcpp::print(thetacpp);
  arma::vec xb    = X * theta.subvec(1, K);
  double lambda   = 1.0/(exp(-theta(0)) + 1);
  double sigma    = exp(theta(K + 1));

  // compute ybar
  fyb(yb, Gyb, G, igroup, ngroup, xb, lambda, sigma, n, tol, maxit);
  
  
  arma::vec logp  = flogp(y, yb, Gyb, xb, h1, lambda, sigma, n);
  double llh      = sum(logp);
  
  if(llh < -1e308) {
    llh           = -1e308;
  }
  return -llh;
}

//// NPL solution

//[[Rcpp::export]]
double foptimREM_NPL(arma::vec& yb,
                     arma::vec& Gyb,
                     const arma::vec& theta,
                     const arma::mat& X,
                     List& G,
                     const arma::mat& igroup,
                     const int& ngroup,
                     const arma::vec& h1,
                     const int& K,
                     const int& n,
                     const arma::vec& y) {
  arma::vec xb    = X * theta.subvec(1, K);
  double lambda   = 1.0/(exp(-theta(0)) + 1);
  double sigma    = exp(theta(K + 1)) + 1e-8;
  
  arma::vec logp  = flogp(y, yb, Gyb, xb, h1, lambda, sigma, n);
  double llh      = sum(logp);
  
  if(llh < -1e293) {
    llh           = -1e293;
  }

  return -llh;
}


//[[Rcpp::export]]
void fL_NPL(arma::vec& u,
            arma::vec& Gu,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::mat& X,
            const arma::vec& theta,
            const int& K,
            const int& n) {
  int nm;
  arma::mat Gm;
  arma::vec xbm, um;
  
  arma::vec xb    = X * theta.subvec(1, K);
  double lambda   = 1.0/(exp(-theta(0)) + 1);
  double sigma    = exp(theta(K + 1));
  
  //loop over group
  for (int m(0); m < ngroup; ++ m) {
    nm                                   = igroup(m,1) - igroup(m,0) + 1;
    um                                   = u.subvec(igroup(m,0), igroup(m,1));
    arma::mat Gm                         = G[m];
    xbm                                  = xb.subvec(igroup(m,0), igroup(m,1));
    
    um = fL(um, lambda, sigma, xbm, Gm, nm);
    
    u.rows(igroup(m,0), igroup(m,1))    = um;
    Gu.rows(igroup(m,0), igroup(m,1))   = Gm * um;
  }
}
