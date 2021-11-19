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

// complet information
//[[Rcpp::export]]
void fySar(arma::vec& y,
           arma::vec& Gy,
           List& G,
           const arma::vec& eps,
           const arma::mat& igroup,
           const int& ngroup,
           const arma::vec& psi,
           const double& lambda) {
  int nm;
  arma::vec Gym, ym, tmp0, tmp1, t(ngroup);
  arma::mat Am;
  //loop over group
  for (int m(0); m < ngroup; ++ m) {
    nm            = igroup(m,1) - igroup(m,0) + 1;
    arma::mat Gm  = G[m];
    Am            = arma::diagmat(arma::ones(nm)) - lambda*Gm;
    ym            = arma::solve(Am, psi.subvec(igroup(m,0), igroup(m,1)) + eps.subvec(igroup(m,0), igroup(m,1)));
    
    y.rows(igroup(m,0), igroup(m,1))    = ym;
    Gy.rows(igroup(m,0), igroup(m,1))   = Gm * ym;
  }
}




// [[Rcpp::export]]
double foptimSAR(const double& alphatilde,
                 const arma::mat& X,
                 const arma::mat& invXX,
                 List& G,
                 List& I,
                 const int& n,
                 const arma::vec& y,
                 const arma::vec& Gy,
                 const int& ngroup){
  double lambda         = 1.0/(1.0 + exp(-alphatilde));
  Rcpp::Rcout<<"---------------"<< endl;
  Rcpp::Rcout<<"Estimate:   "<< lambda << endl;
  arma::vec yrem       = y - lambda*Gy;
  arma::vec beta       = invXX*(X.t()*yrem);
  arma::vec e          = yrem - X*beta;
  double s2            = sum(e%e)/n;
  
  double logdetA(0);
  for (int i(0); i < ngroup; ++ i) {
    double vali, signi;
    arma::mat Gi  = G[i];
    arma::mat Ii  = I[i];
    log_det(vali, signi, Ii - lambda*Gi);
    logdetA      += vali;
    logdetA      += log(signi);
  }
  
  double llh      = - 0.5*n*(log(2*acos(-1)*s2)) + logdetA - 0.5*n;
  Rcpp::Rcout <<"Likelihood: "<< llh << endl;
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}


// [[Rcpp::export]]
double foptimSAR0(const double& alphatilde,
                  const arma::mat& X,
                  const arma::mat& invXX,
                  List& G,
                  List& I,
                  const int& n,
                  const arma::vec& y,
                  const arma::vec& Gy,
                  const int& ngroup){
  double lambda         = 1.0/(1.0 + exp(-alphatilde));
  arma::vec yrem       = y - lambda*Gy;
  arma::vec beta       = invXX*(X.t()*yrem);
  arma::vec e          = yrem - X*beta;
  double s2            = sum(e%e)/n;
  
  double logdetA(0);
  for (int i(0); i < ngroup; ++ i) {
    double vali, signi;
    arma::mat Gi  = G[i];
    arma::mat Ii  = I[i];
    log_det(vali, signi, Ii - lambda*Gi);
    logdetA      += vali;
    logdetA      += log(signi);
  }
  
  double llh      = - 0.5*n*(log(2*acos(-1)*s2)) + logdetA - 0.5*n;
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}



// [[Rcpp::export]]
arma::mat fSARjac(const double& lambda,
                  const double& s2,
                  const arma::mat& X,
                  const arma::mat& XX,
                  const arma::vec& Xbeta,
                  List& G,
                  List& I,
                  const arma::mat igroup,
                  const int& ngroup,
                  const int& n,
                  const int& K){
  
  arma::vec GXbeta(n);
  double trGsG(0);
  double trG(0);
  
  for(int m(0); m < ngroup; ++m) {
    arma::mat Wm    = G[m];
    arma::mat tWm   = Wm.t();
    arma::mat Im    = I[m];
    arma::mat tGm   = arma::solve(Im - lambda*tWm, tWm);
    arma::mat Gm    = tGm.t();
    trGsG          += arma::trace((Gm + tGm)*Gm);
    trG            += arma::trace(Gm);
    GXbeta.subvec(igroup(m,0), igroup(m,1)) = Gm*Xbeta.subvec(igroup(m,0), igroup(m,1));
  }
  
  arma::mat out(K + 2, K + 2, arma::fill::zeros);
  out(0, 0) = sum(GXbeta%GXbeta)/s2 + trGsG;
  out.submat(0, 1, 0, K)  = GXbeta.t()*X/s2;
  out(0, K + 1) = trG/s2;
  out.col(0)  = arma::trans(out.row(0));
  out.submat(1, 1, K, K) = XX/s2;
  out(K + 1, K + 1) = n/(2*pow(s2, 2));
  
  return arma::inv(out);
}

// imcomplet information
// fyb: Takes an initial value of yb and finds the equilibrium
// Gyb is G*yb
//[[Rcpp::export]]
void fybsar(arma::vec& yb,
            arma::vec& Gyb,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::vec& psi,
            const double& lambda) {
  int n1, n2, nm;
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1); 
    nm                 = n2 - n1 +1;
    arma::mat Gm       = G[m];
    yb.subvec(n1, n2)  = solve(arma::eye(nm, nm) - lambda*Gm, psi.subvec(n1, n2));
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
}

// [[Rcpp::export]]
double foptimSAR_RE(const double& alphatilde,
                    const arma::mat& X,
                    List& G,
                    List& I,
                    const arma::vec& y,
                    const arma::vec& Gy,
                    const arma::mat igroup,
                    const int& ngroup,
                    const int& n,
                    const int& K){
  int n1, n2;
  double lambda    = 1.0/(1.0 + exp(-alphatilde));
  Rcpp::Rcout<<"---------------"<< endl;
  Rcpp::Rcout<<"Estimate:   "<< lambda << endl;
  arma::mat M(n, K);
  
  for(int m(0); m < ngroup; ++ m){
    n1             = igroup(m,0);
    n2             = igroup(m,1); 
    arma::mat Gm   = G[m];
    arma::mat Im   = I[m];
    M.rows(n1, n2) = Gm*arma::solve(Im - lambda*Gm, X.rows(n1, n2));
  }
  M                = lambda*M + X;
  arma::mat MpM    = M.t()*M;
  arma::vec Mpy    = M.t()*y;
  arma::vec beta   = arma::solve(MpM, Mpy);
  arma::vec e      = y - M*beta;
  double s2        = sum(e%e)/n;
  double llh       = -0.5*n*log(2*acos(-1)*s2) - 0.5*n;
  Rcpp::Rcout <<"Likelihood: "<< llh << endl;
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}

// [[Rcpp::export]]
double foptimSAR0_RE(const double& alphatilde,
                    const arma::mat& X,
                    List& G,
                    List& I,
                    const arma::vec& y,
                    const arma::vec& Gy,
                    const arma::mat igroup,
                    const int& ngroup,
                    const int& n,
                    const int& K){
  int n1, n2;
  double lambda    = 1.0/(1.0 + exp(-alphatilde));
  arma::mat M(n, K);
  
  for(int m(0); m < ngroup; ++ m){
    n1             = igroup(m,0);
    n2             = igroup(m,1); 
    arma::mat Gm   = G[m];
    arma::mat Im   = I[m];
    M.rows(n1, n2) = arma::solve(Im - lambda*Gm, X.rows(n1, n2));
  }
  arma::mat MpM    = M.t()*M;
  arma::vec Mpy    = M.t()*y;
  arma::vec beta   = arma::solve(MpM, Mpy);
  arma::vec e      = y - M*beta;
  double s2        = sum(e%e)/n;
  double llh       = -0.5*n*log(2*acos(-1)*s2) - 0.5*n;
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}
