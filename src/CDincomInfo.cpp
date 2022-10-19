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

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]

#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <RcppEigen.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// fgamma select gamma
double fgamma(const arma::vec& delta,
              const int& r,
              const int& Rbar){
  if (r == 1)    return 0;
  if (r <= Rbar) return delta(r - 2);
  return delta(Rbar - 1);
}

// compute the log of sum of exponential of each row
arma::vec lsumexp(const arma::mat& x){
  arma::vec xmax = max(x, 1);
  return xmax + log(sum(exp(x.each_col() - xmax), 1));
}

// compute the log of average of exponential of each row
arma::vec laverexp(const arma::mat& x, const int& nc){
  arma::vec xmax = max(x, 1);
  return lsumexp(x) - log(nc);
}

/* LOCAL NOTATIONS
 *  the notation m in front of objects refers to the group m
 *  psi     : is X*beta + cte 
 *  ZtLamba : psi + lambda gi*yb, where psi may be a vector of a matrix (many simulations)
 *  lambda  : is the peer effect
 *  delta   : the vector of deltas as delta2, delta3 ..., deltaRbar.
 */

// computation of fL as the mapping L as in the paper.
// step 1: log fL
// lofLF = log(p_1 + ...)
//       = log(p_1) + log(1 + exp(log(p2) - log(p1)) + exp(log(p3) - log(p1)) + ...)
//       = lF1 + log(1 + sum_r=2^inf exp(lFr - lF1))
NumericVector flogL(const NumericVector& ZtLambda,
                    const arma::vec& delta,
                    const int& Rbar,
                    const int& n) {
  int r     = 1;
  double ar = 0;
  bool next = true;
  NumericVector lFr;
  NumericVector lF1 = Rcpp::pnorm5(ZtLambda, 0, 1, true, true);
  NumericVector lu2 = Rcpp::rep(NumericVector::create(1), n);
  while(next) {
    ++ r;
    ar     += fgamma(delta, r, Rbar);
    lFr     = Rcpp::pnorm5(ZtLambda - ar, 0, 1, true, true); 
    lu2    += exp(lFr - lF1); 
    next    = ((max(lFr) > -1000) || (r <= Rbar));
  }
  return lF1 + log(lu2);
}

// step2 fL
//[[Rcpp::export]]
arma::vec fL(const arma::vec& ZtLambda,
             const arma::vec& delta,
             const int& Rbar,
             const int& n) {
  return exp(flogL(wrap(ZtLambda), delta, Rbar, n));
}

// step2' Sum(fL_u) where u is several draws of explanatory variables
//[[Rcpp::export]]
arma::vec fLncond(const arma::mat& ZtLambda,
                  const arma::vec& delta,
                  const int& Rbar,
                  const int& n,
                  const int& nsimu) {
  NumericMatrix logLs(n, nsimu);
  for(int s(0); s < nsimu; ++s){
    logLs(_, s) = flogL(wrap(ZtLambda.col(s)), delta, Rbar, n);
  }
  return exp(laverexp(as<arma::mat>(logLs), nsimu));
}


// fyb: Takes an initial value of yb and finds the equilibrium
// Gyb is G*yb
//[[Rcpp::export]]
int fyb(arma::vec& yb,
        arma::vec& Gyb,
        List& G,
        const arma::mat& igroup,
        const int& ngroup,
        const arma::vec& psi,
        const double& lambda,
        const arma::vec& delta,
        const int& n, 
        const int& Rbar,
        const double& tol,
        const int& maxit) {
  int n1, n2, t = 0;
  arma::vec ZtLambda(n);
  
  computeL: ++t;
  ZtLambda       = lambda*Gyb + psi;
  arma::vec ybst = fL(ZtLambda, delta, Rbar, n);
  // double dist    = max(arma::abs(ybst/(yb + 1e-50) - 1));
  double dist    = max(arma::abs((ybst - yb)/yb));
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

// nonconditional version
//[[Rcpp::export]]
int fybncond(arma::vec& yb,
             arma::vec& Gyb,
             List& G,
             const arma::mat& igroup,
             const int& ngroup,
             const arma::mat& psi,
             const double& lambda,
             const arma::vec& delta,
             const int& n, 
             const int& nsimu,
             const int& Rbar,
             const double& tol,
             const int& maxit) {
  int n1, n2, t = 0;
  arma::mat ZtLambda;
  
  computeL: ++t;
  ZtLambda       = psi.each_col() + lambda*Gyb;
  arma::vec ybst = fLncond(ZtLambda, delta, Rbar, n, nsimu);
  // double dist    = max(arma::abs(ybst/(yb + 1e-50) - 1));
  double dist    = max(arma::abs((ybst - yb)/yb));
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


// fy returns the vector of y given yst and delta
//[[Rcpp::export]]
arma::vec fy(const arma::vec& yst,
             const double& maxyst,
             const arma::vec& delta,
             const int& n,
             const int& Rbar) {
  bool cont = true;
  arma::vec y(n, arma::fill::zeros);
  int r     = 1;
  double ar = 0;
  while(cont) {
    y.elem(arma::find(yst > ar)) += 1;
    ++r;
    ar     += fgamma(delta, r, Rbar);
    cont    = (maxyst > ar);
  }
  return y;
}

// fmeffect computes the marginal effects
//[[Rcpp::export]]
List fmeffects(const int& n,
               const arma::vec& delta,
               const int& Rbar,
               const NumericVector& ZtLambda,
               const arma::vec& lbeta){
  
  arma::mat lphi(n, 1, arma::fill::zeros);
  
  int Rmax  = 0;
  bool next = true;
  double a  = 0;
  
  // log of f 
  while(next) {
    ++ Rmax;
    a                += fgamma(delta, Rmax, Rbar);
    NumericVector lfr = Rcpp::dnorm4(ZtLambda - a, 0, 1, true);
    lphi    = arma::join_rows(lphi, as<arma::vec>(lfr));
    next    = ((max(lfr) > -1000) || (Rmax <= (Rbar + 2)));
  }
  
  // compute  the marginal effects
  arma::vec meffects = mean(exp(lsumexp(lphi.cols(1, Rmax - 1))))*lbeta;
  return List::create(Named("Rmax") = Rmax, Named("meffects") = meffects);
}

// flogP returns the log-likelihood of each individual
// maxy is max(y)
arma::vec flogp(const arma::uvec& y,
                const arma::vec& ZtLambda,
                const int& maxy,
                const double& lambda,
                const arma::vec& delta,
                const int& Rbar,
                const int& n){
  // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
  // Rbar <= maxy + 1
  arma::vec a(maxy + 2);
  a(0)   = R_NegInf; a(1) = 0;
  for(int r(2); r <= Rbar; ++r) {
    a(r) = a(r - 1) + delta(r - 2);
  }
  for(int r(Rbar + 1); r < (maxy + 2); ++r) {
    a(r) = a(r - 1) + delta(Rbar - 1);
  }
  
  // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
  NumericVector tmp1 = wrap(ZtLambda - a.elem(y));
  NumericVector tmp2 = wrap(ZtLambda - a.elem(y + 1));
  
  // log likelihood
  NumericVector lFleft  = Rcpp::pnorm5(tmp1, 0, 1, true, true);
  NumericVector lFright = Rcpp::pnorm5(tmp2, 0, 1, true, true);
  return lFleft + log(1 - exp(lFright - lFleft));
}

//non conditional version
arma::vec flogpncond(const arma::uvec& y,
                     const arma::mat& ZtLambda,
                     const int& maxy,
                     const double& lambda,
                     const arma::vec& delta,
                     const int& nsimu,
                     const int& Rbar,
                     const int& n){
  // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
  // Rbar <= maxy + 1
  arma::vec a(maxy + 2);
  a(0)   = R_NegInf; a(1) = 0;
  for(int r(2); r <= Rbar; ++r) {
    a(r) = a(r - 1) + delta(r - 2);
  }
  for(int r(Rbar + 1); r < (maxy + 2); ++r) {
    a(r) = a(r - 1) + delta(Rbar - 1);
  }
  
  // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
  NumericMatrix tmp1 = wrap(ZtLambda.each_col() - a.elem(y));
  NumericMatrix tmp2 = wrap(ZtLambda.each_col() - a.elem(y + 1));
  
  // log likelihood
  NumericVector lFlefts, lFrights;
  NumericMatrix logP(n, nsimu);
  for(int s(0); s < nsimu; ++s) {
    lFlefts   = Rcpp::pnorm5(tmp1(_,s), 0, 1, true, true);
    lFrights  = Rcpp::pnorm5(tmp2(_,s), 0, 1, true, true);
    logP(_,s) = lFlefts + log(1 - exp(lFrights - lFlefts));
  }
  
  return laverexp(as<arma::mat>(logP), nsimu);
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
                 const int& K,
                 const int& n,
                 const int& Rbar,
                 const arma::uvec& y,
                 const int& maxy,
                 const double& tol = 1e-13,
                 const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim")   = R_NilValue;
  Rcpp::print(thetacpp);
  
  arma::vec psi   = X * theta.subvec(1, K);
  double lambda   = 1.0/(exp(-theta(0)) + 1);
  arma::vec delta = exp(theta.tail(Rbar)) + lambda;  /*bar delta is in the vector delta. so theta as delta2 to deltaRbar and bardelta*/
  
  // compute ybar
  fyb(yb, Gyb, G, igroup, ngroup, psi, lambda, delta, n, Rbar, tol, maxit);
  
  arma::vec ZtLambda = lambda*Gyb + psi;
  arma::vec logp     = flogp(y, ZtLambda, maxy, lambda, delta, Rbar, n);
  double llh         = sum(logp);
  
  // if(llh < -1e308) {
  //   llh           = -1e308;
  // }
  return -llh;
}

//non conditional version
// one simulation
//[[Rcpp::export]]
double foptimREMncond1(arma::vec& yb,
                       arma::vec& Gyb,
                       const arma::vec& theta,
                       const arma::mat& X,
                       const arma::mat& Simu1,
                       const int& nsimu,
                       List& G,
                       const arma::mat& igroup,
                       const int& ngroup,
                       const int& K,
                       const int& n,
                       const int& Rbar,
                       const arma::uvec& y,
                       const int& maxy,
                       const double& tol = 1e-13,
                       const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim")   = R_NilValue;
  Rcpp::print(thetacpp);
  double lambda   = 1.0/(exp(-theta(0)) + 1);
  arma::mat psi   = Simu1*theta(K+1);
  arma::vec tmp   = X*theta.subvec(1, K);
  psi.each_col() += tmp;
  arma::vec delta = exp(theta.tail(Rbar)) + lambda;
  
  // compute ybar
  fybncond(yb, Gyb, G, igroup, ngroup, psi, lambda, delta, n, nsimu, Rbar, tol, maxit);
  
  arma::vec ZtLambda = psi.each_col() + lambda*Gyb;
  arma::vec logp     = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
  double llh         = sum(logp);
  
  // if(llh < -1e308) {
  //   llh              = -1e308;
  // }
  return -llh;
}


//// NPL solution
//[[Rcpp::export]]
double foptimREM_NPL(const arma::vec& Gyb,
                     const arma::vec& theta,
                     const arma::mat& X,
                     const int& Rbar,
                     const int& maxy,
                     const int& K,
                     const int& n,
                     const arma::uvec& y) {
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  arma::vec ZtLambda = lambda*Gyb + X*theta.subvec(1, K);
  arma::vec delta    = exp(theta.tail(Rbar)) + lambda;
  arma::vec logp     = flogp(y, ZtLambda, maxy, lambda, delta, Rbar, n);
  double llh         = sum(logp);
  //cout<<llh<<endl;
  // if(llh < -1e293) {
  //   llh              = -1e293;
  // }
  
  return -llh;
}

//[[Rcpp::export]]
void fL_NPL(arma::vec& yb,
            arma::vec& Gyb,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::mat& X,
            const arma::vec& theta,
            const int& Rbar,
            const int& K,
            const int& n) {
  int n1, n2;
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  arma::vec ZtLambda = lambda*Gyb + X*theta.subvec(1, K);
  arma::vec delta    = exp(theta.tail(Rbar)) + lambda;
  
  // new yb
  yb.subvec(0, n-1)  = fL(ZtLambda, delta, Rbar, n);
  
  // new Gyb
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1);
    arma::mat Gm       = G[m];
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
}

//[[Rcpp::export]]
void fnewyb(arma::vec& yb,
            arma::vec& Gyb,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::mat& X,
            const arma::vec& theta,
            const int& Rbar,
            const int& K,
            const int& n,
            const double& tol,
            const int& maxit) {
  arma::vec psi = X*theta.subvec(1, K);
  double lambda      = 1.0/(exp(-theta(0)) + 1);
  arma::vec delta    = exp(theta.tail(Rbar)) + lambda;
  fyb(yb, Gyb, G, igroup, ngroup, psi, lambda, delta, n,Rbar, tol, maxit);
}

// non conditional version
// using one factor
//[[Rcpp::export]]
double foptimREM_NPLncond1(const arma::vec& Gyb,
                           const arma::vec& theta,
                           const arma::mat& X,
                           const arma::mat& Simu1,
                           const int& nsimu,
                           const int& Rbar,
                           const int& maxy,
                           const int& K,
                           const int& n,
                           const arma::uvec& y) {
  double lambda        = 1.0/(exp(-theta(0)) + 1);
  arma::mat ZtLambda   = Simu1*theta(K+1);
  arma::vec tmp        = lambda*Gyb + X*theta.subvec(1, K);
  ZtLambda.each_col() += tmp;
  arma::vec delta      = exp(theta.tail(Rbar)) + lambda;
  arma::vec logp       = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
  double llh           = sum(logp);
  // if(llh < -1e293) {
  //   llh                = -1e293;
  // }
  return -llh;
}

//[[Rcpp::export]]
void fL_NPLncond1(arma::vec& yb,
                  arma::vec& Gyb,
                  List& G,
                  const arma::mat& igroup,
                  const int& ngroup,
                  const arma::mat& X,
                  const arma::vec& theta,
                  const arma::mat& Simu1,
                  const int& nsimu,
                  const int& Rbar,
                  const int& K,
                  const int& n) {
  int n1, n2;
  double lambda        = 1.0/(exp(-theta(0)) + 1);
  arma::mat ZtLambda   = Simu1*theta(K+1);
  arma::vec tmp        = lambda*Gyb + X*theta.subvec(1, K);
  ZtLambda.each_col() += tmp;
  arma::vec delta      = exp(theta.tail(Rbar)) + lambda;
  
  // new yb
  yb.subvec(0, n-1)    = fLncond(ZtLambda, delta, Rbar, n, nsimu);
  
  // new Gyb
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1);
    arma::mat Gm       = G[m];
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
}


// using two factors
//[[Rcpp::export]]
double foptimREM_NPLncond2(const arma::vec& Gyb,
                           const arma::vec& theta,
                           const arma::mat& X,
                           const arma::mat& Simu1,
                           const arma::mat& Simu2,
                           const int& nsimu,
                           const int& Rbar,
                           const int& maxy,
                           const int& K,
                           const int& n,
                           const arma::uvec& y) {
  double lambda        = 1.0/(exp(-theta(0)) + 1);
  arma::mat ZtLambda   = Simu1*theta(K+1) + Simu2*theta(K+2);
  arma::vec tmp        = lambda*Gyb + X*theta.subvec(1, K);
  ZtLambda.each_col() += tmp;
  arma::vec delta      = exp(theta.tail(Rbar)) + lambda;
  arma::vec logp       = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
  double llh           = sum(logp);
  // if(llh < -1e293) {
  //   llh                = -1e293;
  // }
  return -llh;
}

//[[Rcpp::export]]
void fL_NPLncond2(arma::vec& yb,
                  arma::vec& Gyb,
                  List& G,
                  const arma::mat& igroup,
                  const int& ngroup,
                  const arma::mat& X,
                  const arma::vec& theta,
                  const arma::mat& Simu1,
                  const arma::mat& Simu2,
                  const int& nsimu,
                  const int& Rbar,
                  const int& K,
                  const int& n) {
  int n1, n2;
  double lambda        = 1.0/(exp(-theta(0)) + 1);
  arma::mat ZtLambda   = Simu1*theta(K+1) + Simu2*theta(K+2);
  arma::vec tmp        = lambda*Gyb + X*theta.subvec(1, K);
  ZtLambda.each_col() += tmp;
  arma::vec delta      = exp(theta.tail(Rbar)) + lambda;
  
  // new yb
  yb.subvec(0, n-1)    = fLncond(ZtLambda, delta, Rbar, n, nsimu);
  
  // new Gyb
  for (int m(0); m < ngroup; ++ m) {
    n1                 = igroup(m,0);
    n2                 = igroup(m,1);
    arma::mat Gm       = G[m];
    Gyb.subvec(n1, n2) = Gm*yb.subvec(n1, n2);
  }
}

//This function computes log of integral of phi(x) from a to b using IS
//[[Rcpp::export]]
arma::vec flogintphi(const int& n,
                     const int& S,
                     const double& a,
                     const double& b,
                     const arma::vec& Mean,
                     const arma::rowvec& simu,
                     const arma::mat& igroup,
                     const int& ngroup){
  arma::vec out(n);
  for(int m(0); m < ngroup; ++m){
    int n1              = igroup(m,0);
    int n2              = igroup(m,1);
    arma::vec Meanm     = Mean.subvec(n1, n2);
    arma::mat tmp       = arma::repmat(Meanm, 1, S);
    arma::mat logphis   = -0.5*log(2*acos(-1)) - 0.5*pow(tmp.each_row() - (simu*(b - a) + a), 2);
    out.subvec(n1, n2)  = laverexp(logphis, S);
  }
  return out + log(b - a);
}


// Numerical optimization using Rcpp

class cdnetreg: public MFuncGrad
{
private:
  const arma::mat& Z;
  const arma::mat& X;
  const int& Rbar;
  const int& maxy;
  const int& K;
  const int& n;
  const arma::uvec& y;
  const double& l2ps2;
  List& lidy;
public:
  cdnetreg(const arma::mat& Z_,
                 const arma::mat& X_,
                 const int& Rbar_,
                 const int& maxy_,
                 const int& K_,
                 const int& n_,
                 const arma::uvec& y_,
                 const double& l2ps2_,
                 List& lidy_) : 
  Z(Z_),
  X(X_),
  Rbar(Rbar_),
  maxy(maxy_),
  K(K_),
  n(n_),
  y(y_),
  l2ps2(l2ps2_),
  lidy(lidy_){}
  
  Eigen::VectorXd Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec beta         = arma::vec(theta0.data(), K + Rbar + 1, false, false); //converte into arma vec
    
    beta(0)                 = 1.0/(exp(-beta(0)) + 1);
    double lambda           = beta(0);
    arma::vec deltat        = exp(beta.tail(Rbar));
    beta.tail(Rbar)         = deltat + lambda;
    arma::vec delta         = beta.tail(Rbar);
    
    // print
    // NumericVector betacpp   = wrap(beta);
    // betacpp.attr("dim")     = R_NilValue;
    // std::printf("beta: \n");
    // Rcpp::print(betacpp);
    
    // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
    // Rbar <= maxy + 1
    arma::vec a(maxy + 2);
    a(0)   = R_NegInf; a(1) = 0;
    for(int r(2); r <= Rbar; ++r) {
      a(r) = a(r - 1) + delta(r - 2);
    }
    for(int r(Rbar + 1); r <= (maxy + 1); ++r) {
      a(r) = a(r - 1) + delta(Rbar - 1);
    }
    
    // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
    arma::vec ZtLambda      = Z*beta.head(K + 1);
    arma::vec Zba1          = ZtLambda - a.elem(y);
    arma::vec Zba2          = ZtLambda - a.elem(y + 1);
    NumericVector Zba1r     = wrap(Zba1);
    NumericVector Zba2r     = wrap(Zba2);
    
    arma::vec lfZba1        = -0.5*Zba1%Zba1 - l2ps2;
    arma::vec lfZba2        = -0.5*Zba2%Zba2 - l2ps2;
    
    NumericVector lFZba1    = Rcpp::pnorm5(Zba1r, 0, 1, true, true);
    NumericVector lFZba2    = Rcpp::pnorm5(Zba2r, 0, 1, true, true);
    
    // log likelihood
    NumericVector ldF       = lFZba1 + log(1 - exp(lFZba2 - lFZba1));
    double f                = sum(ldF);
    
    if(f > 1e293) {
      f                     = 1e293;
    }
    
    // gradient
    arma::vec ldFa         = as<arma::vec>(ldF);
    arma::vec tmp1         = exp(lfZba1 - ldFa);
    arma::vec tmp2         = exp(lfZba2 - ldFa);
    arma::vec gdarm(K + Rbar + 1);
    gdarm.head(K + 1)      = arma::trans(arma::sum(Z.each_col()%(tmp1 - tmp2), 0));
    for(int r(2); r <= Rbar; ++ r){
      arma::uvec id2       = lidy[r - 2]; // greater than r - 1
      arma::uvec id1       = lidy[r - 1]; // greater than r
      double tmp3          = sum(tmp2.elem(id2)) - sum(tmp1.elem(id1));
      gdarm(0)            += tmp3;
      gdarm(K + r - 1)     = tmp3;
    }
    arma::uvec id2         = lidy[Rbar - 1]; // greater than Rbar
    arma::uvec id1         = lidy[Rbar]; // greater than Rbar + 1
    double tmp3            = sum(tmp2.elem(id2)%(y.elem(id2) + 1 - Rbar)) - sum(tmp1.elem(id1)%(y.elem(id1) - Rbar));
    gdarm(0)              += tmp3;
    gdarm(K + Rbar)        = tmp3;
    
    
    gdarm.tail(Rbar)      %= deltat;
    gdarm(0)              *= lambda*(1 - lambda);
    grad                   = -Eigen::Map<Eigen::VectorXd>(gdarm.memptr(), K + Rbar + 1);
    Grad                   = -grad;
    
    // cout<< f <<endl;
    // std::printf("log-likelihood: %f\n", f);
    return -f;
  }
};


class cdnetreg_print: public MFuncGrad
{
private:
  const arma::mat& Z;
  const arma::mat& X;
  const int& Rbar;
  const int& maxy;
  const int& K;
  const int& n;
  const arma::uvec& y;
  const double& l2ps2;
  List& lidy;
public:
  cdnetreg_print(const arma::mat& Z_,
                 const arma::mat& X_,
                 const int& Rbar_,
                 const int& maxy_,
                 const int& K_,
                 const int& n_,
                 const arma::uvec& y_,
                 const double& l2ps2_,
                 List& lidy_) : 
  Z(Z_),
  X(X_),
  Rbar(Rbar_),
  maxy(maxy_),
  K(K_),
  n(n_),
  y(y_),
  l2ps2(l2ps2_),
  lidy(lidy_){}
  
  Eigen::VectorXd Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec beta         = arma::vec(theta0.data(), K + Rbar + 1, false, false); //converte into arma vec
    
    beta(0)                 = 1.0/(exp(-beta(0)) + 1);
    double lambda           = beta(0);
    arma::vec deltat        = exp(beta.tail(Rbar));
    beta.tail(Rbar)         = deltat + lambda;
    arma::vec delta         = beta.tail(Rbar);
    
    // print
    NumericVector betacpp   = wrap(beta);
    betacpp.attr("dim")     = R_NilValue;
    Rcpp::Rcout << "beta: \n";
    Rcpp::print(betacpp);
    
    // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
    // Rbar <= maxy + 1
    arma::vec a(maxy + 2);
    a(0)   = R_NegInf; a(1) = 0;
    for(int r(2); r <= Rbar; ++r) {
      a(r) = a(r - 1) + delta(r - 2);
    }
    for(int r(Rbar + 1); r <= (maxy + 1); ++r) {
      a(r) = a(r - 1) + delta(Rbar - 1);
    }
    
    // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
    arma::vec ZtLambda      = Z*beta.head(K + 1);
    arma::vec Zba1          = ZtLambda - a.elem(y);
    arma::vec Zba2          = ZtLambda - a.elem(y + 1);
    NumericVector Zba1r     = wrap(Zba1);
    NumericVector Zba2r     = wrap(Zba2);
    
    arma::vec lfZba1        = -0.5*Zba1%Zba1 - l2ps2;
    arma::vec lfZba2        = -0.5*Zba2%Zba2 - l2ps2;
    
    NumericVector lFZba1    = Rcpp::pnorm5(Zba1r, 0, 1, true, true);
    NumericVector lFZba2    = Rcpp::pnorm5(Zba2r, 0, 1, true, true);
    
    // log likelihood
    NumericVector ldF       = lFZba1 + log(1 - exp(lFZba2 - lFZba1));
    double f                = sum(ldF);
    
    if(f > 1e293) {
      f                     = 1e293;
    }
    
    // gradient
    arma::vec ldFa         = as<arma::vec>(ldF);
    arma::vec tmp1         = exp(lfZba1 - ldFa);
    arma::vec tmp2         = exp(lfZba2 - ldFa);
    arma::vec gdarm(K + Rbar + 1);
    gdarm.head(K + 1)      = arma::trans(arma::sum(Z.each_col()%(tmp1 - tmp2), 0));
    for(int r(2); r <= Rbar; ++ r){
      arma::uvec id2       = lidy[r - 2]; // greater than r - 1
      arma::uvec id1       = lidy[r - 1]; // greater than r
      double tmp3          = sum(tmp2.elem(id2)) - sum(tmp1.elem(id1));
      gdarm(0)            += tmp3;
      gdarm(K + r - 1)     = tmp3;
    }
    arma::uvec id2         = lidy[Rbar - 1]; // greater than Rbar
    arma::uvec id1         = lidy[Rbar]; // greater than Rbar + 1
    double tmp3            = sum(tmp2.elem(id2)%(y.elem(id2) + 1 - Rbar)) - sum(tmp1.elem(id1)%(y.elem(id1) - Rbar));
    gdarm(0)              += tmp3;
    gdarm(K + Rbar)        = tmp3;
    
    
    gdarm.tail(Rbar)      %= deltat;
    gdarm(0)              *= lambda*(1 - lambda);
    grad                   = -Eigen::Map<Eigen::VectorXd>(gdarm.memptr(), K + Rbar + 1);
    Grad                   = -grad;
    
    // cout<< f <<endl;
    Rcpp::Rcout << "log-likelihood: " << f << "\n";
    return -f;
  }
};


//[[Rcpp::export]]
List cdnetLBFGS(Eigen::VectorXd par,
                const arma::vec& Gyb,
                const arma::mat& X,
                const int& Rbar,
                const int& maxy,
                const int& K,
                const int& n,
                const arma::uvec& y,
                const int& maxit = 300, 
                const double& eps_f = 1e-6, 
                const double& eps_g = 1e-5,
                const bool& print = false) {
  double l2ps2     = 0.5*log(2*acos(-1));
  arma::mat Z      = arma::join_rows(Gyb, X);
  List lidy(Rbar + 1);
  for(int r(0); r <= Rbar; ++ r){
    arma::uvec id  = arma::find(y >= (r + 1));   
    lidy(r)        = id;
  }
  
  double fopt;
  int status;
  Eigen::VectorXd grad;
  
  if(print){
    cdnetreg_print f(Z, X, Rbar, maxy, K, n, y, l2ps2, lidy);
    status = optim_lbfgs(f, par, fopt, maxit, eps_f, eps_g);
    grad  = f.Grad;
  } else {
    cdnetreg f(Z, X, Rbar, maxy, K, n, y, l2ps2, lidy);
    status = optim_lbfgs(f, par, fopt, maxit, eps_f, eps_g);
    grad  = f.Grad;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("par")      = par,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("status")   = status);
}



// variance
//[[Rcpp::export]]
List fcovCDI(const int& n,
             const arma::vec& Gyb,
             const arma::vec& theta,
             const arma::mat& X,
             const int& Rbar,
             const int& K,
             const int& S,
             List& G,
             const arma::mat& igroup,
             const int& ngroup,
             const bool& ccov) {
  double lambda          = 1.0/(exp(-theta(0)) + 1);
  arma::vec beta         = theta.subvec(1, K);
  arma::vec lbeta        = arma::join_cols(arma::ones(1)*lambda, beta);
  arma::vec tdelta       = exp(theta.tail(Rbar));
  arma::vec delta        = tdelta + lambda;
  arma::vec ZtLambdavec  = Gyb*lambda + X*beta;       
  NumericVector ZtLambda = wrap(ZtLambdavec);
  arma::vec Gybpo        = Gyb + 1;
  arma::vec Gybposq      = pow(Gybpo, 2);
  arma::vec tGybpo       = 2*Gyb + 1;
  
  List out;
  
  if(ccov){
    arma::rowvec simu(S, arma::fill::randu);
    int Kz                 = K + 1;
    
    arma::mat Z            = arma::join_rows(Gyb, X);
    arma::mat ldPhi;                          // l(Phir - Phir+1) from r = 0
    {
      NumericVector ldPhi0 = Rcpp::pnorm5(ZtLambda, 0, 1, false, true); 
      ldPhi                = as<arma::vec>(ldPhi0);
    }
    arma::mat lphi         = R_NegInf*arma::ones<arma::vec>(n);  // lphi from r = 0
    arma::rowvec lavec(1, arma::fill::zeros); // ar from r = 0
    // arma::mat zcheck       = Gyb;             // Gye - r from r = 0
    int Rmax  = 0;
    bool next = true;
    double a2 = 0, a1;
    
    // log of f and log(diff F1 F2)
    while(next) {
      ++ Rmax;
      a1                = a2;
      a2               += fgamma(delta, Rmax + 1, Rbar);
      lavec             = arma::join_rows(lavec, arma::ones(1)*a1);
      NumericVector lfr = Rcpp::dnorm4(ZtLambda - a1, 0, 1, true);
      lphi    = arma::join_rows(lphi, as<arma::vec>(lfr));
      ldPhi   = arma::join_rows(ldPhi, flogintphi(n, S, a1, a2, ZtLambda, simu, igroup, ngroup));
      // zcheck  = arma::join_rows(zcheck, Gyb - Rmax);
      next    = ((max(lfr) > -1000) || (Rmax <= (Rbar + 2)));
    }
    lavec     = log(lavec);
    
    // log(0, 1, 2, ...)
    arma::rowvec tmp0       = log(arma::regspace<arma::rowvec>(0, Rmax));
    
    // compute GinvSB, marginal effects, and Qth (used for the marginal effect variance)
    arma::mat GinvSB(n, Kz + Rbar), Qth(Kz, Kz + Rbar);
    arma::vec meffects;
    {
      // Compute d, marg.effect, Qtheta, d is diag(D) of the paper
      arma::mat tmp1(n, Rbar + 2), tmp2(n, Rbar + 2);
      for(int l(0); l <= Rbar; ++ l){
        tmp1.col(l)      = exp(lsumexp(lphi.cols(l, Rmax - 1))); //for sum fr from r = l
        tmp2.col(l)      = exp(lsumexp(lphi.cols(l, Rmax - 1).each_row() + lavec.cols(l, Rmax - 1))); //for sum a_r f_r from r = l, 
      }
      tmp1.col(Rbar + 1) = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1))); //for sum dot a_delta,r fr from r = Rbar + 1
      tmp2.col(Rbar + 1) = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1).each_row() + (lavec.cols(Rbar + 1, Rmax - 1) + tmp0.subvec(1, Rmax - Rbar - 1)))); //for sum (r - Rbar) a_r phi_r from r = Rbar + 1,
      arma::vec tmp3     = exp(lsumexp(lphi.cols(1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - 1))); // sum r*phir from r = 1
      arma::vec tmp4     = exp(lsumexp(lphi.cols(1, Rmax - 1).each_row() + (tmp0.subvec(1, Rmax - 1) + lavec.cols(1, Rmax - 1)))); // sum r*ar*phir from r = 1
      
      // compute d
      arma::vec d   = tmp1.col(1); //for sum fr from r = 1
      
      // marg. effects
      double meand  = mean(d);
      meffects      = meand*lbeta;
      
      // Qtheta
      Qth.col(0)            = arma::eye<arma::mat>(Kz, 1)*lambda*mean(d) - 
        lambda*lbeta*mean(Gyb%(ZtLambdavec%d - tmp2.col(1)) - ZtLambdavec%tmp3 + tmp4);
      Qth.cols(1, Kz - 1)   = arma::join_cols(arma::zeros<arma::rowvec>(K), arma::eye<arma::mat>(K, K)*mean(d)) - 
        (lbeta*mean(X.each_col()%(ZtLambdavec%d - tmp2.col(1)), 0));
      
      // Compute B and complete Qth for column delta 
      // I start by B2
      arma::mat B(n, Rbar);
      for(int k(0); k < Rbar; ++k) {
        B.col(k)            = -tdelta(k)*tmp1.col(2 + k); 
        Qth.col(Kz + k)     = tdelta(k)*lbeta*mean(ZtLambdavec%tmp1.col(2 + k) - tmp2.col(2 + k), 0);
      }
      // I now complete B1 and DZ
      B                     = arma::join_rows(lambda*(Gybpo%d - tmp3), arma::join_rows(X.each_col()%d, B));
      
      // G*inv(S)*B
      for (int m(0); m < ngroup; ++ m) {
        int n1              = igroup(m,0);
        int n2              = igroup(m,1);
        int nm              = n2 - n1 + 1;
        arma::mat Gm        = G[m];
        arma::mat Sm        = Gm.each_col() % d.subvec(n1, n2);
        Sm                  = arma::eye<arma::mat>(nm, nm) - lambda*Sm;
        GinvSB.rows(n1, n2) = Gm * arma::solve(Sm, B.rows(n1, n2));
      }
    }
    // Compute Sigma0 and Omega0
    arma::mat Sigma(Kz + Rbar, Kz + Rbar), Omega(Kz + Rbar, n);
    
    {//Need temporary variables
      arma::mat tmpra(n, Rbar + 1), tmprb(n, Rbar + 1), tmprc(n, Rbar + 1);
      arma::mat tmpr1a(n, Rbar + 1), tmpr1b(n, Rbar + 1), tmpr1c(n, Rbar + 1);
      arma::mat tmprra(n, Rbar + 1), tmprrb(n, Rbar + 1), tmprrc(n, Rbar + 1);
      for(int l(0); l <= Rbar; ++ l){
        tmpra.col(l)   = exp(lsumexp(2*lphi.cols(l, Rmax - 1) - ldPhi.cols(l, Rmax - 1))); //for sum fr^2/(Fr - Fr+1) from r = l
        tmprb.col(l)   = exp(lsumexp(lphi.cols(l, Rmax - 1) + (lphi.cols(l, Rmax - 1).each_row() + tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r*fr^2/(Fr - Fr+1) from r = l
        tmprc.col(l)   = exp(lsumexp(2*(lphi.cols(l, Rmax - 1).each_row() + tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r^2*fr^2/(Fr - Fr+1) from r = l
        tmpr1a.col(l)  = exp(lsumexp(2*lphi.cols(1 + l, Rmax) - ldPhi.cols(l, Rmax - 1))); //for sum fr+1^2/(Fr - Fr+1) from r = l
        tmpr1b.col(l)  = exp(lsumexp(lphi.cols(1 + l, Rmax) + (lphi.cols(1 + l, Rmax).each_row() + tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r*fr+1^2/(Fr - Fr+1) from r = l
        tmpr1c.col(l)  = exp(lsumexp(2*(lphi.cols(1 + l, Rmax).each_row() + tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r^2*fr+1^2/(Fr - Fr+1) from r = l
        tmprra.col(l)  = exp(lsumexp(lphi.cols(l, Rmax - 1) + lphi.cols(1 + l, Rmax) - ldPhi.cols(l, Rmax - 1))); //for sum fr * fr+1/(Fr - Fr+1) from r = l
        tmprrb.col(l)  = exp(lsumexp(lphi.cols(l, Rmax - 1) + (lphi.cols(1 + l, Rmax).each_row() + tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r * fr * fr+1/(Fr - Fr+1) from r = l
        tmprrc.col(l)  = exp(lsumexp(lphi.cols(l, Rmax - 1) + (lphi.cols(1 + l, Rmax).each_row() + 2*tmp0.subvec(l, Rmax - 1)) - ldPhi.cols(l, Rmax - 1))); //for sum r^2 * fr * fr+1/(Fr - Fr+1) from r = l
      }
      
      arma::mat adr1tmprr(n, 2), adr1tmpr1(n, 2), adrtmprr(n, 2), adrtmpr(n, 2);
      adr1tmprr.col(0) = exp(lsumexp(lphi.cols(Rbar, Rmax - 2) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar, Rmax - 2))); //for sum dot a_delta,r+1 fr fr+1/(Fr - Fr+1) from r = Rbar
      adr1tmprr.col(1) = exp(lsumexp(lphi.cols(Rbar, Rmax - 2) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + (tmp0.subvec(1, Rmax - Rbar - 1) + tmp0.subvec(Rbar, Rmax - 2))) - ldPhi.cols(Rbar, Rmax - 2))); //for sum dot a_delta,r+1 r fr fr+1/(Fr - Fr+1) from r = Rbar
      adr1tmpr1.col(0) = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar, Rmax - 2))); //for sum dot a_delta,r+1 fr^2/(Fr - Fr+1) from r = Rbar
      adr1tmpr1.col(1) = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + (tmp0.subvec(1, Rmax - Rbar - 1) + tmp0.subvec(Rbar, Rmax - 2))) - ldPhi.cols(Rbar, Rmax - 2))); //for sum dot a_delta,r+1 r fr^2/(Fr - Fr+1) from r = Rbar
      adrtmprr.col(0)  = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 2, Rmax).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum dot a_delta,r fr fr+1/(Fr - Fr+1) from r = Rbar + 1
      adrtmprr.col(1)  = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 2, Rmax).each_row() + (tmp0.subvec(1, Rmax - Rbar - 1) + tmp0.subvec(Rbar + 1, Rmax - 1))) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum dot a_delta,r r fr fr+1/(Fr - Fr+1) from r = Rbar + 1
      adrtmpr.col(0)   = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum dot a_delta,r fr^2/(Fr - Fr+1) from r = Rbar + 1
      adrtmpr.col(1)   = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 1, Rmax - 1).each_row() + (tmp0.subvec(1, Rmax - Rbar - 1) + tmp0.subvec(Rbar + 1, Rmax - 1))) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum dot a_delta,r r fr^2/(Fr - Fr+1) from r = Rbar + 1
      
      arma::vec adr12tmpr1   = exp(lsumexp(2*(lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar, Rmax - 2))); //for sum (dot a_delta,r+1 fr+1)^2/(Fr - Fr+1) from r = Rbar
      arma::vec adrrtmprr    = exp(lsumexp(lphi.cols(Rbar + 1, Rmax - 1) + (lphi.cols(Rbar + 2, Rmax).each_row() + (tmp0.subvec(1, Rmax - Rbar - 1) + tmp0.subvec(2, Rmax - Rbar))) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum (dot a_delta,r dot a_delta,r+1 fr fr+1)/(Fr - Fr+1) from r = Rbar + 1
      arma::vec adr2tmpr     = exp(lsumexp(2*(lphi.cols(Rbar + 1, Rmax - 1).each_row() + tmp0.subvec(1, Rmax - Rbar - 1)) - ldPhi.cols(Rbar + 1, Rmax - 1))); //for sum (dot a_delta,r fr)^2/(Fr - Fr+1) from r = Rbar + 1
      {// A lambda lambda
        Sigma(0, 0)    = pow(lambda, 2)*arma::accu(pow(Gyb, 2)%tmpr1a.col(0) -2*Gyb%tmpr1b.col(0) + tmpr1c.col(0) + 
          Gybposq%tmpra.col(0) - 2*Gybpo%tmprb.col(0) + tmprc.col(0) -
          2*Gyb%Gybpo%tmprra.col(0) + 2*tGybpo%tmprrb.col(0) - 2*tmprrc.col(0));
      }
      
      {// A lambda Gamma
        arma::vec tmp  = lambda*(Gyb%tmpr1a.col(0) - tmpr1b.col(0) - tGybpo%tmprra.col(0) + 2*tmprrb.col(0) + Gybpo%tmpra.col(0) - tmprb.col(0));
        arma::mat Xtmp           = X.each_col()%tmp; 
        arma::rowvec sXtmp       = arma::sum(Xtmp, 0);
        Sigma.submat(0, 1, 0, K) = sXtmp;
        Sigma.submat(1, 0, K, 0) = sXtmp.t();
        Omega.row(0)             = tmp.t();
      }
      
      {// A lambda deltak
        for(int k(2); k <= Rbar; ++ k){
          Sigma(0, Kz + k - 2)  = lambda*tdelta(k - 2)*arma::accu(Gybpo%tmprra.col(k - 1) - tmprrb.col(k - 1) - Gyb%tmpr1a.col(k - 1) + tmpr1b.col(k - 1) +
            Gyb%tmprra.col(k) - tmprrb.col(k) - Gybpo%tmpra.col(k) + tmprb.col(k));
          Sigma(Kz + k - 2, 0)  = Sigma(0, Kz + k - 2);
        }
      }
      
      {// A lambda deltabar
        Sigma(0, Kz + Rbar - 1)  = lambda*tdelta(Rbar - 1)*arma::accu(Gybpo%adr1tmprr.col(0) - adr1tmprr.col(1) - Gyb%adr1tmpr1.col(0) + adr1tmpr1.col(1) +
          Gyb%adrtmprr.col(0) - adrtmprr.col(1) - Gybpo%adrtmpr.col(0) + adrtmpr.col(1));
        Sigma(Kz + Rbar - 1, 0)  = Sigma(0, Kz + Rbar - 1);
      }
      
      {// A Gamma Gamma
        arma::vec tmp            = tmpr1a.col(0) -2*tmprra.col(0) + tmpra.col(0);
        arma::mat Xtmp           = arma::trans(X.each_col() % tmp); 
        Sigma.submat(1, 1, K, K) = Xtmp*X;
        Omega.rows(1, K)         = Xtmp;
      }
      
      {
        for(int k(2); k <= Rbar; ++ k){
          // A Gamma deltak
          arma::vec tmp         = tdelta(k - 2)*(tmprra.col(k - 1) - tmpr1a.col(k - 1) + tmprra.col(k) - tmpra.col(k));
          Sigma.submat(1, Kz + k - 2, K, Kz + k - 2) = X.t()*tmp;
          Sigma.submat(Kz + k - 2, 1, Kz + k - 2, K) = Sigma.submat(1, Kz + k - 2, K, Kz + k - 2).t();
          Omega.row(Kz + k - 2) = tmp.t();
          // A deltak deltak
          Sigma(Kz + k - 2, Kz + k - 2) = pow(tdelta(k - 2), 2)*arma::accu(tmpr1a.col(k - 1) -2*tmprra.col(k) + tmpra.col(k));
          // A deltak deltal l < k
          for(int l(2); l < k; ++ l){
            Sigma(Kz + k - 2, Kz + l - 2) = -tdelta(l - 2)*arma::accu(tmp);
            Sigma(Kz + l - 2, Kz + k - 2) = Sigma(Kz + k - 2, Kz + l - 2);
          }
        }
      }
      
      {// A Gamma deltabar
        arma::vec tmp            = tdelta(Rbar - 1)*(adr1tmprr.col(0) - adr1tmpr1.col(0) + adrtmprr.col(0) - adrtmpr.col(0));
        Sigma.submat(1, Kz + Rbar - 1, K, Kz + Rbar - 1) = X.t()*tmp;
        Sigma.submat(Kz + Rbar - 1, 1, Kz + Rbar - 1, K) = Sigma.submat(1, Kz + Rbar - 1, K, Kz + Rbar - 1).t();
        Omega.row(Kz + Rbar - 1) = tmp.t();
        // A deltak deltabar
        for(int k(2); k <= Rbar; ++ k){
          Sigma(Kz + k - 2, Kz + Rbar - 1) = -tdelta(k - 2)*arma::accu(tmp);
          Sigma(Kz + Rbar - 1, Kz + k - 2) = Sigma(Kz + k - 2, Kz + Rbar - 1);
        }
      }
      
      {// A deltabar deltabar
        Sigma(Kz + Rbar - 1, Kz + Rbar - 1) = pow(tdelta(Rbar - 1), 2)*arma::accu(adr12tmpr1 - 2*adrrtmprr + adr2tmpr);
      }
      
      Sigma         /= n;
      Omega          = lambda*Omega*GinvSB/n;
      arma::mat covt = arma::inv(Sigma + Omega);
      covt           = covt*Sigma*covt.t()/n;
      arma::mat covm = Qth*covt*Qth.t();
      
      //Finding cov of lambda by delta method
      covt.col(0)   *= lambda;
      covt.row(0)   *= lambda;
      
      out            = List::create(Named("Rmax")        = Rmax, 
                                    Named("meffects")    = meffects,
                                    Named("covtheta")    = covt,
                                    Named("covmeffects") = covm,
                                    Named("var.comp")    = List::create(Named("Sigma") = Sigma, Named("Omega") = -Omega));
    }} else{
      out            = fmeffects(n, delta, Rbar, ZtLambda, lbeta);
    }
    return out;
}
