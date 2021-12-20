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

// fgamma select gamma
double fgamma(const arma::vec& delta,
              const int& r,
              const int& Rbar){
  if (r == 1)    return 0;
  if (r <= Rbar) return delta(r - 2);
  return delta(Rbar - 2);
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

// non conditional version
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
  for(int r(2); r < Rbar; ++r) {
    a(r) = a(r - 1) + delta(r - 2);
  }
  for(int r(Rbar); r < (maxy + 2); ++r) {
    a(r) = a(r - 1) + delta(Rbar - 2);
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
  for(int r(2); r < Rbar; ++r) {
    a(r) = a(r - 1) + delta(r - 2);
  }
  for(int r(Rbar); r < (maxy + 2); ++r) {
    a(r) = a(r - 1) + delta(Rbar - 2);
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
  arma::vec delta = exp(theta.tail(Rbar - 1));
  
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
  arma::vec delta = exp(theta.tail(Rbar - 1));
  
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
  arma::vec delta    = exp(theta.tail(Rbar - 1));
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
  arma::vec delta    = exp(theta.tail(Rbar - 1));
  
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
  arma::vec delta    = exp(theta.tail(Rbar - 1));
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
  arma::vec delta      = exp(theta.tail(Rbar - 1));
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
  arma::vec delta      = exp(theta.tail(Rbar - 1));
  
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
  arma::vec delta      = exp(theta.tail(Rbar - 1));
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
  arma::vec delta      = exp(theta.tail(Rbar - 1));
  
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
    int nm              = n2 - n1 + 1;
    arma::vec Meanm     = Mean.subvec(n1, n2);
    arma::mat logphis   = -0.5*log(2*acos(-1)) - 0.5*pow(arma::repmat(Meanm, 1, S) - arma::repmat(simu*(b - a) + a, nm, 1), 2);
    out.subvec(n1, n2)  = laverexp(logphis, S);
  }
  return out + log(b - a);
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
  arma::vec lbeta        = arma::join_cols(arma::ones(1)*lambda, theta.subvec(1, K));
  arma::vec delta        = exp(theta.tail(Rbar - 1));
  NumericVector ZtLambda = wrap(Gyb*lambda + X*theta.subvec(1, K));
  List out;
  
  if(ccov){
    arma::rowvec simu(S, arma::fill::randu);
    int Kz                 = K + 1;
    
    arma::mat Z            = arma::join_rows(Gyb, X);
    arma::mat ldPhi;
    arma::mat lphi(n, 1, arma::fill::zeros);
    arma::rowvec lavec(1, arma::fill::zeros);
    {
      NumericVector ldPhi0 = Rcpp::pnorm5(ZtLambda, 0, 1, false, true); 
      ldPhi                = as<arma::vec>(ldPhi0);
    }
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
      next    = ((max(lfr) > -1000) || (Rmax <= Rbar));
    }
    lavec     = log(lavec);
    
    // log(1, 2, ...)
    arma::rowvec tmp0       = log(arma::regspace<arma::rowvec>(1, Rmax - Rbar));
    
    // compute GinvSW, marginal effects, and Qtheta
    arma::mat GinvSW(n, Kz + Rbar - 1), Qth(Kz, Kz + Rbar - 1);
    arma::vec meffects;
    {
      // Compute d, marg.effect, Qtheta
      arma::mat tmp1(n, Rbar), tmp2(n, Rbar);
      for(int l(0); l < (Rbar - 1); ++ l){
        tmp1.col(l)      = exp(lsumexp(lphi.cols(1 + l, Rmax - 1))); //for sum fr
        tmp2.col(l)      = exp(lsumexp(lphi.cols(1 + l, Rmax - 1).each_row() + lavec.cols(1 + l, Rmax - 1))); //for sum a_r phi_r, 
      }
      tmp1.col(Rbar - 1) = exp(lsumexp(lphi.cols(Rbar, Rmax - 1).each_row() + tmp0)); //for sum (r - Rbar + 1) fr
      tmp2.col(Rbar - 1) = exp(lsumexp(lphi.cols(Rbar, Rmax - 1).each_row() + (lavec.cols(Rbar, Rmax - 1) + tmp0))); //for sum (r - Rbar + 1) a_r phi_r, 
      // compute d
      arma::vec d   = tmp1.col(0); //for sum fr
      
      // marg. effects
      double meand  = mean(d);
      meffects      = meand*lbeta;
      
      // Qtheta
      Qth.submat(0, 0, Kz - 1, Kz - 1) = arma::eye<arma::mat>(Kz, Kz)*mean(d) - (lbeta*mean(Z.each_col()%(as<arma::vec>(ZtLambda)%tmp1.col(0) - tmp2.col(0)), 0));
      
      // Compute W
      arma::mat W(n, Rbar - 1);
      for(int k(0); k < (Rbar - 1); ++k) {
        W.col(k)            = -delta(k)*tmp1.col(1 + k); 
        Qth.col(Kz + k)     = delta(k)*lbeta*mean(as<arma::vec>(ZtLambda)%tmp1.col(1 + k) - tmp2.col(1 + k), 0);
      }
      W                     = arma::join_rows(Z.each_col() % d, W);
      
      // G*inv(S)*W
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
    // Compute Sigma0 abd Omega0
    arma::mat Sigma(Kz + Rbar - 1, Kz + Rbar - 1), Omega(Kz + Rbar - 1, Kz + Rbar - 1);
    
    {// fill block involving A, B, and C untill Rbar - 2
      arma::mat tmp1(n, Rbar - 1), tmp2(n, Rbar - 1), tmp3(n, Rbar - 1);
      for(int k(0); k < (Rbar - 1); ++ k){
        tmp1.col(k) = exp(lsumexp(2*lphi.cols(1 + k, Rmax - 1) - ldPhi.cols(1 + k, Rmax - 1))); //for sum fr^2/(Fr - Fr+1), 
        tmp2.col(k) = exp(lsumexp(2*lphi.cols(1 + k, Rmax - 1) - ldPhi.cols(k, Rmax - 2)));  //for sum fr^2/(Fr-1 - Fr)
        tmp3.col(k) = exp(lsumexp(lphi.cols(1 + k, Rmax - 1) + lphi.cols(2 + k, Rmax) - ldPhi.cols(1 + k, Rmax - 1))); //for sum fr fr+1/(Fr - Fr+1)
      }
      
      {// fill block of A
        arma::vec tmp                      = tmp1.col(0) + tmp2.col(0) - 2*tmp3.col(0); //A
        arma::mat Ztmp                     = arma::trans(Z.each_col() % tmp); //ZA
        Sigma.submat(0, 0, Kz - 1, Kz - 1) = Ztmp*Z;
        Omega.rows(0, Kz - 1)              = Ztmp*GinvSW;
      }
      
      {
        for(int k(0); k < (Rbar - 2); ++ k){
          {// fill block of B untill Rbar - 2
            arma::vec tmp4                          = delta(k)*(tmp3.col(k) - tmp2.col(k + 1) - tmp1.col(k + 1) + tmp3.col(k + 1));//B_k
            arma::vec Ztmp                          = Z.t()*tmp4;
            Sigma.submat(0, Kz + k, Kz - 1, Kz + k) = Ztmp;
            Sigma.submat(Kz + k, 0, Kz + k, Kz - 1) = Ztmp.t();
            Omega.row(Kz + k)                       = tmp4.t()*GinvSW;
            // fill block of C untill Rbar - 2
            for(int l(0); l <= k; ++ l){
              double tmp5           = -delta(l)*sum(tmp4);//C_kl
              Sigma(Kz + l, Kz + k) = tmp5;
              Sigma(Kz + k, Kz + l) = tmp5;
            }
          }
        }
      }
    }
    
    { // fill block B_Rbar and C_kRbar for all k
      arma::vec tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
      { 
        arma::mat tmp        = 2*lphi.cols(Rbar, Rmax - 1) - ldPhi.cols(Rbar - 1, Rmax - 2);
        tmp1 = exp(lsumexp(tmp.each_row() + tmp0)); //for sum (r - Rbar + 1) fr^2/(Fr-1 - Fr)
        tmp5 = exp(lsumexp(tmp.each_row() + 2*tmp0)); //for sum (r - Rbar + 1)^2 fr^2/(Fr-1 - Fr)
        
        tmp  = 2*lphi.cols(Rbar, Rmax - 1) - ldPhi.cols(Rbar, Rmax - 1);
        tmp2 = exp(lsumexp(tmp.each_row() + tmp0)); //for sum (r - Rbar + 1) fr^2/(Fr - Fr+1)
        tmp6 = exp(lsumexp(tmp.each_row() + 2*tmp0)); //for sum (r - Rbar + 1)^2 fr^2/(Fr - Fr+1) 
        
        tmp  = lphi.cols(Rbar, Rmax - 1) + lphi.cols(Rbar - 1, Rmax - 2) - ldPhi.cols(Rbar - 1, Rmax - 2);
        tmp3 = exp(lsumexp(tmp.each_row() + tmp0)); //for sum (r - Rbar + 1) fr fr-1/(Fr-1 - Fr)
        tmp7 = exp(lsumexp(tmp.each_row() + (tmp0 + log(arma::regspace<arma::rowvec>(0, Rmax - Rbar - 1))))); //for sum (r - Rbar + 1)(r - Rbar) fr fr-1/(Fr-1 - Fr)
        
        tmp  = lphi.cols(Rbar, Rmax - 1) + lphi.cols(Rbar + 1, Rmax) - ldPhi.cols(Rbar, Rmax - 1);
        tmp4 = exp(lsumexp(tmp.each_row() + tmp0)); //for sum (r - Rbar + 1) fr fr+1/(Fr - Fr+1)
        tmp8 = exp(lsumexp(tmp.each_row() + (tmp0 + log(arma::regspace<arma::rowvec>(2, Rmax - Rbar + 1))))); //for sum (r - Rbar + 1)(r - Rbar + 2) fr fr+1/(Fr - Fr+1)
      }
      
      { // fill block B_Rbar 
        arma::vec tmp9 = delta(Rbar - 2)*(tmp3 - tmp1 - tmp2 + tmp4); //B_Rbar
        arma::vec Ztmp = Z.t()*tmp9;
        Sigma.submat(0, Kz + Rbar - 2, Kz - 1, Kz + Rbar - 2) = Ztmp;
        Sigma.submat(Kz + Rbar - 2, 0, Kz + Rbar - 2, Kz - 1) = Ztmp.t();
        Omega.row(Kz + Rbar - 2)                              = tmp9.t()*GinvSW;
        
        // C_kRbar for all k excepted k = Rbar
        for(int k(0); k < (Rbar - 2); ++ k){
          double tmp10                 = -delta(k)*sum(tmp9);//C_lRbar
          Sigma(Kz + k, Kz + Rbar - 2) = tmp10;
          Sigma(Kz + Rbar - 2, Kz + k) = tmp10;
        }
      }
      
      { // fill block C_RbarRbar
        double tmp = pow(delta(Rbar - 2), 2)*sum((tmp5 - tmp7 - tmp8 + tmp6));//C_RbarRbar
        Sigma(Kz + Rbar - 2, Kz + Rbar - 2) = tmp;
      }
    }
    Omega          = lambda*Omega;
    arma::mat covt = arma::inv(Sigma + Omega);
    covt           = covt*Sigma*covt.t();
    arma::mat covm = Qth*covt*Qth.t();
    out            = List::create(Named("Rmax")        = Rmax, 
                                  Named("meffects")    = meffects,
                                  Named("covtheta")    = covt,
                                  Named("covmeffects") = covm,
                                  Named("var.comp")    = List::create(Named("Sigma") = Sigma/n, Named("Omega") = Omega/n));
  } else{
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
      next    = ((max(lfr) > -1000) || (Rmax <= Rbar));
    }
    
    // compute  the marginal effects
    arma::vec meffects = mean(exp(lsumexp(lphi.cols(1, Rmax - 1))))*lbeta;
    out            = List::create(Named("Rmax")       = Rmax, 
                                  Named("meffects")   = meffects);
  }
  return out;
}
