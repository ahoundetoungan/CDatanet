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
 *           igroup[s,] is a 2-dimension vector ans gives the first and the last rows
 *           for the group s.
 * Ca      : Is the group for the cut-points.
 * nCa     : Is the number of groups Ca.
 * Cl      : Is the group for lambda.
 * nCl     + Is the number of groups Cl.
 * ngroup  : is the number of groups.
 * theta   : is the vector of parameters ordered as follow: peer effects, explanatory variables
 *           coefficients, delta2, delta3 ..., deltaRbar.
 * n       : The sample size.
 * tol     : A tolerance value for the iterative method to compute P convergence. 
 * maxit   : The maximum number of iterations of the iterative method to compute P. If this
 *           number is reached, the algorithm stops and the last P is used as the solution. maxit
 *           is important for numerical reasons if tol is too small. For example a tol = 1e-16
 *           may not be reached and the algorithm will stop after maxit iterations.
 * ye      : is the vector of equilibrium outcome expectation.
 * R       : is the upper bound of y.
 */

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]

#include <RcppArmadillo.h>
//#define NDEBUG
#include <RcppNumerical.h>
#include <RcppEigen.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// flambdat converts lambda to lambdatilde given lambda in (a, b)
//[[Rcpp::export]]
arma::vec fcdlambdat(const arma::vec& lambda, 
                     const int& nCa,
                     const double& a, 
                     const double& b){
  arma::vec lt(lambda);
  if(b == R_PosInf){ //a to +Inf 
    for(int ca(0); ca < nCa; ++ ca){
      double sl(sum(lambda.subvec(nCa*ca, nCa*(ca + 1) - 1)));
      lt(nCa*(ca + 1) - 1) = log(sl - a);
    }
    return lt; 
  }
  //a to b
  for(int ca(0); ca < nCa; ++ ca){
    double sl(sum(lambda.subvec(nCa*ca, nCa*(ca + 1) - 1)));
    lt(nCa*(ca + 1) - 1)   = log(sl - a) - log(b - sl);
  }
  return lt; 
}

//fdlambdat computes derivative of lambdatilde wrt lambda
//[[Rcpp::export]]
arma::mat fcddlambdat(const arma::vec& lambda, 
                      const int& nCa,
                      const double& a, 
                      const double& b){
  double sl(sum(lambda)), tp1(1/(sl - a)), tp2(1/(b - sl));
  arma::mat out(nCa, nCa, arma::fill::eye);
  out.row(nCa - 1).zeros();
  if(b == R_PosInf){//a to +Inf
    out.row(nCa - 1) += tp1;
  } else{ //a to b
    out.row(nCa - 1) += (tp1 + tp2);
  }
  return out; 
}


// flambda converts lambdatilde to lambda given lambda in (a, b)
//[[Rcpp::export]]
arma::vec fcdlambda(const arma::vec& lambdat, 
                    const int& nCa,
                    const double& a, 
                    const double& b){
  arma::vec l(lambdat);
  if(b == R_PosInf){//a to +Inf
    for(int ca(0); ca < nCa; ++ ca){
      double sl(exp(lambdat(nCa*(ca + 1) - 1)) + a);
      l(nCa*(ca + 1) - 1) = sl - sum(l.subvec(nCa*ca, nCa*(ca + 1) - 1).head(nCa - 1));
    }
    return l; 
  }
  //a to b
  for(int ca(0); ca < nCa; ++ ca){
    double  sl((b - a)/(1 + exp(-lambdat(nCa*(ca + 1) - 1))) + a);
    l(nCa*(ca + 1) - 1) = sl - sum(l.subvec(nCa*ca, nCa*(ca + 1) - 1).head(nCa - 1));
  }
  return l; 
}

//fdlambda computes derivative of lambda wrt lambdatilde
//[[Rcpp::export]]
arma::mat fcddlambda(const arma::vec& lambda, 
                     const int& nCa,
                     const double& a, 
                     const double& b){
  double sl(sum(lambda));
  arma::mat out(nCa, nCa, arma::fill::eye);
  if(nCa > 1){
    out.submat(nCa - 1, 0, nCa - 1, nCa - 2) = (-arma::ones<arma::rowvec>(nCa - 1));
  }
  if(b == R_PosInf){//a to +Inf
    out(nCa - 1, nCa - 1) = (sl - a);
  } else{ //a to b
    out(nCa - 1, nCa - 1) = (b - sl)*(sl - a)/(b - a);
  }
  return out; 
}

// fgamma select gamma
double fgamma(const arma::vec& delta,
              const int& iida,
              const int& r,
              const int& Rbar,
              const double& R){
  // if (r == 0) return R_NegInf;
  if (r == 1) return 0;
  if (r <= Rbar) return delta(iida + r - 2);
  if (r > R) return R_PosInf;
  return delta(iida + Rbar - 1);
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
 *  ZtLamba : psi + lambda gi*ye, where psi may be a vector of a matrix (many simulations)
 *  lambda  : is the peer effect
 *  delta   : vector of delta for groups, 1, 2, ....
 *  idelta  : matrix of each group delta positions. idelta(s, 0) and idelta(s, 1) give the starting and
 *  ending positions of group s delta.
 */

// computation of fL as the mapping L as in the paper.
// step 1: log fL
// lofLF = log(p_1 + ...)
//       = log(p_1) + log(1 + exp(log(p2) - log(p1)) + exp(log(p3) - log(p1)) + ...)
//       = lF1 + log(1 + sum_r=2^inf exp(lFr - lF1))
arma::vec flogL(const arma::vec& ZtLambda,
                List& lCa,
                const int& nCa,
                const arma::vec& delta,
                const arma::umat& idelta,
                const arma::vec& Rbar,
                const double& R,
                const arma::vec& n,
                const int sumn) {
  arma::vec logL(sumn);
  for(int ca(0); ca < nCa; ++ ca){
    arma::uvec lCas   = lCa[ca];
    double ar         = 0;
    int r             = 1;
    arma::vec ZtLs    = ZtLambda.elem(lCas);
    NumericVector tp0 = wrap(ZtLs);
    NumericVector lF1 = Rcpp::pnorm5(tp0, 0, 1, true, true);
    NumericVector lu2 = Rcpp::rep(NumericVector::create(1), n[ca]);
    NumericVector lFr;
    bool next(r < R);
    while(next) {
      ++ r;
      ar     += fgamma(delta, idelta(ca, 0), r, Rbar(ca), R);
      lFr     = Rcpp::pnorm5(tp0 - ar, 0, 1, true, true); 
      lu2    += exp(lFr - lF1); 
      next    = (((max(lFr) > -1000) || (r <= Rbar(ca))) && (r < R));
    }
    NumericVector tp1 = lF1 + log(lu2);
    logL.elem(lCas)   = as<arma::vec>(tp1);
  }
  return logL;
}

// step2 fL
//[[Rcpp::export]]
arma::vec fL(const arma::vec& ZtLambda,
             List& lCa,
             const int& nCa,
             const arma::vec& delta,
             const arma::umat& idelta,
             const arma::vec& Rbar,
             const double& R,
             const arma::vec& n,
             const int sumn) {
  return exp(flogL(ZtLambda, lCa, nCa, delta, idelta, Rbar, R, n, sumn));
}

// step2' Sum(fL_u) where u is several draws of explanatory variables
//[[Rcpp::export]]
arma::vec fLncond(const arma::mat& ZtLambda,
                  List& lCa,
                  const int& nCa,
                  const arma::mat& delta,
                  const arma::umat& idelta,
                  const arma::vec& Rbar,
                  const double& R,
                  const arma::vec& n,
                  const int sumn,
                  const int& nsimu) {
  arma::mat logLs(sumn, nsimu);
  for(int s(0); s < nsimu; ++s){
    logLs.col(s) = flogL(ZtLambda.col(s), lCa, nCa, delta, idelta, Rbar, R, n, sumn);
  }
  return exp(laverexp(logLs, nsimu));
}


// fye: Takes an initial value of ye and finds the equilibrium
// Gye is (G1*ye, ..., GnCl*ye)
//[[Rcpp::export]]
int fye(arma::vec& ye,
        arma::mat& Gye,
        List& G,
        List& lCa,
        const int& nCa,
        const arma::mat& igroup,
        const int& ngroup,
        const arma::vec& psi,
        const arma::vec& lambda,
        const arma::vec& delta,
        const arma::umat& idelta,
        const arma::vec& n,
        const int sumn,
        const arma::vec& Rbar,
        const double& R,
        const double& tol,
        const int& maxit) {
  int n1, n2, t = 0, nCl(nCa*nCa);
  
  computeL: ++t;
  arma::vec ZtLambda(Gye*lambda + psi);
  arma::vec yest = fL(ZtLambda, lCa, nCa, delta, idelta, Rbar, R, n, sumn);
  
  // double dist    = max(arma::abs(yest/(ye + 1e-50) - 1));
  double dist    = max(arma::abs((yest - ye)/(ye + 1e-50)));
  ye             = yest;
  
  for (int m(0); m < ngroup; ++ m) {
    n1           = igroup(m,0);
    n2           = igroup(m,1);
    List Gm      = G[m];
    for(int cl(0); cl < nCl; ++ cl){
      arma::mat Gms            = Gm[cl];
      Gye.submat(n1, cl, n2, cl) = Gms*ye.subvec(n1, n2);
    }
  }
  if (dist > tol && t < maxit) goto computeL;
  return t; 
}

// nonconditional version
//[[Rcpp::export]]
int fyencond(arma::vec& ye,
             arma::mat& Gye,
             List& G,
             List& lCa,
             const int& nCa,
             const arma::mat& igroup,
             const int& ngroup,
             const arma::mat& psi,
             const double& lambda,
             const arma::vec& delta,
             const arma::umat& idelta,
             const arma::vec& n,
             const int sumn,
             const int& nsimu,
             const arma::vec& Rbar,
             const double& R,
             const double& tol,
             const int& maxit) {
  int n1, n2, t(0), nCl(nCa*nCa);
  
  computeL: ++t;
  arma::vec ZtLambda(psi.each_col() + (Gye*lambda));
  arma::vec yest = fLncond(ZtLambda, lCa, nCa, delta, idelta, Rbar, R, n, sumn, nsimu);
  // double dist    = max(arma::abs(yest/(ye + 1e-50) - 1));
  double dist    = max(arma::abs((yest - ye)/(ye + 1e-50)));
  ye             = yest;
  
  for (int m(0); m < ngroup; ++ m) {
    n1           = igroup(m,0);
    n2           = igroup(m,1);
    List Gm      = G[m];
    for(int cl(0); cl < nCl; ++ cl){
      arma::mat Gms            = Gm[cl];
      Gye.submat(n1, cl, n2, cl) = Gms*ye.subvec(n1, n2);
    }
  }
  
  if (dist > tol && t < maxit) goto computeL;
  return t; 
}


// fy returns the vector of y given yst and delta
//[[Rcpp::export]]
arma::vec fy(const arma::vec& yst,
             const arma::vec& maxyst,
             List& lCa,
             const int& nCa,
             const arma::vec& delta,
             const arma::umat& idelta,
             const arma::vec& n,
             const int& sumn,
             const arma::vec& Rbar,
             const double& R) {
  arma::vec y(sumn);
  for(int ca(0); ca < nCa; ++ ca){
    arma::vec ys(n(ca), arma::fill::zeros);
    arma::uvec lCas = lCa[ca];
    bool cont       = true;
    int r           = 1;
    double ar       = 0;//this is a1
    while(cont) {
      ys.elem(arma::find(yst.elem(lCas) > ar)) += 1;
      ++r;
      ar           += fgamma(delta, idelta(ca, 0), r, Rbar(ca), R);
      cont          = (maxyst(ca) > ar);
    }
    y.elem(lCas)    = ys;
  }
  return y;
}

// fmeffect computes the marginal effects
//[[Rcpp::export]]
List fmeffects(const arma::vec& ZtLambda,
               const arma::vec& lambda,
               const arma::vec& Gamma2,
               List& lCa,
               const int& nCa,
               const arma::vec& delta,
               const arma::umat& idelta,
               const int& sumn,
               const arma::vec& Rbar,
               const double& R){
  arma::mat lphi(sumn, 0);
  int Rmax(0);
  bool next(Rmax < R);
  arma::vec a(nCa, arma::fill::zeros);
  
  // log of f 
  while(next) {
    ++ Rmax;
    arma::vec tp1(sumn);
    for(int ca(0); ca < nCa; ++ ca){
      a(ca)            += fgamma(delta, idelta(ca, 0), Rmax, Rbar(ca), R);
      arma::uvec lCas   = lCa[ca];
      arma::vec ZtLs    = ZtLambda.elem(lCas);
      NumericVector tp0 = wrap(ZtLs);
      NumericVector lfr = Rcpp::dnorm4(tp0 - a(ca), 0, 1, true);
      tp1.elem(lCas)    = as<arma::vec>(lfr);
    }
    lphi    = arma::join_rows(lphi, tp1);
    next    = (((max(tp1) > -1000) || (Rmax <= (Rbar.max() + 2))) && (Rmax < R));
  }
  
  // compute  the marginal effects
  arma::vec tp0(exp(lsumexp(lphi.head_cols(Rmax))));
  // for Gye 
  arma::mat meff1(sumn, nCa*nCa, arma::fill::zeros);
  arma::uvec tp(arma::regspace<arma::uvec>(0, nCa - 1));
  for(int ca(0); ca < nCa; ++ ca){
    arma::uvec lCas  = lCa[ca];
    meff1.submat(lCas, tp) = tp0.elem(lCas)*lambda.elem(tp).t();
    tp += nCa;
  }
  
  // for Z
  arma::mat meff2(tp0*Gamma2.t());
  // for (int m(0); m < ngroup; ++ m) {
  //   int n1(igroup(m,0));
  //   int n2(igroup(m,1));
  //   int nm(n2 - n1 + 1);
  //   List Gm = G[m];
  //   arma::mat Sm(nm, nm, arma::fill::zeros);
  //   for(int cl(0); cl < nCl; ++ cl){
  //     arma::mat Gms = Gm[cl];
  //     Sm           += lambda(cl)*(Gms);//.each_col() % tp0.subvec(n1, n2));
  //   }
  //   meff2.rows(n1, n2)  = arma::solve(arma::eye<arma::mat>(nm, nm) - Sm, meff2.rows(n1, n2));
  // }
  return List::create(Named("imeff") = arma::join_rows(meff1, meff2), Named("Rmax") = Rmax);
}

// flogP returns the log-likelihood of each individual
// maxy is max(y)
arma::vec flogp(const arma::uvec& y,
                const arma::vec& ZtLambda,
                const arma::vec& maxy,
                List& lCa,
                const int& nCa,
                const arma::vec& delta,
                const arma::umat& idelta,
                const arma::vec& n,
                const int& sumn,
                const arma::vec& Rbar,
                const double& R){
  arma::vec out(sumn);
  for(int ca(0); ca < nCa; ++ ca){
    // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
    // Rbar <= maxy + 1 
    arma::vec a(maxy(ca) + 2); 
    a(0)   = R_NegInf; a(1) = 0; 
    for(int r(2); r < (maxy(ca) + 2); ++r) {
      a(r) = a(r - 1) + fgamma(delta, idelta(ca, 0), r, Rbar(ca), R);
    }
    
    // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
    arma::uvec lCas   = lCa[ca];
    NumericVector tp1 = wrap(ZtLambda.elem(lCas) - a.elem(y.elem(lCas)));
    NumericVector tp2 = wrap(ZtLambda.elem(lCas) - a.elem(y.elem(lCas) + 1));
    
    // log likelihood
    NumericVector lFleft  = Rcpp::pnorm5(tp1, 0, 1, true, true);
    NumericVector lFright = Rcpp::pnorm5(tp2, 0, 1, true, true);
    NumericVector tp0     = lFleft + log(1 - exp(lFright - lFleft));
    out.elem(lCas)        = as<arma::vec>(tp0);
  }
  return out;
}

//non conditional version
arma::vec flogpncond(const arma::uvec& y,
                     const arma::mat& ZtLambda,
                     const arma::vec& maxy,
                     List& lCa,
                     const int& nCa,
                     const arma::vec& delta,
                     const arma::umat& idelta,
                     const arma::vec& n,
                     const int& sumn,
                     const arma::vec& Rbar,
                     const double& R,
                     const int& nsimu){
  arma::mat logP(sumn, nsimu);
  for(int ca(0); ca < nCa; ++ ca){
    // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
    // Rbar <= maxy + 1
    arma::vec a(maxy(ca) + 2); 
    a(0)   = R_NegInf; a(1) = 0; 
    for(int r(2); r < (maxy(ca) + 2); ++r) {
      a(r) = a(r - 1) + fgamma(delta, idelta(ca, 0), r, Rbar(ca), R);
    }
    
    // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
    arma::uvec lCas   = lCa[ca];
    arma::mat tp      = ZtLambda.rows(lCas);
    NumericMatrix tp1 = wrap(tp.each_col() - a.elem(y.elem(lCas)));
    NumericMatrix tp2 = wrap(tp.each_col() - a.elem(y.elem(lCas) + 1));
    
    // log likelihood
    arma::mat tp0(n(ca), nsimu);
    NumericVector lFlefts, lFrights, tp3;
    for(int s(0); s < nsimu; ++s) {
      lFlefts       = Rcpp::pnorm5(tp1(_,s), 0, 1, true, true);
      lFrights      = Rcpp::pnorm5(tp2(_,s), 0, 1, true, true);
      tp3           = lFlefts + log(1 - exp(lFrights - lFlefts));
      tp0.col(s)    = as<arma::vec>(tp3);
    }
    logP.rows(lCas) = tp0;
  }
  return laverexp(logP, nsimu);
}

//fdelta computes delta from deltat, lambda and nCa
//[[Rcpp::export]]
arma::vec fdelta(const arma::vec& deltat, 
                 const arma::vec& lambda,
                 const arma::umat& idelta,
                 const arma::uvec& ndelta,
                 const int& nCa){
  arma::vec delta(deltat);
  for(int ca(0); ca < nCa; ++ ca){
    if(ndelta(ca) > 0){
      delta.subvec(idelta(ca, 0), idelta(ca, 1)) += sum(lambda.subvec(nCa*ca, nCa*(ca + 1) - 1));
    }
  }
  return delta;
}

// foptimREM: compute the likelihood given the patameters
//[[Rcpp::export]]
double foptimREM(arma::vec& ye,
                 arma::mat& Gye,
                 const arma::vec& theta,
                 const double& lb_sl,
                 const double& ub_sl,
                 const arma::mat& X,
                 List& G,
                 List& lCa,
                 const int& nCa,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const int& K,
                 const arma::vec& n,
                 const int sumn,
                 const arma::vec& Rbar,
                 const double& R,
                 const arma::umat& idelta,
                 const arma::uvec& ndelta,
                 const arma::uvec& y,
                 const arma::vec& maxy,
                 const double& tol = 1e-13,
                 const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim")   = R_NilValue;
  Rcpp::print(thetacpp);
  int nCl(nCa*nCa);
  
  arma::vec psi(X*theta.subvec(nCl, nCl + K - 1));
  arma::vec lambda(fcdlambda(theta.head(nCl), nCa, lb_sl, ub_sl));
  arma::vec delta(fdelta(exp(theta.tail(sum(ndelta))) + 1e-323, lambda, idelta, ndelta, nCa)); /*bar delta is in the vector delta. so theta as delta2 to deltaRbar and bardelta*/
  
  // compute ye
  fye(ye, Gye, G, lCa, nCa, igroup, ngroup, psi, lambda, delta, idelta, n, sumn, Rbar, R, tol, maxit);
  
  arma::vec ZtLambda = Gye*lambda + psi;
  arma::vec logp     = flogp(y, ZtLambda, maxy, lCa, nCa, delta, idelta, n, sumn, Rbar, R);
  double llh         = sum(logp);
  
  // if(llh < -1e308) {
  //   llh           = -1e308;
  // }
  return -llh;
}

// //non conditional version
// // one simulation
// //[[Rcpp::export]]
// double foptimREMncond1(arma::vec& ye,
//                        arma::vec& Gye,
//                        const arma::vec& theta,
//                        const arma::mat& X,
//                        const arma::mat& Simu1,
//                        const int& nsimu,
//                        List& G,
//                        const arma::mat& igroup,
//                        const int& ngroup,
//                        const int& K,
//                        const int& n,
//                        const int& Rbar,
//                        const arma::uvec& y,
//                        const int& maxy,
//                        const double& tol = 1e-13,
//                        const int& maxit  = 1e3) {
//   NumericVector thetacpp = wrap(theta);
//   thetacpp.attr("dim")   = R_NilValue;
//   Rcpp::print(thetacpp);
//   double lambda   = 1.0/(exp(-theta(0)) + 1);
//   arma::mat psi   = Simu1*theta(K+1);
//   arma::vec tmp   = X*theta.subvec(1, K);
//   psi.each_col() += tmp;
//   arma::vec delta = exp(theta.tail(Rbar)) + lambda + 1e-323;
//   
//   // compute ye
//   fyencond(ye, Gye, G, igroup, ngroup, psi, lambda, delta, n, nsimu, Rbar, tol, maxit);
//   
//   arma::vec ZtLambda = psi.each_col() + lambda*Gye;
//   arma::vec logp     = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
//   double llh         = sum(logp);
//   
//   // if(llh < -1e308) {
//   //   llh              = -1e308;
//   // }
//   return -llh;
// }


//// NPL solution
//[[Rcpp::export]]
double foptimREM_NPL(const arma::mat& Gye,
                     const arma::vec& theta,
                     const double& lb_sl,
                     const double& ub_sl,
                     const arma::mat& X,
                     List& lCa,
                     const int& nCa,
                     const int& K,
                     const arma::vec& n,
                     const int sumn,
                     const arma::umat& idelta,
                     const arma::uvec& ndelta,
                     const arma::vec& Rbar,
                     const double& R,
                     const arma::uvec& y,
                     const arma::vec& maxy,
                     const bool& print = false) {
  // print
  if(print){
    NumericVector thetacpp = wrap(theta);
    thetacpp.attr("dim")   = R_NilValue;
    Rcpp::Rcout << "parms: \n";
    Rcpp::print(thetacpp);
  }
  int nCl(nCa*nCa);
  
  arma::vec lambda(fcdlambda(theta.head(nCl), nCa, lb_sl, ub_sl));
  arma::vec ZtLambda(Gye*lambda + X*theta.subvec(nCl, nCl + K - 1));
  arma::vec delta(fdelta(exp(theta.tail(sum(ndelta))) + 1e-323, lambda, idelta, ndelta, nCa)); /*bar delta is in the vector delta. so theta as delta2 to deltaRbar and bardelta*/
  
  arma::vec logp     = flogp(y, ZtLambda, maxy, lCa, nCa, delta, idelta, n, sumn, Rbar, R);
  double llh         = sum(logp);
  if((llh < -1e250) || R_IsNaN(llh)) {
    llh              = -1e250;
  }
  if(print) Rcpp::Rcout << "log-likelihood: " << llh << "\n";
  return -llh;
}

//[[Rcpp::export]]
void fL_NPL(arma::vec& ye,
            arma::mat& Gye,
            const arma::vec& theta,
            const arma::mat& X,
            List& G,
            List& lCa,
            const int& nCa,
            const arma::mat& igroup,
            const int& ngroup,
            const int& K,
            const arma::vec& n,
            const int sumn,
            const arma::umat& idelta,
            const arma::uvec& ndelta,
            const arma::vec& Rbar,
            const double& R) {
  int n1, n2, nCl(nCa*nCa);
  arma::vec lambda(theta.head(nCl));
  arma::vec ZtLambda(Gye*lambda + X*theta.subvec(nCl, nCl + K - 1));
  arma::vec delta(fdelta(theta.tail(sum(ndelta)) + 1e-323, lambda, idelta, ndelta, nCa)); /*bar delta is in the vector delta. so theta as delta2 to deltaRbar and bardelta*/
  // new ye
  ye.subvec(0, sumn - 1) = fL(ZtLambda, lCa, nCa, delta, idelta, Rbar, R, n, sumn);
  
  // new Gye
  for (int m(0); m < ngroup; ++ m) {
    n1      = igroup(m, 0);
    n2      = igroup(m, 1);
    List Gm = G[m];
    for(int cl(0); cl < nCl; ++ cl){
      arma::mat Gms            = Gm[cl];
      Gye.submat(n1, cl, n2, cl) = Gms*ye.subvec(n1, n2);
    }
  }
}

//[[Rcpp::export]]
void fnewye(arma::vec& ye,
            arma::mat& Gye,
            const arma::vec& theta,
            const arma::mat& X,
            List& G,
            List& lCa,
            const int& nCa,
            const arma::mat& igroup,
            const int& ngroup,
            const int& K,
            const arma::vec& n,
            const int sumn,
            const arma::umat& idelta,
            const arma::uvec& ndelta,
            const arma::vec& Rbar,
            const double& R,
            const double& tol,
            const int& maxit) {
  int nCl(nCa*nCa);
  arma::vec psi(X*theta.subvec(nCl, nCl + K - 1));
  arma::vec lambda(theta.head(nCl));
  arma::vec ZtLambda(Gye*lambda + psi);
  arma::vec delta(fdelta(theta.tail(sum(ndelta)) + 1e-323, lambda, idelta, ndelta, nCa)); /*bar delta is in the vector delta. so theta as delta2 to deltaRbar and bardelta*/
  fye(ye, Gye, G, lCa, nCa, igroup, ngroup, psi, lambda, delta, idelta, n, sumn, Rbar, R, tol, maxit);
}

// // non conditional version
// // using one factor
// //[[Rcpp::export]]
// double foptimREM_NPLncond1(const arma::vec& Gye,
//                            const arma::vec& theta,
//                            const arma::mat& X,
//                            const arma::mat& Simu1,
//                            const int& nsimu,
//                            const int& Rbar,
//                            const int& maxy,
//                            const int& K,
//                            const int& n,
//                            const arma::uvec& y) {
//   double lambda        = 1.0/(exp(-theta(0)) + 1);
//   arma::mat ZtLambda   = Simu1*theta(K+1);
//   arma::vec tmp        = lambda*Gye + X*theta.subvec(1, K);
//   ZtLambda.each_col() += tmp;
//   arma::vec delta      = exp(theta.tail(Rbar)) + lambda + 1e-323;
//   arma::vec logp       = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
//   double llh           = sum(logp);
//   // if(llh < -1e250) {
//   //   llh                = -1e250;
//   // }
//   return -llh;
// }

// //[[Rcpp::export]]
// void fL_NPLncond1(arma::vec& ye,
//                   arma::vec& Gye,
//                   List& G,
//                   const arma::mat& igroup,
//                   const int& ngroup,
//                   const arma::mat& X,
//                   const arma::vec& theta,
//                   const arma::mat& Simu1,
//                   const int& nsimu,
//                   const int& Rbar,
//                   const int& K,
//                   const int& n) {
//   int n1, n2;
//   double lambda        = 1.0/(exp(-theta(0)) + 1);
//   arma::mat ZtLambda   = Simu1*theta(K+1);
//   arma::vec tmp        = lambda*Gye + X*theta.subvec(1, K);
//   ZtLambda.each_col() += tmp;
//   arma::vec delta      = exp(theta.tail(Rbar)) + lambda + 1e-323;
//   
//   // new ye
//   ye.subvec(0, n-1)    = fLncond(ZtLambda, delta, Rbar, n, nsimu);
//   
//   // new Gye
//   for (int m(0); m < ngroup; ++ m) {
//     n1                 = igroup(m,0);
//     n2                 = igroup(m,1);
//     arma::mat Gm       = G[m];
//     Gye.subvec(n1, n2) = Gm*ye.subvec(n1, n2);
//   }
// }


// // using two factors
// //[[Rcpp::export]]
// double foptimREM_NPLncond2(const arma::vec& Gye,
//                            const arma::vec& theta,
//                            const arma::mat& X,
//                            const arma::mat& Simu1,
//                            const arma::mat& Simu2,
//                            const int& nsimu,
//                            const int& Rbar,
//                            const int& maxy,
//                            const int& K,
//                            const int& n,
//                            const arma::uvec& y) {
//   double lambda        = 1.0/(exp(-theta(0)) + 1);
//   arma::mat ZtLambda   = Simu1*theta(K+1) + Simu2*theta(K+2);
//   arma::vec tmp        = lambda*Gye + X*theta.subvec(1, K);
//   ZtLambda.each_col() += tmp;
//   arma::vec delta      = exp(theta.tail(Rbar)) + lambda + 1e-323;
//   arma::vec logp       = flogpncond(y, ZtLambda, maxy, lambda, delta, nsimu, Rbar, n);
//   double llh           = sum(logp);
//   // if(llh < -1e250) {
//   //   llh                = -1e250;
//   // }
//   return -llh;
// }

// //[[Rcpp::export]]
// void fL_NPLncond2(arma::vec& ye,
//                   arma::vec& Gye,
//                   List& G,
//                   const arma::mat& igroup,
//                   const int& ngroup,
//                   const arma::mat& X,
//                   const arma::vec& theta,
//                   const arma::mat& Simu1,
//                   const arma::mat& Simu2,
//                   const int& nsimu,
//                   const int& Rbar,
//                   const int& K,
//                   const int& n) {
//   int n1, n2;
//   double lambda        = 1.0/(exp(-theta(0)) + 1);
//   arma::mat ZtLambda   = Simu1*theta(K+1) + Simu2*theta(K+2);
//   arma::vec tmp        = lambda*Gye + X*theta.subvec(1, K);
//   ZtLambda.each_col() += tmp;
//   arma::vec delta      = exp(theta.tail(Rbar)) + lambda + 1e-323;
//   
//   // new ye
//   ye.subvec(0, n-1)    = fLncond(ZtLambda, delta, Rbar, n, nsimu);
//   
//   // new Gye
//   for (int m(0); m < ngroup; ++ m) {
//     n1                 = igroup(m,0);
//     n2                 = igroup(m,1);
//     arma::mat Gm       = G[m];
//     Gye.subvec(n1, n2) = Gm*ye.subvec(n1, n2);
//   }
// }


// Numerical optimization using Rcpp
class cdnetreg: public MFuncGrad
{
private:
  const double& lb_sl;
  const double& ub_sl;
  const arma::mat& Z;
  List& lCa;
  const int& nCa;
  const arma::umat& idelta;
  const arma::uvec& ndelta;
  const arma::vec& Rbar;
  const double& R;
  const arma::vec& maxy;
  const int& K;
  const arma::uvec& y;
  List& lidy;
  const double print;
public:
  cdnetreg(const double& lb_sl_,
           const double& ub_sl_,
           const arma::mat& Z_,
           List& lCa_,
           const int& nCa_,
           const arma::umat& idelta_,
           const arma::uvec& ndelta_,
           const arma::vec& Rbar_,
           const double& R_,
           const arma::vec& maxy_,
           const int& K_,
           const arma::uvec& y_,
           List& lidy_,
           const double print_) : 
  lb_sl(lb_sl_),
  ub_sl(ub_sl_),
  Z(Z_),
  lCa(lCa_),
  nCa(nCa_),
  idelta(idelta_),
  ndelta(ndelta_),
  Rbar(Rbar_),
  R(R_),
  maxy(maxy_),
  K(K_),
  y(y_),
  lidy(lidy_),
  print(print_){}
  
  Eigen::VectorXd Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0(theta);  //make a copy
    // int ila (0), ilb(nCl - 1), iba(nCl), ibb(nCl + K - 1);
    int nCl(nCa*nCa);
    arma::uvec ida(nCl + K + idelta.col(0)), idb(nCl + K + idelta.col(1));
    int nds(sum(ndelta));
    int parms(nCl + K + nds);
    arma::vec beta(arma::vec(theta0.data(), parms, false, false)); //convert into arma vec
    
    beta.head(nCl)         = fcdlambda(beta.head(nCl), nCa, lb_sl, ub_sl); 
    arma::vec lambda(beta.head(nCl));
    arma::vec deltat(exp(beta.tail(nds)) + 1e-323);
    beta.tail(nds)         = fdelta(deltat, lambda, idelta, ndelta, nCa);
    arma::vec delta(beta.tail(nds));
    
    // print
    if(print){
      NumericVector betacpp   = wrap(beta);
      betacpp.attr("dim")     = R_NilValue;
      Rcpp::Rcout << "Estimate: \n";
      Rcpp::print(betacpp);
    }
    
    arma::vec ZtLambda(Z*beta.head(nCl + K)), gdarm(parms, arma::fill::zeros);
    double f(0);
    for(int ca(0); ca < nCa;  ++ ca){
      arma::uvec lCas   = lCa[ca];
      arma::uvec ys(y.elem(lCas));
      arma::vec ZtLs(ZtLambda.elem(lCas));
      arma::mat Zs(Z.rows(lCas));
      List lidys = lidy[ca];
      bool hdbar(R > Rbar(ca));
      
      // create vector a = [a_0, a_1, ..., a_(maxy + 1)]
      // Rbar <= maxy + 1
      arma::vec a(maxy(ca) + 2); 
      a(0)   = R_NegInf; a(1) = 0; 
      for(int r(2); r < (maxy(ca) + 2); ++r) {
        a(r) = a(r - 1) + fgamma(delta, idelta(ca, 0), r, Rbar(ca), R);
      }
      
      // compute Ztlambda - ar as tmp1 and Ztlambda - ar+1 as tmp2 
      arma::vec Zba1       = ZtLs - a.elem(ys);
      arma::vec Zba2       = ZtLs - a.elem(ys + 1);
      NumericVector Zba1r  = wrap(Zba1);
      NumericVector Zba2r  = wrap(Zba2);
      
      arma::vec lfZba1     = Rcpp::dnorm4(Zba1r, 0, 1, true);
      arma::vec lfZba2     = Rcpp::dnorm4(Zba2r, 0, 1, true);
      
      NumericVector lFZba1 = Rcpp::pnorm5(Zba1r, 0, 1, true, true);
      NumericVector lFZba2 = Rcpp::pnorm5(Zba2r, 0, 1, true, true);
      
      // log likelihood
      NumericVector ldF    = lFZba1 + log(1 - exp(lFZba2 - lFZba1));
      f                   += sum(ldF);
      
      // gradient
      arma::vec ldFa       = as<arma::vec>(ldF);
      arma::vec tmp1       = exp(lfZba1 - ldFa);
      arma::vec tmp2       = exp(lfZba2 - ldFa);
      gdarm.head(nCl + K) += (arma::sum(Zs.each_col()%(tmp1 - tmp2), 0)).t();
      for(int r(2); r <= Rbar(ca); ++ r){
        arma::uvec id2     = lidys[r - 2]; // greater than r - 1
        arma::uvec id1     = lidys[r - 1]; // greater than r
        double tmp3        = sum(tmp2.elem(id2)) - sum(tmp1.elem(id1));
        gdarm.head(nCl)         += tmp3;
        gdarm(ida(ca) + r - 2)  += tmp3;
      }
      if(hdbar){
        arma::uvec id2     = lidys[Rbar(ca) - 1]; // greater than Rbar
        arma::uvec id1     = lidys[Rbar(ca)]; // greater than Rbar + 1
        double tmp3        = sum(tmp2.elem(id2)%(ys.elem(id2) + 1 - Rbar(ca))) - sum(tmp1.elem(id1)%(ys.elem(id1) - Rbar(ca)));
        gdarm.head(nCl)   += tmp3;
        gdarm(idb(ca))    += tmp3;
      }
    }
    
    gdarm.tail(nds) %= deltat;
    for(int ca(0); ca < nCa; ++ ca){
      arma::mat tp(arma::trans(fcddlambda(lambda, nCa, lb_sl, ub_sl)));
      gdarm.subvec(ca*nCa, nCa*(ca + 1) - 1) = tp*gdarm.subvec(ca*nCa, nCa*(ca + 1) - 1);
    }
    grad             = -Eigen::Map<Eigen::VectorXd>(gdarm.memptr(), parms);
    Grad             = -grad;
    
    if((f > 1e250) || R_IsNaN(f)) {
      f              = 1e250;
    }
    
    if(print) Rcpp::Rcout << "log-likelihood: " << f << "\n";
    return -f;
  }
};


//[[Rcpp::export]]
List cdnetLBFGS(Eigen::VectorXd par,
                const double& lb_sl,
                const double& ub_sl,
                const arma::mat& Gye,
                const arma::mat& X,
                List& lCa,
                const int& nCa,
                const arma::vec& n,
                const int sumn,
                const arma::umat& idelta,
                const arma::uvec& ndelta,
                const arma::vec& Rbar,
                const double& R,
                const arma::vec& maxy,
                const int& K,
                const arma::uvec& y,
                const int& maxit = 300, 
                const double& eps_f = 1e-13, 
                const double& eps_g = 1e-13,
                const bool& print = false) {
  arma::mat Z        = arma::join_rows(Gye, X);
  List lidy(nCa);
  for(int ca(0); ca < nCa; ++ ca){
    arma::uvec lCas  = lCa[ca];
    arma::uvec ys(y.elem(lCas));
    List lidys(Rbar(ca) + 1);
    for(int r(0); r <= Rbar(ca); ++ r){
      arma::uvec id  = arma::find(ys >= (r + 1));   
      lidys(r)       = id;
    }
    lidy[ca]         = lidys;
  }
  
  double fopt;
  cdnetreg f(lb_sl, ub_sl, Z, lCa, nCa, idelta, ndelta, Rbar, R, maxy, K, y, lidy, print);
  
  int status(optim_lbfgs(f, par, fopt, maxit, eps_f, eps_g));
  Eigen::VectorXd grad(f.Grad);
  
  return Rcpp::List::create(
    Rcpp::Named("par")      = par,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = grad,
    Rcpp::Named("status")   = status);
}


//This function computes log of integral of phi(x) from a to b using IS
//[[Rcpp::export]]
arma::vec flogintphi(const arma::vec& Mean,
                     List& lCa,
                     const int& nCa,
                     const arma::vec& a,
                     const arma::vec& b,
                     const int& sumn,
                     const int& S,
                     const arma::rowvec& simu){
  arma::vec out(sumn);
  for(int ca(0); ca < nCa; ++ ca){
    arma::uvec lCas = lCa[ca];
    arma::vec Means(Mean.elem(lCas));
    arma::mat tp(arma::repmat(Means, 1, S));
    if(b(ca) == R_PosInf){
      NumericVector tp1(wrap(Means));
      NumericVector ldPhis = Rcpp::pnorm5(tp1 - a(ca), 0, 1, true, true);
      out.elem(lCas)       = as<arma::vec>(ldPhis);
    } else {
      arma::mat logphis(-0.5*arma::square(tp.each_row() - (simu*(b(ca) - a(ca)) + a(ca))));
      laverexp(logphis, S);
      out.elem(lCas)  = laverexp(logphis, S) + log(b(ca) - a(ca)) - 0.5*log(2*acos(-1));
    }
  }
  return out;
}


// variance
//[[Rcpp::export]]
List fcovCDI(const arma::vec& theta,
             const arma::vec& Gamma2,
             const arma::mat& Gye,
             const arma::mat& X,
             const arma::uvec ixWi,//indeces of X without intercept
             List& G,
             List& lCa,
             const int& nCa,
             const arma::mat& igroup,
             const int& ngroup,
             const int& K,
             const arma::vec& n,
             const int sumn,
             const arma::umat& idelta,
             const arma::uvec& ndelta,
             const arma::vec& Rbar,
             const double& R,
             const int& S,
             const bool& ccov) {
  int nCl(nCa*nCa);
  arma::vec lambda(theta.head(nCl));
  arma::vec beta(theta.subvec(nCl, nCl + K - 1));
  arma::vec tdelta(theta.tail(sum(ndelta)) + 1e-323);
  arma::vec delta(fdelta(tdelta, lambda, idelta, ndelta, nCa));
  arma::vec ZtLambda(Gye*lambda + X*beta);
  
  if(ccov){
    arma::rowvec simu(S, arma::fill::randu);
    int Kz(K + nCl), nparms(Kz + sum(ndelta)), Km(Gamma2.n_elem), Kzm(Km + nCl);
    arma::mat Z(arma::join_rows(Gye, X));
    
    // I compute l(Phir - Phir+1) and lphi from r = 0
    arma::mat ldPhi; // l(Phir - Phir+1) from r = 0
    arma::mat lphi(R_NegInf*arma::ones<arma::vec>(sumn)); // lphi from r = 0
    {
      NumericVector ZtL(wrap(ZtLambda));
      NumericVector ldPhi0(Rcpp::pnorm5(ZtL, 0, 1, false, true));
      ldPhi                = as<arma::vec>(ldPhi0);
    }
    arma::mat lamat(nCa, 1, arma::fill::zeros); // ar from r = 0
    // arma::mat zcheck       = Gye;          // Gye - r from r = 0
    int Rmax  = 0;
    bool next = (Rmax < (R + 1));
    arma::vec a2(nCa, arma::fill::zeros), a1(nCa);
    
    // log of f and log(diff F1 F2)
    while(next) {
      ++ Rmax;
      arma::vec lfr(sumn);
      for(int ca(0); ca < nCa; ++ ca){
        arma::uvec lCas = lCa[ca];
        arma::vec ZtL(ZtLambda.elem(lCas));
        NumericVector ZtLs = wrap(ZtL);
        a1(ca)          = a2(ca);
        a2(ca)         += fgamma(delta, idelta(ca, 0), Rmax + 1, Rbar(ca), R);
        NumericVector lfrs(Rcpp::dnorm4(ZtLs - a1(ca), 0, 1, true));
        lfr.elem(lCas)  = as<arma::vec>(lfrs);
      }
      lamat = arma::join_rows(lamat, a1);
      lphi  = arma::join_rows(lphi, lfr);
      ldPhi = arma::join_rows(ldPhi, flogintphi(ZtLambda, lCa, nCa, a1, a2, sumn, S, simu));
      // zcheck  = arma::join_rows(zcheck, Gye - Rmax);
      next  = (((max(lfr) > -1000) || (Rmax <= (Rbar.max() + 2))) && (Rmax < (R + 1)));
    }
    lamat     = log(lamat);
    // log(0, 1, 2, ...)
    arma::rowvec tp0(log(arma::regspace<arma::rowvec>(0, Rmax)));
    
    arma::mat B(sumn, nparms), Sigma(nparms, nparms, arma::fill::zeros);
    arma::mat Omega(nparms, sumn, arma::fill::zeros), Qth(Kzm, nparms, arma::fill::zeros);
    arma::mat Qy(sumn, 1, arma::fill::zeros);
    arma::vec d(sumn);
    arma::uvec il(arma::regspace<arma::uvec>(0, nCa - 1));
    for(int ca(0); ca < nCa; ++ ca){
      arma::uvec lCas = lCa[ca];
      arma::vec ZtL(ZtLambda.elem(lCas)), lbeta(arma::join_cols(arma::zeros<arma::vec>(nCl), Gamma2));
      lbeta.elem(il) = lambda.elem(il);
      arma::mat Gyes(Gye.submat(lCas, il)), lphis(lphi.rows(lCas)), ldPhis(ldPhi.rows(lCas)), Xs(X.rows(lCas)), 
      Gyepo(Gyes + 1), tGyepo(2*Gyes + 1);
      int iida(idelta(ca, 0)), ida(Kz + iida);
      bool hdbar(R > Rbar(ca));
      
      {
        // compute B and d, d is diag(D) of the paper
        arma::mat tp1(n(ca), Rbar(ca) + 1 + hdbar), tp2(n(ca), Rbar(ca) + 1 + hdbar);
        for(int l(0); l <= Rbar(ca); ++ l){
          tp1.col(l) = exp(lsumexp(lphis.cols(l, Rmax - 1))); //for sum fr from r = l
          tp2.col(l) = exp(lsumexp(lphis.cols(l, Rmax - 1).each_row() + lamat.submat(ca, l, ca, Rmax - 1))); //for sum a_r f_r from r = l, 
        }
        if(hdbar){
          tp1.col(Rbar(ca) + 1) = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1))); //for sum (r - Rbar) fr from r = Rbar + 1
          tp2.col(Rbar(ca) + 1) = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + (lamat.submat(ca, Rbar(ca) + 1, ca, Rmax - 1) + tp0.subvec(1, Rmax - Rbar(ca) - 1)))); //for sum (r - Rbar) a_r phi_r from r = Rbar + 1,
        }
        arma::vec tp3 = exp(lsumexp(lphis.cols(1, Rmax - 1).each_row() + tp0.subvec(0, Rmax - 2))); // sum (r - 1)*phir from r = 1
        arma::vec tp4 = exp(lsumexp(lphis.cols(1, Rmax - 1).each_row() + (tp0.subvec(0, Rmax - 2) + lamat.submat(ca, 1, ca, Rmax - 1)))); // sum (r-1)*ar*phir from r = 1
        if(R == 1){
          tp2.col(1).zeros();
          tp3.zeros();
          tp4.zeros();
        }
        
        // compute d
        arma::vec ds(tp1.col(1)); //for sum fr from r = 1
        
        // Qtheta (used to estimate cov(marginal effects))
        arma::vec tp5(tp2.col(1) - ZtL%ds);
        Qth.cols(il)          += (lbeta*sum((Gyes.each_col()%tp5).each_col() + (ZtL%tp3  - tp4), 0));//missinf Eye(kz)*...
        Qth.cols(nCl, Kz - 1) += (lbeta*sum(Xs.each_col()%tp5, 0));
        Qy.rows(lCas)          = tp5;
        // Compute B
        // I start with B2
        arma::mat Bs(n(ca), nparms, arma::fill::zeros);
        for(int k(0); k < Rbar(ca) - 1 + hdbar; ++ k) {
          Bs.col(ida + k)    = -tp1.col(2 + k);
          Qth.col(ida + k)  += lbeta*sum(ZtL%tp1.col(2 + k) - tp2.col(2 + k), 0);
        }
        // add B1
        Bs.cols(il)          = (Gyes.each_col()%ds).each_col() - tp3;
        
        // add DZ
        Bs.cols(nCl, Kz - 1) = Xs.each_col()%ds;
        B.rows(lCas)         = Bs;
        d.elem(lCas)         = ds;
      }
      
      // Fill Sigma0 and Omega0
      arma::mat Omegas(nparms, n(ca), arma::fill::zeros);
      {//Need temporary variables
        arma::mat tpra(n(ca), Rbar(ca) + 1), tprb(n(ca), Rbar(ca) + 1), tprc(n(ca), Rbar(ca) + 1);
        arma::mat tpr1a(n(ca), Rbar(ca) + 1), tpr1b(n(ca), Rbar(ca) + 1), tpr1c(n(ca), Rbar(ca) + 1);
        arma::mat tprra(n(ca), Rbar(ca) + 1), tprrb(n(ca), Rbar(ca) + 1), tprrc(n(ca), Rbar(ca) + 1);
        for(int l(0); l <= Rbar(ca); ++ l){
          tpra.col(l)   = exp(lsumexp(2*lphis.cols(l, Rmax - 1) - ldPhis.cols(l, Rmax - 1))); //for sum fr^2/(Fr - Fr+1) from r = l
          tprb.col(l)   = exp(lsumexp(lphis.cols(l, Rmax - 1) + (lphis.cols(l, Rmax - 1).each_row() + tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r*fr^2/(Fr - Fr+1) from r = l
          tprc.col(l)   = exp(lsumexp(2*(lphis.cols(l, Rmax - 1).each_row() + tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r^2*fr^2/(Fr - Fr+1) from r = l
          tpr1a.col(l)  = exp(lsumexp(2*lphis.cols(1 + l, Rmax) - ldPhis.cols(l, Rmax - 1))); //for sum fr+1^2/(Fr - Fr+1) from r = l
          tpr1b.col(l)  = exp(lsumexp(lphis.cols(1 + l, Rmax) + (lphis.cols(1 + l, Rmax).each_row() + tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r*fr+1^2/(Fr - Fr+1) from r = l
          tpr1c.col(l)  = exp(lsumexp(2*(lphis.cols(1 + l, Rmax).each_row() + tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r^2*fr+1^2/(Fr - Fr+1) from r = l
          tprra.col(l)  = exp(lsumexp(lphis.cols(l, Rmax - 1) + lphis.cols(1 + l, Rmax) - ldPhis.cols(l, Rmax - 1))); //for sum fr * fr+1/(Fr - Fr+1) from r = l
          tprrb.col(l)  = exp(lsumexp(lphis.cols(l, Rmax - 1) + (lphis.cols(1 + l, Rmax).each_row() + tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r * fr * fr+1/(Fr - Fr+1) from r = l
          tprrc.col(l)  = exp(lsumexp(lphis.cols(l, Rmax - 1) + (lphis.cols(1 + l, Rmax).each_row() + 2*tp0.subvec(l, Rmax - 1)) - ldPhis.cols(l, Rmax - 1))); //for sum r^2 * fr * fr+1/(Fr - Fr+1) from r = l
        }
        if(R == Rbar(ca)){
          tpr1a.col(Rbar(ca)).zeros();
          tpr1b.col(Rbar(ca)).zeros();
          tpr1c.col(Rbar(ca)).zeros();
          tprra.col(Rbar(ca)).zeros();
          tprrb.col(Rbar(ca)).zeros();
          tprrc.col(Rbar(ca)).zeros();
        }
        if(R == 1){
          tpr1b.col(0).zeros();
          tpr1c.col(0).zeros();
          tprra.col(0).zeros();
          tprrb.col(0).zeros();
          tprrc.col(0).zeros();
        }
        
        arma::mat adr1tprr(n(ca), 2), adr1tpr1(n(ca), 2), adrtprr(n(ca), 2), adrtpr(n(ca), 2);
        arma::vec adr12tpr1, adrrtprr, adr2tpr;
        if(hdbar){
          adr1tprr.col(0) = exp(lsumexp(lphis.cols(Rbar(ca), Rmax - 2) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca), Rmax - 2))); //for sum (r + 1 - Rbar) fr fr+1/(Fr - Fr+1) from r = Rbar
          adr1tprr.col(1) = exp(lsumexp(lphis.cols(Rbar(ca), Rmax - 2) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + (tp0.subvec(1, Rmax - Rbar(ca) - 1) + tp0.subvec(Rbar(ca), Rmax - 2))) - ldPhis.cols(Rbar(ca), Rmax - 2))); //for sum (r + 1 - Rbar) r fr fr+1/(Fr - Fr+1) from r = Rbar
          adr1tpr1.col(0) = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca), Rmax - 2))); //for sum (r + 1 - Rbar) fr+1^2/(Fr - Fr+1) from r = Rbar
          adr1tpr1.col(1) = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + (tp0.subvec(1, Rmax - Rbar(ca) - 1) + tp0.subvec(Rbar(ca), Rmax - 2))) - ldPhis.cols(Rbar(ca), Rmax - 2))); //for sum (r + 1 - Rbar) r fr+1^2/(Fr - Fr+1) from r = Rbar
          adrtprr.col(0)  = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 2, Rmax).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum (r - Rbar) fr fr+1/(Fr - Fr+1) from r = Rbar + 1
          adrtprr.col(1)  = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 2, Rmax).each_row() + (tp0.subvec(1, Rmax - Rbar(ca) - 1) + tp0.subvec(Rbar(ca) + 1, Rmax - 1))) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum (r - Rbar) r fr fr+1/(Fr - Fr+1) from r = Rbar + 1
          adrtpr.col(0)   = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum (r - Rbar) fr^2/(Fr - Fr+1) from r = Rbar + 1
          adrtpr.col(1)   = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + (tp0.subvec(1, Rmax - Rbar(ca) - 1) + tp0.subvec(Rbar(ca) + 1, Rmax - 1))) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum (r - Rbar) r fr^2/(Fr - Fr+1) from r = Rbar + 1
          
          adr12tpr1       = exp(lsumexp(2*(lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca), Rmax - 2))); //for sum ((r + 1 - Rbar) fr+1)^2/(Fr - Fr+1) from r = Rbar
          adrrtprr        = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1) + (lphis.cols(Rbar(ca) + 2, Rmax).each_row() + (tp0.subvec(1, Rmax - Rbar(ca) - 1) + tp0.subvec(2, Rmax - Rbar(ca)))) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum ((r - Rbar) (r + 1 - Rbar) fr fr+1)/(Fr - Fr+1) from r = Rbar + 1
          adr2tpr         = exp(lsumexp(2*(lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1)) - ldPhis.cols(Rbar(ca) + 1, Rmax - 1))); //for sum ((r - Rbar) fr)^2/(Fr - Fr+1) from r = Rbar + 1
          
          if((Rbar(ca) + 1) == R){
            adrtprr.zeros();
            adrrtprr.zeros();
          }
        }
        {// A lambda lambda
          for(int k1(0); k1 < nCa; ++ k1){
            Sigma(il(k1), il(k1))   += sum(arma::square(Gyes.col(k1))%tpr1a.col(0) -2*Gyes.col(k1)%tpr1b.col(0) + tpr1c.col(0) +
              arma::square(Gyepo.col(k1))%tpra.col(0) - 2*Gyepo.col(k1)%tprb.col(0) + tprc.col(0) -
              (2*Gyes.col(k1)%Gyepo.col(k1))%tprra.col(0) + 2*tGyepo.col(k1)%tprrb.col(0) - 2*tprrc.col(0));
            for(int k2(0); k2 < k1; ++ k2){
              Sigma(il(k1), il(k2)) += sum(Gyes.col(k1)%Gyes.col(k2)%tpr1a.col(0) - (Gyes.col(k1) + Gyes.col(k2))%tpr1b.col(0) + tpr1c.col(0) +
                Gyepo.col(k1)%Gyepo.col(k2)%tpra.col(0) - (Gyepo.col(k1) + Gyepo.col(k2))%tprb.col(0) + tprc.col(0) -
                (2*Gyes.col(k1)%Gyes.col(k2) + Gyes.col(k1) + Gyes.col(k2))%tprra.col(0) + 2*(Gyepo.col(k1) + Gyes.col(k2))%tprrb.col(0) - 2*tprrc.col(0));
            }
          }
        }
        
        {// A lambda Gamma
          arma::mat tp  = ((Gyes.each_col()%tpr1a.col(0)).each_col() - tpr1b.col(0)) - ((tGyepo.each_col()%tprra.col(0)).each_col() - 2*tprrb.col(0)) + ((Gyepo.each_col()%tpra.col(0)).each_col() - tprb.col(0));
          Sigma.submat(nCl, ca*nCa, Kz - 1, (ca + 1)*nCa - 1) += Xs.t()*tp;
          Omegas.rows(il) = tp.t();
        }
        
        {// A lambda deltak
          for(int k(2); k <= Rbar(ca); ++ k){
            arma::rowvec tp = sum(((Gyepo.each_col()%tprra.col(k - 1)).each_col() - tprrb.col(k - 1)) - ((Gyes.each_col()%tpr1a.col(k - 1)).each_col() - tpr1b.col(k - 1)) +
              ((Gyes.each_col()%tprra.col(k)).each_col() - tprrb.col(k)) - ((Gyepo.each_col()%tpra.col(k)).each_col() - tprb.col(k)), 0);
            Sigma.submat(ida + k - 2, ca*nCa, ida + k - 2, (ca + 1)*nCa - 1)  += tp;
          }
        }
        
        if(hdbar){// A lambda deltabar
          arma::rowvec tp = sum(((Gyepo.each_col()%adr1tprr.col(0)).each_col() - adr1tprr.col(1)) - ((Gyes.each_col()%adr1tpr1.col(0)).each_col() - adr1tpr1.col(1)) +
            ((Gyes.each_col()%adrtprr.col(0)).each_col() - adrtprr.col(1)) - ((Gyepo.each_col()%adrtpr.col(0)).each_col() - adrtpr.col(1)), 0);
          Sigma.submat(ida + Rbar(ca) - 1, ca*nCa, ida + Rbar(ca) - 1, (ca + 1)*nCa - 1)  += tp;
        }
        
        {// A Gamma Gamma
          arma::vec tp             = tpr1a.col(0) -2*tprra.col(0) + tpra.col(0);
          arma::mat Xtp            = arma::trans(Xs.each_col() % tp);
          Sigma.submat(nCl, nCl, Kz - 1, Kz - 1) += Xtp*Xs;
          Omegas.rows(nCl, Kz - 1)                = Xtp;
        }
        
        {
          for(int k(2); k <= Rbar(ca); ++ k){
            // A Gamma deltak
            arma::rowvec tp = (tprra.col(k - 1) - tpr1a.col(k - 1) + tprra.col(k) - tpra.col(k)).t();
            Sigma.submat(ida + k - 2, nCl, ida + k - 2, Kz - 1) += tp*Xs;
            Omegas.row(ida + k - 2) = tp;
            // A deltak deltak
            Sigma(ida + k - 2, ida + k - 2) += sum(tpr1a.col(k - 1) -2*tprra.col(k) + tpra.col(k));
            // A deltak deltal l < k
            for(int l(2); l < k; ++ l){
              Sigma(ida + k - 2, ida + l - 2) += (-sum(tp));
            }
          }
        }
        
        if(hdbar){
          // A Gamma deltabar
          arma::rowvec tp = (adr1tprr.col(0) - adr1tpr1.col(0) + adrtprr.col(0) - adrtpr.col(0)).t();
          Sigma.submat(ida + Rbar(ca) - 1, nCl, ida + Rbar(ca) - 1, Kz - 1) += tp*Xs;
          Omegas.row(ida + Rbar(ca) - 1) = tp;
          // A deltak deltabar
          for(int k(2); k <= Rbar(ca); ++ k){
            Sigma(ida + Rbar(ca) - 1, ida + k - 2) += (-sum(tp));
          }
        }
        
        if(hdbar){// A deltabar deltabar
          Sigma(ida + Rbar(ca) - 1, ida + Rbar(ca) - 1) += sum(adr12tpr1 - 2*adrrtprr + adr2tpr);
        }
      }
      Omega.cols(lCas) = Omegas;
      il              += nCa;
    }
    
    // Fill upper triangular Sigma
    for(int k(0); k < (nparms - 1); ++ k){
      Sigma.submat(k, k + 1, k, nparms - 1) = arma::trans(Sigma.submat(k + 1, k, nparms - 1, k));
    }
    
    // Compute G*inv(S)*B
    arma::mat GinvSB(sumn, nparms, arma::fill::zeros);
    for (int m(0); m < ngroup; ++ m) {
      int n1(igroup(m,0));
      int n2(igroup(m,1));
      int nm(n2 - n1 + 1);
      List Gm = G[m];
      arma::mat Sm(nm, nm, arma::fill::zeros), SlG(nm, nm, arma::fill::zeros);
      for(int cl(0); cl < nCl; ++ cl){
        arma::mat Gms = Gm[cl];
        SlG          += lambda(cl)*Gms;
        Sm           += lambda(cl)*(Gms.each_col() % d.subvec(n1, n2));
      }
      Sm                  = arma::eye<arma::mat>(nm, nm) - Sm;
      GinvSB.rows(n1, n2) = SlG*arma::solve(Sm, B.rows(n1, n2));
    }
    
    // Covariances
    Sigma         /= sumn;
    Omega          = Omega*GinvSB/sumn;
    arma::mat covt(arma::inv(Sigma + Omega));
    covt           = covt*Sigma*covt.t()/sumn;
    
    // Marginal effects
    arma::mat meff1(sumn, nCl, arma::fill::zeros);
    arma::uvec tp1(arma::regspace<arma::uvec>(0, nCa - 1));
    arma::mat tp2(nCl, nparms, arma::fill::eye); 
    arma::mat QyGiSB(Kzm, nparms, arma::fill::zeros);
    QyGiSB.tail_rows(Km) = Gamma2*Qy.t()*GinvSB;
    for(int ca(0); ca < nCa; ++ ca){
      arma::uvec lCas  = lCa[ca];
      meff1.submat(lCas, tp1) = d.elem(lCas)*lambda.elem(tp1).t();
      tp2.submat(tp1, tp1)   *= sum(d.elem(lCas));
      QyGiSB.rows(tp1)       += lambda.elem(tp1)*Qy.elem(lCas).t()*GinvSB.rows(lCas);
      tp1 += nCa;
    }
    arma::mat meff2(d*Gamma2.t());
    arma::mat meffects(arma::join_rows(meff1, meff2));
    
    arma::mat tp3(Kzm - nCl, nparms, arma::fill::zeros); 
    tp3.cols(nCl + ixWi) = arma::eye(Kzm - nCl, Kzm - nCl)*sum(d);
    Qth                 += arma::join_cols(tp2, tp3);
    arma::mat Q((Qth + QyGiSB)/sumn);
    arma::mat covm(Q*covt*Q.t());
    return List::create(Named("Rmax")  = Rmax,
                        Named("imeff") = meffects, 
                        Named("covt")  = covt,
                        Named("covm")  = covm);
  }
  return fmeffects(ZtLambda, lambda, Gamma2, lCa, nCa, delta, idelta, sumn, Rbar, R);
}

// Compute the derivative of the expected outcome wrt theta
//[[Rcpp::export]]
arma::mat fcddEy(const arma::vec& theta,
                 const arma::mat& Gye,
                 const arma::mat& X,
                 const arma::vec& psi,
                 List& G,
                 List& lCa,
                 const int& nCa,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const int& K,
                 const arma::vec& n,
                 const int sumn,
                 const arma::umat& idelta,
                 const arma::uvec& ndelta,
                 const arma::vec& Rbar,
                 const double& R,
                 const int& S) {
  int nCl(nCa*nCa);
  arma::vec lambda(theta.head(nCl));
  arma::vec beta(theta.subvec(nCl, nCl + K - 1));
  arma::vec tdelta(theta.tail(sum(ndelta)) + 1e-323);
  arma::vec delta(fdelta(tdelta, lambda, idelta, ndelta, nCa));
  arma::vec ZtLambda(Gye*lambda + psi);
  
  arma::rowvec simu(S, arma::fill::randu);
  int Kz(K + nCl), nparms(Kz + sum(ndelta));
  
  // I compute l(Phir - Phir+1) and lphi from r = 0
  arma::mat ldPhi; // l(Phir - Phir+1) from r = 0
  arma::mat lphi(R_NegInf*arma::ones<arma::vec>(sumn)); // lphi from r = 0
  {
    NumericVector ZtL(wrap(ZtLambda));
    NumericVector ldPhi0(Rcpp::pnorm5(ZtL, 0, 1, false, true));
    ldPhi                = as<arma::vec>(ldPhi0);
  }
  // arma::mat zcheck       = Gye;          // Gye - r from r = 0
  int Rmax  = 0;
  bool next = (Rmax < (R + 1));
  arma::vec a2(nCa, arma::fill::zeros), a1(nCa);
  
  // log of f and log(diff F1 F2)
  while(next) {
    ++ Rmax;
    arma::vec lfr(sumn);
    for(int ca(0); ca < nCa; ++ ca){
      arma::uvec lCas = lCa[ca];
      arma::vec ZtL(ZtLambda.elem(lCas));
      NumericVector ZtLs = wrap(ZtL);
      a1(ca)          = a2(ca);
      a2(ca)         += fgamma(delta, idelta(ca, 0), Rmax + 1, Rbar(ca), R);
      NumericVector lfrs(Rcpp::dnorm4(ZtLs - a1(ca), 0, 1, true));
      lfr.elem(lCas)  = as<arma::vec>(lfrs);
    }
    lphi  = arma::join_rows(lphi, lfr);
    ldPhi = arma::join_rows(ldPhi, flogintphi(ZtLambda, lCa, nCa, a1, a2, sumn, S, simu));
    // zcheck  = arma::join_rows(zcheck, Gye - Rmax);
    next  = (((max(lfr) > -1000) || (Rmax <= (Rbar.max() + 2))) && (Rmax < (R + 1)));
  }
  // log(0, 1, 2, ...)
  arma::rowvec tp0(log(arma::regspace<arma::rowvec>(0, Rmax)));
  
  arma::mat B(sumn, nparms);
  arma::vec d(sumn);
  arma::uvec il(arma::regspace<arma::uvec>(0, nCa - 1));
  for(int ca(0); ca < nCa; ++ ca){
    arma::uvec lCas = lCa[ca];
    arma::vec ZtL(ZtLambda.elem(lCas));
    arma::mat Gyes(Gye.submat(lCas, il)), lphis(lphi.rows(lCas)), ldPhis(ldPhi.rows(lCas)), Xs(X.rows(lCas));
    int iida(idelta(ca, 0)), ida(Kz + iida);
    bool hdbar(R > Rbar(ca));
    
    {
      // compute B and d, d is diag(D) of the paper
      arma::mat tp1(n(ca), Rbar(ca) + 1 + hdbar);
      for(int l(0); l <= Rbar(ca); ++ l){
        tp1.col(l) = exp(lsumexp(lphis.cols(l, Rmax - 1))); //for sum fr from r = l
      }
      if(hdbar){
        tp1.col(Rbar(ca) + 1) = exp(lsumexp(lphis.cols(Rbar(ca) + 1, Rmax - 1).each_row() + tp0.subvec(1, Rmax - Rbar(ca) - 1))); //for sum (r - Rbar) fr from r = Rbar + 1
      }
      arma::vec tp3 = exp(lsumexp(lphis.cols(1, Rmax - 1).each_row() + tp0.subvec(0, Rmax - 2))); // sum (r - 1)*phir from r = 1
      if(R == 1){
        tp3.zeros();
      }
      
      // compute d
      arma::vec ds(tp1.col(1)); //for sum fr from r = 1
      
      // Compute B
      // I start with B2
      arma::mat Bs(n(ca), nparms, arma::fill::zeros);
      for(int k(0); k < Rbar(ca) - 1 + hdbar; ++ k) {
        Bs.col(ida + k)    = -tp1.col(2 + k);
      }
      // add B1
      Bs.cols(il)          = (Gyes.each_col()%ds).each_col() - tp3;
      
      // add DZ
      Bs.cols(nCl, Kz - 1) = Xs.each_col()%ds;
      B.rows(lCas)         = Bs;
      d.elem(lCas)         = ds;
    }
    il                    += nCa;
  }
  
  // Compute G*inv(S)*B
  arma::mat invSB(sumn, nparms, arma::fill::zeros);
  for (int m(0); m < ngroup; ++ m) {
    int n1(igroup(m,0));
    int n2(igroup(m,1));
    int nm(n2 - n1 + 1);
    List Gm = G[m];
    arma::mat Sm(nm, nm, arma::fill::zeros);
    for(int cl(0); cl < nCl; ++ cl){
      arma::mat Gms = Gm[cl];
      Sm           += lambda(cl)*(Gms.each_col() % d.subvec(n1, n2));
    }
    Sm                  = arma::eye<arma::mat>(nm, nm) - Sm;
    invSB.rows(n1, n2)  = arma::solve(Sm, B.rows(n1, n2));
  }
  
  return invSB;
}