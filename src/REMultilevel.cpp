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
 *           igroup[s,] gives the firt and the last rows number for the group s
 * ngroup  : is the number of groups.
 * R       : is an integer value such that y can only take R + 1 values.
 * a       : is the vector of boundaries a_0, a_1, a_2, ..., a_(R+1).
 * w       : is the vector of y possibilities.
 * A       : is the matrix of boundaries. For the same reason that I removed the first colum
 *           of G, I remove the fist column of A. Moreother, the matrix A1 and A2 from t
 *           paper are fitted together in A. A is then an identical row matrix where each row
 *           is a_1, a_2, ..., a_(R+1).
 * theta   : is the vector of parameters ordered as follow: peer effects, explanatory variables
 *           coefficients, sigma (not sigma^2).
 * N       : The sample size.
 * tol     : A tolerance value for the iterative method to compute P convergence. 
 * maxit   : The maximum number of iterations of the iterative method to compute P. If this
 *           number is reached, the algorithm stops and the last P is used as the solution. maxit
 *           is important for numerical reasons if tol is too small. For example a tol = 1e-16
 *           may not be reached and the algorithm will stop after maxit iterations.
 * P       : is the matrix of probabilities. P(i,r) is the probability of y(i) = r
 *           For numerical reasons, P does not contain the colunm of probability of y = 0
 *           P is then N * R matrix and rumrow of P is  <= 0. 
 *           The zero probabilities can be computed by 1 - sumrow of P.
 * GPw     : A machine vector variable to save G*P*w 
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

/* LOCAL NOTATIONS
 *  xb     : is X*beta vector
 *  alpha  : is the peer effect
 *  sigma  : is sigma as in the paper
 */

// H: computes Pnew from a previous P.
//[[Rcpp::export]]
arma::mat H(List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::mat& P,
            const arma::vec& w,
            arma::vec& GPw,
            const arma::vec xb,
            const arma::mat& A,
            const double& alpha,
            const double& sigma,
            const int& N,
            const int& R) {
  for (int i(0); i < ngroup; ++ i) {
    arma::mat Gi                         = G[i];
    GPw.subvec(igroup(i,0), igroup(i,1)) = Gi*(P.rows(igroup(i,0), igroup(i,1))*w.tail(R));
  }
  arma::vec tmp   =  alpha*GPw + xb;
  
  arma::mat out(N, R + 1);
  for (int i(0); i < N; ++ i) {
    for (int r(0); r <= R; ++r) {
      out(i, r) = R::pnorm(tmp(i) + A(i,r), 0, sigma, true, false);
    }
  }
  return out.head_cols(R) - out.tail_cols(R);
}


// updateP: Takes an initial value of P and finds the equilibrium
//[[Rcpp::export]]
int updateP(arma::mat& P,
            List& G,
            const arma::mat& igroup,
            const int& ngroup,
            const arma::vec& w,
            arma::vec& GPw,
            const arma::vec xb,
            const arma::mat& A,
            const double& alpha,
            const double& sigma,
            const int& N,
            const int& R,
            const double& tol = 1e-13,
            const int& maxit = 1e3) {
  //Rcpp::Rcout<<"****** start updating P ******"<<endl;
  int t = 0;
  newP: arma::mat tmp1 = H(G, igroup, ngroup, P, w, GPw, xb, A, alpha, sigma, N, R);
  arma::mat       tmp2 = arma::abs(tmp1 - P);
  //arma::rowvec    tmp4 = arma::pow(arma::sum(arma::pow(tmp3, 2), 0), 0.5);
  double          tmp3 = arma::accu(tmp2);
  ++t;
  P                    = tmp1;
  //cout<<t<<" "<<tmp4.max()<<endl;
  if (tmp3 > tol && t < maxit) goto newP;
  //if (tmp4 > tol) goto newP;
  //cout<<t<<endl;
  //if (tmp3.max() > tol && t < 1000) goto newP;
  //Rcpp::Rcout<<"          iteration: "<<t<<endl;
  //Rcpp::Rcout<<"****** stop  updating P ******"<<endl;
  return t;
}

// flopn: compute numerical log.
//        lopn(x) = log(1e-323) if x <= 1e-323
//        logn(x) = log(x) if x > 1e-323
arma::mat flogn(const arma::mat& x) {
  //Rcpp::NumericVector out  = log(x);
  //Rcpp::LogicalVector corr = (x <= 1e-323);
  //out[corr]                = log(1e-323);
  return log(x);
}

// flogP: compute log of probabilities
//        Unlike P, log of probabilities includes zero probabilities
//        The ouput is then a N * (R+1) matrix
//[[Rcpp::export]]
arma::mat flogP(const arma::mat& P,
                const arma::vec& w,
                const arma::vec& GPw,
                const arma::vec xb,
                const arma::mat& A,
                const double& alpha,
                const double& sigma,
                const int& N,
                const int& R) {
  arma::vec tmp   =  alpha*GPw + xb;
  arma::mat out(N, R + 1);
  for (int i(0); i < N; ++ i) {
    for (int r(0); r <= R; ++r) {
      out(i, r) = R::pnorm(tmp(i) + A(i,r), 0, sigma, true, true);
    }
  }
  return out.head_cols(R) + log(1 - exp(out.tail_cols(R) - out.head_cols(R)));
  
  
}

// floglih: computes the likelihood given log of probabilities
//[[Rcpp::export]]
double floglih(const int& N,
               const arma::vec& y,
               const arma::vec& a,
               const arma::mat& logP) {
  double out = 0;
  int    g;
  for (int i(0); i < N; ++ i) {
    g    = sum(a < y(i)) - 1;
    out += logP(i,g);
  }
  return out;
}

// foptimREM: compute the likelihood given the patameters
//[[Rcpp::export]]
double foptimREM(arma::mat& P,
                 const arma::vec& theta,
                 const arma::mat& X,
                 List& G,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const arma::vec& w,
                 const arma::vec& a,
                 const arma::mat& A,
                 const int& K,
                 const int& N,
                 const int& R,
                 const arma::vec& y,
                 const double& tol = 1e-13,
                 const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  Rcpp::print(thetacpp);
  arma::vec xb    = X * theta.subvec(1, K);
  double alpha    = 1.0/(exp(-theta(0)) + 1);
  double sigma    = exp(theta(K + 1)) + 1e-8;
  arma::vec GPw(N);
  updateP(P, G, igroup, ngroup, w, GPw, xb, A, alpha, sigma, N, R, tol,  maxit);
  arma::mat logP  = flogP(P, w, GPw, xb, A, alpha, sigma, N, R);
  double llh      = floglih(N, y, a, logP);
  Rcpp::Rcout << -llh << endl;
  return -llh;
}

//' foptimREMvec: vectorized form of foptimREM
//'               returns the likelihood for each
//'               individuals given the patameters
//[[Rcpp::export]]
NumericVector foptimREMvec(arma::mat& P,
                           const arma::vec& theta,
                           const arma::mat& X,
                           List& G,
                           const arma::mat& igroup,
                           const int& ngroup,
                           const arma::vec& w,
                           const arma::vec& a,
                           const arma::mat& A,
                           const int& K,
                           const int& N,
                           const int& R,
                           const arma::vec& y,
                           const double& tol = 1e-13,
                           const int& maxit  = 1e3) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  //Rcpp::print(thetacpp);
  arma::vec xb    = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1);
  double sigma         = exp(theta(K + 1)) + 1e-8;
  arma::vec GPw(N);
  updateP(P, G, igroup, ngroup, w, GPw, xb, A, alpha, sigma, N, R, tol,  maxit);
  arma::mat logP  = flogP(P, w, GPw, xb, A, alpha, sigma, N, R);
  NumericVector outllh(N);
  int    g;
  for (int i(0); i < N; ++ i) {
    g      = sum(a < y(i)) - 1;
    outllh(i) = logP(i,g);
  }
  
  return outllh;
}

//' foptimREM_nopeer: computes the likelihood with the constraint alpha = 0
//'                   theta does not contain alpha here: just the explanatory
//'                   variables and sigma.
//[[Rcpp::export]]
double foptimREM_nopeer(const arma::vec& theta,
                        const arma::mat& X,
                        const arma::vec& w,
                        const arma::vec& a,
                        const arma::mat& A,
                        const int& K,
                        const int& N,
                        const int& R,
                        const arma::vec& y) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  Rcpp::print(thetacpp);
  arma::vec xb    = X * theta.subvec(0, K - 1);
  double sigma    = exp(theta(K)) + 1e-8;
  // log P
  arma::mat out(N, R + 1);
  for (int i(0); i < N; ++ i) {
    for (int r(0); r <= R; ++r) {
      out(i, r) = R::pnorm(xb(i) + A(i,r), 0, sigma, true, true);
    }
  }
  
  arma::mat logP = out.head_cols(R) + log(1 - exp(out.tail_cols(R) - out.head_cols(R)));
  
  double llh     = floglih(N, y, a, logP);
  Rcpp::Rcout << -llh << endl;
  return -llh;
}

//' foptimREM_nopeervec: the vectorized form of foptimREM_nopeer
//'                      returns the likelihood for each individuals
//'                       given the patameters
//[[Rcpp::export]]
NumericVector foptimREM_nopeervec(const arma::vec& theta,
                                  const arma::mat& X,
                                  const arma::vec& w,
                                  const arma::vec& a,
                                  const arma::mat& A,
                                  const int& K,
                                  const int& N,
                                  const int& R,
                                  const arma::vec& y) {
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  Rcpp::print(thetacpp);
  arma::vec xb    = X * theta.subvec(0, K - 1);
  double sigma    = exp(theta(K)) + 1e-8;
  // log P
  arma::mat out(N, R + 1);
  for (int i(0); i < N; ++ i) {
    for (int r(0); r <= R; ++r) {
      out(i, r) = R::pnorm(xb(i) + A(i,r), 0, sigma, true, true);
    }
  }
  
  arma::mat logP = out.head_cols(R) + log(1 - exp(out.tail_cols(R) - out.head_cols(R)));
  
  NumericVector outllh(N);
  int    g;
  for (int i(0); i < N; ++ i) {
    g      = sum(a < y(i)) - 1;
    outllh(i) = logP(i,g);
  }
  return outllh;
}

/*
 * Fmvnorm samples one vector fastly from numtivariate normal
 * dim   : is the vector dimension
 * u     : is the mean vector
 * sigma : is the covariance matrix
 */
arma::vec Fmvnorm(const double& dim, arma::vec u, arma::mat sigma) {
  arma::vec x = arma::randn(dim,1);
  return arma::chol(sigma).t()*x + u;
}

// END OF THE MACHINE FUNCTIONS 

/* NOTATIONS (CONTINUED)
 * yst        : is a vector of the latent variable y^*
 * XtX        : is a matrice to save X.t()*X
 * P0         : Initial value of P
 * theta0     : Initial values of theta
 * Prior distribution of (beta sigma) is NIV(as/2, bs/2, mub, inverse of(invVb))  
 * Prior distribution of the peer effect is N(mua, inverse of (invVa)) 
 * truncated in (alphamin, alphamax)
 * NOTE       : invVb is  the inverse of Vb and invVa is the inverse of Va
 * niteration : is the number of the Gibbs iterations
 */

// updateyst updates yst
void updateyst(arma::vec& yst,
               const arma::vec& y,
               const arma::vec& a,
               const arma::vec& GPw,
               const arma::vec& xb,
               const double& alpha,
               const double& sigma,
               const int& N,
               Function rtnormR){
  int    g;
  NumericVector ysttmp;
  arma::vec ztheta1 = alpha*GPw + xb;
  for (int i(0); i < N; ++ i) {
    g      = sum(a < y(i)) - 1;
    ysttmp = rtnormR(1, ztheta1(i),  sigma, a(g), a(g+1));
    yst(i) = ysttmp(0);
  }
}




//update alpha by MCMC
void updateaMH(double& alpha,
               double& zeta,
               double& jumpzeta,
               double& alphaaccept,
               arma::mat& P,
               arma::vec& GPw,
               List& G,
               const arma::mat& igroup,
               const int& ngroup,
               const arma::vec& w,
               const arma::mat& A,
               const arma::vec& yst,
               const arma::vec& xb,
               const double& sigma,
               const int& N,
               const int& R,
               const double& tol,
               const int& maxit) {
  double zetast     = R::rnorm(zeta, jumpzeta);
  double alphast    = 1/(1 + exp(-zetast));
  arma::mat Pst     = P;
  arma::vec GPwst   = GPw;
  
  // Declaration of some variable 
  double logalphazeta, logalpha2zeta;
  
  //Compute logalpha2
  updateP(Pst, G, igroup, ngroup, w, GPwst, xb, A, alphast, sigma, N, R, tol, maxit);
  arma::vec tmp     = yst - GPw*alpha - xb, tmpst = yst - GPwst*alphast - xb;
  
  logalpha2zeta     = 0.5*(arma::dot(tmp, tmp) - arma::dot(tmpst, tmpst))/pow(sigma, 2);
  logalpha2zeta    += zeta - zetast +2*log(1 + exp(-zeta)) -2*log(1 + exp(-zetast));
  
  //Compute logalpha
  logalphazeta      = min(Rcpp::NumericVector::create(0, logalpha2zeta));
  
  if(unif_rand()<exp(logalphazeta)){
    zeta            = zetast;
    alpha           = alphast;
    P               = Pst;
    GPw             = GPwst;
    alphaaccept    += 1;     
  }
}


//update sigma by MCMC
void updatesMH(double& sigma,
               double& lsigma2,
               double& jumplsigma2,
               double& sigmaaccept,
               arma::mat& P,
               arma::vec& GPw,
               List& G,
               const arma::mat& igroup,
               const int& ngroup,
               const arma::vec& w,
               const arma::mat& A,
               const arma::vec& yst,
               const arma::vec& xb,
               const double& alpha,
               const arma::vec& beta,
               const arma::vec& mub,
               const arma::mat& invVb,
               const double& as,
               const double& bs,
               const int& N,
               const int& K,
               const int& R,
               const double& tol,
               const int& maxit) {
  double sigma2       = pow(sigma, 2);
  arma::mat Pst       = P;
  arma::vec GPwst     = GPw;
  
  double lsigma2st    = R::rnorm(lsigma2, jumplsigma2);
  double sigma2st     = exp(lsigma2st);
  double sigmast      = sqrt(sigma2st);
  
  // Declaration of some variable 
  double logalphalsigma2, logalpha2lsigma2;
  
  //Compute logalpha2
  updateP(Pst, G, igroup, ngroup, w, GPwst, xb, A, alpha, sigmast, N, R, tol, maxit);
  arma::vec tmp       = yst - GPw*alpha - xb, tmpst = yst - GPwst*alpha - xb;
  
  logalpha2lsigma2    = 0.5*(bs + arma::dot(tmp, tmp) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2;
  logalpha2lsigma2   -= 0.5*(bs + arma::dot(tmpst, tmpst) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2st;
  logalpha2lsigma2   += 0.5*(N + K + as)*(lsigma2 - lsigma2st);
  
  //Compute logalpha
  logalphalsigma2     = min(Rcpp::NumericVector::create(0, logalpha2lsigma2));
  
  if(unif_rand()<exp(logalphalsigma2)){
    sigma             = sigmast;
    lsigma2           = lsigma2st,
      P                 = Pst;
    GPw               = GPwst;
    sigmaaccept      += 1;     
  }
}


//update sigma by MCMC
void updatebMH(arma::vec& beta,
               NumericVector& jumpbeta,
               NumericVector& betaaccept,
               arma::mat& P,
               arma::vec& GPw,
               arma::vec& xb,
               List& G,
               const arma::mat& igroup,
               const int& ngroup,
               const arma::vec& w,
               const arma::mat& A,
               const arma::vec& yst,
               const arma::mat& X,
               const double& alpha,
               const arma::vec& mub,
               const arma::mat& invVb,
               const double& sigma,
               const int& N,
               const int& K,
               const int& R,
               const double& tol,
               const int& maxit) {
  // Declaration of some variable 
  double logalphabeta, logalpha2beta;
  arma::vec betast, xbst, GPwst;
  arma::mat Pst;
  
  for (int k(0); k < K; ++ k) {
    Pst               = P;
    GPwst             = GPw;
    betast            = beta;
    betast(k)         = R::rnorm(beta(k), jumpbeta(k));
    xbst              = X*betast;
    
    //Compute logalpha2
    updateP(Pst, G, igroup, ngroup, w, GPwst, xbst, A, alpha, sigma, N, R, tol, maxit);
    arma::vec tmp     = yst - GPw*alpha - xb, tmpst = yst - GPwst*alpha - xbst;
    
    logalpha2beta     = arma::dot(tmp, tmp) - arma::dot(tmpst, tmpst);
    logalpha2beta    += arma::dot(beta - mub, invVb*(beta - mub));
    logalpha2beta    -= arma::dot(betast - mub, invVb*(betast - mub));
    logalpha2beta    *= (0.5/pow(sigma, 2));
    
    //Compute logalpha
    logalphabeta      = min(Rcpp::NumericVector::create(0, logalpha2beta));
    if(unif_rand()<exp(logalphabeta)){
      beta            = betast;
      xb              = xbst;
      P               = Pst;
      GPw             = GPwst;
      betaaccept(k)  += 1;     
    }
  }
}


void updatebsMH(arma::vec& beta,
                double& lsigma2,
                double& sigma,
                double& jump,
                double& accept,
                arma::mat& P,
                arma::vec& GPw,
                arma::vec& xb,
                const arma::mat& covmcmc,
                List& G,
                const arma::mat& igroup,
                const int& ngroup,
                const arma::vec& w,
                const arma::mat& A,
                const arma::vec& yst,
                const arma::mat& X,
                const double& alpha,
                const arma::vec& mub,
                const arma::mat& invVb,
                const double& as,
                const double& bs,
                const int& N,
                const int& K,
                const int& R,
                const double& tol,
                const int& maxit) {
  // Declaration of some variable 
  double logalpha, logalpha2, lsigma2st, sigmast, sigma2st;
  double sigma2 = exp(lsigma2);
  arma::vec betast, xbst, GPwst = GPw, thetatmp(K + 1);
  arma::mat Pst =  P;
  
  thetatmp.head(K)     = beta;
  thetatmp(K)          = lsigma2;
  
  
  arma::vec thetatmpst = Fmvnorm(K + 1, thetatmp, pow(jump, 2)*covmcmc);
  betast               = thetatmpst.head(K);
  lsigma2st            = thetatmpst(K);
  sigma2st             = exp(lsigma2st);
  sigmast              = pow(sigma2st, 0.5);
  xbst                 = X*betast;
  
  //Compute logalpha2
  updateP(Pst, G, igroup, ngroup, w, GPwst, xbst, A, alpha, sigmast, N, R, tol, maxit);
  arma::vec tmp        = yst - GPw*alpha - xb, tmpst = yst - GPwst*alpha - xbst;
  
  logalpha2            = 0.5*(bs + arma::dot(tmp, tmp) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2;
  logalpha2           -= 0.5*(bs + arma::dot(tmpst, tmpst) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2st;
  logalpha2           += 0.5*(N + K + as)*(lsigma2 - lsigma2st);
  
  //Compute logalpha
  logalpha            = min(Rcpp::NumericVector::create(0, logalpha2));
  if(unif_rand()<exp(logalpha)){
    beta              = betast;
    lsigma2           = lsigma2st;
    sigma             = sigmast;
    xb                = xbst;
    P                 = Pst;
    GPw               = GPwst;
    accept            += 1;     
  }
}



void updateallMH(double& alpha,
                 double& zeta,
                 arma::vec& beta,
                 double& lsigma2,
                 double& sigma,
                 double& jump,
                 double& accept,
                 arma::mat& P,
                 arma::vec& GPw,
                 arma::vec& xb,
                 const arma::mat& covmcmc,
                 List& G,
                 const arma::mat& igroup,
                 const int& ngroup,
                 const arma::vec& w,
                 const arma::mat& A,
                 const arma::vec& yst,
                 const arma::mat& X,
                 const arma::vec& mub,
                 const arma::mat& invVb,
                 const double& as,
                 const double& bs,
                 const int& N,
                 const int& K,
                 const int& R,
                 const double& tol,
                 const int& maxit) {
  // Declaration of some variable 
  double logalpha, logalpha2, zetast, alphast, lsigma2st, sigmast, sigma2st;
  double sigma2 = exp(lsigma2);
  arma::vec betast, xbst, GPwst = GPw, thetatmp(K + 2);
  arma::mat Pst =  P;
  
  thetatmp(0)          = zeta;
  thetatmp.subvec(1, K)= beta;
  thetatmp(K + 1)      = lsigma2;
  
  
  arma::vec thetatmpst = Fmvnorm(K + 2, thetatmp, pow(jump, 2)*covmcmc);
  zetast               = thetatmpst(0);
  alphast              = 1/(exp(-zetast) + 1.0);
  betast               = thetatmpst.subvec(1, K);
  lsigma2st            = thetatmpst(K + 1);
  sigma2st             = exp(lsigma2st);
  sigmast              = pow(sigma2st, 0.5);
  xbst                 = X*betast;
  
  //Compute logalpha2
  updateP(Pst, G, igroup, ngroup, w, GPwst, xbst, A, alphast, sigmast, N, R, tol, maxit);
  arma::vec tmp        = yst - GPw*alpha - xb, tmpst = yst - GPwst*alphast - xbst;
  
  logalpha2            = 0.5*(bs + arma::dot(tmp, tmp) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2;
  logalpha2           -= 0.5*(bs + arma::dot(tmpst, tmpst) + 
    arma::dot(beta - mub, invVb*(beta - mub)))/sigma2st;
  logalpha2           += 0.5*(N + K + as)*(lsigma2 - lsigma2st);
  logalpha2           += zeta - zetast +2*log(1 + exp(-zeta)) -2*log(1 + exp(-zetast));
  
  //Compute logalpha
  logalpha            = min(Rcpp::NumericVector::create(0, logalpha2));
  if(unif_rand()<exp(logalpha)){
    alpha             = alphast;
    zeta              = zetast;
    beta              = betast;
    lsigma2           = lsigma2st;
    sigma             = sigmast;
    xb                = xbst;
    P                 = Pst;
    GPw               = GPwst;
    accept            += 1;     
  }
}


// burnin     : Burn-in iteration
// nitcov     : number of iterations to compute the covariance. This willbe ignored if the covariance is given
// niteration : Number of iteration of MCMC

// Without covariance
//[[Rcpp::export]]
arma::mat bREMCMC0(const arma::mat& P0,
                   const arma::vec& theta0,
                   const arma::mat& X,
                   List& G,
                   const arma::mat& igroup,
                   const arma::vec& w,
                   const arma::vec& a,
                   const arma::mat& A,
                   const int& K,
                   const int& N,
                   const int& R,
                   const arma::vec& y,
                   const arma::vec& mub,
                   const arma::mat& invVb,
                   const double& as,
                   const double& bs,
                   const double& tol,
                   const int& maxit, 
                   const int& burnin,
                   const int& nitcov,
                   const int& niteration,
                   const double& target,
                   const double& jumpmin,
                   const double& jumpmax,
                   const double& c,
                   Function rtnormR) {
  // variables declaration
  int iter;
  int burn1     = burnin;
  int burn2     = burnin + nitcov;
  int totsim    = burn2 + niteration;
  
  arma::vec theta(theta0), thetatmp(theta0), beta(theta.subvec(1, K));
  arma::vec yst(N), GPw(N), xb(N);
  double alpha(theta(0)), zeta(log(alpha/(1 - alpha)));
  double sigma(theta(K + 1)), jumpzeta(1), alphaaccept(0), jump(1), accept(0);
  double jumplsigma2(1), sigmaaccept(0), lsigma2(log(pow(theta(K + 1),2)));
  arma::mat P(P0), out(totsim + 1, K + 2), outtmp(nitcov, K + 2);
  NumericVector thetacpp, jumpbeta = wrap(arma::ones(K)), betaaccept(K);
  
  xb                   = X * beta;  // AH Error
  
  int ngroup           = G.size();
  updateP(P, G, igroup, ngroup, w, GPw, xb, A, alpha, sigma, N, R, tol,  maxit);
  arma::rowvec thetam;     // temp
  NumericVector thetamcpp; // temp
  
  thetatmp(0)          = zeta;
  thetatmp(K + 1)      = lsigma2;
  out.row(0)           = theta.t();
  outtmp.row(0)        = thetatmp.t();
  
  
  // FIRST PHASE ** BURN-IN ***************************************
  for (iter = 1; iter <= burn1; ++ iter) {
    Rcpp::Rcout<<"Burn-in-1    : "<<iter<<"/"<<burn1<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateaMH(alpha, zeta, jumpzeta, alphaaccept, P, GPw, G, igroup, ngroup, w, A, yst,
              xb, sigma, N, R, tol, maxit);
    updatebMH(beta, jumpbeta, betaaccept, P, GPw, xb, G, igroup, ngroup, w, A, yst, X,
              alpha, mub, invVb, sigma, N, K, R, tol, maxit);
    updatesMH(sigma, lsigma2, jumplsigma2, sigmaaccept, P, GPw, G, igroup, ngroup, w, A, yst, xb,
              alpha, beta, mub, invVb, as, bs, N, K, R, tol, maxit);
    
    theta.subvec(1, K)    = beta; 
    theta(0)              = alpha; 
    theta(K + 1)          = sigma; 
    out.row(iter)         = theta.t(); 
    
    // update jumping scale
    double jumpzetast = jumpzeta + (alphaaccept/iter - target)/pow(iter, c);
    if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
    for (int k(0); k < K; ++ k) {
      double jumpbetast = jumpbeta(k) + (betaaccept(k)/iter - target)/pow(iter, c);
      if((jumpbetast > jumpmin) & (jumpbetast < jumpmax)){jumpbeta(k) = jumpbetast;}
    }
    double jumplsigma2st = jumplsigma2 + (sigmaaccept/iter - target)/pow(iter, c);
    if((jumplsigma2st > jumpmin) & (jumplsigma2st < jumpmax)){jumplsigma2 = jumplsigma2st;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
  }
  
  // SECOND STAGE ** COMPUTE COVARIANCE ******************
  for (iter = burn1 + 1; iter <= burn2; ++ iter) {
    Rcpp::Rcout<<"Burn-in-2    : "<<iter - burn1<<"/"<<nitcov<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateaMH(alpha, zeta, jumpzeta, alphaaccept, P, GPw, G, igroup, ngroup, w, A, yst,
              xb, sigma, N, R, tol, maxit);
    updatebMH(beta, jumpbeta, betaaccept, P, GPw, xb, G, igroup, ngroup, w, A, yst, X,
              alpha, mub, invVb, sigma, N, K, R, tol, maxit);
    updatesMH(sigma, lsigma2, jumplsigma2, sigmaaccept, P, GPw, G, igroup, ngroup, w, A, yst, xb,
              alpha, beta, mub, invVb, as, bs, N, K, R, tol, maxit);
    
    theta.subvec(1, K)    = beta; thetatmp.subvec(1, K)             = beta;
    theta(0)              = alpha; thetatmp(0)                      = zeta;
    theta(K + 1)          = sigma; thetatmp(K + 1)                  = lsigma2;
    out.row(iter)         = theta.t(); outtmp.row(iter - burn1 - 1) = thetatmp.t();
    
    // update jumping scale
    double jumpzetast = jumpzeta + (alphaaccept/iter - target)/pow(iter, c);
    if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
    for (int k(0); k < K; ++ k) {
      double jumpbetast = jumpbeta(k) + (betaaccept(k)/iter - target)/pow(iter, c);
      if((jumpbetast > jumpmin) & (jumpbetast < jumpmax)){jumpbeta(k) = jumpbetast;}
    }
    double jumplsigma2st = jumplsigma2 + (sigmaaccept/iter - target)/pow(iter, c);
    if((jumplsigma2st > jumpmin) & (jumplsigma2st < jumpmax)){jumplsigma2 = jumplsigma2st;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
  }
  
  // FIRD STAGE ** MV MCMC ******************
  arma::mat covmcmc = arma::cov(outtmp);
  accept           = (alphaaccept + sum(betaaccept) + sigmaaccept)/(K + 2.0);
  
  for (iter = burn2 + 1; iter <= totsim; ++ iter) {
    Rcpp::Rcout<<"Iteration    : "<<iter - burn2<<"/"<<totsim - burn2<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateallMH(alpha, zeta, beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc,
                G, igroup, ngroup, w, A, yst, X, mub, invVb, as, bs, N, K, R, tol,
                maxit);
    
    theta.subvec(1, K) = beta;
    theta(0)           = alpha;
    theta(K + 1)       = sigma;
    out.row(iter)      = theta.t();
    
    // update jumping scale
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
    
  }
  
  // acceptance rate
  betaaccept             = betaaccept / burn2;
  betaaccept.attr("dim") = R_NilValue;
  Rcpp::Rcout<<"Acceptance Burn-in ** "<<endl;
  Rcpp::Rcout<<"                   alpha: "<< alphaaccept/burn2<<endl;
  Rcpp::Rcout<<"                   beta : "<< endl;
  Rcpp::print(betaaccept);
  Rcpp::Rcout<<"                   sigma: "<< sigmaaccept/burn2<<endl;
  Rcpp::Rcout<<"         --------         "<<endl;
  
  Rcpp::Rcout<<"Acceptance Overall ** "<< accept/(double)(totsim)<<endl;
  return out;
}

// Covariance for b and s only
//[[Rcpp::export]]
arma::mat bREMCMC1(const arma::mat& P0,
                   const arma::vec& theta0,
                   const arma::mat& X,
                   List& G,
                   const arma::mat& igroup,
                   const arma::vec& w,
                   const arma::vec& a,
                   const arma::mat& A,
                   const int& K,
                   const int& N,
                   const int& R,
                   const arma::vec& y,
                   const arma::vec& mub,
                   const arma::mat& invVb,
                   const double& as,
                   const double& bs,
                   const double& tol,
                   const int& maxit, 
                   const int& burnin,
                   const int& nitcov,
                   const int& niteration,
                   const arma::mat& covmcmc,
                   const double& target,
                   const double& jumpmin,
                   const double& jumpmax,
                   const double& c,
                   Function rtnormR) {
  // variables declaration
  int iter, burn1, burn2, totsim;
  burn1     = burnin;
  burn2     = burnin + nitcov;
  totsim    = burn2 + niteration;
  
  arma::vec theta(theta0), thetatmp(theta0), beta(theta.subvec(1, K));
  arma::vec yst(N), GPw(N), xb(N);
  double alpha(theta(0)), zeta(log(alpha/(1 - alpha)));
  double sigma(theta(K + 1)), jumpzeta(1), alphaaccept(0), jump(1), accept(0);
  double lsigma2(log(pow(theta(K + 1),2)));
  arma::mat P(P0), out(totsim + 1, K + 2), outtmp(nitcov, K + 2);
  NumericVector thetacpp;
  
  xb                   = X * beta;  // AH Error
  
  int ngroup           = G.size();
  updateP(P, G, igroup, ngroup, w, GPw, xb, A, alpha, sigma, N, R, tol,  maxit);
  arma::rowvec thetam;     // temp
  NumericVector thetamcpp; // temp
  
  thetatmp(0)          = zeta;
  thetatmp(K + 1)      = lsigma2;
  out.row(0)           = theta.t();
  outtmp.row(0)        = thetatmp.t();
  
  
  // FIRST PHASE ** BURN-IN ***************************************
  for (iter = 1; iter <= burn1; ++ iter) {
    Rcpp::Rcout<<"Burn-in-1    : "<<iter<<"/"<<burn1<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateaMH(alpha, zeta, jumpzeta, alphaaccept, P, GPw, G, igroup, ngroup, w, A, yst,
              xb, sigma, N, R, tol, maxit);
    
    updatebsMH(beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc, G, igroup, ngroup,
               w, A, yst, X, alpha, mub, invVb, as, bs, N, K, R, tol, maxit);
    
    theta.subvec(1, K)    = beta; 
    theta(0)              = alpha; 
    theta(K + 1)          = sigma; 
    out.row(iter)         = theta.t(); 
    
    // update jumping scale
    double jumpzetast = jumpzeta + (alphaaccept/iter - target)/pow(iter, c);
    if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
    
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
  }
  
  // SECOND STAGE ** COMPUTE COVARIANCE ******************
  for (iter = burn1 + 1; iter <= burn2; ++ iter) {
    Rcpp::Rcout<<"Burn-in-2    : "<<iter - burn1<<"/"<<nitcov<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateaMH(alpha, zeta, jumpzeta, alphaaccept, P, GPw, G, igroup, ngroup, w, A, yst,
              xb, sigma, N, R, tol, maxit);
    
    updatebsMH(beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc, G, igroup, ngroup,
               w, A, yst, X, alpha, mub, invVb, as, bs, N, K, R, tol, maxit);
    
    theta.subvec(1, K)    = beta; thetatmp.subvec(1, K)             = beta;
    theta(0)              = alpha; thetatmp(0)                      = zeta;
    theta(K + 1)          = sigma; thetatmp(K + 1)                  = lsigma2;
    out.row(iter)         = theta.t(); outtmp.row(iter - burn1 - 1) = thetatmp.t();
    
    // update jumping scale
    double jumpzetast = jumpzeta + (alphaaccept/iter - target)/pow(iter, c);
    if((jumpzetast > jumpmin) & (jumpzetast < jumpmax)){jumpzeta = jumpzetast;}
    
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
  }
  
  
  // FIRD STAGE ** MV MCMC ******************
  arma::mat covmcmcnew = arma::cov(outtmp);
  accept = 0.5*(accept + alphaaccept);
  
  for (iter = burn2 + 1; iter <= totsim; ++ iter) {
    Rcpp::Rcout<<"Iteration    : "<<iter - burn2<<"/"<<totsim - burn2<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateallMH(alpha, zeta, beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc,
                G, igroup, ngroup, w, A, yst, X, mub, invVb, as, bs, N, K, R, tol,
                maxit);
    
    theta.subvec(1, K) = beta;
    theta(0)           = alpha;
    theta(K + 1)       = sigma;
    out.row(iter)      = theta.t();
    
    // update jumping scale
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
    
  }
  
  // acceptance rate
  Rcpp::Rcout<<"Acceptance Overall ** "<< accept/(double)(totsim)<<endl;
  return out;
}



//[[Rcpp::export]]
arma::mat bREMCMC2(const arma::mat& P0,
                   const arma::vec& theta0,
                   const arma::mat& X,
                   List& G,
                   const arma::mat& igroup,
                   const arma::vec& w,
                   const arma::vec& a,
                   const arma::mat& A,
                   const int& K,
                   const int& N,
                   const int& R,
                   const arma::vec& y,
                   const arma::vec& mub,
                   const arma::mat& invVb,
                   const double& as,
                   const double& bs,
                   const double& tol,
                   const int& maxit, 
                   const int& burnin,
                   const int& niteration,
                   const arma::mat& covmcmc,
                   const double& target,
                   const double& jumpmin,
                   const double& jumpmax,
                   const double& c,
                   Function rtnormR) {
  // variables declaration
  int iter;
  int burn1  = burnin;
  int totsim = burn1 + niteration;
  
  arma::vec theta(theta0), beta(theta.subvec(1, K));
  arma::vec yst(N), GPw(N), xb(N);
  double alpha(theta(0)), zeta(log(alpha/(1 - alpha)));
  double sigma(theta(K + 1)), jump(1), accept(0);
  double lsigma2(log(pow(theta(K + 1),2)));
  arma::mat P(P0), out(totsim + 1, K + 2);
  NumericVector thetacpp;
  
  xb                   = X * beta;  // AH Error
  
  int ngroup           = G.size();
  updateP(P, G, igroup, ngroup, w, GPw, xb, A, alpha, sigma, N, R, tol,  maxit);
  arma::rowvec thetam;     // temp
  NumericVector thetamcpp; // temp
  
  
  out.row(0)           = theta.t();
  
  
  // FIRST PHASE ** BURN-IN ***************************************
  for (iter = 1; iter <= burn1; ++ iter) {
    Rcpp::Rcout<<"Iteration    : "<<iter<<"/"<<burn1<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateallMH(alpha, zeta, beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc,
                G, igroup, ngroup, w, A, yst, X, mub, invVb, as, bs, N, K, R, tol,
                maxit);
    
    theta.subvec(1, K) = beta;
    theta(0)           = alpha;
    theta(K + 1)       = sigma;
    out.row(iter)      = theta.t();
    
    // update jumping scale
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
    
  }
  
  // FIRD STAGE ** MV MCMC ******************
  for (iter = burn1 + 1; iter <= totsim; ++ iter) {
    Rcpp::Rcout<<"Iteration    : "<<iter - burn1<<"/"<<totsim - burn1<<endl;
    thetacpp             = wrap(theta);
    thetacpp.attr("dim") = R_NilValue;
    Rcpp::print(thetacpp);
    
    updateyst(yst, y, a, GPw, xb, alpha, sigma, N, rtnormR);
    updateallMH(alpha, zeta, beta, lsigma2, sigma, jump, accept, P, GPw, xb, covmcmc,
                G, igroup, ngroup, w, A, yst, X, mub, invVb, as, bs, N, K, R, tol,
                maxit);
    
    theta.subvec(1, K) = beta;
    theta(0)           = alpha;
    theta(K + 1)       = sigma;
    out.row(iter)      = theta.t();
    
    // update jumping scale
    double jumpst = jump + (accept/iter - target)/pow(iter, c);
    if((jumpst > jumpmin) & (jumpst < jumpmax)){jump = jumpst;}
    
    // start temp
    if (iter < 1000) {
      thetam           = arma::mean(out.rows(1, iter), 0);
    } else {
      thetam           = arma::mean(out.rows(iter - 999, iter), 0);
    }
    thetamcpp             = wrap(thetam);
    thetamcpp.attr("dim") = R_NilValue;
    Rcpp::print(thetamcpp);
    // end temp
    Rcpp::Rcout<<"************************"<<endl;
    
  }
  
  // acceptance rate
  Rcpp::Rcout<<"Acceptance Overall ** "<< accept/(double)(totsim)<<endl;
  return out;
}




//[[Rcpp::export]]
arma::mat bREMCMC(const arma::mat& P0,
                  const arma::vec& theta0,
                  const arma::mat& X,
                  List& G,
                  const arma::mat& igroup,
                  const arma::vec& w,
                  const arma::vec& a,
                  const arma::mat& A,
                  const int& K,
                  const int& N,
                  const int& R,
                  const arma::vec& y,
                  const arma::vec& mub,
                  const arma::mat& invVb,
                  const double& as,
                  const double& bs,
                  const double& tol     = 1e-13,
                  const int& maxit      = 1e3, 
                  const int& burnin     = 5e2,
                  const int& nitcov     = 1e3,
                  const int& niteration = 2e3,
                  const SEXP& covMCMC   = R_NilValue,
                  const double& target  = 0.44,
                  const double& jumpmin = 1e-12,
                  const double& jumpmax = 1e2,
                  const double& c       = 0.6,
                  const int& type       = 0){
  
  // call rtnorm from R
  Rcpp::Environment base("package:MCMCglmm"); 
  Rcpp::Function rtnormR = base["rtnorm"];
  
  arma::mat covmcmc, out;
  if (covMCMC != R_NilValue) {
    covmcmc   = Rcpp::as<arma::mat>(covMCMC);
  }
  
  if (type == 0) {
    out  = bREMCMC0(P0, theta0, X, G, igroup, w, a, A, K, N, R, y, mub, invVb,
                    as, bs, tol, maxit, burnin, nitcov, niteration, target, 
                    jumpmin, jumpmax, c, rtnormR);
  }
  if (type == 1) {
    out  = bREMCMC1(P0, theta0, X, G, igroup, w, a, A, K, N, R, y, mub, invVb,
                    as, bs, tol, maxit, burnin, nitcov, niteration, covmcmc, target, 
                    jumpmin, jumpmax, c, rtnormR);
  }
  if (type == 2) {
    out  = bREMCMC2(P0, theta0, X, G, igroup, w, a, A, K, N, R, y, mub, invVb,
                    as, bs, tol, maxit, burnin, niteration, covmcmc, target, 
                    jumpmin, jumpmax, c, rtnormR);
  }
  
  return out;
  
}



