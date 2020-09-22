// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
arma::vec fyTobit(arma::vec& yst,
                  arma::vec& y,
                  arma::vec& Gy,
                  List& G,
                  const arma::vec& eps,
                  const arma::mat& igroup,
                  const int& ngroup,
                  const arma::vec& psi,
                  const double& lambda,
                  const double& tol,
                  const int& maxit) {
  int tm, nm;
  arma::vec ystm, Gym, ym, tmp0, tmp1, t(ngroup); 
  //loop over group
  for (int m(0); m < ngroup; ++ m) {
    tm            = 0;
    nm            = igroup(m,1) - igroup(m,0) + 1;
    tmp0          = arma::zeros(nm);
    ym            = y.subvec(igroup(m,0), igroup(m,1));
    arma::mat Gm  = G[m];
    newstep: Gym  = Gm*ym;
    
    ystm          = lambda*Gym + psi.subvec(igroup(m,0), igroup(m,1)) + eps.subvec(igroup(m,0), igroup(m,1));
    
    tmp1          = arma::max(arma::join_rows(tmp0, ystm), 1);
    
    double tmp2   = arma::accu(arma::abs(tmp1 - ym));
    
    ++tm;
    ym            = tmp1;
    if (tmp2 > tol && tm < maxit) goto newstep;
    
    y.rows(igroup(m,0), igroup(m,1))    = ym;
    Gy.rows(igroup(m,0), igroup(m,1))   = Gm * ym;
    
    t(m)                                = tm;
  }
  return t;  // the number of iteration if needed
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
                   const arma::vec& Npos,
                   const int& ngroup,
                   List& I,
                   List& W,
                   const int& N, 
                   const arma::vec& indzero,
                   const arma::vec& indpos,
                   const arma::mat igroup){
  NumericVector thetacpp = wrap(theta);
  thetacpp.attr("dim") = R_NilValue;
  Rcpp::Rcout<<"---------------"<< endl;
  Rcpp::Rcout<<"Estimate: "<< endl;
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1.0);
  double sigma         = exp(theta(K + 1));
  arma::vec tmp        = alpha*Gy + xb;
  arma::vec tmpzer     = tmp.elem(idzero);
  NumericVector tzcpp  = wrap(tmpzer);
  tmp                  = y - tmp;
  thetacpp(0)          = alpha;
  thetacpp(K + 1)      = sigma;
  Rcpp::print(thetacpp);


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
    sum(Npos)*(0.5*log(2*acos(-1)) + log(sigma)) + logdetA2(0) - 0.5*sum(pow(tmp.elem(idpos)/sigma, 2));
  Rcpp::Rcout <<"Likelihood: "<< llh << endl;
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}


// [[Rcpp::export]]
double foptimTobit0(const arma::vec& theta,
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
                    const arma::vec& Npos,
                    const int& ngroup,
                    List& I,
                    List& W,
                    const int& N, 
                    const arma::vec& indzero,
                    const arma::vec& indpos,
                    const arma::mat igroup){
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1.0);
  double sigma         = exp(theta(K + 1)) + 1e-8;
  arma::vec tmp        = alpha*Gy + xb;
  arma::vec tmpzer     = tmp.elem(idzero);
  NumericVector tzcpp  = wrap(tmpzer);
  tmp                  = y - tmp;
  
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
    sum(Npos)*(0.5*log(2*acos(-1)) + log(sigma)) + logdetA2(0) - 0.5*sum(pow(tmp.elem(idpos)/sigma, 2));
  
  if(llh < -1e293) {
    llh           = -1e293;
  }
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
                        const arma::vec& Npos,
                        const int& ngroup,
                        List& I,
                        List& W,
                        const int& N, 
                        const arma::vec& indzero,
                        const arma::vec& indpos,
                        const arma::mat igroup) {
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = 1.0/(exp(-theta(0)) + 1.0);
  
  double sigma         = exp(theta(K + 1)) + 1e-8;
  arma::vec tmp        = y - alpha*Gy - xb;
  NumericVector tmpcpp = wrap(tmp);
  
  NumericVector irmcpp = Rcpp::dnorm(tmpcpp/sigma, 0, 1, false)/Rcpp::pnorm(tmpcpp/sigma, 0, 1, true, false);
  arma::vec     irm    = as<arma::vec>(irmcpp);
  
  arma::vec rvec(N);
  
  for (int i(0); i < ngroup; ++ i) {
    arma::mat Wi   = W[i];
    arma::mat Ii   = I[i];
    arma::mat Rmat = arma::inv(Ii - alpha*Wi)*Wi;
    rvec.subvec(igroup(i,0), igroup(i,1)) = arma::diagvec(Rmat);
  }
  
  arma::mat out(N, K + 2);
  
  // derivation with respect alpha
  out.col(0)       = ((indzero%(irm)/sigma - indpos%tmp/pow(sigma, 2))%Gy + rvec)*exp(theta(0))/pow(1 + exp(theta(0)), 2);
  out.cols(1, K)   = arma::repmat(indzero%irm/sigma - indpos%tmp/pow(sigma,2), 1, K)%X;
  out.col(K + 1)   = (-indzero%irm%(alpha*Gy + xb)/pow(sigma, 2) + indpos%(1 - pow(tmp/sigma, 2))/sigma)*sigma;
  
  NumericVector outcpp = wrap(arma::trans(sum(out, 0)));
  outcpp.attr("dim") = R_NilValue;
  
  return outcpp;
}


//[[Rcpp::export]]
arma::mat fqTobit(const arma::vec& theta,
                  const arma::mat& X,
                  List& G2,
                  List& I,
                  List& W,
                  const int& K,
                  const int& N,
                  const arma::vec& y,
                  const arma::vec& Gy,
                  const arma::vec& indzero,
                  const arma::vec& indpos,
                  const arma::mat igroup,
                  const int& ngroup) {
  arma::vec xb         = X * theta.subvec(1, K);
  double alpha         = theta(0);
  double sigma         = theta(K + 1);
  arma::vec tmp        = y - alpha*Gy - xb;
  NumericVector tmpcpp = wrap(tmp);
  
  NumericVector irmcpp = Rcpp::dnorm(tmpcpp/sigma, 0, 1, false)/Rcpp::pnorm(tmpcpp/sigma, 0, 1, true, false);
  arma::vec     irm    = as<arma::vec>(irmcpp);
  
  arma::vec rvec(N);
  
  for (int i(0); i < ngroup; ++ i) {
    arma::mat Wi   = W[i];
    arma::mat Ii   = I[i];
    arma::mat Rmat = arma::inv(Ii - alpha*Wi)*Wi;
    rvec.subvec(igroup(i,0), igroup(i,1)) = arma::diagvec(Rmat);
  }
  
  arma::mat out(N, K + 2);
  
  // derivation with respect alpha
  out.col(0)       = ((-indzero%(irm)/sigma + indpos%tmp/pow(sigma, 2))%Gy - rvec);
  out.cols(1, K)   = arma::repmat(-indzero%irm/sigma + indpos%tmp/pow(sigma, 2), 1, K)%X;
  out.col(K + 1)   = (indzero%irm%(alpha*Gy + xb)/pow(sigma, 2) - indpos%(1 - pow(tmp/sigma, 2))/sigma);
  
  return out;
}


//------------------------------ SAR -----------------------------
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
                 const int& N,
                 const arma::vec& y,
                 const arma::vec& Gy,
                 const int& ngroup){
  double alpha         = 1.0/(1.0 + exp(-alphatilde));
  Rcpp::Rcout<<"---------------"<< endl;
  Rcpp::Rcout<<"Estimate:   "<< alpha << endl;
  arma::vec yrem       = y - alpha*Gy;
  arma::vec beta       = invXX*(X.t()*yrem);
  arma::vec e          = yrem - X*beta;
  double s2            = sum(e%e)/N;
  
  double logdetA(0);
  for (int i(0); i < ngroup; ++ i) {
    double vali, signi;
    arma::mat Gi  = G[i];
    arma::mat Ii  = I[i];
    log_det(vali, signi, Ii - alpha*Gi);
    logdetA      += vali;
    logdetA      += log(signi);
  }
  
  double llh      = - 0.5*N*(log(2*acos(-1)*s2)) + logdetA - 0.5*N;
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
                  const int& N,
                  const arma::vec& y,
                  const arma::vec& Gy,
                  const int& ngroup){
  double alpha         = 1.0/(1.0 + exp(-alphatilde));
  arma::vec yrem       = y - alpha*Gy;
  arma::vec beta       = invXX*(X.t()*yrem);
  arma::vec e          = yrem - X*beta;
  double s2            = sum(e%e)/N;
  
  double logdetA(0);
  for (int i(0); i < ngroup; ++ i) {
    double vali, signi;
    arma::mat Gi  = G[i];
    arma::mat Ii  = I[i];
    log_det(vali, signi, Ii - alpha*Gi);
    logdetA      += vali;
    logdetA      += log(signi);
  }
  
  double llh      = - 0.5*N*(log(2*acos(-1)*s2)) + logdetA - 0.5*N;
  if(llh < -1e293) {
    llh           = -1e293;
  }
  return -llh;
}



// [[Rcpp::export]]
arma::mat fSARjac(const double& alpha,
                  const double& s2,
                  const arma::mat& X,
                  const arma::mat& XX,
                  const arma::vec& Xbeta,
                  List& G,
                  List& I,
                  const arma::mat igroup,
                  const int& ngroup,
                  const int& N,
                  const int& K){
  
  arma::vec GXbeta(N);
  double trGsG(0);
  double trG(0);
  
  for(int m(0); m < ngroup; ++m) {
    arma::mat Wm    = G[m];
    arma::mat tWm   = Wm.t();
    arma::mat Im    = I[m];
    arma::mat tGm   = arma::solve(Im - alpha*tWm, tWm);
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
  out(K + 1, K + 1) = N/(2*pow(s2, 2));
  
  return arma::inv(out);
}

