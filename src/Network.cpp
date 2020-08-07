// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;

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
//[[Rcpp::export]]
arma::vec fmusum(const arma::vec& mu,
                 const arma::mat& index,
                 const arma::mat& indexgr,
                 const int& M,
                 const int& N) {
  int igr1, igr2, nm, j = 0, j1, j2;
  double  mui;
  arma::vec mum, mumj; 
  
  arma::vec musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0); 
    igr2 = indexgr(m, 1);
    nm   = igr2 - igr1 + 1;
    mum  = mu.subvec(igr1, igr2);
    
    for(int i(0); i < nm; ++ i){
      mui    = mum(i);
      mumj   = mum;
      mumj.shed_row(i);
      j1     = index(j, 0);
      j2     = index(j, 1);
      musum.subvec(j1, j2) = mui + mumj;
      ++j;
    }
  }
  return musum; 
}

void updatebeta1(arma::vec& beta,
                 arma::vec& dxbeta,
                 arma::vec& dxbpmu,
                 double& llh,
                 const arma::vec A,
                 const arma::mat& dX,
                 const arma::vec& musum,
                 const arma::vec& jsbeta,
                 arma::vec& betaaccept,
                 const int& K){
  
  arma::vec betast, dxbetast; 
  double llhst, lalpha2;
  arma::vec dxbpmust;
  
  for(int k(0); k < K; ++ k) {
    betast           = beta;
    betast(k)        = R::rnorm(betast(k), jsbeta(k));
    dxbetast         = dX*betast;
    dxbpmust         = dxbetast + musum;
    llhst            = arma::sum(A%(dxbpmust) - log(1 + exp(dxbpmust)));
    lalpha2          = llhst - llh;
    if(unif_rand() < exp(lalpha2)){
      beta(k)        = betast(k);    
      dxbeta         = dxbetast; 
      dxbpmu         = dxbpmust;
      llh            = llhst;
      betaaccept(k) +=1.0;
    }
  }
}



void updatebeta2(arma::vec& beta,
                 arma::vec& dxbeta,
                 arma::vec& dxbpmu,
                 double& llh,
                 const arma::vec A,
                 const arma::mat& dX,
                 const arma::vec& musum,
                 const arma::mat& jscovbeta,
                 double& betaaccept,
                 const int& K){
  
  arma::vec betast   = Fmvnorm(K, beta, jscovbeta);
  arma::vec dxbetast = dX*betast;
  arma::vec dxbpmust = dxbetast + musum;
  double llhst       = arma::sum(A%(dxbpmust) - log(1 + exp(dxbpmust)));
  double lalpha2     = llhst - llh;
  if(unif_rand() < exp(lalpha2)){
    beta             = betast; 
    dxbeta           = dxbetast; 
    dxbpmu           = dxbpmust;
    llh              = llhst;
    betaaccept      += 1.0;
  }
}


void updatemu(arma::vec& mu,
              const arma::vec A,
              const arma::vec& dxbeta,
              const arma::vec& sigmau2,
              const arma::vec& uu,
              const arma::uvec& possigma,
              const int& M,
              const arma::vec& nvec,
              const arma::mat& index,
              const arma::mat& indexgr,
              const arma::vec& jsmu,
              arma::vec& muaccept){
  int i1, i2, igr1, igr2, nm, j = 0;
  double muist, mui, lalpha2;
  arma::vec tmp, tmp1, mum, mumj, muveci, jsmum; 
  arma::uvec indexsup;
  arma::mat indexm;
  
  for (int m(0); m < M; ++ m) {
    igr1                    = indexgr(m, 0); 
    igr2                    = indexgr(m, 1);
    
    nm                      = nvec(m);
    mum                     = mu.subvec(igr1, igr2);
    indexm                  = index.rows(igr1, igr2); 
    
    for (int i(0); i < nm; ++ i) { 
      i1                    = indexm(i, 0); 
      i2                    = indexm(i, 1);
      indexsup              = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
      indexsup.head(i + 1) -= 1;
      indexsup.shed_row(i);
      muveci                = mum;
      muveci.shed_row(i);
      mui                   = mum(i);
      muist                 = R::rnorm(mui, jsmu(j));
      
      tmp                   = dxbeta.subvec(i1, i2) + muveci;
      tmp1                  = dxbeta.elem(indexsup) + muveci;
      
      lalpha2               = arma::sum((muist - mui)*A.subvec(i1, i2) - log(1 + exp(tmp + muist)) +  log(1 + exp(tmp + mui))) +
        arma::sum((muist - mui)*A.elem(indexsup) - log(1 + exp(tmp1  + muist)) +  log(1 + exp(tmp1  + mui))) -
        0.5*(pow(muist - uu(possigma(j)), 2) - pow(mui - uu(possigma(j)), 2))/sigmau2(possigma(j));
      if(unif_rand() < exp(lalpha2)){
        mum(i)              = muist;    
        muaccept(j)        += 1.0;
      }
      ++ j;
    }
    mu.subvec(igr1, igr2)   = mum;
  }
}


void updateusigmamu2(arma::vec& sigmau2,
                     arma::vec& uu,
                    const arma::vec& mu, 
                    const double& Msigma,
                    const arma::mat& indexsigma) {
  for (int m(0); m < Msigma; ++ m){
    int n1                 = indexsigma(m, 0); 
    int n2                 = indexsigma(m, 1); 

    int ahat               = n2 - n1 + 1;
    arma::vec mum          = mu.subvec(n1, n2);
    double mmum            = mean(mum);
    uu(m)                  = R::rnorm(mmum, sqrt(sigmau2(m)/ahat));
      
    double bhat            = arma::accu(pow(mum - uu(m), 2));
    sigmau2(m)             = 1./R::rgamma(0.5 * (ahat - 1), 2/bhat);
  }
}



//[[Rcpp::export]]
List updategparms1(const arma::vec A,
                   const arma::mat& dX,
                   const arma::vec& beta0,
                   const arma::vec& mu0,
                   const arma::vec& sigmau20,
                   const arma::vec& uu0,
                   arma::vec& jsbeta,
                   arma::vec& jsmu,
                   const arma::mat& index,
                   const arma::mat& indexgr,
                   const arma::mat& indexsigma,
                   const arma::uvec& possigma,
                   const int& N,
                   const int& M,
                   const int& K,
                   const int& Msigma,
                   arma::vec& nvec, 
                   const int& iteration1,
                   const int& iteration2,
                   const double tbeta,
                   const double tmu){
  arma::vec beta       = beta0;
  arma::vec mu         = mu0;
  arma::vec sigmau2    = sigmau20;
  arma::vec uu         = uu0; 
  
  int iteration        = iteration1 + iteration2;
  int n                = sum(nvec);
  arma::vec musum      = fmusum(mu, index, indexgr, M, N);
  arma::vec dxbeta     = dX*beta;
  arma::vec dxbpmu     = dxbeta + musum;
  double llh           = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
  arma::vec betaaccept = arma::zeros(K);
  double    betaallac  = 0;
  arma::vec muaccept   = arma::zeros(n);
  arma::mat smu(n, iteration);
  arma::mat sbeta(K, iteration);
  arma::mat suu(Msigma, iteration);
  arma::mat ssigmamu2(Msigma, iteration);
  NumericVector spdis(iteration);
  spdis.attr("dim") = R_NilValue;
  
  double jsbetatmp, jsmutmp, jsbetaall = 1;
  Progress p(iteration, true);
  
  for(int t(0); t < iteration1; ++ t) {
    p.increment(); 
    updatebeta1(beta, dxbeta, dxbpmu, llh, A, dX, musum, jsbeta, betaaccept, K);
    
    updatemu(mu, A, dxbeta, sigmau2, uu, possigma, M, nvec, index, indexgr, jsmu, muaccept);
    musum     = fmusum(mu, index, indexgr, M, N);
    dxbpmu    = dxbeta + musum;
    llh       = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
    
    updateusigmamu2(sigmau2, uu, mu, Msigma, indexsigma);
    
    for(int k(0); k < K; ++ k) {
      jsbetatmp       = jsbeta[k] + (betaaccept[k]/t - tbeta)/pow(t,0.6);
      if((jsbetatmp > 1e-8) & (jsbetatmp < 100)){jsbeta[k] = jsbetatmp;}
    }
    for(int i(0); i < n; ++ i) {
      jsmutmp         = jsmu[i] + (muaccept[i]/t - tmu)/pow(t,0.6);
      if((jsmutmp > 1e-8) & (jsmutmp < 100)){jsmu[i] = jsmutmp;}
    }
    
    
    smu.col(t)        = mu;
    sbeta.col(t)      = beta;
    suu.col(t)        = uu;
    ssigmamu2.col(t)  = sigmau2;
    spdis(t)          = llh - arma::sum(pow(mu - uu.elem(possigma), 2)/sigmau2.elem(possigma)) - 
      arma::sum(0.5*(indexsigma.col(1) - indexsigma.col(0) + 2)%log(sigmau2));
  }
  
  betaallac           = arma::sum(betaaccept)/K;
  arma::mat jscovbeta = arma::cov(arma::trans(sbeta.cols(0.5*iteration1, iteration1 - 1)));
  for(int t(iteration1); t < iteration; ++ t) {
    p.increment(); 
    updatebeta2(beta, dxbeta, dxbpmu, llh, A, dX, musum, pow(jsbetaall, 2)*jscovbeta, betaallac, K);
    
    updatemu(mu, A, dxbeta, sigmau2, uu, possigma, M, nvec, index, indexgr, jsmu, muaccept);
    musum     = fmusum(mu, index, indexgr, M, N);
    dxbpmu    = dxbeta + musum;
    llh       = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
    
    updateusigmamu2(sigmau2, uu, mu, Msigma, indexsigma);
    
    double jsbetatmp  = jsbetaall + (betaallac/t - tbeta)/pow(t,0.6);
    if((jsbetatmp > 1e-8) & (jsbetatmp < 100)){jsbetaall = jsbetatmp;}
    for(int i(0); i < n; ++ i) {
      double jsmutmp  = jsmu[i] + (muaccept[i]/t - tmu)/pow(t,0.6);
      if((jsmutmp > 1e-8) & (jsmutmp < 100)){jsmu[i] = jsmutmp;}
    }
    
    
    smu.col(t)        = mu;
    sbeta.col(t)      = beta;
    suu.col(t)        = uu;
    ssigmamu2.col(t)  = sigmau2;
    spdis(t)          = llh - arma::sum(pow(mu - uu.elem(possigma), 2)/sigmau2.elem(possigma)) - 
      arma::sum(0.5*(indexsigma.col(1) - indexsigma.col(0) + 2)%log(sigmau2));
  }
  
  List posterior      = List::create(Named("beta")        = sbeta.t(),
                                     Named("mu")          = smu.t(),
                                     Named("uu")          = suu.t(),
                                     Named("sigmamu2")    = ssigmamu2.t(),
                                     Named("log.density") = spdis);
  List accept         = List::create(Named("beta")        = betaallac/iteration,
                                     Named("mu")          = muaccept/iteration);
  return List::create(Named("posterior")                  = posterior, 
                      Named("acceptance.rate")            = accept);
  
}


//[[Rcpp::export]]
List updategparms2(const arma::vec A,
                   const arma::mat& dX,
                   const arma::vec& beta0,
                   const arma::vec& mu0,
                   const arma::vec& sigmau20,
                   const arma::vec& uu0,
                   arma::vec& jsbeta,
                   arma::vec& jsmu,
                   const arma::mat& index,
                   const arma::mat& indexgr,
                   const arma::mat& indexsigma,
                   const arma::uvec& possigma, 
                   const int& N,
                   const int& M,
                   const int& K,
                   const int& Msigma,
                   arma::vec& nvec, 
                   const int& iteration1,
                   const int& iteration2,
                   const double tbeta,
                   const double tmu){
  arma::vec beta       = beta0;
  arma::vec mu         = mu0;
  arma::vec sigmau2    = sigmau20; 
  arma::vec uu         = uu0; 
  
  int iteration        = iteration1 + iteration2;
  int n                = sum(nvec);
  arma::vec musum      = fmusum(mu, index, indexgr, M, N);
  arma::vec dxbeta     = dX*beta;
  arma::vec dxbpmu     = dxbeta + musum;
  double llh           = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
  arma::vec betaaccept = arma::zeros(K);
  double    betaallac  = 0;
  arma::vec muaccept   = arma::zeros(n);
  arma::mat smu(n, iteration);
  arma::mat sbeta(K, iteration);
  arma::mat suu(Msigma, iteration);
  arma::mat ssigmamu2(Msigma, iteration);
  NumericVector spdis(iteration);
  spdis.attr("dim") = R_NilValue;
  
  double jsbetatmp, jsmutmp, jsbetaall = 1;
  
  for(int t(0); t < iteration1; ++ t) {
    updatebeta1(beta, dxbeta, dxbpmu, llh, A, dX, musum, jsbeta, betaaccept, K);
    
    updatemu(mu, A, dxbeta, sigmau2, uu, possigma, M, nvec, index, indexgr, jsmu, muaccept);
    musum     = fmusum(mu, index, indexgr, M, N);
    dxbpmu    = dxbeta + musum;
    llh       = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
    
    updateusigmamu2(sigmau2, uu, mu, Msigma, indexsigma);
    
    for(int k(0); k < K; ++ k) {
      jsbetatmp       = jsbeta[k] + (betaaccept[k]/t - tbeta)/pow(t,0.6);
      if((jsbetatmp > 1e-8) & (jsbetatmp < 100)){jsbeta[k] = jsbetatmp;}
    }
    for(int i(0); i < n; ++ i) {
      jsmutmp         = jsmu[i] + (muaccept[i]/t - tmu)/pow(t,0.6);
      if((jsmutmp > 1e-8) & (jsmutmp < 100)){jsmu[i] = jsmutmp;}
    }
    
    
    smu.col(t)        = mu;
    sbeta.col(t)      = beta;
    suu.col(t)        = uu;
    ssigmamu2.col(t)  = sigmau2;
    spdis(t)          = llh - arma::sum(pow(mu - uu.elem(possigma), 2)/sigmau2.elem(possigma)) - 
      arma::sum(0.5*(indexsigma.col(1) - indexsigma.col(0) + 2)%log(sigmau2));
  }
  
  betaallac           = arma::sum(betaaccept)/K;
  arma::mat jscovbeta = arma::cov(arma::trans(sbeta.cols(0.5*iteration1, iteration1 - 1)));

  for(int t(iteration1); t < iteration; ++ t) {
    updatebeta2(beta, dxbeta, dxbpmu, llh, A, dX, musum, pow(jsbetaall, 2)*jscovbeta, betaallac, K);
    
    updatemu(mu, A, dxbeta, sigmau2, uu, possigma, M, nvec, index, indexgr, jsmu, muaccept);
    musum     = fmusum(mu, index, indexgr, M, N);
    dxbpmu    = dxbeta + musum;
    llh       = arma::sum(A%(dxbpmu) - log(1 + exp(dxbpmu)));
    
    updateusigmamu2(sigmau2, uu, mu, Msigma, indexsigma);
    
    double jsbetatmp  = jsbetaall + (betaallac/t - tbeta)/pow(t,0.6);
    if((jsbetatmp > 1e-8) & (jsbetatmp < 100)){jsbetaall = jsbetatmp;}
    for(int i(0); i < n; ++ i) {
      double jsmutmp  = jsmu[i] + (muaccept[i]/t - tmu)/pow(t,0.6);
      if((jsmutmp > 1e-8) & (jsmutmp < 100)){jsmu[i] = jsmutmp;}
    }
    
    smu.col(t)        = mu;
    sbeta.col(t)      = beta;
    suu.col(t)        = uu;
    ssigmamu2.col(t)  = sigmau2;
    spdis(t)          = llh - arma::sum(pow(mu - uu.elem(possigma), 2)/sigmau2.elem(possigma)) - 
      arma::sum(0.5*(indexsigma.col(1) - indexsigma.col(0) + 2)%log(sigmau2));
  }
  
  List posterior      = List::create(Named("beta")        = sbeta.t(),
                                     Named("mu")          = smu.t(),
                                     Named("uu")          = suu.t(),
                                     Named("sigmamu2")    = ssigmamu2.t(),
                                     Named("log.density") = spdis);
  List accept         = List::create(Named("beta")        = betaallac/iteration,
                                     Named("mu")          = muaccept/iteration);
  return List::create(Named("posterior")                  = posterior, 
                      Named("acceptance.rate")            = accept);
  
}
