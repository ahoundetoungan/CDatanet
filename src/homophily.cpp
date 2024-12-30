// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppDist, RcppEigen, RcppNumerical)]]
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//#define NDEBUG
#include <progress.hpp>
#include <progress_bar.hpp>
#include <mvnorm.h>
#include <wishart.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// creates dummies
//[[Rcpp::export]]
void fdummies(arma::mat& out,
              const arma::mat& limit,
              const int& M,
              const int& n){
  for(int m(0); m < M; ++ m){
    out.submat(limit(m, 0), m, limit(m, 1), m) += 1;
  }
}

// this function organizes the data in the suitable shape
// We first compute intermediate function
// Use i characteristics
arma::mat fdatai(const arma::vec& Xk, const int& n){return arma::repmat(Xk, 1, n);}

// Use j characteristics
arma::mat fdataj(const arma::vec& Xk, const int& n){return arma::repmat(Xk.t(), n, 1);}

// sum
arma::mat fdatasum(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1) + Xk(i));
  }
  return out + out.t();
}

// prod
arma::mat fdataprod(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1)*Xk(i));
  }
  return out + out.t();
}

// Use Xi = Xj
arma::umat fdatasame(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = (Xk.subvec(i + 1, n - 1) == Xk(i));
  }
  return out + out.t();
}

// abs(Xi - Xj)
arma::mat fdatadiff(const arma::vec& Xk, const int& n){
  arma::mat out(n, n, arma::fill::zeros);
  for(int i(0); i < (n - 1); ++ i){
    out.submat(i + 1, i, n - 1, i) = abs(Xk.subvec(i + 1, n - 1) - Xk(i));
  }
  return out + out.t();
}

// Xi < Xj
arma::umat fdatalower(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    out.col(i) = (Xk < Xk(i));
  }
  return out;
}

// Xi > Xj
arma::umat fdatagreater(const arma::vec& Xk, const int& n){
  arma::umat out(n, n, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    out.col(i) = (Xk > Xk(i));
  }
  return out;
}


//[[Rcpp::export]]
arma::cube fdatar(const arma::mat X, List ftovar, const int& nvar, const int& K){
  int n(X.n_rows);
  arma::cube out(n, n, nvar);
  int s(0);
  for(int k(0); k < K; ++ k){
    arma::uvec ftv = ftovar[k];
    int L          = ftv.n_elem;
    for(int l(0); l < L; ++ l){
      switch(ftv(l)) {
      case 1: //Xi
        out.slice(s) = fdatai(X.col(k), n);
        ++ s;
        break;
      case 2: //Xj
        out.slice(s) = fdataj(X.col(k), n);
        ++ s;
        break;
      case 3: //Xi + Xj
        out.slice(s) = fdatasum(X.col(k), n);
        ++ s;
        break;
      case 4: //Xi*Xj
        out.slice(s) = fdataprod(X.col(k), n);
        ++ s;
        break;
      case 5: //Same
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatasame(X.col(k), n));
        ++ s;
        break;
      case 6: //Adiff
        out.slice(s) = fdatadiff(X.col(k), n);
        ++ s;
        break;
      case 7: //1{Xi < Xj}
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatalower(X.col(k), n));
        ++ s;
        break;
      case 8: //1{Xi > Xj}
        out.slice(s) = arma::conv_to<arma::mat>::from(fdatagreater(X.col(k), n));
        ++ s;
      }
    }
  }
  return out;
}

// /*
//  * Fmvnorm samples one vector fastly from numtivariate normal
//  * dim   : is the vector dimension
//  * u     : is the mean vector
//  * sigma : is the covariance matrix
//  */
// arma::vec Fmvnorm(const double& dim, arma::vec u, arma::mat sigma) {
//   arma::vec x = arma::randn(dim,1);
//   return arma::chol(sigma).t()*x + u;
// }
// 
/* Simulation from truncated normal [low, +inf]
 * following Robert, C. P. (1995). Simulation of truncated normal variables.
 Statistics and computing, 5(2), 121-125.*/
// when low <= 0
inline double tnormneg(const double& low) {
  bool repeat = true;
  double z = 0;
  
  while (repeat) {
    z = Rf_rnorm(0.0, 1.0) ;
    repeat = (z <= low);
  }
  return z;
}

// when low > 0

double tnormpos(const double& low){
  // Init Values
  const double alpha = 0.5*(low + sqrt(pow(low, 2.0) + 4.0));
  double e = 0 ;
  double z = 0 ;
  double rho = 0 ;
  double u = 0 ;
  //
  
  // Init Valid Flag
  bool repeat = true;
  //
  
  // Loop Until Valid Draw
  while (repeat) {
    e = Rf_rexp(1.0) ;
    z = low + e / alpha ;
    
    rho = exp(-pow(alpha - z, 2.0) / 2) ;
    u = Rf_runif(0, 1) ;
    if (u <= rho) {
      repeat = false;
    }
  }
  return z ;
}

double tnorm(const double& low) {
  if(low <= 0){
    return tnormneg(low);
  }
  return tnormpos(low);
}

// Convert matrix to vec
// [[Rcpp::export]]
Eigen::VectorXd frMtoVbyCOL(List& u,
                            const Rcpp::IntegerVector& N,
                            const double& M) {
  int r(0), n, outn = sum(N*N - N);
  Eigen::VectorXd out(outn);
  for(int m(0); m < M; ++m) {
    Eigen::MatrixXd um = u[m];
    n                  = N(m) - 1;
    out.segment(r, n)  = um.block(1, 0, n, 1);    
    r                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)         = um.block(0, i, i, 1);
      out.segment(r + i, n - i) = um.block(i + 1, i, n - i, 1);
      r                        += n;
    }
    
    out.segment(r, n) = um.block(0, n, n, 1);
    r                += n;
  }
  return out;
}


// Convert symmetric matrix to vec
// [[Rcpp::export]]
Eigen::VectorXd frMtoVbyCOLsym(List& u,
                               const Rcpp::IntegerVector& N,
                               const double& M) {
  int r(0), n, outn = sum(N*N - N)/2;
  Eigen::VectorXd out(outn);
  for(int m(0); m < M; ++m) {
    Eigen::MatrixXd um = u[m];
    n                  = N(m) - 1;
    out.segment(r, n)  = um.block(1, 0, n, 1);    
    r                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, n - i)     = um.block(i + 1, i, n - i, 1);
      r                        += n - i;
    }
  }
  return out;
}


// compute mu_i + nu_j for i and j in the same group
// export is a vector as (2, 1), (3, 1), ... , (n, 1), (1, 2) (3, 2), ....
//[[Rcpp::export]]
arma::vec fmusum(const arma::vec& mu,
                 const arma::vec& nu,
                 const arma::umat& index, // where each j starts and ends in a
                 const arma::umat& indexgr, // where each group starts and ends in mu
                 const int& M,
                 const int& N) {
  int igr1, igr2, nm, j = 0, j1, j2;
  arma::vec mum, num, mumj; 
  
  arma::vec musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0); 
    igr2 = indexgr(m, 1);
    nm   = igr2 - igr1 + 1;
    mum  = mu.subvec(igr1, igr2);
    num  = nu.subvec(igr1, igr2);
    
    for(int i(0); i < nm; ++ i){
      mumj   = mum;   // all mu in his group
      mumj.shed_row(i); // all nu excepted nu(i)
      j1     = index(j, 0); //starting individual i observation in dx
      j2     = index(j, 1); //ending individual i observation in dx
      musum.subvec(j1, j2) = num(i) + mumj;
      ++j;
    }
  }
  return musum; 
}

// for symmetric models
// export is a vector as (2, 1), (3, 1), ... , (n, 1), (3, 2) (4, 2), ....
//[[Rcpp::export]]
arma::vec fmusumsym(const arma::vec& mu,
                    const arma::umat& index,  // where each j starts and ends in a
                    const arma::umat& indexgr,
                    const int& M,
                    const int& N) {
  int igr1, igr2, nm, j = 0, j1, j2;
  arma::vec mum; 
  
  arma::vec musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0); 
    igr2 = indexgr(m, 1);
    nm   = igr2 - igr1 + 1;
    mum  = mu.subvec(igr1, igr2);
    
    for(int i(0); i < (nm - 1); ++ i){
      j1     = index(j, 0); //starting individual i observation in dx
      j2     = index(j, 1); //ending individual i observation in dx
      musum.subvec(j1, j2) = mum(i) + mum.tail(nm - 1 - i);
      ++j;
    }
    ++ j;
  }
  return musum; 
}

void updateast(arma::vec& ast,
               const arma::vec& dxbeta,
               const arma::vec& mupnu,
               const arma::vec& a,
               const int& N) {
  
  arma::vec mean = dxbeta + mupnu;
  for (int i(0); i < N; ++ i){
    if(a(i) == 1){
      ast(i) = tnorm(-mean(i)) + mean(i);
    } else {
      ast(i) = -tnorm(mean(i)) + mean(i);
    }
  }
}

void updatebeta(arma::vec& beta,
                arma::vec& dxbeta,
                const arma::umat& INDEXgr,
                const int& nfix,
                const int& Kx,
                const arma::mat& dx,
                const arma::mat& invdxdx,
                const arma::vec& mupnu,
                const arma::vec& ast) {
  arma::vec tmp    = ast - mupnu;
  arma::vec mubeta = dx.t()*tmp;
  if(nfix > 1) {
    mubeta         = arma::join_cols(arma::ones(nfix), mubeta);
    for(int m(0); m < nfix; ++ m){
      mubeta(m)    = sum(tmp.subvec(INDEXgr(m, 0), INDEXgr(m, 1)));
    }
  }
  mubeta           = invdxdx*mubeta;
  beta             = rmvnorm(1, mubeta, invdxdx).t();
  dxbeta           = dx*beta.tail(Kx);
  if(nfix > 1) {
    for(int m(0); m < nfix; ++m){
      dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
    }
  }
}



void updatemunu(arma::vec& mu,
                arma::vec& nu,
                arma::vec& mupnu,
                arma::vec& beta,
                arma::vec& dxbeta,
                const arma::mat& dx,
                const arma::vec& ast,
                const int& Kx,
                const arma::umat& INDEXgr,
                const int& M,
                const int& N,
                const int& n,
                const int& nfix,
                const arma::vec& nvec,
                const arma::umat& index,  // where each j starts and ends in a
                const arma::umat& indexgr,
                const double& smu2,
                const double& snu2,
                const double& rho){
  int j1, j2, igr1, igr2, nm, j = 0;
  double tmpmu, tmpnu;
  arma::vec num, mum, numj, mumj, munub, munu; 
  arma::uvec indexj;
  arma::umat indexm;
  arma::mat smunu(2, 2), vmunu(2, 2), invvmunu;
  // ast minus dxbeta
  arma::vec astmdxbeta = ast - dxbeta;
  
  // prior variance
  vmunu(0, 0) = smu2;
  vmunu(1, 1) = snu2;
  vmunu(0, 1) = rho*sqrt(smu2*snu2);
  vmunu(1, 0) = vmunu(0, 1);
  invvmunu    = arma::inv(vmunu);
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0); 
    igr2   = indexgr(m, 1);
    nm     = nvec(m);
    indexm = index.rows(igr1, igr2);
    // posterior variance
    smunu  = arma::inv((nm - 1)*arma::eye(2, 2) + invvmunu);
    
    for(int i(0); i < nm; ++ i){
      mum     = mu.subvec(igr1, igr2);
      num     = nu.subvec(igr1, igr2);
      
      // elements that vary when nuj is fixed
      mumj    = mum;
      mumj.shed_row(i);
      j1      = index(j, 0);
      j2      = index(j, 1);
      tmpnu   = sum(astmdxbeta.subvec(j1, j2) - mumj);
      
      // elements that vary when mui is fixed
      indexj  = indexm.col(0) + i;
      indexj.head(i + 1) -= 1;
      indexj.shed_row(i);  
      numj    = num;
      numj.shed_row(i);
      tmpmu   = sum(astmdxbeta.elem(indexj) - numj);
      
      // posterior mean
      arma::vec tmp = {tmpmu, tmpnu};
      munub         = smunu*tmp;
      
      // simulation
      munu    = rmvnorm(1, munub, smunu).t();
      mu(j)   = munu(0);
      nu(j)   = munu(1);
      ++ j;
    }
  }
  
  // normalize mu and nu to zero mean
  if(nfix == 1){
    double mub       = sum(mu)/n;
    double nub       = sum(nu)/n;
    mu              -= mub;
    nu              -= nub;
    beta(0)         += (mub + nub);
    dxbeta           = dx*beta;
  } else{
    if(nfix > 1){
      for(int m(0); m < M; ++ m) {
        igr1         = indexgr(m, 0);
        igr2         = indexgr(m, 1);
        double mub   = sum(mu.subvec(igr1, igr2))/nvec(m);
        double nub   = sum(nu.subvec(igr1, igr2))/nvec(m);
        mu.subvec(igr1, igr2) -= mub;
        nu.subvec(igr1, igr2) -= nub;
        beta(m)     += (mub + nub);
      }
      dxbeta         = dx*beta.tail(Kx);
      for(int m(0); m < nfix; ++m){
        dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
      }
    }
  }
  mupnu              = fmusum(mu, nu, index, indexgr, M, N);
}


void updatemu(arma::vec& mu,
              arma::vec& mupmu,
              arma::vec& beta,
              arma::vec& dxbeta,
              const arma::mat& dx,
              const arma::vec& ast,
              const int& Kx,
              const arma::umat& INDEXgr, 
              const int& M,
              const int& N,
              const int& n,
              const int& nfix,
              const arma::vec& nvec,
              const arma::umat& index,
              const arma::umat& indexgr, // where each j starts and ends in a
              const double& smu2){
  int j1, j2, igr1, igr2, nm, j = 0;
  double tmp, smu, mub, invvmu;
  arma::vec mum, mumj; 
  arma::uvec indexj;
  arma::umat indexm;
  // ast minus dxbeta
  arma::vec astmdxbeta = ast - dxbeta;
  
  // prior variance
  invvmu   = 1.0/smu2;
  
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0); 
    igr2   = indexgr(m, 1);
    nm     = nvec(m);
    indexm = index.rows(igr1, igr2);
    // posterior variance
    smu    = sqrt(1.0/(2.0*(nm - 1) + invvmu));
    
    for(int i(0); i < nm; ++ i){
      mum     = mu.subvec(igr1, igr2);
      
      // elements that vary when muj is fixed
      mumj    = mum;
      mumj.shed_row(i);
      j1      = index(j, 0);
      j2      = index(j, 1);
      tmp     = sum(astmdxbeta.subvec(j1, j2) - mumj);
      
      // elements that vary when mui is fixed
      indexj  = indexm.col(0) + i;
      indexj.head(i + 1) -= 1;
      indexj.shed_row(i);  
      tmp    += sum(astmdxbeta.elem(indexj) - mumj);
      
      // posterior mean
      mub     = pow(smu, 2)*tmp;
      
      // simulation
      mu(j)   = R::rnorm(mub, smu);
      ++ j;
    }
  }
  
  // normalize mu to zero mean
  if(nfix == 1){
    double mub       = sum(mu)/n;
    mu              -= mub;
    beta(0)         += (2*mub);
    dxbeta           = dx*beta;
  } else{
    if(nfix > 1){
      for(int m(0); m < M; ++ m) {
        igr1         = indexgr(m, 0);
        igr2         = indexgr(m, 1);
        double mub   = sum(mu.subvec(igr1, igr2))/nvec(m);
        mu.subvec(igr1, igr2) -= mub;
        beta(m)     += (2*mub);
      }
      dxbeta         = dx*beta.tail(Kx);
      for(int m(0); m < nfix; ++m){
        dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
      }
    }
  }
  mupmu              = fmusum(mu, mu, index, indexgr, M, N);
}

void updatemusym(arma::vec& mu,
                 arma::vec& mupmu,
                 arma::vec& beta,
                 arma::vec& dxbeta,
                 const arma::mat& dx,
                 const arma::vec& ast,
                 const int& Kx,
                 const arma::umat& INDEXgr, 
                 const int& M,
                 const int& N,
                 const int& n,
                 const int& nfix,
                 const arma::vec& nvec,
                 const arma::umat& index,
                 const arma::umat& indexgr, // where each j starts and ends in a
                 const double& smu2){
  int j1, j2, igr1, igr2, nm, j = 0;
  double smu, mub, invvmu;
  arma::vec mum, mumj; 
  arma::uvec indexj;
  arma::umat indexm;
  // ast minus dxbeta
  arma::vec astmdxbeta = ast - dxbeta;
  
  // prior variance
  invvmu   = 1.0/smu2;
  
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0); 
    igr2   = indexgr(m, 1);
    nm     = nvec(m);
    indexm = index.rows(igr1, igr2);
    // posterior variance
    smu    = sqrt(1.0/(nm - 1 + invvmu));
    mum    = mu.subvec(igr1, igr2);
    
    for(int i(0); i < nm; ++ i){
      // elements that vary when muj is fixed
      double tmp(0);
      if(i < (nm - 1)){
        j1  = index(j, 0);
        j2  = index(j, 1);
        tmp = sum(astmdxbeta.subvec(j1, j2) - mum.tail(nm - 1 - i));
      }
      
      // elements that vary when mui is fixed
      if(i > 0){
        indexj  = indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i);
        tmp    += sum(astmdxbeta.elem(indexj) - mum.head(i));
      }
      
      // posterior mean
      mub     = pow(smu, 2)*tmp;
      
      // simulation
      mu(j)   = R::rnorm(mub, smu);
      ++ j;
    }
  }
  
  // normalize mu to zero mean
  if(nfix == 1){
    double mub       = sum(mu)/n;
    mu              -= mub;
    beta(0)         += (2*mub);
    dxbeta           = dx*beta;
  } else{
    if(nfix > 1){
      for(int m(0); m < M; ++ m) {
        igr1         = indexgr(m, 0);
        igr2         = indexgr(m, 1);
        double mub   = sum(mu.subvec(igr1, igr2))/nvec(m);
        mu.subvec(igr1, igr2) -= mub;
        beta(m)     += (2*mub);
      }
      dxbeta         = dx*beta.tail(Kx);
      for(int m(0); m < nfix; ++m){
        dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
      }
    }
  }
  mupmu              = fmusumsym(mu, index, indexgr, M, N);
}


void updatesigma(double& smu2,
                 double& snu2,
                 double& rho,
                 arma::mat& Sigma,
                 const int& n,
                 const arma::vec& mu,
                 const arma::vec& nu) {
  
  arma::mat d   = arma::join_rows(mu, nu);
  Sigma = riwish(n, d.t()*d);
  smu2  = Sigma(0, 0);
  snu2  = Sigma(1, 1);
  rho   = Sigma(0, 1)/sqrt(smu2*snu2);
}


void updatesmu2(double& smu2,
                const int& n,
                const arma::vec& mu){
  smu2 = 1.0/R::rgamma(0.5*n, 1.0/(0.5*sum(mu%mu)));
}

void updatellhmunu(double& llh,
                   const arma::vec& mu,
                   const arma::vec& nu,
                   const arma::vec& a,
                   const arma::vec& dxbeta,
                   const arma::vec& mupnu,
                   const arma::mat& Sigma,
                   const int& n){
  arma::vec tmp  = dxbeta + mupnu;
  arma::vec tmpp = tmp.elem(arma::find(a == 0));
  NumericVector tmp0 = wrap(tmpp);
  tmpp = tmp.elem(arma::find(a == 1));
  NumericVector tmp1 = wrap(tmpp);
  arma::mat tmp2 = arma::join_rows(mu, nu);
  llh = sum(Rcpp::pnorm5(tmp0, 0, 1, false, true)) + sum(Rcpp::pnorm5(tmp1, 0, 1, true, true)) - 
    0.5*arma::accu((tmp2*arma::inv(Sigma))%tmp2) -  0.5*n*log(arma::det(Sigma));
}

//[[Rcpp::export]]
List bayesmunu(const arma::vec& a,
               const arma::mat& dx,
               const arma::mat& invdxdx,
               const arma::vec& beta0,
               const arma::vec& mu0,
               const arma::vec& nu0,
               const double& smu20,
               const double& snu20,
               const double& rho0,
               const arma::umat& index,
               const arma::umat& indexgr,
               const arma::umat& INDEXgr,
               const int& nfix,
               const int& N,
               const int& M,
               const int& K,
               const int& Kx,
               const arma::vec& nvec, 
               const int& n,
               const int& iteration,
               const bool& Print){
  // init
  arma::vec beta       = beta0;
  arma::vec mu         = mu0;
  arma::vec nu         = nu0;
  double smu2          = smu20;
  double snu2          = snu20;
  double rho           = rho0;
  arma::vec mupnu      = fmusum(mu, nu, index, indexgr, M, N);
  arma::vec dxbeta     = dx*beta.tail(Kx);
  if(nfix > 1) {
    for(int m(0); m < nfix; ++m){
      dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
    }
  }
  arma::vec ast(N, arma::fill::zeros);
  arma::mat Sigma;
  // output
  arma::mat Smu(n, iteration), Snu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration), Ssnu2(iteration), Srho(iteration); //, Sllh(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  Ssnu2.attr("dim") = R_NilValue;
  Srho.attr("dim")  = R_NilValue;
  //Sllh.attr("dim")  = R_NilValue;
  
  if(Print){
    Progress p(iteration, true);
    for(int t(0); t < iteration; ++ t) {
      p.increment(); 
      
      //update ast
      updateast(ast, dxbeta, mupnu, a, N);
      
      //update beta
      updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupnu, ast);
      
      //update mu and nu
      updatemunu(mu, nu, mupnu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2, snu2, rho);
      
      //update sigmas
      updatesigma(smu2, snu2, rho, Sigma, n, mu, nu);
      
      //likelihood
      //updatellhmunu(llh, mu, nu, a, dxbeta, mupnu, Sigma, n);
      
      //save
      Smu.col(t)    = mu;
      Snu.col(t)    = nu;
      Sbeta.col(t)  = beta;
      Ssmu2(t)      = smu2;
      Ssnu2(t)      = snu2;
      Srho(t)       = rho;
      //Sllh(t)       = llh;
    }
  } else {
    for(int t(0); t < iteration; ++ t) {
      
      //update ast
      updateast(ast, dxbeta, mupnu, a, N);
      
      //update beta
      updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupnu, ast);
      
      //update mu and nu
      updatemunu(mu, nu, mupnu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2, snu2, rho);
      
      //update sigmas
      updatesigma(smu2, snu2, rho, Sigma, n, mu, nu);
      
      //likelihood
      //updatellhmunu(llh, mu, nu, a, dxbeta, mupnu, Sigma, n);
      
      //save
      Smu.col(t)    = mu;
      Snu.col(t)    = nu;
      Sbeta.col(t)  = beta;
      Ssmu2(t)      = smu2;
      Ssnu2(t)      = snu2;
      Srho(t)       = rho;
      //Sllh(t)       = llh;
    }
  }
  return List::create(Named("beta")      = Sbeta.t(),
                      Named("mu")        = Smu.t(),
                      Named("nu")        = Snu.t(),
                      Named("sigma2_mu") = Ssmu2,
                      Named("sigma2_nu") = Ssnu2,
                      Named("rho")       = Srho);//,
  //Named("loglike")   = Sllh);
}


//[[Rcpp::export]]
List bayesmu(const arma::vec& a,
             const arma::mat& dx,
             const arma::mat& invdxdx,
             const arma::vec& beta0,
             const arma::vec& mu0,
             const double& smu20,
             const arma::umat& index,
             const arma::umat& indexgr,
             const arma::umat& INDEXgr,
             const int& nfix,
             const int& N,
             const int& M,
             const int& K,
             const int& Kx,
             const arma::vec& nvec, 
             const int& n,
             const int& iteration,
             const bool& sym,
             const bool& Print){
  // init
  arma::vec beta       = beta0;
  arma::vec mu         = mu0;
  double smu2          = smu20;
  arma::vec mupmu;
  if(sym){
    mupmu              = fmusumsym(mu, index, indexgr, M, N);
  } else {
    mupmu              = fmusum(mu, mu, index, indexgr, M, N);
  }
  
  arma::vec dxbeta     = dx*beta.tail(Kx);
  if(nfix > 1) {
    for(int m(0); m < nfix; ++m){
      dxbeta.subvec(INDEXgr(m, 0), INDEXgr(m, 1)) += beta(m);
    }
  }
  arma::vec ast(N, arma::fill::zeros);
  
  // output
  arma::mat Smu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration); //, Sllh(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  
  if(Print){
    Progress p(iteration, true);
    for(int t(0); t < iteration; ++ t) {
      p.increment(); 
      
      //update ast
      updateast(ast, dxbeta, mupmu, a, N);
      
      //update beta
      updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupmu, ast);
      
      //update mu 
      if(sym){
        updatemusym(mu, mupmu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2);
      } else {
        updatemu(mu, mupmu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2);
      }
      
      //update smu2
      updatesmu2(smu2, n, mu);
      
      //save
      Smu.col(t)    = mu;
      Sbeta.col(t)  = beta;
      Ssmu2(t)      = smu2;
      //Sllh(t)       = llh;
    }
  } else {
    for(int t(0); t < iteration; ++ t) {
      //update ast
      updateast(ast, dxbeta, mupmu, a, N);
      
      //update beta
      updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupmu, ast);
      
      //update mu and nu
      if(sym){
        updatemusym(mu, mupmu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2);
      } else {
        updatemu(mu, mupmu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2);
      }
      
      //update smu2
      updatesmu2(smu2, n, mu);
      
      //save
      Smu.col(t)    = mu;
      Sbeta.col(t)  = beta;
      Ssmu2(t)      = smu2;
      //Sllh(t)       = llh;
    }
  }
  return List::create(Named("beta")      = Sbeta.t(),
                      Named("mu")        = Smu.t(),
                      Named("sigma2_mu") = Ssmu2);//,
  //Named("loglike")   = Sllh);
}


/////////////////// Logit model
// Estimation using fixed effects two sides
class llhhomo2f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& index;
  const arma::umat& indexgr;
  const arma::uvec& nvec;
  const int& M;
  const int& n;
  const int& Kx;
  const int& nparms;
  const bool& hasX;
  const bool& Print;
public:
  llhhomo2f(const arma::vec& a_,
            const arma::mat& dx_,
            const arma::mat& adx_,
            const arma::vec& d_,
            const arma::vec& b_,
            const arma::umat& index_,
            const arma::umat& indexgr_,
            const arma::uvec& nvec_,
            const int& M_,
            const int& n_,
            const int& Kx_,
            const int& nparms_,
            const bool& hasX_,
            const bool& Print_) : 
  a(a_),
  dx(dx_),
  adx(adx_),
  d(d_),
  b(b_),
  index(index_),
  indexgr(indexgr_),
  nvec(nvec_),
  M(M_),
  n(n_),
  Kx(Kx_),
  nparms(nparms_),
  hasX(hasX_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec thetaa       = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    
    // int b1                 = 0;
    int b2                 = Kx - 1;
    int m1                 = b2 + 1;
    int m2                 = m1 + n - 1;
    int n1                 = m2 + 1;
    int n2                 = n1 + n - M - 1;
    
    arma::vec mu           = thetaa.subvec(m1, m2);
    arma::vec nu           = thetaa.subvec(n1, n2);
    
    arma::vec gd(nparms, arma::fill::zeros);
    
    arma::vec beta, dxb, adxb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dxb                  = dx*beta;
      adxb                 = a%dxb;
      llh                  = sum(adxb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, num, numj, mumj, tmp, ai, exbmn, smunu;
    arma::mat dxi;
    arma::uvec indexi;
    arma::umat indexm;
    
    for (int m(0); m < M; ++ m) {
      igr1                 = indexgr(m, 0); 
      igr2                 = indexgr(m, 1); 
      nm                   = nvec(m);
      indexm               = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
      mum                  = mu.subvec(igr1, igr2);
      num                  = arma::zeros<arma::vec>(nm); num.head(nm - 1) = nu.subvec(igr1 - m, igr2 - m - 1);
      for(int i(0); i < nm; ++ i){
        j1                 = index(j, 0);
        j2                 = index(j, 1);
        ai                 = a.subvec(j1, j2);
        
        // nuj when mui is fixed
        numj               = num;
        numj.shed_row(i);
        
        // muj when nui is fixed
        mumj               = mum;
        mumj.shed_row(i);
        
        // rows on which nui is used
        indexi             = indexm.col(0) + i;
        indexi.head(i + 1)-= 1;
        indexi.shed_row(i);  
        
        // sum of muj + num(i)
        smunu              = mumj + num(i);
        
        if(hasX){
          exbmn            = exp(dxb.subvec(j1, j2) + smunu);
          tmp              = exbmn/(1 + exbmn);
          llh             += sum(ai%smunu - log(1 + exbmn));
          
          // grad X
          dxi              = dx.rows(j1, j2);
          gd.head(Kx)     -= arma::trans(arma::sum(dxi.each_col()%tmp, 0));
          
          // grad nui
          if(i < (nm - 1)){
            gd(n1 + j - m) = b(j) - sum(tmp);
          }
          
          // grad mui
          tmp              = exp(dxb.elem(indexi) + mum(i) + numj);
          gd(m1 + j)       = d(j) - sum(tmp/(1 + tmp));
        } else {
          exbmn            = exp(smunu);
          tmp              = exbmn/(1 + exbmn);
          llh             += sum(ai%smunu - log(1 + exbmn));
          
          // grad nui
          if(i < (nm - 1)){
            gd(n1 + j - m) = b(j) - sum(tmp);
          }
          
          // grad mui
          tmp              = exp(mum(i) + numj);
          gd(m1 + j)       = d(j) - sum(tmp/(1 + tmp));
        }
        ++ j;
      }
    }
    
    grad                   = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                   = gd;
    
    if(Print){
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "\nbeta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
    }
    return -llh;
  }
};


//[[Rcpp::export]]
List fhomobeta2f(Eigen::VectorXd theta,
                 const arma::vec& a,
                 const arma::mat& dx,
                 const arma::uvec& nvec,
                 const arma::umat& index,
                 const arma::umat& indexgr,
                 const int& M,       
                 const int maxit = 300, 
                 const double& eps_f = 1e-6, 
                 const double& eps_g = 1e-5,
                 const bool& hasX = true,
                 const bool& Print = true){
  int n(sum(nvec)), Kx(0);
  arma::mat adx;
  if(hasX){
    Kx          = dx.n_cols;
    adx         = dx.each_col()%a;
  } 
  
  int nparms    = Kx + 2*n - M;
  arma::vec d(n), b(n);
  int j(0);
  for (int m(0); m < M; ++ m) {
    int igr1              = indexgr(m, 0); // group m starts from igr1 in X
    int igr2              = indexgr(m, 1); // group m ends at igr2 in X
    int nm                = nvec(m);
    arma::umat indexm     = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = indexm.col(0) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      d(j)                = sum(a.elem(indexi));
      ++ j;
    }
  }
  
  llhhomo2f f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms, hasX, Print);
  
  double fopt;
  int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}

// optimize only mu, nu
class llhhomomunu2f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::vec& dxb;
  const double& sadxb;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& indexm;
  arma::vec& llh;
  const int& m;
  const int& nm;
  const int& nparms;
  const bool& Print;
public:
  llhhomomunu2f(const arma::vec& a_,
                const arma::vec& dxb_,
                const double& sadxb_,
                const arma::vec& d_,
                const arma::vec& b_,
                const arma::umat& indexm_,
                arma::vec& llh_,
                const int& m_,
                const int& nm_,
                const int& nparms_,
                const bool& Print_) : 
  a(a_),
  dxb(dxb_),
  sadxb(sadxb_),
  d(d_),
  b(b_),
  indexm(indexm_),
  llh(llh_),
  m(m_),
  nm(nm_),
  nparms(nparms_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec thetaa       = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    arma::vec gd(nparms, arma::fill::zeros);
    llh(m)                 = sadxb;
    int j1, j2;
    arma::vec numj, mumj, tmp, ai, exbmn, smunu;
    arma::uvec indexi;
    
    arma::vec mum          = thetaa.head(nm);
    arma::vec num          = arma::zeros<arma::vec>(nm); num.head(nm - 1) = thetaa.tail(nm - 1);
    for(int i(0); i < nm; ++ i){
      j1                   = indexm(i, 0);
      j2                   = indexm(i, 1);
      ai                   = a.subvec(j1, j2);
      
      // nuj when mui is fixed
      numj                 = num;
      numj.shed_row(i);
      
      // muj when nui is fixed
      mumj                 = mum;
      mumj.shed_row(i);
      
      // rows on which nui is used
      indexi               = indexm.col(0) + i;
      indexi.head(i + 1)  -= 1;
      indexi.shed_row(i);  
      
      // sum of muj + num(i)
      smunu                = mumj + num(i);
      
      exbmn                = exp(dxb.subvec(j1, j2) + smunu);
      tmp                  = exbmn/(1 + exbmn);
      llh(m)              += sum(ai%smunu - log(1 + exbmn));
      
      // grad nui
      if(i < (nm - 1)){
        gd(nm + i)         = b(i) - sum(tmp);
      }
      
      // grad mui
      tmp                  = exp(dxb.elem(indexi) + mum(i) + numj);
      gd(i)                = d(i) - sum(tmp/(1 + tmp));
    }
    
    grad                   = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                   = gd;
    return -llh(m);
  }
};

//[[Rcpp::export]]
double fhomobetamunu2f(arma::vec& theta,
                       const arma::vec& a,
                       const arma::vec& dxb,
                       const arma::uvec& nvec,
                       const arma::umat& index,
                       const arma::umat& indexgr,
                       const int& N,
                       const int& M,   
                       const int& Kx,
                       const int maxit = 300, 
                       const double& eps_f = 1e-6, 
                       const double& eps_g = 1e-5,
                       const bool& hasX = true,
                       const bool& Print = true){
  int n(sum(nvec)), j(0);
  arma::vec adxb(N, arma::fill::zeros), mu(theta.subvec(Kx, Kx + n - 1)), nu(theta.subvec(Kx + n, Kx + 2*n - M - 1)),
  d(n), b(n), mumj, llh(M), musum(arma::zeros<arma::vec>(N));
  if(hasX){
    adxb        = a%dxb;
  } 
  
  for (int m(0); m < M; ++ m) {
    int nm(nvec(m)), igr1(indexgr(m, 0)), igr2(indexgr(m, 1));
    arma::umat indexm(index.rows(igr1, igr2)); // ith row is the row at each i interacts with others, where the link goes from i
    arma::vec mum(mu.subvec(igr1, igr2));
    arma::vec num(arma::zeros(nm)); num.head(nm - 1) = nu.subvec(igr1 - m, igr2 - m - 1);
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = indexm.col(0) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      d(j)                = sum(a.elem(indexi));
      
      // sum mu nu
      mumj                = mum;   // all mu in his group
      mumj.shed_row(i);            // all mu excepted mu(i)
      musum.subvec(index(j, 0), index(j, 1)) = num(i) + mumj;
      ++ j;
    }
    arma::vec tp(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum.subvec(index(igr1, 0), index(igr2, 1)));
    llh(m)                = sum(a.subvec(index(igr1, 0), index(igr2, 1))%tp - log(1 + exp(tp)));
  }
  
  
  
  for (int m(0); m < M; ++ m) {
    int igr1(indexgr(m, 0)), igr2(indexgr(m, 1)), N1(index(igr1, 0)), N2(index(igr2, 1)), nm(nvec(m)), nparms(2*nm - 1);
    arma::umat indexm(index.rows(igr1, igr2) - index(igr1, 0));
    arma::vec mum(mu.subvec(igr1, igr2)), num(nu.subvec(igr1 - m, igr2 - m - 1)), thetacpp(arma::join_cols(mum, num)),
    am(a.subvec(N1, N2)), dxbm(dxb.subvec(N1, N2)), dm(d.subvec(igr1, igr2)), bm(b.subvec(igr1, igr2));
    Eigen::VectorXd thetaEi(Eigen::Map<Eigen::VectorXd>(thetacpp.memptr(), nparms));
    double fopt, adxbm(sum(adxb.subvec(N1, N2)));
    llhhomomunu2f f(am, dxbm, adxbm, dm, bm, indexm, llh, m, nm, nparms, Print);
    optim_lbfgs(f, thetaEi, fopt, maxit, eps_f, eps_g);
    llh(m)     = -fopt;
    mu.subvec(igr1, igr2)             = arma::vec(thetaEi.segment(0, nm).data(), nm, false, false); 
    nu.subvec(igr1 - m, igr2 - m - 1) = arma::vec(thetaEi.segment(nm, nm - 1).data(), nm - 1, false, false); 
    if(Print){
      Rcpp::Rcout << "group: " << m + 1 << " -- log-likelihood: " << sum(llh) << "\n";
    }
  }
  
  theta.subvec(Kx, Kx + n - 1)           = mu;
  theta.subvec(Kx + n, Kx + 2*n - M - 1) = nu;
  return sum(llh);
}


// Estimation using fixed effects only mu, one way
class llhhomo1f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& index;
  const arma::umat& indexgr;
  const arma::uvec& nvec;
  const int& M;
  const int& n;
  const int& Kx;
  const int& nparms;
  const bool& hasX;
  const bool& Print;
public:
  llhhomo1f(const arma::vec& a_,
            const arma::mat& dx_,
            const arma::mat& adx_,
            const arma::vec& d_,
            const arma::vec& b_,
            const arma::umat& index_,
            const arma::umat& indexgr_,
            const arma::uvec& nvec_,
            const int& M_,
            const int& n_,
            const int& Kx_,
            const int& nparms_,
            const bool& hasX_,
            const bool& Print_) : 
  a(a_),
  dx(dx_),
  adx(adx_),
  d(d_),
  b(b_),
  index(index_),
  indexgr(indexgr_),
  nvec(nvec_),
  M(M_),
  n(n_),
  Kx(Kx_),
  nparms(nparms_),
  hasX(hasX_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec thetaa       = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    
    // int b1                 = 0;
    int b2                 = Kx - 1;
    int m1                 = b2 + 1;
    int m2                 = m1 + n - 1;
    
    arma::vec mu           = thetaa.subvec(m1, m2);
    arma::vec gd(nparms, arma::fill::zeros);
    arma::vec beta, dxb, adxb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dxb                  = dx*beta;
      adxb                 = a%dxb;
      llh                  = sum(adxb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, mumj, tmp, ai, exbmn, smunu;
    arma::mat dxi;
    arma::uvec indexi;
    arma::umat indexm;
    
    for (int m(0); m < M; ++ m) {
      igr1                 = indexgr(m, 0); 
      igr2                 = indexgr(m, 1); 
      nm                   = nvec(m);
      indexm               = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
      mum                  = mu.subvec(igr1, igr2);
      
      for(int i(0); i < nm; ++ i){
        j1                 = index(j, 0);
        j2                 = index(j, 1);
        ai                 = a.subvec(j1, j2);
        
        // muj when nui is fixed
        mumj               = mum;
        mumj.shed_row(i);
        
        // sum of mu(i) + numj
        smunu              = mum(i) + mumj;
        
        if(hasX){
          dxi                = dx.rows(j1, j2);
          exbmn              = exp(dxb.subvec(j1, j2) + smunu);
          tmp                = exbmn/(1 + exbmn);
          llh               += sum(ai%smunu - log(1 + exbmn));
          
          // grad X
          gd.head(Kx)       -= arma::trans(arma::sum(dxi.each_col()%tmp, 0));
          
          // grad mui
          gd(m1 + j)         = b(j) - sum(tmp);
          indexi             = indexm.col(0) + i;
          indexi.head(i + 1)-= 1;
          indexi.shed_row(i); 
          tmp                = exp(dxb.elem(indexi) + mumj + mum(i));
          gd(m1 + j)        += d(j) - sum(tmp/(1 + tmp));
        } else {
          exbmn              = exp(smunu);
          tmp                = exbmn/(1 + exbmn);
          llh               += sum(ai%smunu - log(1 + exbmn));
          
          // grad mui
          gd(m1 + j)         = b(j) - sum(tmp);
          indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
          indexi.head(i + 1)-= 1;
          indexi.shed_row(i); 
          tmp                = exp(mumj + mum(i));
          gd(m1 + j)        += d(j) - sum(tmp/(1 + tmp));
        }
        ++ j;
      }
    }
    
    grad                   = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                   = gd;
    
    if(Print){
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "\nbeta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
    }
    return -llh;
  }
};

//[[Rcpp::export]]
List fhomobeta1f(Eigen::VectorXd theta,
                 const arma::vec& a,
                 const arma::mat& dx,
                 const arma::uvec& nvec,
                 const arma::umat& index,
                 const arma::umat& indexgr,
                 const int& M,       
                 const int maxit = 300, 
                 const double& eps_f = 1e-6, 
                 const double& eps_g = 1e-5,
                 const bool& hasX = true,
                 const bool& Print = true){
  int n(sum(nvec)), Kx(0);
  arma::mat adx;
  if(hasX){
    Kx          = dx.n_cols;
    adx         = dx.each_col()%a;
  } 
  
  int nparms    = Kx + n;
  arma::vec d(n), b(n);
  int j(0);
  for (int m(0); m < M; ++ m) {
    int igr1              = indexgr(m, 0); // group m starts from igr1 in X
    int igr2              = indexgr(m, 1); // group m ends at igr2 in X
    int nm                = nvec(m);
    arma::umat indexm     = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = indexm.col(0) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      d(j)                = sum(a.elem(indexi));
      ++ j;
    }
  }
  
  llhhomo1f f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms, hasX, Print);
  
  double fopt;
  int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}

// optimize only mu
class llhhomomu1f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::vec& dxb;
  const double& sadxb;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& indexm;
  arma::vec& llh;
  const int& m;
  const int& nm;
  const int& nparms;
  const bool& Print;
public:
  llhhomomu1f(const arma::vec& a_,
              const arma::vec& dxb_,
              const double& sadxb_,
              const arma::vec& d_,
              const arma::vec& b_,
              const arma::umat& indexm_,
              arma::vec& llh_,
              const int& m_,
              const int& nm_,
              const int& nparms_,
              const bool& Print_) : 
  a(a_),
  dxb(dxb_),
  sadxb(sadxb_),
  d(d_),
  b(b_),
  indexm(indexm_),
  llh(llh_),
  m(m_),
  nm(nm_),
  nparms(nparms_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec mum          = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    
    arma::vec gd(nparms, arma::fill::zeros);
    llh(m)                 = sadxb;
    
    int j1, j2;
    arma::vec mumj, tmp, ai, exbmn, smunu;
    arma::uvec indexi;
    for(int i(0); i < nm; ++ i){
      j1                 = indexm(i, 0);
      j2                 = indexm(i, 1);
      ai                 = a.subvec(j1, j2);
      
      // muj when nui is fixed
      mumj               = mum;
      mumj.shed_row(i);
      
      // sum of mu(i) + numj
      smunu              = mum(i) + mumj;
      exbmn              = exp(dxb.subvec(j1, j2) + smunu);
      tmp                = exbmn/(1 + exbmn);
      
      // grad mui
      gd(i)              = b(i) - sum(tmp);
      indexi             = indexm.col(0) + i;
      indexi.head(i + 1)-= 1;
      indexi.shed_row(i); 
      tmp                = exp(dxb.elem(indexi) + mumj + mum(i));
      gd(i)             += d(i) - sum(tmp/(1 + tmp));
      llh(m)            += sum(ai%smunu - log(1 + exbmn));
    }
    
    grad                 = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                 = gd;
    return -llh(m);
  }
};

//[[Rcpp::export]]
double fhomobetamu1f(arma::vec& theta,
                     const arma::vec& a,
                     const arma::vec& dxb,
                     const arma::uvec& nvec,
                     const arma::umat& index,
                     const arma::umat& indexgr,
                     const int& N,
                     const int& M,     
                     const int& Kx,
                     const int maxit = 300, 
                     const double& eps_f = 1e-6, 
                     const double& eps_g = 1e-5,
                     const bool& hasX = true,
                     const bool& Print = true){
  int n(sum(nvec)), j(0);
  arma::vec adxb(N, arma::fill::zeros), mu(theta.subvec(Kx, Kx + n - 1)), d(n), b(n), mumj, 
  llh(M), musum(arma::zeros<arma::vec>(N));
  if(hasX){
    adxb        = a%dxb;
  } 

  for (int m(0); m < M; ++ m) {
    int nm(nvec(m)), igr1(indexgr(m, 0)), igr2(indexgr(m, 1));
    arma::umat indexm(index.rows(igr1, igr2)); // ith row is the row at each i interacts with others, where the link goes from i
    arma::vec mum(mu.subvec(igr1, igr2));
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = indexm.col(0) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      d(j)                = sum(a.elem(indexi));
      
      // sum mu nu
      mumj                = mum;   // all mu in his group
      mumj.shed_row(i);            // all mu excepted mu(i)
      musum.subvec(index(j, 0), index(j, 1)) = mum(i) + mumj;
      ++ j;
    }
    arma::vec tp(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum.subvec(index(igr1, 0), index(igr2, 1)));
    llh(m)                = sum(a.subvec(index(igr1, 0), index(igr2, 1))%tp - log(1 + exp(tp)));
  }
  
  
  for (int m(0); m < M; ++ m) {
    int igr1(indexgr(m, 0)), igr2(indexgr(m, 1)), N1(index(igr1, 0)), N2(index(igr2, 1)), nm(nvec(m)), nparms(nm);
    arma::umat indexm(index.rows(igr1, igr2) - index(igr1, 0));
    arma::vec mum(mu.subvec(igr1, igr2)), am(a.subvec(N1, N2)), dxbm(dxb.subvec(N1, N2)), dm(d.subvec(igr1, igr2)), 
    bm(b.subvec(igr1, igr2));
    Eigen::VectorXd thetaEi(Eigen::Map<Eigen::VectorXd>(mum.memptr(), nparms));
    double fopt, adxbm(sum(adxb.subvec(N1, N2)));
    llhhomomu1f f(am, dxbm, adxbm, dm, bm, indexm, llh, m, nm, nparms, Print);
    optim_lbfgs(f, thetaEi, fopt, maxit, eps_f, eps_g);
    llh(m)     = -fopt;
    mu.subvec(igr1, igr2) = arma::vec(thetaEi.data(), nm, false, false); 
    if(Print){
      Rcpp::Rcout << "group: " << m + 1 << " -- log-likelihood: " << sum(llh) << "\n";
    }
  }
  
  theta.subvec(Kx, Kx + n - 1)        = mu;
  return sum(llh);
}


// Estimation using fixed effects with symmetric networks
class llhhomosym: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& index;
  const arma::umat& indexgr;
  const arma::uvec& nvec;
  const int& M;
  const int& n;
  const int& Kx;
  const int& nparms;
  const bool& hasX;
  const bool& Print;
public:
  llhhomosym(const arma::vec& a_,
             const arma::mat& dx_,
             const arma::mat& adx_,
             const arma::vec& d_,
             const arma::vec& b_,
             const arma::umat& index_,
             const arma::umat& indexgr_,
             const arma::uvec& nvec_,
             const int& M_,
             const int& n_,
             const int& Kx_,
             const int& nparms_,
             const bool& hasX_,
             const bool& Print_) : 
  a(a_),
  dx(dx_),
  adx(adx_),
  d(d_),
  b(b_),
  index(index_),
  indexgr(indexgr_),
  nvec(nvec_),
  M(M_),
  n(n_),
  Kx(Kx_),
  nparms(nparms_),
  hasX(hasX_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec thetaa       = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    
    // int b1                 = 0;
    int b2                 = Kx - 1;
    int m1                 = b2 + 1;
    int m2                 = m1 + n - 1;
    
    arma::vec mu           = thetaa.subvec(m1, m2);
    arma::vec gd(nparms, arma::fill::zeros);
    arma::vec beta, dxb, adxb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dxb                  = dx*beta;
      adxb                 = a%dxb;
      llh                  = sum(adxb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, tmp, ai, exbmn, smunu;
    arma::mat dxi;
    arma::uvec indexi;
    arma::umat indexm;
    
    for (int m(0); m < M; ++ m) {
      igr1   = indexgr(m, 0); 
      igr2   = indexgr(m, 1); 
      nm     = nvec(m);
      indexm = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
      mum    = mu.subvec(igr1, igr2);
      
      for(int i(0); i < nm; ++ i){
        if(i < (nm - 1)){
          j1           = index(j, 0);
          j2           = index(j, 1);
          ai           = a.subvec(j1, j2);
          
          // sum of mu(i) + numj
          smunu        = mum(i) + mum.tail(nm - 1 - i);
          
          if(hasX){
            dxi          = dx.rows(j1, j2);
            exbmn        = exp(dxb.subvec(j1, j2) + smunu);
            tmp          = exbmn/(1 + exbmn);
            llh         += sum(ai%smunu - log(1 + exbmn));
            
            // grad X
            gd.head(Kx) -= arma::trans(arma::sum(dxi.each_col()%tmp, 0));
            
            // grad mui
            gd(m1 + j)   = b(j) - sum(tmp);
          } else {
            exbmn        = exp(smunu);
            tmp          = exbmn/(1 + exbmn);
            llh         += sum(ai%smunu - log(1 + exbmn));
            
            // grad mui
            gd(m1 + j)   = b(j) - sum(tmp);
          }
        }
        
        if(i > 0){
          indexi       = indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i);
          if(hasX){
            tmp        = exp(dxb.elem(indexi) + mum.head(i) + mum(i));
          } else {
            tmp        = exp(mum.head(i) + mum(i));
          }
          gd(m1 + j)  += d(j) - sum(tmp/(1 + tmp));
        }
        ++ j;
      }
    }
    
    grad               = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad               = gd;
    
    if(Print){
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "\nbeta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
    }
    
    return -llh;
  }
};

//[[Rcpp::export]]
List fhomobetasym(Eigen::VectorXd theta,
                  const arma::vec& a,
                  const arma::mat& dx,
                  const arma::uvec& nvec,
                  const arma::umat& index,
                  const arma::umat& indexgr,
                  const int& M,       
                  const int maxit = 300, 
                  const double& eps_f = 1e-6, 
                  const double& eps_g = 1e-5,
                  const bool& hasX = true,
                  const bool& Print = true){
  int n(sum(nvec)), Kx(0);
  arma::mat adx;
  if(hasX){
    Kx          = dx.n_cols;
    adx         = dx.each_col()%a;
  } 
  
  int nparms    = Kx + n;
  arma::vec d(n), b(n);
  int j(0);
  for (int m(0); m < M; ++ m) {
    int igr1          = indexgr(m, 0); // group m starts from igr1 in X
    int igr2          = indexgr(m, 1); // group m ends at igr2 in X
    int nm            = nvec(m);
    arma::umat indexm = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        int j1       = index(j, 0);
        int j2       = index(j, 1);
        b(j)         = sum(a.subvec(j1, j2));
      }
      
      // rows on which nui is used
      if(i > 0){
        arma::uvec indexi = indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i);
        d(j)              = sum(a.elem(indexi));
      }
      
      ++ j;
    }
  }
  
  llhhomosym f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms, hasX, Print);
  
  double fopt;
  int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}

// optimize only mu
class llhhomomusym: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::vec& dxb;
  const double& sadxb;
  const arma::vec& d;
  const arma::vec& b;
  const arma::umat& indexm;
  arma::vec& llh;
  const int& m;
  const int& nm;
  const int& nparms;
  const bool& Print;
public:
  llhhomomusym(const arma::vec& a_,
               const arma::vec& dxb_,
               const double& sadxb_,
               const arma::vec& d_,
               const arma::vec& b_,
               const arma::umat& indexm_,
               arma::vec& llh_,
               const int& m_,
               const int& nm_,
               const int& nparms_,
               const bool& Print_) : 
  a(a_),
  dxb(dxb_),
  sadxb(sadxb_),
  d(d_),
  b(b_),
  indexm(indexm_),
  llh(llh_),
  m(m_),
  nm(nm_),
  nparms(nparms_),
  Print(Print_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& theta, Refvec grad)
  {
    Eigen::VectorXd theta0 = theta;  //make a copy
    arma::vec mum          = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
    
    arma::vec gd(nparms, arma::fill::zeros);
    llh(m)                 = sadxb;
    
    int j1, j2;
    arma::vec mumj, tmp, ai, exbmn, smunu;
    arma::uvec indexi;
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        j1           = indexm(i, 0);
        j2           = indexm(i, 1);
        ai           = a.subvec(j1, j2);
        
        // sum of mu(i) + numj
        smunu        = mum(i) + mum.tail(nm - 1 - i);
        
        exbmn        = exp(dxb.subvec(j1, j2) + smunu);
        tmp          = exbmn/(1 + exbmn);
        llh(m)      += sum(ai%smunu - log(1 + exbmn));
        // grad mui
        gd(i)        = b(i) - sum(tmp);
      }
      
      if(i > 0){
        indexi       = indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i);
        tmp          = exp(dxb.elem(indexi) + mum.head(i) + mum(i));
        gd(i)       += d(i) - sum(tmp/(1 + tmp));
      }
    }
    
    grad             = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad             = gd;
    return -llh(m);
  }
};

//[[Rcpp::export]]
double fhomobetamusym(arma::vec& theta,
                      const arma::vec& a,
                      const arma::vec& dxb,
                      const arma::uvec& nvec,
                      const arma::umat& index,
                      const arma::umat& indexgr,
                      const int& N,
                      const int& M,       
                      const int& Kx,
                      const int maxit = 300, 
                      const double& eps_f = 1e-6, 
                      const double& eps_g = 1e-5,
                      const bool& hasX = true,
                      const bool& Print = true){
   int n(sum(nvec)), j(0);
  arma::vec adxb(N, arma::fill::zeros), mu(theta.subvec(Kx, Kx + n - 1)), d(n), b(n), mumj, 
  llh(M), musum(arma::zeros<arma::vec>(N));
  if(hasX){
    adxb        = a%dxb;
  } 
  
  for (int m(0); m < M; ++ m) {
    int nm(nvec(m)), igr1(indexgr(m, 0)), igr2(indexgr(m, 1));
    arma::umat indexm(index.rows(igr1, igr2)); // ith row is the row at each i interacts with others, where the link goes from i
    arma::vec mum(mu.subvec(igr1, igr2));
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        int j1       = index(j, 0);
        int j2       = index(j, 1);
        b(j)         = sum(a.subvec(j1, j2));
        // sum of mu(i) + numj
        musum.subvec(index(j, 0), index(j, 1)) = mum(i) + mum.tail(nm - 1 - i);
      }
      // rows on which nui is used
      if(i > 0){
        arma::uvec indexi = indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i);
        d(j)              = sum(a.elem(indexi));
      }
      ++ j;
    }
    arma::vec tp(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum.subvec(index(igr1, 0), index(igr2, 1)));
    llh(m)                = sum(a.subvec(index(igr1, 0), index(igr2, 1))%tp - log(1 + exp(tp)));
  }
  
  
  for (int m(0); m < M; ++ m) {
    int igr1(indexgr(m, 0)), igr2(indexgr(m, 1)), N1(index(igr1, 0)), N2(index(igr2, 1)), nm(nvec(m)), nparms(nm);
    arma::umat indexm(index.rows(igr1, igr2) - index(igr1, 0));
    arma::vec mum(mu.subvec(igr1, igr2)), am(a.subvec(N1, N2)), dxbm(dxb.subvec(N1, N2)), dm(d.subvec(igr1, igr2)), 
    bm(b.subvec(igr1, igr2));
    Eigen::VectorXd thetaEi(Eigen::Map<Eigen::VectorXd>(mum.memptr(), nparms));
    double fopt, adxbm(sum(adxb.subvec(N1, N2)));
    llhhomomusym f(am, dxbm, adxbm, dm, bm, indexm, llh, m, nm, nparms, Print);
    optim_lbfgs(f, thetaEi, fopt, maxit, eps_f, eps_g);
    llh(m)     = -fopt;
    mu.subvec(igr1, igr2) = arma::vec(thetaEi.data(), nm, false, false); 
    if(Print){
      Rcpp::Rcout << "group: " << m + 1<< " -- log-likelihood: " << sum(llh) << "\n";
    }
  }
  
  theta.subvec(Kx, Kx + n - 1)        = mu;
  return sum(llh);
}

// These functions compute likelihoods
double fllh2f(const arma::vec& a,
              const arma::vec& dxb,
              const arma::vec& mu,
              const arma::vec& nu,
              const arma::umat& index,
              const arma::umat& indexgr,
              const arma::uvec& nvec,
              const arma::uvec& Nvec,
              const int& M) {
  double llh(0);
  int igr1, igr2, nm;
  arma::umat indexm;
  arma::vec num, mum, mumj;
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0);
    igr2   = indexgr(m, 1);
    indexm = index.rows(igr1, igr2) - index(igr1, 0);
    nm     = nvec(m);
    mum    = mu.subvec(igr1, igr2);
    num    = arma::zeros(nm); num.head(nm - 1) = nu.subvec(igr1 - m, igr2 - m - 1);
    arma::vec musum(Nvec(m), arma::fill::zeros);
    for(int i(0); i < nm; ++ i){
      mumj = mum;   // all mu in his group
      mumj.shed_row(i); // all mu excepted mu(i)
      musum.subvec(indexm(i, 0), indexm(i, 1)) = num(i) + mumj;
    }
    arma::vec axb(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum);
    llh   += sum(a.subvec(index(igr1, 0), index(igr2, 1))%axb - log(1 + exp(axb)));
  }
  return llh;
}


double fllh1f(const arma::vec& a,
              const arma::vec& dxb,
              const arma::vec& mu,
              const arma::umat& index,
              const arma::umat& indexgr,
              const arma::uvec& nvec,
              const arma::uvec& Nvec,
              const int& M) {
  double llh(0);
  int igr1, igr2, nm;
  arma::umat indexm;
  arma::vec mum, mumj;
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0);
    igr2   = indexgr(m, 1);
    indexm = index.rows(igr1, igr2) - index(igr1, 0);
    nm     = nvec(m);
    mum    = mu.subvec(igr1, igr2);
    arma::vec musum(Nvec(m), arma::fill::zeros);
    for(int i(0); i < nm; ++ i){
      mumj = mum;   // all mu in his group
      mumj.shed_row(i); // all mu excepted mu(i)
      musum.subvec(indexm(i, 0), indexm(i, 1)) = mum(i) + mumj;
    }
    arma::vec axb(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum);
    llh   += sum(a.subvec(index(igr1, 0), index(igr2, 1))%axb - log(1 + exp(axb)));
  }
  return llh;
}


double fllhsym(const arma::vec& a,
              const arma::vec& dxb,
              const arma::vec& mu,
              const arma::umat& index,
              const arma::umat& indexgr,
              const arma::uvec& nvec,
              const arma::uvec& Nvec,
              const int& M) {
  double llh(0);
  int igr1, igr2, nm;
  arma::umat indexm;
  arma::vec mum, mumj;
  for (int m(0); m < M; ++ m) {
    igr1   = indexgr(m, 0);
    igr2   = indexgr(m, 1);
    indexm = index.rows(igr1, igr2) - index(igr1, 0);
    nm     = nvec(m);
    mum    = mu.subvec(igr1, igr2);
    arma::vec musum(Nvec(m), arma::fill::zeros);
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        musum.subvec(indexm(i, 0), indexm(i, 1)) = mum(i) + mum.tail(nm - 1 - i);
      }
    }
    arma::vec axb(dxb.subvec(index(igr1, 0), index(igr2, 1)) + musum);
    llh   += sum(a.subvec(index(igr1, 0), index(igr2, 1))%axb - log(1 + exp(axb)));
  }
  return llh;
}

// Block of Newton Raphson
// These functions compute Hessian and grandients for beta
void fHGbeta2f(arma::vec& grad,
               arma::mat& Hess,
               const arma::mat& dx,
               const arma::vec& a,
               const arma::vec& beta,
               const arma::vec& mu,
               const arma::vec& nu,
               const arma::uvec& nvec,
               const int& M,
               const int& N,
               const arma::umat& index,
               const arma::umat& indexgr){
  int igr1, igr2, nm, j(0);
  arma::vec mum, num, mumj, musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0);
    igr2 = indexgr(m, 1);
    nm   = nvec(m);
    mum  = mu.subvec(igr1, igr2);
    num  = arma::zeros(nm); num.head(nm - 1) = nu.subvec(igr1 - m, igr2 - m - 1);
    for(int i(0); i < nm; ++ i){
      mumj   = mum;   // all mu in his group
      mumj.shed_row(i); // all mu excepted mu(i)
      musum.subvec(index(j, 0), index(j, 1)) = num(i) + mumj;
      ++j;
    }
  }
  
  arma::vec p(1/(1 + exp(-dx*beta - musum)));
  arma::vec pq(p%(p - 1));
  grad = dx.t()*(a - p);
  Hess = arma::trans(dx.each_col()%pq)*dx;
}

void fHGbeta1f(arma::vec& grad,
               arma::mat& Hess,
               const arma::mat& dx,
               const arma::vec& a,
               const arma::vec& beta,
               const arma::vec& mu,
               const arma::uvec& nvec,
               const int& M,
               const int& N,
               const arma::umat& index,
               const arma::umat& indexgr){
  int igr1, igr2, nm, j(0);
  arma::vec mum, mumj, musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0);
    igr2 = indexgr(m, 1);
    nm   = nvec(m);
    mum  = mu.subvec(igr1, igr2);
    for(int i(0); i < nm; ++ i){
      mumj   = mum;   // all mu in his group
      mumj.shed_row(i); // all mu excepted mu(i)
      musum.subvec(index(j, 0), index(j, 1)) = mum(i) + mumj;
      ++j;
    }
  }
  
  arma::vec p(1/(1 + exp(-dx*beta - musum)));
  arma::vec pq(p%(p - 1));
  grad = dx.t()*(a - p);
  Hess = arma::trans(dx.each_col()%pq)*dx;
}

void fHGbetasym(arma::vec& grad,
                arma::mat& Hess,
                const arma::mat& dx,
                const arma::vec& a,
                const arma::vec& beta,
                const arma::vec& mu,
                const arma::uvec& nvec,
                const int& M,
                const int& N,
                const arma::umat& index,
                const arma::umat& indexgr){
  int igr1, igr2, nm, j(0);
  arma::vec mum, mumj, musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0);
    igr2 = indexgr(m, 1);
    nm   = nvec(m);
    mum  = mu.subvec(igr1, igr2);
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        musum.subvec(index(j, 0), index(j, 1)) = mum(i) + mum.tail(nm - 1 - i);
      }
      ++j;
    }
  }
  arma::vec p(1/(1 + exp(-dx*beta - musum)));
  arma::vec pq(p%(p - 1));
  grad = dx.t()*(a - p);
  Hess = arma::trans(dx.each_col()%pq)*dx;
}

void fGHmunu2fm(arma::vec& grad,
                arma::mat& Hess, 
                const arma::umat& indexm,
                const arma::vec& ampm,
                const arma::vec& pqm,
                const int& nm) {
  grad = arma::zeros<arma::vec>(2*nm - 1);
  Hess = arma::zeros<arma::mat>(2*nm - 1, 2*nm - 1);
  arma::uvec idj;
  for (int i(0); i < nm; ++ i) {
    int id1(indexm(i, 0)), id2(indexm(i, 1));
    idj              = indexm.col(0) + i;
    idj.head(i + 1) -= 1;
    idj.shed_row(i);
    
    if (i < (nm - 1)) {
      grad(i)    = sum(ampm.subvec(id1, id2));
      arma::uvec tp(arma::linspace<arma::uvec>(0, nm - 1, nm));
      tp.shed_row(i);
      arma::vec Hi(Hess.col(i));
      Hi.elem(tp + nm - 1) = pqm.subvec(id1, id2);
      Hess.col(i)          = Hi;
      Hess.row(i)          = Hi.t();
      Hess(i, i)           = sum(pqm.subvec(id1, id2));
    }
    grad(nm - 1 + i)             = sum(ampm.elem(idj));
    Hess(nm - 1 + i, nm - 1 + i) = sum(pqm.elem(idj));
  }
}

void fGHmu1fm(arma::vec& grad,
              arma::mat& Hess, 
              const arma::umat& indexm,
              const arma::vec& ampm,
              const arma::vec& pqm,
              const int& nm) {
  grad = arma::zeros<arma::vec>(nm);
  Hess = arma::zeros<arma::mat>(nm, nm);
  arma::uvec idj;
  for (int i(0); i < nm; ++ i) {
    int id1(indexm(i, 0)), id2(indexm(i, 1));
    idj              = indexm.col(0) + i;
    idj.head(i + 1) -= 1;
    idj.shed_row(i);
    
    grad(i)    = sum(ampm.subvec(id1, id2)) + sum(ampm.elem(idj));
    
    arma::uvec tp(arma::linspace<arma::uvec>(0, nm - 1, nm));
    tp.shed_row(i);
    arma::vec Hi(Hess.col(i));
    Hi.elem(tp) = pqm.subvec(id1, id2) + pqm.elem(idj);
    Hess.col(i) = Hi;
    Hess.row(i) = Hi.t();
    Hess(i, i)  = sum(pqm.subvec(id1, id2)) + sum(pqm.elem(idj));
  }
}

void fGHmusym(arma::vec& grad,
              arma::mat& Hess, 
              const arma::umat& indexm,
              const arma::vec& ampm,
              const arma::vec& pqm,
              const int& nm) {
  grad = arma::zeros<arma::vec>(nm);
  Hess = arma::zeros<arma::mat>(nm, nm);
  arma::uvec idj;
  double gi, Hii;
  
  for (int i(0); i < nm; ++ i) {
    arma::vec Hi(Hess.col(i));
    if (i < (nm - 1)) {
      int id1(indexm(i, 0)), id2(indexm(i, 1));
      gi           = sum(ampm.subvec(id1, id2));
      Hii          = sum(pqm.subvec(id1, id2));
      arma::uvec tp(arma::linspace<arma::uvec>(i + 1, nm - 1, nm - i - 1));
      Hi.elem(tp)  = pqm.subvec(id1, id2);
    }
    
    if (i > 0) {
      arma::uvec idj(indexm.col(0).head(i) + arma::linspace<arma::uvec>(i - 1, 0, i));
      gi          += sum(ampm.elem(idj));
      Hii         += sum(pqm.elem(idj));
      arma::uvec tp(arma::linspace<arma::uvec>(0, i - 1, i));
      Hi.elem(tp) += pqm.elem(idj);
    }
    grad(i)      = gi;
    Hess.col(i)  = Hi;
    Hess.row(i)  = Hi.t();
    Hess(i, i)   = Hii;
  }
}

// Newton Raphson
//[[Rcpp::export]]
Rcpp::List NewRaph2f(arma::vec& theta,
                     const arma::vec& a,
                     const arma::mat& dx,
                     const arma::uvec& nvec,
                     const arma::uvec& Nvec,
                     const arma::umat& index,
                     const arma::umat& indexgr,
                     const int& M,
                     const int& N,      
                     const bool& hasX = true,
                     const bool& Print = true,
                     const double& tol = 1e-4,
                     const int& maxit = 50) {
  int Kx(0), n = sum(nvec);
  if (hasX) {
    Kx = dx.n_cols;
  }
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1), n1(m2 + 1), n2(n1 + n - M - 1);
  arma::vec grad, beta(theta.head(Kx)), mu(theta.subvec(m1, m2)), nu(theta.subvec(n1, n2));
  arma::mat Hess;
  
  int igr1, igr2, nm, k(0);
  arma::vec mum, num, mumj, pm, pqm, ampm, dxb(N, arma::fill::zeros);
  arma::umat indexm;
  double llh, dist;
  bool cont(true);
  
  if (hasX) {
    dxb   = dx*beta;
  }
  
  while (cont) {
    arma::vec betap(beta), mup(mu), nup(nu);
    
    // mu, nu
    for (int m(0); m < M; ++ m) {
      igr1   = indexgr(m, 0);
      igr2   = indexgr(m, 1);
      indexm = index.rows(igr1, igr2) - index(igr1, 0);
      nm     = nvec(m);
      mum    = mu.subvec(igr1, igr2);
      num    = arma::zeros(nm); num.head(nm - 1) = nu.subvec(igr1 - m, igr2 - m - 1);
      arma::vec musum = arma::zeros<arma::vec>(Nvec(m));
      for(int i(0); i < nm; ++ i){
        mumj = mum;   // all mu in his group
        mumj.shed_row(i); // all mu excepted mu(i)
        musum.subvec(indexm(i, 0), indexm(i, 1)) = num(i) + mumj;
      }
      pm     = 1/(1 + exp(-dxb.subvec(index(igr1, 0), index(igr2, 1)) - musum));
      pqm    = pm%(pm - 1);
      ampm   = a.subvec(index(igr1, 0), index(igr2, 1)) - pm;
      fGHmunu2fm(grad, Hess, indexm, ampm, pqm, nm);
      arma::vec tp = arma::join_cols(nu.subvec(igr1 - m, igr2 - m - 1), mu.subvec(igr1, igr2)) - arma::solve(Hess, grad);
      nu.subvec(igr1 - m, igr2 - m - 1)  = tp.head(nm - 1);
      mu.subvec(igr1, igr2)              = tp.tail(nm);
    }
    
    // beta
    if (hasX) {
      fHGbeta2f(grad, Hess, dx, a, beta, mu, nu, nvec, M, N, index, indexgr);
      beta = betap - arma::solve(Hess, grad);
      dxb  = dx*beta;
      arma::vec tpdis{max(abs(beta - betap)), max(abs(mu - mup)), max(abs(nu - nup))};
      dist   = max(tpdis);
    } else {
      arma::vec tpdis{max(abs(mu - mup)), max(abs(nu - nup))};
      dist   = max(tpdis);
    }
    ++ k;
    cont   = ((dist >= tol) & (k < maxit));
    
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      llh      = fllh2f(a, dxb, mu, nu, index, indexgr, nvec, Nvec, M);
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  if (!Print) {
    llh        = fllh2f(a, dxb, mu, nu, index, indexgr, nvec, Nvec, M);
  }
  theta.head(Kx)       = beta; 
  theta.subvec(m1, m2) = mu; 
  theta.subvec(n1, n2) = nu;
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}


//[[Rcpp::export]]
Rcpp::List NewRaph1f(arma::vec& theta,
                     const arma::vec& a,
                     const arma::mat& dx,
                     const arma::uvec& nvec,
                     const arma::uvec& Nvec,
                     const arma::umat& index,
                     const arma::umat& indexgr,
                     const int& M,
                     const int& N,      
                     const bool& hasX = true,
                     const bool& Print = true,
                     const double& tol = 1e-4,
                     const int& maxit = 50) {
  int Kx(0), n = sum(nvec);
  if (hasX) {
    Kx = dx.n_cols;
  }
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1);
  arma::vec grad, beta(theta.head(Kx)), mu(theta.subvec(m1, m2));
  arma::mat Hess;
  
  int igr1, igr2, nm, k(0);
  arma::vec mum, mumj, pm, pqm, ampm, dxb(N, arma::fill::zeros);
  arma::umat indexm;
  double llh, dist;
  bool cont(true);
  
  if (hasX) {
    dxb   = dx*beta;
  }
  
  while (cont) {
    arma::vec betap(beta), mup(mu);
    
    // mu
    for (int m(0); m < M; ++ m) {
      igr1   = indexgr(m, 0);
      igr2   = indexgr(m, 1);
      indexm = index.rows(igr1, igr2) - index(igr1, 0);
      nm     = nvec(m);
      mum    = mu.subvec(igr1, igr2);
      arma::vec musum = arma::zeros<arma::vec>(Nvec(m));
      for(int i(0); i < nm; ++ i){
        mumj = mum;   // all mu in his group
        mumj.shed_row(i); // all mu excepted mu(i)
        musum.subvec(indexm(i, 0), indexm(i, 1)) = mum(i) + mumj;
      }
      pm     = 1/(1 + exp(-dxb.subvec(index(igr1, 0), index(igr2, 1)) - musum));
      pqm    = pm%(pm - 1);
      ampm   = a.subvec(index(igr1, 0), index(igr2, 1)) - pm;
      fGHmu1fm(grad, Hess, indexm, ampm, pqm, nm);
      mu.subvec(igr1, igr2) = mu.subvec(igr1, igr2) - arma::solve(Hess, grad);
    }
    
    // beta
    if (hasX) {
      fHGbeta1f(grad, Hess, dx, a, beta, mu, nvec, M, N, index, indexgr);
      beta = betap - arma::solve(Hess, grad);
      dxb  = dx*beta;
      arma::vec tpdis{max(abs(beta - betap)), max(abs(mu - mup))};
      dist   = max(tpdis);
    } else {
      dist   = max(abs(mu - mup));
    }
    ++ k;
    cont   = ((dist >= tol) & (k < maxit));
    
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      llh  = fllh1f(a, dxb, mu, index, indexgr, nvec, Nvec, M);
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  if (!Print) {
    llh    = fllh1f(a, dxb, mu, index, indexgr, nvec, Nvec, M);
  }
  theta.head(Kx)       = beta; 
  theta.subvec(m1, m2) = mu; 
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}

//[[Rcpp::export]]
Rcpp::List NewRaphsym(arma::vec& theta,
                      const arma::vec& a,
                      const arma::mat& dx,
                      const arma::uvec& nvec,
                      const arma::uvec& Nvec,
                      const arma::umat& index,
                      const arma::umat& indexgr,
                      const int& M,
                      const int& N,      
                      const bool& hasX = true,
                      const bool& Print = true,
                      const double& tol = 1e-4,
                      const int& maxit = 50) {
  int Kx(0), n = sum(nvec);
  if (hasX) {
    Kx = dx.n_cols;
  }
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1);
  arma::vec grad, beta(theta.head(Kx)), mu(theta.subvec(m1, m2));
  arma::mat Hess;
  
  int igr1, igr2, nm, k(0);
  arma::vec mum, mumj, pm, pqm, ampm, dxb(N, arma::fill::zeros);
  arma::umat indexm;
  double llh, dist;
  bool cont(true);
  
  if (hasX) {
    dxb   = dx*beta;
  }
  
  while (cont) {
    arma::vec betap(beta), mup(mu);

    // mu
    for (int m(0); m < M; ++ m) {
      igr1   = indexgr(m, 0);
      igr2   = indexgr(m, 1);
      indexm = index.rows(igr1, igr2) - index(igr1, 0);
      nm     = nvec(m);
      mum    = mu.subvec(igr1, igr2);
      arma::vec musum = arma::zeros<arma::vec>(Nvec(m));
      for(int i(0); i < nm; ++ i){
        if(i < (nm - 1)){
          musum.subvec(indexm(i, 0), indexm(i, 1)) = mum(i) + mum.tail(nm - 1 - i);
        }
      }
      pm     = 1/(1 + exp(-dxb.subvec(index(igr1, 0), index(igr2, 1)) - musum));
      pqm    = pm%(pm - 1);
      ampm   = a.subvec(index(igr1, 0), index(igr2, 1)) - pm;
      fGHmusym(grad, Hess, indexm, ampm, pqm, nm);
      mu.subvec(igr1, igr2) = mu.subvec(igr1, igr2) - arma::solve(Hess, grad);
    }
    
    // beta
    if (hasX) {
      fHGbetasym(grad, Hess, dx, a, beta, mu, nvec, M, N, index, indexgr);
      beta = betap - arma::solve(Hess, grad);
      dxb  = dx*beta;
      arma::vec tpdis{max(abs(beta - betap)), max(abs(mu - mup))};
      dist   = max(tpdis);
    } else {
      dist   = max(abs(mu - mup));
    }
    ++ k;
    cont   = ((dist >= tol) & (k < maxit));
    
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
      if(hasX){
        NumericVector betacpp  = wrap(beta);
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      llh  = fllhsym(a, dxb, mu, index, indexgr, nvec, Nvec, M);
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  if (!Print) {
    llh    = fllhsym(a, dxb, mu, index, indexgr, nvec, Nvec, M);
  }
  theta.head(Kx)       = beta; 
  theta.subvec(m1, m2) = mu; 
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}


// Mix of Newton Raphson and LBFGS
//[[Rcpp::export]]
Rcpp::List NewRaphLBFGS2f(arma::vec& theta,
                          const arma::vec& a,
                          const arma::mat& dx,
                          const arma::uvec& nvec,
                          const arma::uvec& Nvec,
                          const arma::umat& index,
                          const arma::umat& indexgr,
                          const int& M,
                          const int& N,      
                          const bool& hasX = true,
                          const bool& Print = true,
                          const double& tol = 1e-4,
                          const int& maxitNR = 50,
                          const int& maxitopt = 1e9,
                          const double& eps_f = 1e-6, 
                          const double& eps_g = 1e-5) {
  int Kx(0), n = sum(nvec);
  arma::vec grad, dxb(N, arma::fill::zeros);
  arma::mat Hess;
  if (hasX) {
    Kx  = dx.n_cols;
    dxb = dx*theta.head(Kx);
  }
  
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1), n1(m2 + 1), n2(n1 + n - M - 1), k(0);
  double llh, dist;
  bool cont(true);

  while (cont) {
    arma::vec thetap(theta);
    ++ k;
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
    }
    
    // mu, nu
    llh    = fhomobetamunu2f(theta, a, dxb, nvec, index, indexgr, N, M, Kx, maxitopt, eps_f, eps_g, hasX, Print);
    
    // beta
    arma::vec mu(theta.subvec(m1, m2)), nu(theta.subvec(n1, n2));
    if (hasX) {
      fHGbeta2f(grad, Hess, dx, a, theta.head(Kx), mu, nu, nvec, M, N, index, indexgr);
      theta.head(Kx) -= arma::solve(Hess, grad);
      dxb             = dx*theta.head(Kx);
    }
    
    dist   = max(abs(theta - thetap));
    cont   = ((dist >= tol) & (k < maxitNR));
    llh    = fllh2f(a, dxb, mu, nu, index, indexgr, nvec, Nvec, M);
    
    if (Print) {
      if(hasX){
        NumericVector betacpp  = wrap(theta.head(Kx));
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}


//[[Rcpp::export]]
Rcpp::List NewRaphLBFGS1f(arma::vec& theta,
                          const arma::vec& a,
                          const arma::mat& dx,
                          const arma::uvec& nvec,
                          const arma::uvec& Nvec,
                          const arma::umat& index,
                          const arma::umat& indexgr,
                          const int& M,
                          const int& N,      
                          const bool& hasX = true,
                          const bool& Print = true,
                          const double& tol = 1e-4,
                          const int& maxitNR = 50,
                          const int& maxitopt = 1e9,
                          const double& eps_f = 1e-6, 
                          const double& eps_g = 1e-5) {
  int Kx(0), n = sum(nvec);
  arma::vec grad, dxb(N, arma::fill::zeros);
  arma::mat Hess;
  if (hasX) {
    Kx = dx.n_cols;
    dxb   = dx*theta.head(Kx);
  }
  
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1), k(0);
  double llh, dist;
  bool cont(true);
  
  while (cont) {
    arma::vec thetap(theta);
    ++ k;
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
    }
    
    // mu
    llh    = fhomobetamu1f(theta, a, dxb, nvec, index, indexgr, N, M, Kx, maxitopt, eps_f, eps_g, hasX, Print);
    
    // beta
    arma::vec mu(theta.subvec(m1, m2));
    if (hasX) {
      fHGbeta1f(grad, Hess, dx, a, theta.head(Kx), mu, nvec, M, N, index, indexgr);
      theta.head(Kx) -= arma::solve(Hess, grad);
      dxb             = dx*theta.head(Kx);
    }

    dist   = max(abs(theta - thetap));
    cont   = ((dist >= tol) & (k < maxitNR));
    llh    = fllh1f(a, dxb, mu, index, indexgr, nvec, Nvec, M);
    
    if (Print) {
      if(hasX){
        NumericVector betacpp  = wrap(theta.head(Kx));
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}


//[[Rcpp::export]]
Rcpp::List NewRaphLBFGSsym(arma::vec& theta,
                           const arma::vec& a,
                           const arma::mat& dx,
                           const arma::uvec& nvec,
                           const arma::uvec& Nvec,
                           const arma::umat& index,
                           const arma::umat& indexgr,
                           const int& M,
                           const int& N,      
                           const bool& hasX = true,
                           const bool& Print = true,
                           const double& tol = 1e-4,
                           const int& maxitNR = 50,
                           const int& maxitopt = 1e9,
                           const double& eps_f = 1e-6, 
                           const double& eps_g = 1e-5) {
  int Kx(0), n = sum(nvec);
  arma::vec grad, dxb(N, arma::fill::zeros);
  arma::mat Hess;
  if (hasX) {
    Kx  = dx.n_cols;
    dxb = dx*theta.head(Kx);
  }
  
  int b2(Kx - 1), m1(b2 + 1), m2(m1 + n - 1), k(0);
  double llh, dist;
  bool cont(true);
  
  while (cont) {
    arma::vec thetap(theta);
    ++ k;
    if (Print) {
      Rcpp::Rcout << "\nIteration: " << k << "\n";
    }
    
    // mu 
    llh    = fhomobetamusym(theta, a, dxb, nvec, index, indexgr, N, M, Kx, maxitopt, eps_f, eps_g, hasX, Print);
    
    // beta
    arma::vec mu(theta.subvec(m1, m2));
    if (hasX) {
      fHGbetasym(grad, Hess, dx, a, theta.head(Kx), mu, nvec, M, N, index, indexgr);
      theta.head(Kx) -= arma::solve(Hess, grad);
      dxb             = dx*theta.head(Kx);
    }

    dist   = max(abs(theta - thetap));
    cont   = ((dist >= tol) & (k < maxitNR));
    llh    = fllhsym(a, dxb, mu, index, indexgr, nvec, Nvec, M);
    
    if (Print) {
      if(hasX){
        NumericVector betacpp  = wrap(theta.head(Kx));
        betacpp.attr("dim")    = R_NilValue;
        Rcpp::Rcout << "beta: \n";
        Rcpp::print(betacpp);
      }
      Rcpp::Rcout << "log-likelihood: " << llh << "\n";
      Rcpp::Rcout << "Distance: " << dist << "\n";
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("estimate")  = theta,
    Rcpp::Named("value")     = llh,
    Rcpp::Named("iteration") = k,
    Rcpp::Named("distance")  = dist);
}

// Change X from Full to Lower
//[[Rcpp::export]]
arma::mat hdataF2L(const arma::mat& data,
                   const arma::vec& nvec,
                   const arma::mat& index,
                   const int& M){
  int N(sum(nvec%(nvec - 1))/2), j(0), r(0);
  arma::mat out(N, data.n_cols);
  for (int m(0); m < M; ++ m) {
    int nm(nvec(m));
    for(int i(0); i < (nm - 1); ++ i){
      int j1 = index(j, 0) + i;
      int j2 = index(j, 1);
      out.rows(r, r + nm - 2 - i) = data.rows(j1, j2);
      ++ j;
      r      = r + nm - 1 - i;
    }
    ++ j;
  }
  return out;
}


// Change X from Full to Upper
//[[Rcpp::export]]
arma::mat hdataF2U(const arma::mat& data,
                   const arma::vec& nvec,
                   const arma::mat& index,
                   const arma::mat& indexgr,
                   const int& M){
  int N(sum(nvec%(nvec - 1))/2), j(0), r(0);
  arma::mat out(N, data.n_cols);
  for (int m(0); m < M; ++ m) {
    int nm(nvec(m)), igr1(indexgr(m, 0)), igr2(indexgr(m, 1));
    arma::mat indexm(index.rows(igr1, igr2)); 
    for(int i(0); i < (nm - 1); ++ i){
      arma::uvec indexi(arma::conv_to<arma::uvec>::from(indexm.col(0).tail(nm - 1 - i)) + i);
      out.rows(r, r + nm - 2 - i) = data.rows(indexi);
      ++ j;
      r      = r + nm - 1 - i;
    }
    ++ j;
  }
  return out;
}


// Change X from lower to Full
//[[Rcpp::export]]
arma::mat hdata2S(const arma::mat& data,
                  const arma::vec& nvec,
                  const arma::mat& index,
                  const arma::mat& indexgr,
                  const int& M){
  int N(sum(nvec%(nvec - 1))), j(0), r(0);
  arma::mat out(N, data.n_cols);
  for (int m(0); m < M; ++ m) {
    int igr1(indexgr(m, 0)), igr2(indexgr(m, 1)), nm(nvec(m));
    arma::mat indexm(index.rows(igr1, igr2)); 
    out.rows(r, r + nm - 2) = data.rows(index(j, 0), index(j, 1));
    ++ j;
    r   += nm - 1;
    for(int i(1); i < (nm - 1); ++ i){
      arma::uvec indexi      = arma::conv_to<arma::uvec>::from(indexm.col(0).head(i) + arma::linspace(i - 1, 0, i));
      out.rows(r, r + i - 1) = data.rows(indexi);
      r                               += i;
      out.rows(r, r + nm - i - 2) = data.rows(index(j, 0), index(j, 1));
      r                               += nm - i - 1;
      ++ j;
    }
    arma::uvec indexi            = arma::conv_to<arma::uvec>::from(indexm.col(0).head(nm - 1) + arma::linspace(nm - 2, 0, nm - 1));
    out.rows(r, r + nm - 2) = data.rows(indexi);
    r                               += nm - 1;
    ++ j;
  }
  return out;
}

// Functions to compute initial guess
//[[Rcpp::export]]
Rcpp::List finithasX(const arma::mat& dx,
                     const arma::vec& theta,
                     const arma::vec& a,
                     const bool& updateb = false,
                     const double& tol   = 1e-9,
                     const int& maxit    = 100) {
  int Kx(dx.n_cols), K(0);
  bool cont(true);
  double dist, llh;
  arma::vec xb, p, pq, amp, beta, betap;
  
  if (updateb) {
    beta    = arma::join_cols(arma::zeros(1), theta);
    arma::vec grad(Kx + 1, arma::fill::zeros);
    arma::mat hess(Kx + 1, Kx + 1, arma::fill::zeros);
    while (cont) {
      ++ K;
      betap = beta;
      xb    = betap(0) + dx*beta.tail(Kx);
      p     = 1/(1 + exp(-xb));
      pq    = p%(p - 1);
      amp   = a - p;
      grad.tail(Kx) = dx.t()*amp;
      grad(0)       = sum(amp);
      hess(0, 0)    = sum(pq);
      hess.submat(1, 0, Kx, 0)  = dx.t()*pq;
      hess.submat(0, 1, 0, Kx)  = arma::trans(hess.submat(1, 0, Kx, 0)); 
      hess.submat(1, 1, Kx, Kx) = arma::trans(dx.each_col()%pq)*dx;
      beta -= arma::solve(hess, grad);
      dist  = arma::max(abs(betap - beta));
      cont  = (dist >= tol) && (K < maxit);
    }
    xb      = betap(0) + dx*beta.tail(Kx);
  } else {
    beta    = arma::zeros(1);
    arma::vec xb0(dx*theta);
    while (cont) {
      ++ K;
      betap = beta;
      xb    = betap(0) + xb0;
      p     = 1/(1 + exp(-xb));
      pq    = p%(p - 1);
      amp   = a - p;
      beta -= sum(amp)/sum(pq);
      dist  = arma::max(abs(betap - beta));
      cont  = (dist >= tol) && (K < maxit);
    }
    xb      = betap(0) + xb0;
    beta    = arma::join_cols(beta, theta);
  }
  llh       = sum(a%xb - log(1 + exp(xb)));
  return Rcpp::List::create(_["theta"] = beta, _["llh"] = llh);
}

//[[Rcpp::export]]
Rcpp::List finit(const arma::vec& a,
                 const double& tol   = 1e-9,
                 const int& maxit    = 100) {
  int K(0), n(a.n_elem);
  bool cont(true);
  double dist, llh,  p, pq, beta(0), betap;
  arma::vec amp;
  
  while (cont) {
    ++ K;
    betap = beta;
    p     = 1/(1 + exp(-beta));
    pq    = p*(p - 1);
    amp   = a - p;
    beta -= sum(amp)/(n*pq);
    dist  = max(abs(betap - beta));
    cont  = (dist >= tol) && (K < maxit);
  }
  
  llh     = sum(beta*a - log(1 + exp(beta)));
  return Rcpp::List::create(_["theta"] = beta, _["llh"] = llh);
}
// 
// // same function where beta is printer
// // Estimation using fixed effects
// class llhhomop: public MFuncGrad
// {
// private:
//   const arma::vec& a;
//   const arma::mat& dx;
//   const arma::mat& adx;
//   const arma::vec& d;
//   const arma::vec& b;
//   const arma::mat& index;
//   const arma::mat& indexgr;
//   const arma::vec& nvec;
//   const int& M;
//   const int& n;
//   const int& Kx;
//   const int& nparms;
// public:
//   llhhomop(const arma::vec& a_,
//           const arma::mat& dx_,
//           const arma::mat& adx_,
//           const arma::vec& d_,
//           const arma::vec& b_,
//           const arma::mat& index_,
//           const arma::mat& indexgr_,
//           const arma::vec& nvec_,
//           const int& M_,
//           const int& n_,
//           const int& Kx_,
//           const int& nparms_) : 
//   a(a_),
//   dx(dx_),
//   adx(adx_),
//   d(d_),
//   b(b_),
//   index(index_),
//   indexgr(indexgr_),
//   nvec(nvec_),
//   M(M_),
//   n(n_),
//   Kx(Kx_),
//   nparms(nparms_){}
//   
//   arma::vec Grad;
//   
//   double f_grad(Constvec& theta, Refvec grad)
//   {
//     Eigen::VectorXd theta0 = theta;  //make a copy
//     arma::vec thetaa       = arma::vec(theta0.data(), theta0.size(), false, false); //converte into arma vec
//     
//     // int b1                 = 0;
//     int b2                 = Kx - 1;
//     int m1                 = b2 + 1;
//     int m2                 = m1 + n - 1;
//     int n1                 = m2 + 1;
//     int n2                 = n1 + n - M - 1;
//     
//     arma::vec beta         = thetaa.head(Kx);
//     arma::vec mu           = thetaa.subvec(m1, m2);
//     arma::vec nu           = thetaa.subvec(n1, n2);
//     // cout<<nu.t()<<endl;
//     NumericVector betacpp  = wrap(beta);
//     betacpp.attr("dim")    = R_NilValue;
//     Rcpp::Rcout << "beta: \n";
//     Rcpp::print(betacpp);
//     
//     arma::vec dxb          = dx*beta;
//     arma::vec adxb         = a%dxb;
//     double llh             = sum(adxb);
//     
//     arma::vec gd(nparms, arma::fill::zeros);
//     gd.head(Kx)            = arma::trans(sum(adx, 0));
//     
//     int igr1, igr2, nm, j(0), j1, j2;
//     arma::vec mum, num, numj, mumj, tmp, ai, exbmn, smunu;
//     arma::mat indexm, dxi;
//     arma::uvec indexi;
//     
//     for (int m(0); m < M; ++ m) {
//       igr1                 = indexgr(m, 0); 
//       igr2                 = indexgr(m, 1); 
//       nm                   = nvec(m);
//       indexm               = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
//       mum                  = mu.subvec(igr1, igr2);
//       num                  = nu.subvec(igr1 - m, igr2 - m - 1);
//       num                  = arma::join_cols(num, arma::zeros(1));
//       
//       for(int i(0); i < nm; ++ i){
//         j1                 = index(j, 0);
//         j2                 = index(j, 1);
//         ai                 = a.subvec(j1, j2);
//         dxi                = dx.rows(j1, j2);
//         
//         // nuj when mui is fixed
//         numj               = num;
//         numj.shed_row(i);
//         
//         // muj when nui is fixed
//         mumj               = mum;
//         mumj.shed_row(i);
//         
//         // rows on which nui is used
//         indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
//         indexi.head(i + 1)-= 1;
//         indexi.shed_row(i);  
//         
//         // sum of mu(i) + numj
//         smunu              = mum(i) + numj;
//         
//         exbmn              = exp(dxb.subvec(j1, j2) + smunu);
//         tmp                = exbmn/(1 + exbmn);
//         llh               += sum(ai%smunu - log(1 + exbmn));
//         
//         // grad X
//         gd.head(Kx)       -= arma::trans(arma::sum(dxi.each_col()%tmp, 0));
//         
//         // grad mui
//         gd(m1 + j)         = d(j) - sum(tmp);
//         // gd(m1 + j)         = sum(ai - tmp);
//         // grad nui
//         if(i < (nm - 1)){
//           tmp              = exp(dxb.elem(indexi) + mumj + num(i));
//           gd(n1 + j - m)   = b(j) - sum(tmp/(1 + tmp));
//           // gd(n1 + j - m)   = sum(a.elem(indexi)- tmp/(1 + tmp));
//         }
//         ++ j;
//       }
//     }
//     
//     grad                    = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
//     Grad                    = gd;
//     
//     Rcpp::Rcout << "log-likelihood: " << llh << "\n";
//     return -llh;
//   }
// };
// 
// //[[Rcpp::export]]
// List fhomobetap(Eigen::VectorXd theta,
//                const arma::vec& a,
//                const arma::mat& dx,
//                const arma::vec& nvec,
//                const arma::mat& index,
//                const arma::mat& indexgr,
//                const int& M,       
//                const int& maxit = 300, 
//                const double& eps_f = 1e-15, 
//                const double& eps_g = 1e-15){
// 
//   int n         = sum(nvec);
//   int Kx        = dx.n_cols;
//   
//   int nparms    = Kx + 2*n - M;
//   arma::mat adx = dx.each_col()%a;
//   
//   arma::vec d(n), b(n);
//   int j(0);
//   for (int m(0); m < M; ++ m) {
//     int igr1              = indexgr(m, 0); // group m starts from igr1 in X
//     int igr2              = indexgr(m, 1); // group m ends at igr2 in X
//     int nm                = nvec(m);
//     arma::mat indexm      = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
//     
//     for(int i(0); i < nm; ++ i){
//       int j1              = index(j, 0);
//       int j2              = index(j, 1);
//       d(j)                = sum(a.subvec(j1, j2));
//       
//       // rows on which nui is used
//       arma::uvec indexi   = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
//       indexi.head(i + 1) -= 1;
//       indexi.shed_row(i);
//       b(j)                = sum(a.elem(indexi));
//       ++ j;
//     }
//   }
// 
//   llhhomop f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms);
//   
//   double fopt;
//   int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
//   
//   return Rcpp::List::create(
//     Rcpp::Named("estimate") = theta,
//     Rcpp::Named("value")    = fopt,
//     Rcpp::Named("gradient") = f.Grad,
//     Rcpp::Named("status")   = status);
// }

// /////////////////////////////
// // Linear model compute X´X with fixed effects
// List Object_LinearFE(const arma::vec& a,
//                      const arma::mat& dx,
//                      const arma::vec& nvec,
//                      const arma::mat& index,
//                      const arma::mat& indexgr,
//                      const int& Kx,
//                      const int& n,
//                      const int& N,
//                      const int& M){
// 
//   int b2                 = Kx - 1;
//   int m1                 = b2 + 1;
//   int m2                 = m1 + n - 1;
//   int n1                 = m2 + 1;
//   int n2                 = n1 + n - M - 1;
//   
//   
//   int igr1, igr2, nm, j(0), j1, j2;
//   arma::vec ai, Xpa(N, arma::fill::zeros);
//   arma::mat indexm, dxi, XpX(n2 + 1, n2 + 1, arma::fill::zeros);
//   arma::uvec indexi;
//   
//   XpX.submat(0, 0, b2, b2) = dx.t()*dx;
//   Xpa.subvec(0, b2)        = arma::trans(sum(dx, 0))
//   for (int m(0); m < M; ++ m) {
//     igr1                   = indexgr(m, 0); 
//     igr2                   = indexgr(m, 1); 
//     nm                     = nvec(m);
//     indexm                 = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
// 
//     for(int i(0); i < nm; ++ i){
//       j1                   = index(j, 0);
//       j2                   = index(j, 1);
//       ai                   = a.subvec(j1, j2);
//       dxi                  = dx.rows(j1, j2);
//       XpX(m1 + j, m1 + j)  =         
//       // rows on which nui is used
//       indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
//       indexi.head(i + 1)-= 1;
//       indexi.shed_row(i);  
//       
//       // sum of mu(i) + numj
//       smunu              = mum(i) + numj;
//       
//       exbmn              = exp(dxb.subvec(j1, j2) + smunu);
//       tmp                = exbmn/(1 + exbmn);
//       llh               += sum(ai%smunu - log(1 + exbmn));
//       
//       // grad X
//       gd.head(Kx)       -= arma::trans(arma::sum(dxi.each_col()%tmp, 0));
//       
//       // grad mui
//       gd(m1 + j)         = d(j) - sum(tmp);
//       // gd(m1 + j)         = sum(ai - tmp);
//       // grad nui
//       if(i < (nm - 1)){
//         tmp              = exp(dxb.elem(indexi) + mumj + num(i));
//         gd(n1 + j - m)   = b(j) - sum(tmp/(1 + tmp));
//         // gd(n1 + j - m)   = sum(a.elem(indexi)- tmp/(1 + tmp));
//       }
//       ++ j;
//     }
//   }
//   
//   grad                    = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
//   Grad                    = gd;
//   
//   Rcpp::Rcout << "log-likelihood: " << llh << "\n";
//   return -llh;
// }
