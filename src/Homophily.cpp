// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppDist, RcppEigen, RcppNumerical)]]
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
                 const arma::mat& index, // where each j starts and ends in a
                 const arma::mat& indexgr, // where each group starts and ends in mu
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
      j1     = index(j, 0); //starting individual i observation in dX
      j2     = index(j, 1); //ending individual i observation in dX
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
                 const arma::mat& index,  // where each j starts and ends in a
                 const arma::mat& indexgr,
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
      j1     = index(j, 0); //starting individual i observation in dX
      j2     = index(j, 1); //ending individual i observation in dX
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
                const arma::mat& INDEXgr,
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
                const arma::mat& INDEXgr,
                const int& M,
                const int& N,
                const int& n,
                const int& nfix,
                const arma::vec& nvec,
                const arma::mat& index,  // where each j starts and ends in a
                const arma::mat& indexgr,
                const double& smu2,
                const double& snu2,
                const double& rho){
  int j1, j2, igr1, igr2, nm, j = 0;
  double tmpmu, tmpnu;
  arma::vec num, mum, numj, mumj, munub, munu; 
  arma::uvec indexj;
  arma::mat indexm, smunu(2, 2), vmunu(2, 2), invvmunu;
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
      indexj  = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
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
                const arma::mat& INDEXgr, 
                const int& M,
                const int& N,
                const int& n,
                const int& nfix,
                const arma::vec& nvec,
                const arma::mat& index,
                const arma::mat& indexgr, // where each j starts and ends in a
                const double& smu2){
  int j1, j2, igr1, igr2, nm, j = 0;
  double tmp, smu, mub, invvmu;
  arma::vec mum, mumj; 
  arma::uvec indexj;
  arma::mat indexm;
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
      indexj  = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
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
              const arma::mat& INDEXgr, 
              const int& M,
              const int& N,
              const int& n,
              const int& nfix,
              const arma::vec& nvec,
              const arma::mat& index,
              const arma::mat& indexgr, // where each j starts and ends in a
              const double& smu2){
  int j1, j2, igr1, igr2, nm, j = 0;
  double smu, mub, invvmu;
  arma::vec mum, mumj; 
  arma::uvec indexj;
  arma::mat indexm;
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
        indexj  = arma::conv_to<arma::uvec>::from(indexm.col(0).head(i) + arma::linspace(i - 1, 0, i));
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
               const arma::mat& index,
               const arma::mat& indexgr,
               const arma::mat& INDEXgr,
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
               const arma::mat& index,
               const arma::mat& indexgr,
               const arma::mat& INDEXgr,
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


// Estimation using fixed effects two sides
class llhhomo2f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::mat& index;
  const arma::mat& indexgr;
  const arma::vec& nvec;
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
          const arma::mat& index_,
          const arma::mat& indexgr_,
          const arma::vec& nvec_,
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
    
    arma::vec beta, dXb, adXb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dXb                  = dx*beta;
      adXb                 = a%dXb;
      llh                  = sum(adXb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, num, numj, mumj, tmp, ai, exbmn, smunu;
    arma::mat indexm, dXi;
    arma::uvec indexi;
    
    for (int m(0); m < M; ++ m) {
      igr1                 = indexgr(m, 0); 
      igr2                 = indexgr(m, 1); 
      nm                   = nvec(m);
      indexm               = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
      mum                  = mu.subvec(igr1, igr2);
      num                  = nu.subvec(igr1 - m, igr2 - m - 1);
      num                  = arma::join_cols(num, arma::zeros(1));
      
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
        indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
        indexi.head(i + 1)-= 1;
        indexi.shed_row(i);  
        
        // sum of muj + num(i)
        smunu              = mumj + num(i);
        
        if(hasX){
          exbmn            = exp(dXb.subvec(j1, j2) + smunu);
          tmp              = exbmn/(1 + exbmn);
          llh             += sum(ai%smunu - log(1 + exbmn));
          
          // grad X
          dXi              = dx.rows(j1, j2);
          gd.head(Kx)     -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
          
          // grad nui
          if(i < (nm - 1)){
            gd(n1 + j - m) = b(j) - sum(tmp);
          }
          
          // grad mui
          tmp              = exp(dXb.elem(indexi) + mum(i) + numj);
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
        Rcpp::Rcout << "beta: \n";
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
               const arma::vec& nvec,
               const arma::mat& index,
               const arma::mat& indexgr,
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
    arma::mat indexm      = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
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

// Estimation using fixed effects only mu, one way
class llhhomo1f: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::mat& index;
  const arma::mat& indexgr;
  const arma::vec& nvec;
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
          const arma::mat& index_,
          const arma::mat& indexgr_,
          const arma::vec& nvec_,
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
    arma::vec beta, dXb, adXb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dXb                  = dx*beta;
      adXb                 = a%dXb;
      llh                  = sum(adXb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, mumj, tmp, ai, exbmn, smunu;
    arma::mat indexm, dXi;
    arma::uvec indexi;
    
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
          dXi                = dx.rows(j1, j2);
          exbmn              = exp(dXb.subvec(j1, j2) + smunu);
          tmp                = exbmn/(1 + exbmn);
          llh               += sum(ai%smunu - log(1 + exbmn));
          
          // grad X
          gd.head(Kx)       -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
          
          // grad mui
          gd(m1 + j)         = b(j) - sum(tmp);
          indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
          indexi.head(i + 1)-= 1;
          indexi.shed_row(i); 
          tmp                = exp(dXb.elem(indexi) + mumj + mum(i));
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
        Rcpp::Rcout << "beta: \n";
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
               const arma::vec& nvec,
               const arma::mat& index,
               const arma::mat& indexgr,
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
    arma::mat indexm      = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      int j1              = index(j, 0);
      int j2              = index(j, 1);
      b(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
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

// Estimation using fixed effects with symmetric networks
class llhhomosym: public MFuncGrad
{
private:
  const arma::vec& a;
  const arma::mat& dx;
  const arma::mat& adx;
  const arma::vec& d;
  const arma::vec& b;
  const arma::mat& index;
  const arma::mat& indexgr;
  const arma::vec& nvec;
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
            const arma::mat& index_,
            const arma::mat& indexgr_,
            const arma::vec& nvec_,
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
    arma::vec beta, dXb, adXb;
    double llh(0);
    if(hasX){
      beta                 = thetaa.head(Kx);
      dXb                  = dx*beta;
      adXb                 = a%dXb;
      llh                  = sum(adXb);
      gd.head(Kx)          = arma::trans(sum(adx, 0));
    }
    
    int igr1, igr2, nm, j(0), j1, j2;
    arma::vec mum, tmp, ai, exbmn, smunu;
    arma::mat indexm, dXi;
    arma::uvec indexi;
    
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
            dXi          = dx.rows(j1, j2);
            exbmn        = exp(dXb.subvec(j1, j2) + smunu);
            tmp          = exbmn/(1 + exbmn);
            llh         += sum(ai%smunu - log(1 + exbmn));
            
            // grad X
            gd.head(Kx) -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
            
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
          indexi       = arma::conv_to<arma::uvec>::from(indexm.col(0).head(i) + arma::linspace(i - 1, 0, i));
          if(hasX){
            tmp        = exp(dXb.elem(indexi) + mum.head(i) + mum(i));
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
        Rcpp::Rcout << "beta: \n";
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
                 const arma::vec& nvec,
                 const arma::mat& index,
                 const arma::mat& indexgr,
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
    int igr1         = indexgr(m, 0); // group m starts from igr1 in X
    int igr2         = indexgr(m, 1); // group m ends at igr2 in X
    int nm           = nvec(m);
    arma::mat indexm = index.rows(igr1, igr2); // ith row is the row at each i interacts with others, where the link goes from i
    
    for(int i(0); i < nm; ++ i){
      if(i < (nm - 1)){
        int j1       = index(j, 0);
        int j2       = index(j, 1);
        b(j)         = sum(a.subvec(j1, j2));
      }
      
      // rows on which nui is used
      if(i > 0){
        arma::uvec indexi = arma::conv_to<arma::uvec>::from(indexm.col(0).head(i) + arma::linspace(i - 1, 0, i));
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
//     arma::vec dXb          = dx*beta;
//     arma::vec adXb         = a%dXb;
//     double llh             = sum(adXb);
//     
//     arma::vec gd(nparms, arma::fill::zeros);
//     gd.head(Kx)            = arma::trans(sum(adx, 0));
//     
//     int igr1, igr2, nm, j(0), j1, j2;
//     arma::vec mum, num, numj, mumj, tmp, ai, exbmn, smunu;
//     arma::mat indexm, dXi;
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
//         dXi                = dx.rows(j1, j2);
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
//         exbmn              = exp(dXb.subvec(j1, j2) + smunu);
//         tmp                = exbmn/(1 + exbmn);
//         llh               += sum(ai%smunu - log(1 + exbmn));
//         
//         // grad X
//         gd.head(Kx)       -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
//         
//         // grad mui
//         gd(m1 + j)         = d(j) - sum(tmp);
//         // gd(m1 + j)         = sum(ai - tmp);
//         // grad nui
//         if(i < (nm - 1)){
//           tmp              = exp(dXb.elem(indexi) + mumj + num(i));
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
//   arma::mat indexm, dXi, XpX(n2 + 1, n2 + 1, arma::fill::zeros);
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
//       dXi                  = dx.rows(j1, j2);
//       XpX(m1 + j, m1 + j)  =         
//       // rows on which nui is used
//       indexi             = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
//       indexi.head(i + 1)-= 1;
//       indexi.shed_row(i);  
//       
//       // sum of mu(i) + numj
//       smunu              = mum(i) + numj;
//       
//       exbmn              = exp(dXb.subvec(j1, j2) + smunu);
//       tmp                = exbmn/(1 + exbmn);
//       llh               += sum(ai%smunu - log(1 + exbmn));
//       
//       // grad X
//       gd.head(Kx)       -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
//       
//       // grad mui
//       gd(m1 + j)         = d(j) - sum(tmp);
//       // gd(m1 + j)         = sum(ai - tmp);
//       // grad nui
//       if(i < (nm - 1)){
//         tmp              = exp(dXb.elem(indexi) + mumj + num(i));
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
