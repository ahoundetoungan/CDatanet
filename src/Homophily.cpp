// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppDist)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <mvnorm.h>
#include <wishart.h>

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


// compute mu_i + nu_j for i and j in the same group
// export is a vector as (1, 2), (1, 3), ... , (1, n), (2, 1) (2, 3), ....
//[[Rcpp::export]]
arma::vec fmusum(const arma::vec& mu,
                 const arma::vec& nu,
                 const arma::mat& index,
                 const arma::mat& indexgr,
                 const int& M,
                 const int& N) {
  int igr1, igr2, nm, j = 0, j1, j2;
  double  nui;
  arma::vec mum, num, mumj; 
  
  arma::vec musum(N, arma::fill::zeros);
  for (int m(0); m < M; ++ m) {
    igr1 = indexgr(m, 0); 
    igr2 = indexgr(m, 1);
    nm   = igr2 - igr1 + 1;
    mum  = mu.subvec(igr1, igr2);
    num  = nu.subvec(igr1, igr2);
    
    for(int i(0); i < nm; ++ i){
      nui    = num(i);// individual i, nu
      mumj   = mum;   // all mu in his group
      mumj.shed_row(i); // all nu excepted nu(i)
      j1     = index(j, 0); //starting individual i observation in dX
      j2     = index(j, 1); //ending individual i observation in dX
      musum.subvec(j1, j2) = nui + mumj;
      ++j;
    }
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
                const arma::mat& index,
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
      
      // elements which vary when nuj is fixed
      mumj    = mum;
      mumj.shed_row(i);
      j1      = index(j, 0);
      j2      = index(j, 1);
      tmpnu   = sum(astmdxbeta.subvec(j1, j2) - mumj);
      
      // elements which vary when mui is fixed
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
  
  // normalize mu and nu to mean zero
  if(nfix == 1){
    double mub       = sum(mu)/n;
    double nub       = sum(nu)/n;
    mu              -= mub;
    nu              -= nub;
    beta(0)         += (mub + nub);
    dxbeta           = dx*beta.tail(Kx);
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
  mupnu                       = fmusum(mu, nu, index, indexgr, M, N);
}


void updateusigma(double& smu2,
                  double& snu2,
                  double& rho,
                  const int& n,
                  const arma::vec& mu,
                  const arma::vec& nu) {
  
  arma::mat d   = arma::join_rows(mu, nu);
  
  arma::mat tmp = riwish(n, d.t()*d);
  smu2          = tmp(0, 0);
  snu2          = tmp(1, 1);
  rho           = tmp(0, 1)/sqrt(smu2*snu2);
}



//[[Rcpp::export]]
List updategparms1(const arma::vec& a,
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
                   const int& iteration){
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
  
  // output
  arma::mat Smu(n, iteration), Snu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration), Ssnu2(iteration), Srho(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  Ssnu2.attr("dim") = R_NilValue;
  Srho.attr("dim")  = R_NilValue;
  
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
    updateusigma(smu2, snu2, rho, n, mu, nu);
    
    //save
    Smu.col(t)    = mu;
    Snu.col(t)    = nu;
    Sbeta.col(t)  = beta;
    Ssmu2(t)      = smu2;
    Ssnu2(t)      = snu2;
    Srho(t)       = rho;
  }
  return List::create(Named("beta")      = Sbeta.t(),
                      Named("mu")        = Smu.t(),
                      Named("nu")        = Snu.t(),
                      Named("sigma2_mu") = Ssmu2,
                      Named("sigma2_nu") = Ssnu2,
                      Named("rho")       = Srho);
}


//[[Rcpp::export]]
List updategparms2(const arma::vec& a,
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
                   const int& iteration){
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
  arma::vec ast(N, arma::fill::zeros), astmdxbeta;
  
  // output
  arma::mat Smu(n, iteration), Snu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration), Ssnu2(iteration), Srho(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  Ssnu2.attr("dim") = R_NilValue;
  Srho.attr("dim")  = R_NilValue;

  for(int t(0); t < iteration; ++ t) {
    //update ast
    updateast(ast, dxbeta, mupnu, a, N);
    
    //update beta
    updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupnu, ast);
    
    //update mu and nu
    updatemunu(mu, nu, mupnu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2, snu2, rho);
    
    //update sigmas
    updateusigma(smu2, snu2, rho, n, mu, nu);
    
    //save
    Smu.col(t)    = mu;
    Snu.col(t)    = nu;
    Sbeta.col(t)  = beta;
    Ssmu2(t)      = smu2;
    Ssnu2(t)      = snu2;
    Srho(t)       = rho;
  }
  return List::create(Named("beta")      = Sbeta.t(),
                      Named("mu")        = Smu.t(),
                      Named("nu")        = Snu.t(),
                      Named("sigma2_mu") = Ssmu2,
                      Named("sigma2_nu") = Ssnu2,
                      Named("rho")       = Srho);
}
