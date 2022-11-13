// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppDist, RcppEigen, RcppNumerical)]]
#include <RcppArmadillo.h>
#define NDEBUG
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

void updatellh(double& llh,
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
  arma::mat Sigma;
  // output
  arma::mat Smu(n, iteration), Snu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration), Ssnu2(iteration), Srho(iteration); //, Sllh(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  Ssnu2.attr("dim") = R_NilValue;
  Srho.attr("dim")  = R_NilValue;
  //Sllh.attr("dim")  = R_NilValue;
  
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
    updateusigma(smu2, snu2, rho, Sigma, n, mu, nu);
    
    //likelihood
    //updatellh(llh, mu, nu, a, dxbeta, mupnu, Sigma, n);
    
    //save
    Smu.col(t)    = mu;
    Snu.col(t)    = nu;
    Sbeta.col(t)  = beta;
    Ssmu2(t)      = smu2;
    Ssnu2(t)      = snu2;
    Srho(t)       = rho;
    //Sllh(t)       = llh;
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
  arma::vec ast(N, arma::fill::zeros);
  arma::mat Sigma;
  // output
  arma::mat Smu(n, iteration), Snu(n, iteration), Sbeta(K, iteration);
  NumericVector Ssmu2(iteration), Ssnu2(iteration), Srho(iteration); //, Sllh(iteration);
  Ssmu2.attr("dim") = R_NilValue;
  Ssnu2.attr("dim") = R_NilValue;
  Srho.attr("dim")  = R_NilValue;
  //Sllh.attr("dim")  = R_NilValue;
  
  //Progress p(iteration, true);
  for(int t(0); t < iteration; ++ t) {
    //p.increment(); 
    
    //update ast
    updateast(ast, dxbeta, mupnu, a, N);
    
    //update beta
    updatebeta(beta, dxbeta, INDEXgr, nfix, Kx, dx, invdxdx, mupnu, ast);
    
    //update mu and nu
    updatemunu(mu, nu, mupnu, beta, dxbeta, dx, ast, Kx, INDEXgr, M, N, n, nfix, nvec, index, indexgr, smu2, snu2, rho);
    
    //update sigmas
    updateusigma(smu2, snu2, rho, Sigma, n, mu, nu);
    
    //likelihood
    //updatellh(llh, mu, nu, a, dxbeta, mupnu, Sigma, n);
    
    //save
    Smu.col(t)    = mu;
    Snu.col(t)    = nu;
    Sbeta.col(t)  = beta;
    Ssmu2(t)      = smu2;
    Ssnu2(t)      = snu2;
    Srho(t)       = rho;
    //Sllh(t)       = llh;
  }
  return List::create(Named("beta")      = Sbeta.t(),
                      Named("mu")        = Smu.t(),
                      Named("nu")        = Snu.t(),
                      Named("sigma2_mu") = Ssmu2,
                      Named("sigma2_nu") = Ssnu2,
                      Named("rho")       = Srho);//,
  //Named("loglike")   = Sllh);
}



// Estimation using fixed effects
class llhhomo: public MFuncGrad
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
public:
  llhhomo(const arma::vec& a_,
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
          const int& nparms_) : 
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
  nparms(nparms_){}
  
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
    
    arma::vec beta         = thetaa.head(Kx);
    arma::vec mu           = thetaa.subvec(m1, m2);
    arma::vec nu           = thetaa.subvec(n1, n2);
    // cout<< "Beta: "<<endl;
    // cout<< beta.t() <<endl;
    
    arma::vec dXb          = dx*beta;
    arma::vec adXb         = a%dXb;
    double llh             = sum(adXb);
    
    arma::vec gd(nparms, arma::fill::zeros);
    gd.head(Kx)            = arma::trans(sum(adx, 0));
    
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
        dXi                = dx.rows(j1, j2);
        
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
        
        // sum of mu(i) + numj
        smunu              = mum(i) + numj;
        
        exbmn              = exp(dXb.subvec(j1, j2) + smunu);
        tmp                = exbmn/(1 + exbmn);
        llh               += sum(ai%smunu - log(1 + exbmn));
        
        // grad X
        gd.head(Kx)       -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
        
        // grad mui
        gd(m1 + j)         = d(j) - sum(tmp);
        
        // grad nui
        if(i < (nm - 1)){
          tmp              = exp(dXb.elem(indexi) + mumj + num(i));
          gd(n1 + j - m)   = b(j) - sum(tmp/(1 + tmp));
        }
        ++ j;
      }
    }
    
    grad                   = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                   = gd;
    
    // cout<< llh <<endl;
    return -llh;
  }
};

//[[Rcpp::export]]
List fhomobeta(Eigen::VectorXd theta,
               const arma::vec& a,
               const arma::mat& dx,
               const arma::vec& nvec,
               const arma::mat& index,
               const arma::mat& indexgr,
               const int& M,       
               const int maxit = 300, 
               const double& eps_f = 1e-6, 
               const double& eps_g = 1e-5){
  int n         = sum(nvec);
  int Kx        = dx.n_cols;
  
  int nparms    = Kx + 2*n - M;
  arma::mat adx = dx.each_col()%a;
  
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
      d(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      b(j)                = sum(a.elem(indexi));
      ++ j;
    }
  }
  
  llhhomo f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms);
  
  double fopt;
  int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}


// same function where beta is printer
// Estimation using fixed effects
class llhhomop: public MFuncGrad
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
public:
  llhhomop(const arma::vec& a_,
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
          const int& nparms_) : 
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
  nparms(nparms_){}
  
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
    
    arma::vec beta         = thetaa.head(Kx);
    arma::vec mu           = thetaa.subvec(m1, m2);
    arma::vec nu           = thetaa.subvec(n1, n2);
    // cout<<nu.t()<<endl;
    NumericVector betacpp  = wrap(beta);
    betacpp.attr("dim")    = R_NilValue;
    Rcpp::Rcout << "beta: \n";
    Rcpp::print(betacpp);
    
    arma::vec dXb          = dx*beta;
    arma::vec adXb         = a%dXb;
    double llh             = sum(adXb);
    
    arma::vec gd(nparms, arma::fill::zeros);
    gd.head(Kx)            = arma::trans(sum(adx, 0));
    
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
        dXi                = dx.rows(j1, j2);
        
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
        
        // sum of mu(i) + numj
        smunu              = mum(i) + numj;
        
        exbmn              = exp(dXb.subvec(j1, j2) + smunu);
        tmp                = exbmn/(1 + exbmn);
        llh               += sum(ai%smunu - log(1 + exbmn));
        
        // grad X
        gd.head(Kx)       -= arma::trans(arma::sum(dXi.each_col()%tmp, 0));
        
        // grad mui
        gd(m1 + j)         = d(j) - sum(tmp);
        // gd(m1 + j)         = sum(ai - tmp);
        // grad nui
        if(i < (nm - 1)){
          tmp              = exp(dXb.elem(indexi) + mumj + num(i));
          gd(n1 + j - m)   = b(j) - sum(tmp/(1 + tmp));
          // gd(n1 + j - m)   = sum(a.elem(indexi)- tmp/(1 + tmp));
        }
        ++ j;
      }
    }
    
    grad                    = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad                    = gd;
    
    Rcpp::Rcout << "log-likelihood: " << llh << "\n";
    return -llh;
  }
};

//[[Rcpp::export]]
List fhomobetap(Eigen::VectorXd theta,
               const arma::vec& a,
               const arma::mat& dx,
               const arma::vec& nvec,
               const arma::mat& index,
               const arma::mat& indexgr,
               const int& M,       
               const int& maxit = 300, 
               const double& eps_f = 1e-15, 
               const double& eps_g = 1e-15){

  int n         = sum(nvec);
  int Kx        = dx.n_cols;
  
  int nparms    = Kx + 2*n - M;
  arma::mat adx = dx.each_col()%a;
  
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
      d(j)                = sum(a.subvec(j1, j2));
      
      // rows on which nui is used
      arma::uvec indexi   = arma::conv_to<arma::uvec>::from(indexm.col(0)) + i;
      indexi.head(i + 1) -= 1;
      indexi.shed_row(i);
      b(j)                = sum(a.elem(indexi));
      ++ j;
    }
  }

  llhhomop f(a, dx, adx, d, b, index, indexgr, nvec, M, n, Kx, nparms);
  
  double fopt;
  int status = optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = theta,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}
