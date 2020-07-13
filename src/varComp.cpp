// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List cABC(const int& n,
          const arma::mat& miq,
          const double& sigma,
          const int& R,
          const int& S,
          const arma::vec& tu,
          const arma::vec& tutu,
          const arma::mat& cpdf,
          const arma::mat& ccdf) {
  
  arma::vec A(arma::zeros(n)), B(arma::zeros(n)), C(arma::zeros(n)), d(arma::zeros(n)), b(arma::zeros(n)), m2d(arma::zeros(n));
  double sigma2 = sigma*sigma;
  
  for(int i(0); i < n; ++ i) {
    // positive index or negative index
    // in index p and index n xf(x) - (x -1)f(x -1) are negative
    // but in index p the num of C is positive and negetive in index n
    // in index 0 it is positive and the computation can be done without using log transform
    // moreover in index the first is treated a part
    bool cont = true;
    int r1(1), r2(1), indp(0), ind0, indn, lenr;
    double x1, x2, ltm, lf1, lf2, f1, f2, lfm, lm2fm, lF1, lF2, F1, F2, lFm, ldFir, lSAm, lSBm, lSCm, lSbm;
    arma::vec lt;
    
    while(cont) {
      if (miq(i, r1) > sigma) {
        r1         += 1;
        r2         += 1;
      } else{
        if (abs(miq(i, r2)) < sigma) {
          r2       += 1;
        } else {
          cont      = false;
        }
      }
    }
    if(r1 >= 2) {
      indp           = r1 - 1;
      ind0           = r2;
    } 
    ind0             = r2;
  
    
    // compute A, B, C, d, b for positive index p1
    lenr             = indp;
    if (lenr > 0) {
      arma::vec lSAi(lenr), lSBi(lenr), lSCi(lenr), lSbi(lenr - 1);
      
      // ABC
      int j          = 0;
      
      for (int r(0); r < indp; ++ r) {
        x1           = miq(i, r);
        x2           = miq(i,1 + r);
           
        lt           = -0.5*(x2*x2 + 2*x2*tu + tutu)/sigma2 - 0.5*log(2*acos(-1)*sigma2);
        
        ltm          = max(lt);
        
        ldFir        = ltm + log(sum(exp(lt - ltm))) - log(S);
        
        lf1          = cpdf(i, r);
        lf2          = cpdf(i, r + 1);
        lfm          = lf2;
        
        double tmp1  = lfm + log(1 - exp(lf1 - lf2));
        double tmp2  = lfm + log(x2 -  x1*exp(lf1 - lf2));
        
        lSAi(j)      = 2*tmp1 - ldFir;
        lSBi(j)      = 2*tmp2 - ldFir;
        lSCi(j)      = tmp1 + tmp2 - ldFir;
        j           += 1;
        
      }
      
      lSAm           = max(lSAi);
      lSBm           = max(lSBi);
      lSCm           = max(lSCi);
      
      A(i)          += exp(lSAm + log(sum(exp(lSAi - lSAm))));
      B(i)          += exp(lSBm + log(sum(exp(lSBi - lSBm))));
      C(i)          += exp(lSCm + log(sum(exp(lSCi - lSCm))));
      
      // b
      if (lenr > 1) {
        int j        = 0;
        
        for (int r(1); r < indp; ++ r) {
          x1         = miq(i, r);
          lf1        = cpdf(i, r);
          lSbi(j)    = log(x1) + lf1;
          j         += 1;
        }
        lSbm         = max(lSbi);
        b(i)        += exp(lSbm + log(sum(exp(lSbi - lSbm))));
      }
    }
    
    // compute A, B, C, d, b for positive index 0
    lenr             = ind0 - indp;
    if (lenr > 0){
      for (int r(indp); r < ind0; ++ r) {
        x1           = miq(i, r); 
        x2           = miq(i, r + 1);
    
        f1           = exp(cpdf(i, r)); 
        f2           = exp(cpdf(i, 1 + r)); 
        F1           = exp(ccdf(i, r)); 
        F2           = exp(ccdf(i, 1 + r)); 

        A(i)        += pow(f1 - f2, 2)/(F1 - F2); 
        B(i)        += pow(x1*f1 - x2*f2, 2)/(F1 - F2); 
        C(i)        += pow(x1*f1 - x2*f2, 2)/(F1 - F2); 
        b(i)        += x1*f1;
      }
    }
    
    // compute A, B, C, d, b for positive index n
    lenr             = R + 1 - ind0;
    if (lenr > 0){
      arma::vec lSAi(lenr), lSBi(lenr), lSCi(lenr), lSbi(lenr);
      int j = 0;
      
      for (int r(ind0); r <= R; ++ r) {
        x1           = miq(i, r); 
        x2           = miq(i, r + 1); 
        
        lf1          = cpdf(i, r); 
        lf2          = cpdf(i, r + 1); 
        lfm          = lf1; 
        
        lF1          = ccdf(i, r); 
        lF2          = ccdf(i, r + 1);
        lFm          = lF1; 
        
        double tmp1  = 2*lfm - lFm; 
        double tmp2  = log(1 - exp(lf2 - lfm)); 
        double tmp3  = log(1 - exp(lF2 - lFm)); 
        double tmp4  = log(x2*exp(lf2 - lfm) - x1); 
        
        lSAi(j)      = tmp1 + 2*tmp2 - tmp3; 
        lSBi(j)      = tmp1 + 2*tmp4 - tmp3; 
        lSCi(j)      = tmp1 + tmp4 + tmp2 - tmp3; 
        lSbi(j)      = log(-x1) + lf1; 
        j           += 1;
      }
      
      lSAm           = max(lSAi);
      lSBm           = max(lSBi);
      lSCm           = max(lSCi);
      lSbm           = max(lSbi);
      
      A(i)          += exp(lSAm + log(sum(exp(lSAi - lSAm))));
      B(i)          += exp(lSBm + log(sum(exp(lSBi - lSBm))));
      C(i)          -= exp(lSCm + log(sum(exp(lSCi - lSCm))));
      b(i)          -= exp(lSbm + log(sum(exp(lSbi - lSbm))));
    }
    
    //compute d m2d
    arma::rowvec tp1 = cpdf.submat(i, 0, i, R); 
    arma::rowvec tp2 = log(pow(miq.submat(i, 0, i, R), 2));
    lfm              = max(tp1); 
    lm2fm            = max(tp1 + tp2);
    d(i)             = exp(lfm + log(sum(exp(tp1 - lfm)))); 
    m2d(i)           = exp(lm2fm + log(sum(exp(tp1 + tp2 - lm2fm)))); 
  }
  
  
  B                 /=  sigma2; 
  C                 /= -sigma;       
  b                 /= -sigma;              

  return List::create(Named("A") = A, Named("B") = B, 
                      Named("C") = C, Named("d") = d, 
                      Named("b") = b, Named("m2d") = m2d);
}