// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#define NDEBUG 1
#include <RcppEigen.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// This function simulate a G from a given dnetwork
//[[Rcpp::export]]
List simG (List& dnetwork,
           const arma::vec& N,
           const int& M) {
  List out(M);
    for (int m(0); m < M; ++ m) {
      int Nm        = N(m);
      arma::mat dnm = dnetwork(m);
      arma::mat matunif(Nm, Nm, arma::fill::randu);
      //out[m]        = arma::normalise(conv_to<mat>::from((matunif < dnm)),1,1);
      dnm           = arma::conv_to<mat>::from((matunif < dnm));

      dnm.diag()    = arma::zeros(Nm);
      out[m]        = dnm;
    }

    return out;
}


//[[Rcpp::export]]
List simGnorm (List& dnetwork,
              const arma::vec& N,
              const int& M) {
  List out(M);
    for (int m(0); m < M; ++ m) {
      int Nm        = N(m);
      arma::mat dnm = dnetwork(m);
      arma::mat matunif(Nm, Nm, arma::fill::randu);
      dnm           = arma::normalise(arma::conv_to<mat>::from((matunif < dnm)),1,1);
      //out[m]        = arma::conv_to<mat>::from((matunif < dnm));

      dnm.diag()    = arma::zeros(Nm);
      out[m]        = dnm;
    }

    return out;
}

// This function normalizes network
// [[Rcpp::export]]
List fGnormalise(List& u, const double& M) {
  List out(M);

  for(int m(0); m < M; ++m) {
    arma::mat um = u[m];
    um           = arma::normalise(um, 1, 1);
    out[m]       = um;
  }
  return out;
}


// Remove NA from network data
//[[Rcpp::export]]
List rem_non_fin (const arma::mat& net) {
  arma::mat out = net;
  int n         = net.n_rows;
  arma::uvec id = arma::regspace<arma::uvec>(0, n - 1);
  arma::umat R  = arma::repmat(id, 1, n);
  arma::umat C  = arma::repmat(id.t(), n, 1);
  arma::uvec tp = arma::find_nan(out);
  while(tp.n_elem > 0){
    R             = R.elem(tp);
    C             = C.elem(tp);
    arma::uvec z  = arma::join_cols(R, C);
    arma::uvec zu = arma::sort(arma::unique(z));
    int nzu       = zu.n_elem;
    arma::vec czu(nzu);
    for(int s(0); s < nzu; ++ s){
      czu(s)      = sum(z == zu(s));
    }
    arma::uvec x  = zu.elem(arma::find(czu == czu.max()));
    id.shed_rows(x);
    out           = net.rows(id);
    out           = out.cols(id);
    tp            = arma::find_nan(out);
    n             = out.n_rows;
    arma::umat y  = arma::regspace<arma::uvec>(0, n - 1);
    R             = arma::repmat(y, 1, n);
    C             = arma::repmat(y.t(), n, 1);
  }
  return List::create(Named("net") = out, Named("id") = id + 1);
}

// Create a list a square matrixes from a given vector and sizes.
// The size of the m-th matrice in the list is N[m] with zeros on the diagonal
// The elements in the generated matrix are placed column-wise (ie. the first column is filled up before filling the second column)
// [[Rcpp::export]]
List frVtoM(const Eigen::VectorXd& u,
            const Rcpp::IntegerVector& N,
            const double& M) {
  List out(M);

  int r                                = 0;
  int n;

  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);

    n                                  = Nm - 1;

    Eigen::MatrixXd outm(Eigen::MatrixXd::Zero(Nm, Nm));
    outm.block(1, 0, n, 1)             = u.segment(r, n);

    r                                 += n;
    for(int i(1); i < n; ++i) {
      outm.block(0, i, i, 1)          = u.segment(r, i);
      outm.block(i + 1, i, n - i, 1)  = u.segment(r + i, n - i);
      r                              += n;
    }

    outm.block(0, n, n, 1)            = u.segment(r, n);
    r                                += n;

    out[m]                            = outm;
  }
  return out;
}


// Same function but the returned matrix are normalized
// [[Rcpp::export]]
List frVtoMnorm(const arma::vec& u,
                const IntegerVector& N,
                const double& M) {
  List out(M);

  int r2                               = -1;
  int r1;

  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);

    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;

    arma::mat outm(Nm, Nm, arma::fill::zeros);
    outm.submat(1, 0, Nm - 1, 0)       = u.subvec(r1, r2);

    for(int i(1); i < (Nm - 1); ++i) {
      r2                              += Nm - 1;
      r1                               = r2 - Nm + 2;
      outm.submat(0, i, i - 1, i)      = u.subvec(r1, r1 + i - 1);
      outm.submat(i + 1, i, Nm - 1, i) = u.subvec(r1 + i, r2);
    }

    r2                                += Nm - 1;
    r1                                 = r2 - Nm + 2;
    outm.submat(0, Nm - 1, Nm - 2, Nm - 1) = u.subvec(r1, r2);

    outm                                   = arma::normalise(outm, 1, 1);

    out[m]                                 = outm;
  }
  return out;
}


// Create a vector from a given list a square matrixes
// The size of the length of the vector is the sum(N), where N is the vector of matrice sizes
// The elements in the generated vector are taken from column-wise (ie. the first column is filled up before filling the second column)
// and from the first matrix of the list to the last matrix of the list.
// [[Rcpp::export]]
Eigen::VectorXd frMtoV(List& u,
                       const Rcpp::IntegerVector& N,
                       const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);

  int r                                = 0;
  int n;

  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    //um                                 = um.array().ceil();

    n                                  = Nm - 1;

    out.segment(r, n)                  = um.block(1, 0, n, 1);
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }

    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}

// same function but the matrixes are ceiled first
// [[Rcpp::export]]
Eigen::VectorXd frMceiltoV(List& u,
                           const Rcpp::IntegerVector& N,
                           const double& M) {
  int sN                               = sum(N*N - N);
  Eigen::VectorXd out(sN);

  int r                                = 0;
  int n;

  for(int m(0); m < M; ++m) {
    int Nm                             = N(m);
    Eigen::MatrixXd um                 = u[m];
    um                                 = um.array().ceil();

    n                                  = Nm - 1;

    out.segment(r, n)                  = um.block(1, 0, n, 1);
    r                                 += n;
    for(int i(1); i < n; ++i) {
      out.segment(r, i)                = um.block(0, i, i, 1);
      out.segment(r + i, n - i)        = um.block(i + 1, i, n - i, 1);
      r                               += n;
    }

    out.segment(r, n)                  = um.block(0, n, n, 1);
    r                                 += n;
  }
  return out;
}
