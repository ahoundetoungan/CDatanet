#' @title Simulating data from linear-in-mean models with social interactions
#' @param formula a class object \link[stats]{formula}: a symbolic description of the model. `formula` must be as, for example, \code{y ~ x1 + x2 + gx1 + gx2}
#' where `y` is the endogenous vector and `x1`, `x2`, `gx1` and `gx2` are control variables, which can include contextual variables, i.e. averages among the peers.
#' Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist The network matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' @param theta a vector defining the true value of \eqn{\theta = (\lambda, \Gamma, \sigma)} (see the model specification in details). 
#' @param cinfo a Boolean indicating whether information is complete (`cinfo = TRUE`) or incomplete (`cinfo = FALSE`). In the case of incomplete information, the model is defined under rational expectations. 
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `simsar` is called.
#' @description
#' `simsar` simulates continuous variables with social interactions (see Lee, 2004 and Lee et al., 2010). 
#' @references  
#' Lee, L. F. (2004). Asymptotic distributions of quasi-maximum likelihood estimators for spatial autoregressive models. \emph{Econometrica}, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}.
#' @references  
#' Lee, L. F., Liu, X., & Lin, X. (2010). Specification and estimation of social interaction models with network structures. The Econometrics Journal, 13(2), 145-176, \doi{10.1111/j.1368-423X.2010.00310.x}
#' @details 
#' For a complete information model, the outcome \eqn{y_i} is defined as:
#' \deqn{y_i = \lambda \bar{y}_i + \mathbf{z}_i'\Gamma + \epsilon_i,}
#' where \eqn{\bar{y}_i} is the average of \eqn{y} among peers, 
#' \eqn{\mathbf{z}_i} is a vector of control variables, 
#' and \eqn{\epsilon_i \sim N(0, \sigma^2)}. 
#' In the case of incomplete information models with rational expectations, \eqn{y_i} is defined as:
#' \deqn{y_i = \lambda E(\bar{y}_i) + \mathbf{z}_i'\Gamma + \epsilon_i.}
#' @seealso \code{\link{sar}}, \code{\link{simsart}}, \code{\link{simcdnet}}.
#' @return A list consisting of:
#'     \item{y}{the observed count data.}
#'     \item{Gy}{the average of y among friends.}
#' @examples 
#' \donttest{
#' # Groups' size
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 1000))
#' n      <- sum(nvec)
#' 
#' # Parameters
#' lambda <- 0.4
#' Gamma  <- c(2, -1.9, 0.8, 1.5, -1.2)
#' sigma  <- 1.5
#' theta  <- c(lambda, Gamma, sigma)
#' 
#' # X
#' X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))
#' 
#' # Network
#' G      <- list()
#' 
#' for (m in 1:M) {
#'   nm           <- nvec[m]
#'   Gm           <- matrix(0, nm, nm)
#'   max_d        <- 30
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
#'     Gm[i, tmp] <- 1
#'   }
#'   rs           <- rowSums(Gm); rs[rs == 0] <- 1
#'   Gm           <- Gm/rs
#'   G[[m]]       <- Gm
#' }
#' 
#' # data
#' data   <- data.frame(X, peer.avg(G, cbind(x1 = X[,1], x2 =  X[,2])))
#' colnames(data) <- c("x1", "x2", "gx1", "gx2")
#' 
#' ytmp    <- simsar(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, 
#'                   theta = theta, data = data) 
#' y       <- ytmp$y
#' }
#' @importFrom Rcpp sourceCpp
#' @export
simsar   <- function(formula,
                     Glist,
                     theta,
                     cinfo = TRUE,
                     data) {
  stopifnot(length(cinfo) == 1)
  stopifnot(cinfo %in% c(TRUE, FALSE))
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  
  f.t.data <- formula.to.data(formula = formula, contextual = FALSE, Glist = Glist, M = M, igr = igr, 
                              data = data, type = "sim", theta0  = 0)
  X        <- f.t.data$X
  # print(colnames(X))
  K        <- length(theta)
  if(K != (ncol(X) + 2)) {
    stop("Length of theta is not suited.")
  }
  lambda   <- theta[1]
  b        <- theta[2:(K - 1)]
  sigma    <- theta[K ]
  
  
  
  xb       <- c(X %*% b)
  eps      <- rnorm(n, 0, sigma)
  
  out      <- NULL
  if(cinfo){
    y      <- rep(0, n)
    Gy     <- numeric(n)
    fySar(y, Gy, Glist, eps, igr, M, xb, lambda)
    out    <- list("y"  = y,
                   "Gy" = Gy)
  } else {
    y      <- rep(0, n)
    ye     <- rep(0, n)
    Gye    <- numeric(n)
    fySarRE(y, Gye, ye, Glist, eps, igr, M, xb, lambda)
    out    <- list("y"  = y,
                   "Gy" = Gy)
  }
  out
}


#' @title Estimating linear-in-mean models with social interactions
#' @param formula a class object \link[stats]{formula}: a symbolic description of the model. `formula` must be as, for example, \code{y ~ x1 + x2 + gx1 + gx2}
#' where `y` is the endogenous vector and `x1`, `x2`, `gx1` and `gx2` are control variables, which can include contextual variables, i.e. averages among the peers.
#' Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist The network matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' @param lambda0 an optional starting value of \eqn{\lambda}.
#' @param fixed.effects a Boolean indicating whether group heterogeneity must be included as fixed effects.
#' @param optimizer is either `nlm` (referring to the function \link[stats]{nlm}) or `optim` (referring to the function \link[stats]{optim}). 
#' Arguments for these functions such as, `control` and `method` can be set via the argument `opt.ctr`.
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as `control`, `method`, etc.
#' @param print a Boolean indicating if the estimate should be printed at each step.
#' @param cov a Boolean indicating if the covariance should be computed.
#' @param cinfo a Boolean indicating whether information is complete (`cinfo = TRUE`) or incomplete (`cinfo = FALSE`). In the case of incomplete information, the model is defined under rational expectations. 
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `sar` is called.
#' @description
#' `sar` computes quasi-maximum likelihood estimators for linear-in-mean models with social interactions (see Lee, 2004 and Lee et al., 2010). 
#' @details 
#' For a complete information model, the outcome \eqn{y_i} is defined as:
#' \deqn{y_i = \lambda \bar{y}_i + \mathbf{z}_i'\Gamma + \epsilon_i,}
#' where \eqn{\bar{y}_i} is the average of \eqn{y} among peers, 
#' \eqn{\mathbf{z}_i} is a vector of control variables, 
#' and \eqn{\epsilon_i \sim N(0, \sigma^2)}. 
#' In the case of incomplete information models with rational expectations, \eqn{y_i} is defined as:
#' \deqn{y_i = \lambda E(\bar{y}_i) + \mathbf{z}_i'\Gamma + \epsilon_i.}
#' @references  
#' Lee, L. F. (2004). Asymptotic distributions of quasi-maximum likelihood estimators for spatial autoregressive models. \emph{Econometrica}, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}.
#' @references  
#' Lee, L. F., Liu, X., & Lin, X. (2010). Specification and estimation of social interaction models with network structures. The Econometrics Journal, 13(2), 145-176, \doi{10.1111/j.1368-423X.2010.00310.x}
#' @seealso \code{\link{sart}}, \code{\link{cdnet}}, \code{\link{simsar}}.
#' @examples 
#' \donttest{
#' # Groups' size
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 1000))
#' n      <- sum(nvec)
#' 
#' # Parameters
#' lambda <- 0.4
#' Gamma  <- c(2, -1.9, 0.8, 1.5, -1.2)
#' sigma  <- 1.5
#' theta  <- c(lambda, Gamma, sigma)
#' 
#' # X
#' X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))
#' 
#' # Network
#' G      <- list()
#' 
#' for (m in 1:M) {
#'   nm           <- nvec[m]
#'   Gm           <- matrix(0, nm, nm)
#'   max_d        <- 30
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
#'     Gm[i, tmp] <- 1
#'   }
#'   rs           <- rowSums(Gm); rs[rs == 0] <- 1
#'   Gm           <- Gm/rs
#'   G[[m]]       <- Gm
#' }
#' 
#' # data
#' data   <- data.frame(X, peer.avg(G, cbind(x1 = X[,1], x2 =  X[,2])))
#' colnames(data) <- c("x1", "x2", "gx1", "gx2")
#' 
#' ytmp    <- simsar(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, 
#'                   theta = theta, data = data) 
#' data$y  <- ytmp$y
#' 
#' out     <- sar(formula = y ~ x1 + x2 + + gx1 + gx2, Glist = G, 
#'                optimizer = "optim", data = data)
#' summary(out)
#' }
#' @return A list consisting of:
#'     \item{info}{list of general information on the model.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{details}{outputs as returned by the optimizer.}
#' @export
sar <- function(formula,
                Glist, 
                lambda0       = NULL, 
                fixed.effects = FALSE,
                optimizer     = "optim",
                opt.ctr       = list(), 
                print         = TRUE, 
                cov           = TRUE,
                cinfo         = TRUE,
                data) {
  stopifnot(length(cinfo) == 1)
  stopifnot(cinfo %in% c(TRUE, FALSE))
  if(!cinfo) stop("Incomplete information is not supported in this version.")
  contextual  <- FALSE
  stopifnot(optimizer %in% c("optim", "nlm"))
  env.formula <- environment(formula)
  #size 
  if (!is.list(Glist)) {
    Glist    <- list(Glist)
  }
  M          <- length(Glist)
  nvec       <- unlist(lapply(Glist, nrow))
  n          <- sum(nvec)
  igr        <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)

  f.t.data   <- formula.to.data(formula, contextual, Glist, M, igr, data, 
                                theta0 = NULL, fixed.effects = fixed.effects)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  Gy         <- f.t.data$Gy
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X), "sigma")
  
  K          <- ncol(X)

  # variables
  Nvec       <- sapply(Glist, nrow)
  
  XX         <- t(X)%*%X
  invXX      <- solve(XX)
  
  Ilist      <- lapply(1:M, function(w) diag(Nvec[w]))
  
  lambdat    <- NULL
  if (!is.null(lambda0)) {
    lambdat  <- log(lambda0/(1- lambda0))
  } else {
    Xtmp     <- cbind(f.t.data$Gy, X)
    b        <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    lambdat  <- log(max(b[1]/(1 - b[1]), 0.01))
  }
  
  
  # arguments
  if ((length(opt.ctr) == 0) & optimizer == "optim") {
    opt.ctr  <- list("method" = "Brent",
                     "upper"  = 37,
                     "lower"  = -710)
  }
  ctr        <- c(list(X = X,invXX = invXX, G = Glist, I = Ilist, n = n, y = y, Gy = Gy, 
                       ngroup = M, FE = fixed.effects, print = print), opt.ctr)
  
  if (optimizer == "optim") {
    ctr      <- c(ctr, list(par = lambdat, fn = foptimSAR))
    par0     <- "par"
    par1     <- "par"
    like     <- "value"
  } else {
    ctr     <- c(ctr, list(p = lambdat, f = foptimSAR))
    par0    <- "p"
    par1    <- "estimate"
    like    <- "minimum"
  }
  
  resSAR    <- do.call(get(optimizer), ctr)
  lambdat   <- resSAR[[par1]]
  llh       <- -resSAR[[like]]
  
  lambda    <- 1/(1 + exp(-lambdat))
  
  hessian   <- jacobSAR(lambda, X, invXX, XX, y, n, Glist, Ilist, Gy, M, igr, cov, fixed.effects)
  
  
  beta      <- hessian$beta
  sigma2    <- hessian$sigma2
  
  covout    <- NULL
  if(cov) {
    covout           <- hessian$cov
    colnames(covout) <- coln
    rownames(covout) <- coln
  }
  
  theta             <- c(lambda, beta, sqrt(sigma2))
  names(theta)      <- coln
  
  environment(formula) <- env.formula
  sdata                <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist)),
    "nfriends"      = unlist(lapply(Glist, function(u) sum(u > 0))) 
  )
  if (!missing(data)) {
    sdata              <- c(sdata, list("data" = deparse(substitute(data))))
  }  
  
  
  INFO                 <- list("M"             = M,
                               "n"             = n,
                               "fixed.effects" = fixed.effects,
                               "nlinks"        = unlist(lapply(Glist, function(u) sum(u > 0))),
                               "formula"       = formula,
                               "log.like"      = llh)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = theta, 
                               "cov"        = covout,
                               "details"    = resSAR)
  class(out)           <- "sar"
  out
}

jacobSAR <- function(alpha, X, invXX, XX, y, N, G, I, Gy, ngroup, igroup, 
                     cov = TRUE, fixed.effects = FALSE) {
  Ay           <- (y - alpha * Gy)
  beta         <- invXX %*% (t(X) %*% Ay)
  Xbeta        <- X %*% beta
  s2           <- sum((Ay - Xbeta) ^ 2) / ifelse(fixed.effects, N - ngroup, N)
  nbeta        <- length(beta)
  
  covout       <- NULL
  if (cov) {
    covout   <- fSARjac(alpha, s2, X, XX, Xbeta, G, I, igroup, ngroup, N, 
                        nbeta, fixed.effects)
  }
  
  
  list("alpha"   = alpha,
       "beta"    = beta,
       "sigma2"  = s2,
       "cov"     = covout)
}


#' @title Summary for the estimation of linear-in-mean models with social interactions
#' @description Summary and print methods for the class `sar` as returned by the function \link{sar}.
#' @param object an object of class `sar`, output of the function \code{\link{sar}}.
#' @param x an object of class `summary.sar`, output of the function \code{\link{summary.sar}} or 
#' class `sar`, output of the function \code{\link{sar}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
"summary.sar" <- function(object,
                          ...) {
  stopifnot(class(object) == "sar")
  out           <- c(object, list("..."       = ...)) 
  if(is.null(object$cov)){
    stop("Covariance was not computed")
  }
  class(out)    <- "summary.sar"
  out
}


#' @rdname summary.sar
#' @export
"print.summary.sar"  <- function(x, ...) {
  stopifnot(class(x) == "summary.sar")
  
  M                    <- x$info$M
  n                    <- x$info$n
  estimate             <- x$estimate
  formula              <- x$info$formula
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  std                  <- sqrt(diag(x$cov[-K, -K, drop = FALSE]))
  sigma                <- estimate[K]
  llh                  <- x$info$log.like
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:4)], list(...))
  
  
  nfr                  <- x$info$nlinks
  cat("SAR Model\n\n")
  cat("Call:\n")
  cat(paste0(formula, ", fixed.effects = ", x$info$fixed.effects), "\n")
  # print(formula)
  cat("\nMethod: Quasi-Maximum Likelihood (ML)", "\n\n")
  
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: ", llh, "\n")
  
  invisible(x)
}

#' @rdname summary.sar
#' @export
"print.sar" <- function(x, ...) {
  stopifnot(class(x) == "sar")
  print(summary(x, ...))
}
