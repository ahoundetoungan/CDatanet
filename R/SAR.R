#' @title Estimate SAR model
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param lambda0 (optional) starting value of \eqn{\lambda}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param optimizer is either `nlm` (referring to the function \link[stats]{nlm}) or `optim` (referring to the function \link[stats]{optim}). 
#' Other arguments 
#' of these functions such as, the control values and the method can be defined through the argument `opt.ctr`.
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param print a boolean indicating if the estimate should be printed at each step.
#' @param cov a boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @details 
#' ## Model
#' The variable \eqn{\mathbf{y}}{y} is given for all i as
#' \deqn{y_i = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{y_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.
#' ## `codedata`
#' The \link[base]{class} of the output of this function is \code{SAR}. This class has a \link[base]{summary} 
#' and \link[base]{print} \link[utils]{methods} to summarize and print the results. 
#' In order to save 
#' memory, the function does not return neither the adjacency matrix nor the data. Instead, it returns `codedata` which contains among others, the `formula` 
#' and the name of the adjacency matrix passed through the argument `Glist`.
#' @seealso \code{\link{CDnetNPL}} and \code{\link{SARTML}}.
#' @examples 
#' \donttest{
#' # Groups' size
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 1000))
#' n      <- sum(nvec)
#' 
#' # Parameters
#' lambda <- 0.4
#' beta   <- c(2, -1.9, 0.8)
#' gamma  <- c(1.5, -1.2)
#' sigma  <- 1.5
#' theta  <- c(lambda, beta, gamma, sigma)
#' 
#' # X
#' X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))
#' 
#' # Network
#' Glist  <- list()
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
#'   Glist[[m]]   <- Gm
#' }
#' 
#' 
#' # data
#' data    <- data.frame(x1 = X[,1], x2 =  X[,2])
#' 
#' rm(list = ls()[!(ls() %in% c("Glist", "data", "theta"))])
#' 
#' ytmp    <- simSARnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                      theta = theta, data = data) 
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
#' 
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out     <- SARML(formula = yt ~ x1 + x2, contextual = TRUE, 
#'                  Glist = Glist, optimizer = "optim", data = data)
#' summary(out)
#' }
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network and name of the object `data` (see details).}
#' @export
SARML <- function(formula,
                  contextual,
                  Glist, 
                  lambda0   = NULL, 
                  optimizer = "optim",
                  opt.ctr   = list(), 
                  print     = TRUE, 
                  cov       = TRUE,
                  data) {
  stopifnot(optimizer %in% c("optim", "nlm"))
  env.formula <- environment(formula)
  #size 
  if (missing(contextual)) {
    contextual <- FALSE
  }
  if (!is.list(Glist)) {
    Glist    <- list(Glist)
  }
  M          <- length(Glist)
  nvec       <- unlist(lapply(Glist, nrow))
  n          <- sum(nvec)
  igr        <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  f.t.data   <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = lambda0)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X), "sigma")
  
  K          <- ncol(X)
  
  # variables
  ylist      <- lapply(1:M, function(x) y[(igr[x,1]:igr[x,2]) + 1])
  
  Nvec       <- unlist(lapply(ylist, length))
  
  Gy         <- unlist(lapply(1:M, function(w) Glist[[w]] %*% ylist[[w]]))
  
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
  ctr        <- c(list(X = X,invXX = invXX, G = Glist, I = Ilist, N = n,
                       y = y, Gy = Gy, ngroup = M), opt.ctr)
  
  if (optimizer == "optim") {
    ctr      <- c(ctr, list(par = lambdat))
    par0     <- "par"
    par1     <- "par"
    like     <- "value"
    if (print) {
      ctr    <- c(ctr, list(fn = foptimSAR))  
    } else {
      ctr    <- c(ctr, list(fn = foptimSAR0))  
    }
  } else {
    ctr     <- c(ctr, list(p = lambdat))
    par0    <- "p"
    par1    <- "estimate"
    like    <- "minimum"
    if (print) {
      ctr   <- c(ctr, list(f = foptimSAR))  
    } else {
      ctr   <- c(ctr, list(f = foptimSAR0))  
    }
  }
  
  resSAR    <- do.call(get(optimizer), ctr)
  lambdat   <- resSAR[[par1]]
  llh       <- -resSAR[[like]]
  
  lambda    <- 1/(1 + exp(-lambdat))
  
  hessian   <- jacobSAR(lambda, X, invXX, XX, y, n, Glist, Ilist, Gy, M, igr, cov)
  
  
  beta      <- hessian$beta
  sigma2    <- hessian$sigma2
  
  covout    <- NULL
  if(cov) {
    covout           <- hessian$cov
    colnames(covout) <- coln
    rownames(covout) <- coln
  }
  
  theta           <- c(lambda, beta, sqrt(sigma2))
  names(theta)    <- coln
  
  environment(formula) <- env.formula
  sdata                <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist)),
    "nfriends"      = unlist(lapply(Glist, function(u) sum(u > 0))) 
  )
  if (!missing(data)) {
    sdata              <- c(sdata, list("data" = deparse(substitute(data))))
  }  
  
  out                  <- list("M"             = M,
                          "n"             = n,
                          "estimate"      = theta, 
                          "likelihood"    = llh, 
                          "cov"           = covout,
                          "optimization"  = resSAR,
                          "codedata"      = sdata)
  class(out)           <- "SARML"
  out
}




#' @title Summarize SAR Model
#' @description Summary and print methods for the class `SARML` as returned by the function \link{SARML}.
#' @param object an object of class `SARML`, output of the function \code{\link{SARML}}.
#' @param x an object of class `summary.SARML`, output of the function \code{\link{summary.SARML}} or 
#' class `SARML`, output of the function \code{\link{SARML}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network and name of the object `data`.}
#' @param ... further arguments passed to or from other methods.
#' @export 
"summary.SARML" <- function(object,
                            ...) {
  stopifnot(class(object) == "SARML")
  out           <- c(object, list("..."       = ...)) 
  
  
  if(is.null(object$cov)){
    stop("Covariance was not computed")
  }
  class(out)    <- "summary.SARML"
  out
}


#' @rdname summary.SARML
#' @export
"print.summary.SARML"  <- function(x, ...) {
  stopifnot(class(x) == "summary.SARML")
  
  M                    <- x$M
  n                    <- x$n
  estimate             <- x$estimate
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  std                  <- sqrt(diag(x$cov[-K, -K, drop = FALSE]))
  sigma                <- estimate[K]
  llh                  <- x$likelihood
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:7)], list(...))
  

  nfr                  <- x$codedata$nfriends
  cat("SAR Model\n\n")
  cat("Method: Maximum Likelihood (ML)", "\n\n")
  
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: ", llh, "\n")
  
  out                  <- c(x[1:4], list("coefficients" = out), x[-(1:4)])
  class(out)           <- "print.summary.SARML"
  invisible(out)
}

#' @rdname summary.SARML
#' @export
"print.SARML" <- function(x, ...) {
  stopifnot(class(x) == "SARML")
  print(summary(x, ...))
}

#' @rdname summary.SARML
#' @export
"print.summary.SARMLs" <- function(x, ...) {
  stopifnot(class(x) %in% c("list", "summary.SARMLs", "print.summary.SARMLs")) 
  
  type2               <- (class(x) == "print.summary.SARMLs")
  nsim                <- NULL
  estimate            <- NULL
  vcoef               <- NULL
  llh                 <- NULL
  n                   <- NULL
  M                   <- NULL
  
  if (type2) {
    nsim              <- x$simulation
    estimate          <- x$estimate
    vcoef             <- x$cov
    llh               <- x$likelihood
    n                 <- x$n
    M                 <- x$M
  } else {
    lclass            <- unique(unlist(lapply(x, class)))
    if (!all(lclass %in% "summary.SARML")) {
      stop("All the components in `x` should be from `summary.SARML` class")
    }
    
    nsim              <- length(x)
    coef              <- do.call("rbind", lapply(x, function(z) t(z$estimate)))
    estimate          <- colSums(coef)/nsim
    
    vcoef2            <- Reduce("+", lapply(x, function(z) z$cov))/nsim
    
    vcoef1            <- cov(coef)
    
    vcoef             <- vcoef1 + vcoef2
    
    
    llh               <- unlist(lapply(x, function(z) z$likelihood))
    llh               <- c("min" = min(llh), "mean" = mean(llh), "max" = max(llh))
    
    M                 <- x[[1]]$M
    n                 <- x[[1]]$n
  }
  
  
  
  K                   <- length(estimate)
  coef                <- estimate[-K]
  std                 <- sqrt(diag(vcoef)[-K])
  sigma               <- estimate[K]
  
  tmp                 <- fcoefficients(coef, std)
  out_print           <- tmp$out_print
  out                 <- tmp$out
  
  
  if (type2) {
    out_print         <- c(list(out_print), x[-(1:6)], list(...))
  } else {
    out_print         <- c(list(out_print), x[[1]][-(1:7)], list(...))
  }
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Replication of SAR-ML \nReplication: ", nsim, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: ", "\n")
  print(llh)
  
  out                  <- list("M"          = M,
                               "n"          = n,
                               "simulation" = nsim, 
                               "estimate"   = estimate, 
                               "likelihood" = llh, 
                               "cov"        = vcoef, 
                               ...          = ...)
  class(out)           <- "print.summary.SARMLs"
  invisible(out)
} 