#' @title Estimate SAR model
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @return A list consisting of:
#'     \item{theta0}{starting values.}
#'     \item{formula}{input value of `formula`.}
#'     \item{contextual}{input value of `contextual`.}
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
  coln       <- c("lambda", colnames(X))
  
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
    covout           <- solve(hessian$hessian)
    colnames(covout) <- coln
    rownames(covout) <- coln
  }
  
  theta           <- c(lambda, beta, sqrt(sigma2))
  names(theta)    <- c(coln, "sigma")
  
  sdata <- c(as.character(formula), deparse(substitute(Glist)))
  if (!missing(data)) {
    sdata         <- c(sdata, deparse(substitute(data)))
  }
  
  out             <- list("M"             = M,
                          "n"             = n,
                          "estimate"      = theta, 
                          "likelihood"    = llh, 
                          "cov"           = covout,
                          "optimization"  = resSAR,
                          "codedata"      = sdata)
  class(out)      <- "SARML"
  out
}




#' @title Summarize SAR Model
#' @description Summary and print methods for the class `SARML` as returned by the function \link{SARML}.
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
  std                  <- sqrt(diag(x$cov))
  sigma                <- estimate[K]
  llh                  <- x$likelihood
  Glist                <- get(x$codedata[2])
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:7)])
  
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  nfr                  <- unlist(lapply(Glist, function(u) sum(u > 0)))
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