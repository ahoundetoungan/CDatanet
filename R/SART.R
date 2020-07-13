#' @title Estimate SART model
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
SARTML <- function(formula,
                   contextual,
                   Glist,
                   theta0 = NULL,
                   optimizer = "optim",
                   opt.ctr = list(),
                   print = TRUE,
                   cov = TRUE,
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
  
  f.t.data   <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = theta0)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X))
  
  K          <- ncol(X)
  
  # variables
  
  ylist      <- lapply(1:M, function(x) y[(igr[x,1]:igr[x,2]) + 1])
  
  idpos      <- which(y > 0) - 1
  idzero     <- which(!(y > 0)) - 1
  
  idposlis   <- lapply(ylist, function(w) which(w > 0))
  Npos       <- unlist(lapply(idposlis, length))
  
  G2list     <- lapply(1:M, function(w) Glist[[w]][idposlis[[w]], idposlis[[w]]])
  Gy         <- unlist(lapply(1:M, function(w) Glist[[w]] %*% ylist[[w]]))
  I2list     <- lapply(Npos, function(w) diag(w))
  
  
  alphatde   <- Inf
  logdetA2   <- 0
  
  
  theta     <- NULL
  if (!is.null(theta0)) {
    if(length(theta0) != (K + 2)) {
      stop("Length of theta0 is not suited.")
    }
    theta    <- c(log(theta0[1]/(1 -theta0[1])), theta0[2:(K+1)], log(theta0[K+2]))
  } else {
    Xtmp     <- cbind(f.t.data$Gy, X)
    b        <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    s        <- sqrt(sum((y - Xtmp%*%b)^2)/n)
    theta    <- c(log(max(b[1]/(1 - b[1]), 0.01)), b[-1], log(s))
  }
  
  alphatilde <- theta[1]
  logdetA2   <- 0
  
  
  # arguments
  ctr        <- c(list(X = X, G2 = G2list, I2 = I2list, K = K, y = y, Gy = Gy,
                       idpos = idpos, idzero = idzero, Npos = Npos, ngroup = M,
                       alphatilde = alphatde, logdetA2 = logdetA2,
                       hessian = cov), opt.ctr)
  if (optimizer == "optim") {
    ctr    <- c(ctr, list(par = theta))
    par1   <- "par"
    like   <- "value"
    if (print) {
      ctr    <- c(ctr, list(fn = foptimTobit))  
    } else {
      ctr    <- c(ctr, list(fn = foptimTobit0))  
    }
  } else {
    ctr    <- c(ctr, list(p = theta))
    par1   <- "estimate"
    like   <- "minimum"
    if (print) {
      ctr    <- c(ctr, list(f = foptimTobit))  
    } else {
      ctr    <- c(ctr, list(f = foptimTobit0))  
    }
  }
  
  
  
  resTO    <- do.call(get(optimizer), ctr)
  theta    <- resTO[[par1]]
  theta    <- c(1/(1 + exp(-theta[1])), theta[2:(K + 1)], exp(theta[K + 2]))
  llh      <- -resTO[[like]]
  
  covout             <- NULL
  if (cov) {
    covtmp           <- solve(resTO$hessian)[-(K + 2), -(K + 2)]
    Rmat             <- diag(K + 1)
    Rmat[1,1]        <- theta[1]*(1 - theta[1])
    covout           <- Rmat %*% covtmp %*% t(Rmat)
    colnames(covout) <- coln
    rownames(covout) <- coln
  }
  
  names(theta)       <- c(coln, "sigma")
  
  sdata <- c(as.character(formula), deparse(substitute(Glist)))
  if (!missing(data)) {
    sdata         <- c(sdata, deparse(substitute(data)))
  }
  
  out                <- list("M"             = M,
                             "n"             = n,
                             "estimate"      = theta, 
                             "likelihood"    = llh, 
                             "cov"           = covout,
                             "optimization"  = resTO,
                             "codedata"      = sdata)
  class(out)         <- "SARTML"
  out
  
}



#' @title Summarize SART Model
#' @description Summary and print methods for the class `SARTML` as returned by the function \link{SARTML}.
#' @param ... further arguments passed to or from other methods.
#' @export 
"summary.SARTML" <- function(object,
                            ...) {
  stopifnot(class(object) == "SARTML")
  out           <- c(object, list("..."       = ...)) 
  
  if(is.null(object$cov)){
    stop("Covariance was not computed")
  }
  
  class(out)    <- "summary.SARTML"
  out
}


#' @rdname summary.SARTML
#' @export
"print.summary.SARTML"  <- function(x, ...) {
  stopifnot(class(x) == "summary.SARTML")
  
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
  cat("SART Model\n\n")
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
  class(out)           <- "print.summary.SARTML"
  invisible(out)
}

#' @rdname summary.SARTML
#' @export
"print.SARTML" <- function(x, ...) {
  stopifnot(class(x) == "SARTML")
  print(summary(x, ...))
}