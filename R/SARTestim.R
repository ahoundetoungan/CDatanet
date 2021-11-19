#' @title Estimate SART model
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param theta0 (optional) starting value of \eqn{\theta = (\lambda, \beta, \gamma, \sigma)}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param yb0 (optional) expectation of y.
#' @param optimizer is either `nlm` (referring to the function \link[stats]{nlm}) or `optim` (referring to the function \link[stats]{optim}). 
#' Other arguments 
#' of these functions such as, the control values and the method can be defined through the argument `opt.ctr`.
#' @param npl.ctr list of controls for the NPL method (see \code{\link{CDnetNPL}}).
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param cov a boolean indicating if the covariance should be computed.
#' @param RE a boolean which indicates if the model if under rational expectation of not.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `SARTML` is called.
#' @return A list consisting of:
#'     \item{info}{list of general information on the model.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{yb}{ybar (see details), expectation of y.}
#'     \item{Gyb}{average of the expectation of y among friends.}
#'     \item{cov}{List of covariances.}
#'     \item{details}{outputs as returned by the optimizer.}
#' @details 
#' ## Model
#' The left-censored variable \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.\cr
#' The count variable \eqn{y_i} is then define that is \eqn{y_i = 0} if  
#' \eqn{y_i^* \leq 0}{ys_i \le 0} and \eqn{y_i = y_i^*}{y_i = ys_i} otherwise.
#' @seealso \code{\link{CDnetNPL}} and \code{\link{SARML}}.
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
#' ytmp    <- simTobitnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                        theta = theta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y)
#' 
#' opt.ctr <- list(method  = "Nelder-Mead", 
#'                 control = list(abstol = 1e-16, abstol = 1e-11, maxit = 5e3))
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out     <- SARTML(formula = yt ~ x1 + x2, optimizer = "nlm",
#'                   contextual = TRUE, Glist = Glist, data = data)
#' summary(out)
#' }
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @importFrom numDeriv hessian
#' @export
SARTML <- function(formula,
                   contextual,
                   Glist,
                   theta0 = NULL,
                   yb0  = NULL,
                   optimizer = "optim",
                   npl.ctr  = list(), 
                   opt.ctr = list(),
                   print = TRUE,
                   cov = TRUE,
                   RE = FALSE,
                   data) {
  stopifnot(optimizer %in% c("optim", "nlm"))
  env.formula <- environment(formula)
  # controls
  npl.print   <- print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  
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
  
  f.t.data   <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = NULL)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X), "sigma")
  
  K          <- ncol(X)
  
  # variables
  indpos     <- (y > 1e-323)
  indzero    <- !indpos
  idpos      <- which(indpos) - 1
  idzero     <- which(indzero) - 1
  ylist      <- lapply(1:M, function(x) y[(igr[x,1]:igr[x,2]) + 1])
  idposlis   <- lapply(ylist, function(w) which(w > 0))
  npos       <- unlist(lapply(idposlis, length))  
  
  thetat     <- NULL
  if (!is.null(theta0)) {
    if(length(theta0) != (K + 2)) {
      stop("Length of theta0 is not suited.")
    }
    thetat   <- c(log(theta0[1]/(1 -theta0[1])), theta0[2:(K+1)], log(theta0[K+2]))
  } else {
    Xtmp     <- cbind(f.t.data$Gy, X)
    b        <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    s        <- sqrt(sum((y - Xtmp%*%b)^2)/n)
    thetat   <- c(log(max(b[1]/(1 - b[1]), 0.01)), b[-1], log(s))
  }
  
  llh        <- NULL
  resTO      <- list()
  covt       <- NULL
  covm       <- NULL
  t          <- NULL
  ybt        <- NULL
  Gybt       <- NULL
  theta      <- NULL
  Ztlambda   <- NULL
  
  
  if(RE){
    # yb 
    ybt      <- rep(0, n)
    if (!is.null(yb0)) {
      ybt    <- yb0
    }
    # Gybt
    Gybt     <- unlist(lapply(1:M, function(x) Glist[[x]] %*% ybt[(igr[x,1]:igr[x,2])+1]))
    if (is.null(npl.print)) {
      npl.print <- TRUE
    }
    if (is.null(npl.tol)) {
      npl.tol   <- 1e-4
    }
    if (is.null(npl.maxit)) {
      npl.maxit <- 500L
    }
    # other variables
    cont     <- TRUE
    t        <- 0
    REt      <- NULL
    llh      <- 0
    par0     <- NULL
    par1     <- NULL
    like     <- NULL
    yidpos   <- y[indpos]
    ctr      <- c(list("yidpos" = yidpos, "Gyb" = Gybt, "X" = X, "npos" = sum(npos), 
                       "idpos" = idpos, "idzero" = idzero, "K" = K), opt.ctr)
    
    if (optimizer == "optim") {
      ctr    <- c(ctr, list(fn = foptimTBT_NPL, par = thetat))
      
      par0   <- "par"
      par1   <- "par"
      like   <- "value"
    } else {
      ctr    <- c(ctr, list(f = foptimTBT_NPL,  p = thetat))
      par0   <- "p"
      par1   <- "estimate"
      like   <- "minimum"
    }
    
    if(npl.print) {
      while(cont) {
        ybt0        <- ybt + 0    #copy in different memory
        
        # compute theta
        REt         <- do.call(get(optimizer), ctr)
        thetat      <- REt[[par1]]
        llh         <- -REt[[like]]
        
        theta       <- c(1/(1 + exp(-thetat[1])), thetat[2:(K + 1)], exp(thetat[K + 2]))
        
        # compute y
        fLTBT_NPL(ybt, Gybt, Glist, X, thetat, igr, M, n, K)
        
        # distance
        dist        <- sum(abs(ctr[[par0]] - thetat)) + sum(abs(ybt0 - ybt))
        cont        <- (dist > npl.tol & t < (npl.maxit - 1))
        t           <- t + 1
        REt$dist    <- dist
        ctr[[par0]] <- thetat
        resTO[[t]]  <- REt
        
        cat("---------------\n")
        cat(paste0("Step          : ", t), "\n")
        cat(paste0("Distance      : ", round(dist,3)), "\n")
        cat(paste0("Likelihood    : ", round(llh,3)), "\n")
        cat("Estimate:", "\n")
        print(theta)
        
      }
    } else {
      while(cont) {
        ybt0        <- ybt + 0    #copy in different memory
        
        # compute theta
        REt         <- do.call(get(optimizer), ctr)
        thetat      <- REt[[par1]]
        
        # compute y
        fLTBT_NPL(ybt, Gybt, Glist, X, thetat, igr, M, n, K)
        
        # distance
        dist        <- sum(abs(ctr[[par0]] - thetat)) + sum(abs(ybt0 - ybt))
        cont        <- (dist > npl.tol & t < (npl.maxit - 1))
        t           <- t + 1
        REt$dist    <- dist
        
        ctr[[par0]] <- thetat
        resTO[[t]]  <- REt
      }
      llh           <- -REt[[like]]
      theta         <- c(1/(1 + exp(-thetat[1])), thetat[2:(K +1)], exp(thetat[K + 2]))
    }
    
    if (npl.maxit == t) {
      warning("The maximum number of iterations of the NPL algorithm has been reached.")
    }
    
    if (cov) {
      covt           <- solve(hessian(foptimRE_TBT, thetat, method="Richardson",
                                      yb = ybt, Gyb = Gybt, yidpos = yidpos, X = X, G = Glist, 
                                      igroup = igr, ngroup = M, npos = sum(npos), 
                                      idpos = idpos, idzero = idzero, K = K, n = n, 
                                      tol = npl.tol, maxit  = npl.maxit))
      
      Rm               <- diag(K + 2)
      Rm[1,1]          <- exp(thetat[1])/(1 + exp(thetat[1]))^2
      Rm[K + 2, K + 2] <- exp(thetat[K + 2])
      covt             <- Rm %*% covt %*% Rm
      colnames(covt)   <- coln
      rownames(covt)   <- coln
    }
    
    Ztlambda <- Gybt*theta[1] + X%*%theta[2:(K +1)]
  } else{
    G2list     <- lapply(1:M, function(w) Glist[[w]][idposlis[[w]], idposlis[[w]]])
    Gy         <- unlist(lapply(1:M, function(w) Glist[[w]] %*% ylist[[w]]))
    I2list     <- lapply(npos, diag)
    Ilist      <- lapply(nvec, diag)  
    Wlist      <- lapply(1:M, function(x) (indpos[(igr[x,1]:igr[x,2]) + 1] %*% t(indpos[(igr[x,1]:igr[x,2]) + 1]))*Glist[[x]])
    if(exists("alphatde")) rm("alphatde")
    if(exists("logdetA2")) rm("logdetA2")
    alphatde   <- Inf
    logdetA2   <- 0
    
    # arguments
    ctr        <- c(list("X" = X, "G2" = G2list, "I2" = I2list, "K" = K, "y" = y, "Gy" = Gy,
                         "idpos" = idpos, "idzero" = idzero, "npos" = sum(npos), "ngroup" = M,
                         "alphatilde" = alphatde, "logdetA2" = logdetA2, "N" = n, 
                         "I" = Ilist,  "W" = Wlist, "igroup" = igr), opt.ctr)
    if (optimizer == "optim") {
      ctr    <- c(ctr, list(par = thetat))
      par1   <- "par"
      like   <- "value"
      if (npl.print) {
        ctr  <- c(ctr, list(fn = foptimTobit))  
      } else {
        ctr  <- c(ctr, list(fn = foptimTobit0))  
      }
    } else {
      ctr    <- c(ctr, list(p = thetat))
      par1   <- "estimate"
      like   <- "minimum"
      if (npl.print) {
        ctr  <- c(ctr, list(f = foptimTobit)) 
      } else {
        ctr  <- c(ctr, list(f = foptimTobit0)) 
      }
    }
    
    resTO    <- do.call(get(optimizer), ctr)
    theta    <- resTO[[par1]]
    theta    <- c(1/(1 + exp(-theta[1])), theta[2:(K + 1)], exp(theta[K + 2]))
    llh      <- -resTO[[like]]
    rm(list = c("alphatde", "logdetA2"))
    if (cov) {
      covtmp           <- fqTobit(theta, X, G2list, Ilist, Wlist, K, n, y, Gy,  indzero,
                                  indpos, igr, M)
      covt           <- solve(stats::cov(covtmp))/n
      colnames(covt) <- coln
      rownames(covt) <- coln
    }
    Ztlambda <- Gy*theta[1] + X%*%theta[2:(K +1)]
  }
  
  names(theta)         <- coln
  
  # Marginal effects
  thetaME              <- head(theta, K + 1)
  colnME               <- head(coln, K + 1)
  if("(Intercept)" %in% coln) {
    thetaME            <- thetaME[-2]
    colnME             <- colnME[-2]
  }
  
  ZLst                 <- Ztlambda/theta[K + 2]
  meanPhiZLst          <- mean(pnorm(ZLst))
  meffects             <- thetaME*meanPhiZLst
  names(meffects)      <- colnME
  
  if(cov) {
    phiZLst            <- c(dnorm(ZLst))
    meanphiZLstz       <- NULL 
    if(RE){
      meanphiZLstz     <- colSums(phiZLst*as.matrix(cbind(Gybt, X)))/n
    } else {
      meanphiZLstz     <- colSums(phiZLst*as.matrix(cbind(f.t.data$Gy, X)))/n 
    }
    tmp1               <- diag(K + 1)*meanPhiZLst + theta[-(K + 2)] %*% matrix(meanphiZLstz, nrow = 1)/theta[K + 2]
    tmp2               <- - mean(ZLst*phiZLst)*theta[-(K + 2)]/theta[K + 2]
    tmp3               <- cbind(tmp1, tmp2)
    
    covm               <- tmp3 %*% covt %*% t(tmp3)
    
    if("(Intercept)" %in% coln) {
      covm             <- covm[-2, -2]
    }
    colnames(covm)     <- colnME
    rownames(covm)     <- colnME
  }
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "formula"    = formula,
                               "nlinks"     = unlist(lapply(Glist, function(u) sum(u > 0))),
                               "censured"   = sum(indzero),
                               "uncensured" = n - sum(indzero),
                               "log.like"   = llh, 
                               "npl.iter"   = t)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(theta = theta, marg.effects = meffects),
                               "yb"         = ybt, 
                               "Gyb"        = Gybt,
                               "cov"        = list(parms = covt, marg.effects = covm),
                               "details"    = resTO)
  class(out)            <- "SARTML"
  out
}



#' @title Summarize SART Model
#' @description Summary and print methods for the class `SARTML` as returned by the function \link{SARTML}.
#' @param object an object of class `SARTML`, output of the function \code{\link{SARTML}}.
#' @param x an object of class `summary.SARTML`, output of the function \code{\link{summary.SARTML}} 
#' or class `SARTML`, output of the function \code{\link{SARTML}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
#' @param ... further arguments passed to or from other methods.
#' @export 
"summary.SARTML" <- function(object,
                             ...) {
  stopifnot(class(object) == "SARTML")
  out           <- c(object, list("..." = ...)) 
  if(is.null(object$cov$parms)){
    stop("Covariance was not computed")
  }
  class(out)    <- "summary.SARTML"
  out
}


#' @rdname summary.SARTML
#' @export
"print.summary.SARTML"  <- function(x, ...) {
  stopifnot(class(x) == "summary.SARTML")
  
  M                    <- x$info$M
  n                    <- x$info$n
  estimate             <- x$estimate$theta
  iteration            <- x$info$npl.iter
  RE                   <- !is.null(iteration)
  formula              <- x$info$formula
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  meff                 <- x$estimate$marg.effects
  std                  <- sqrt(diag(x$cov$parms)[-K])
  std.meff             <- sqrt(diag(x$cov$marg.effects))
  sigma                <- estimate[K]
  llh                  <- x$info$log.like
  censored             <- x$info$censured
  uncensored           <- x$info$uncensured
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:6)], list(...))
  
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:6)], list(...))
  
  
  nfr                  <- x$info$nlinks
  cat("SART Model", ifelse(RE, "with Rational Expectation", ""), "\n\n")
  cat("Call:\n")
  print(formula)
  if(RE){
    cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, "\n\n")
  } else{
    cat("\nMethod: Maximum Likelihood (ML)", "\n\n")
  }
  
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("No. censored             : ", paste0(censored, "(", round(censored/n*100, 0), "%)"), "\n")
  cat("No. uncensored           : ", paste0(uncensored, "(", round(uncensored/n*100, 0), "%)"), "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: ", llh, "\n")
  
  invisible(x)
}

#' @rdname summary.SARTML
#' @export
"print.SARTML" <- function(x, ...) {
  stopifnot(class(x) == "SARTML")
  print(summary(x, ...))
}


#' @rdname summary.SARTML
#' @export
"print.summary.SARTMLs" <- function(x, ...) {
  stopifnot(class(x) %in% c("list", "SARTMLs", "summary.SARTMLs")) 
  
  lclass          <- unique(unlist(lapply(x, class)))
  if (!all(lclass %in% c("SARTML", "summary.SARTML"))) {
    stop("All the components in `x` should be from `SARTML` or `summary.SARTML` class")
  }
  
  nsim            <- length(x)
  coef            <- do.call("rbind", lapply(x, function(z) t(z$estimate$theta)))
  meff            <- do.call("rbind", lapply(x, function(z) t(z$estimate$marg.effects)))
  estimate        <- colSums(coef)/nsim
  meffects        <- colSums(meff)/nsim
  
  vcoef2          <- Reduce("+", lapply(x, function(z) z$cov$parms))/nsim
  vmeff2          <- Reduce("+", lapply(x, function(z) z$cov$marg.effects))/nsim
  
  vcoef1          <- cov(coef)
  vmeff1          <- cov(meff)
  
  vcoef           <- vcoef1 + vcoef2
  vmeff           <- vmeff1 + vmeff2
  
  
  llh             <- sapply(x, function(z) z$info$log.like)
  llh             <- c("min" = min(llh), "mean" = mean(llh), "max" = max(llh))
  
  M               <- x[[1]]$info$M
  n               <- x[[1]]$info$n
  
  
  K               <- length(estimate)
  coef            <- estimate[-K]
  std             <- sqrt(diag(vcoef)[-K])
  std.meff        <- sqrt(diag(vmeff))
  sigma           <- estimate[K]
  
  tmp             <- fcoefficients(coef, std)
  out_print       <- tmp$out_print
  out             <- tmp$out
  
  tmp.meff        <- fcoefficients(meffects, std.meff)
  out_print.meff  <- tmp.meff$out_print
  out.meff        <- tmp.meff$out
  
  out_print       <- c(list(out_print), x[[1]][-(1:6)], list(...))
  out_print.meff  <- c(list(out_print.meff), x[[1]][-(1:6)], list(...))
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Replication of SART-ML \nReplication: ", nsim, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log pseudo-likelihood: ", "\n")
  print(llh)
  
  out                  <- list("M"          = M,
                               "n"          = n,
                               "simulation" = nsim, 
                               "estimate"   = list(theta = estimate, marg.effects = meffects), 
                               "likelihood" = llh, 
                               "cov"        = list(parms = vcoef, marg.effects = vmeff), 
                               ...          = ...)
  class(out)           <- "print.summary.SARTMLs"
  invisible(out)
} 