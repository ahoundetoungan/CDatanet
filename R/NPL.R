#' @title Estimate Count Data Model With Social Interactions
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
CDnetNPL    <- function(formula,
                        contextual, 
                        Glist, 
                        theta0    = NULL, 
                        yb0       = NULL,
                        optimizer = "optim", 
                        fix.ctr   = list(),
                        npl.ctr   = list(), 
                        opt.ctr   = list(), 
                        data) {
  
  stopifnot(optimizer %in% c("optim", "nlm"))
  
  # controls
  fix.tol     <- fix.ctr$tol 
  fix.maxit   <- fix.ctr$maxit
  
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  
  if (is.null(fix.tol)) {
    fix.tol   <- 1e-13
  }
  if (is.null(fix.maxit)) {
    fix.maxit <- 1000L
  }
  if (is.null(npl.print)) {
    npl.print <- TRUE
  }
  if (is.null(npl.tol)) {
    npl.tol   <- 1e-4
  }
  if (is.null(npl.maxit)) {
    npl.maxit <- 500L
  }
  
  # data
  if (missing(contextual)) {
    contextual <- FALSE
  }
  
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = theta0)
  formula  <- f.t.data$formula
  y        <- f.t.data$y
  X        <- f.t.data$X
  coln     <- c("lambda", colnames(X))
  
  h1       <- c(-Inf, seq(0, max(y)))
  K        <- ncol(X)
  
  # yb Gyb theta
  ybt      <- rep(0, n)
  if (!is.null(yb0)) {
    ybt    <- yb0
  }
  thetat   <- NULL
  if (!is.null(theta0)) {
    if(length(theta0) != (K + 2)) {
      stop("Length of theta0 is not suited.")
    }
    thetat <- c(log(theta0[1]/(1 -theta0[1])), theta0[2:(K+1)], log(theta0[K+2]))
  } else {
    Xtmp   <- cbind(f.t.data$Gy, X)
    b      <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    s      <- sqrt(sum((y - Xtmp%*%b)^2)/n)
    thetat <- c(log(max(b[1]/(1 - b[1]), 0.01)), b[-1], log(s))
  }
  
  Gybt     <- unlist(lapply(1:M, function(x) Glist[[x]] %*% ybt[(igr[x,1]:igr[x,2])+1]))
  
  # other variables
  cont     <- TRUE
  t        <- 0
  theta    <- NULL
  REt      <- NULL
  llht     <- 0
  par0     <- NULL
  par1     <- NULL
  like     <- NULL
  steps    <- list()
  
  # arguments
  ctr      <- c(list(yb = ybt, Gyb = Gybt, X = X, G = Glist, igroup = igr, ngroup = M,
                     h1 = h1, K = K, n = n, y = y), opt.ctr)
  if (optimizer == "optim") {
    ctr    <- c(ctr, list(fn = foptimREM_NPL, par = thetat))
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else {
    ctr    <- c(ctr, list(f = foptimREM_NPL, p = thetat))
    par0   <- "p"
    par1   <- "estimate"
    like   <- "minimum"
  }
  
  
  if(npl.print) {
    while(cont) {
      ybt0        <- ctr$yb + 0
      REt         <- do.call(get(optimizer), ctr)
      thetat      <- REt[[par1]]
      llht        <- -REt[[like]]
      
      theta       <- c(1/(1 + exp(-thetat[1])), thetat[2:(K+1)], exp(thetat[K+2]))
      
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, K, n, fix.tol, fix.maxit)
      
      dist        <- sum(abs(ctr[[par0]] - thetat)) + sum(abs(ybt0 - ybt))
      cont        <- (dist > npl.tol & t < (npl.maxit - 1))
      t           <- t + 1
      REt$dist    <- dist
      
      ctr[[par0]] <- thetat
      steps[[t]]  <- REt
      
      cat("---------------\n")
      cat(paste0("Step          : ", t), "\n")
      cat(paste0("Distance      : ", round(dist,3)), "\n")
      cat(paste0("Likelihood    : ", round(llht,3)), "\n")
      cat("Estimate", "\n")
      print(theta)
      
    }
  } else {
    while(cont) {
      ybt0        <- ctr$yb + 0
      REt         <- do.call(get(optimizer), ctr)
      thetat      <- REt[[par1]]
      
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, K, n, fix.tol, fix.maxit)
      
      dist        <- sum(abs(ctr[[par0]] - thetat)) + sum(abs(ybt0 - ybt))
      cont        <- (dist > npl.tol & t < (npl.maxit - 1))
      t           <- t + 1
      REt$dist    <- dist
      
      ctr[[par0]] <- thetat
      steps[[t]]  <- REt
    }
    llht          <- -REt[[like]]
    theta         <- c(1/(1 + exp(-thetat[1])), thetat[2:(K+1)], exp(thetat[K+2]))
  }
  
  names(theta)    <- c(coln, "sigma")
  
  
  sdata <- c(deparse(substitute(formula)), deparse(substitute(contextual)), deparse(substitute(Glist)))
  if (!missing(data)) {
    sdata         <- c(sdata, deparse(substitute(data)))
  }
  out             <- list("M"          = M,
                          "n"          = n,
                          "iteration"  = t, 
                          "estimate"   = theta, 
                          "likelihood" = llht, 
                          "yb"         = ybt, 
                          "Gyb"        = Gybt,
                          "steps"      = steps,
                          "codedata"   = sdata)
  class(out)      <- "CDnetNPL"
  out
}


#' @title Summarize Count Data Model With Social Interactions
#' @description Summary and print methods for the class `CDnetNPL` as returned by the function \link{CDnetNPL}.
#' @param ... further arguments passed to or from other methods.
#' @export 
"summary.CDnetNPL" <- function(object,
                               cov.ctr   = list(),
                               ...) {
  stopifnot(class(object) == "CDnetNPL")
  codedata      <- object$codedata
  formula       <- as.formula(codedata[1])
  contextual    <- as.logical(codedata[2])
  Glist         <- get(codedata[3])
  data          <- environment(formula)
  if (length(codedata) == 4) {
    data        <- get(codedata[4])
  }
  theta         <- object$estimate
  Gyb           <- object$Gyb
  
  J             <- length(theta)
  
  lambda        <- theta[1]
  b             <- theta[2:(J - 1)]
  sigma         <- theta[J]
  
  
  #data
  if (!is.list(Glist)) {
    Glist       <- list(Glist)
  }
  
  M             <- length(Glist)
  nvec          <- unlist(lapply(Glist, nrow))
  n             <- sum(nvec)
  igr           <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  f.t.data      <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = 0)
  X             <- f.t.data$X
  coln          <- c("lambda", colnames(X))
  
  
  Z             <- cbind(Gyb, X)
  Zd            <- Z %*% c(lambda, b)
  
  # controls
  R0            <- cov.ctr$R
  S0            <- cov.ctr$S
  
  if(is.null(R0)) {
    R0          <- ceiling(max(Zd) + 10*sigma) + 1
  }
  if(is.null(S0)) {
    S0          <- 1000
  }
  
  # compute ABCbd
  tu            <- runif(S0)
  tutu          <- tu*tu
  
  
  Rvec          <- matrix(c(-Inf, 0:R0), nrow = 1)
  
  miq           <- kronecker(Zd, Rvec, "-")
  
  cpdf          <- dnorm(miq, mean = 0, sd = sigma, log = TRUE)
  ccdf          <- pnorm(miq, mean = 0, sd = sigma, lower.tail = TRUE, log.p = TRUE)
  
  
  if(!all(miq[,R0 + 1]/sigma < -10)){
    stop("R is not sufficieltly large")
  }
  
  miq[,1]       <- 0
  tmp           <- cABC(n, miq, sigma, R0, S0, tu, tutu, cpdf, ccdf)
  A             <- c(tmp[[1]])
  B             <- c(tmp[[2]])
  C             <- c(tmp[[3]])
  d             <- c(tmp[[4]])
  b             <- c(tmp[[5]])
  
  # Sigma and Omega
  AZ            <- A*Z
  Sigma         <- matrix(NA, J, J)     
  Sigma[-J, -J] <- t(AZ)%*%Z
  Sigma[-J, J]  <- colSums(C*Z)
  Sigma[J, -J]  <- Sigma[-J, J] 
  Sigma[J, J]   <- sum(B)
  
  
  Omega         <- matrix(0, J, J)  
  for (m in 1:M) {
    r1          <- igr[m, 1] + 1
    r2          <- igr[m, 2] + 1
    dm          <- d[r1:r2]
    Zm          <- Z[r1:r2,]
    Gm          <- Glist[[m]]
    bm          <- b[r1:r2]
    AZm         <- AZ[r1:r2,]
    Bm          <- B[r1:r2]
    
    S           <- diag(r2 - r1 + 1) - lambda*dm*Gm
    M           <- cbind(dm*Zm, bm)
    SinvM       <- solve(S, M)
    
    Omega       <- Omega + (t(cbind(AZm, Bm)) %*% Gm) %*% SinvM
  }
  
  Omega         <- Omega*lambda
  
  
  tmp1          <- solve(Sigma + Omega)
  tmp2          <- tmp1 %*% Sigma %*% t(tmp1)
  
  covout          <- tmp2[-J, -J]
  
  colnames(covout)<- coln
  rownames(covout)<- coln
  
  cov.ctr         <- list("R0" = R0, "S0" = S0)
  
  out             <- c(object[-9], 
                       list("cov"       = covout, 
                            "cov.ctr"   = cov.ctr, 
                            "codedata"  = codedata,
                            "..."       = ...)) 
  
  class(out)      <- "summary.CDnetNPL"
  out
}


#' @rdname summary.CDnetNPL
#' @export
"print.summary.CDnetNPL"  <- function(x, ...) {
  stopifnot(class(x) == "summary.CDnetNPL")
  
  M                    <- x$M
  n                    <- x$n
  iteration            <- x$iteration
  estimate             <- x$estimate
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  std                  <- sqrt(diag(x$cov))
  sigma                <- estimate[K]
  llh                  <- x$likelihood
  Glist                <- get(x$codedata[3])
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:11)])
  
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  nfr                  <- unlist(lapply(Glist, function(u) sum(u > 0)))
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log pseudo-likelihood: ", llh, "\n")
  
  out                  <- c(x[1:4], list("coefficients" = out), x[-(1:4)])
  class(out)           <- "print.summary.CDnetNPL"
  invisible(out)
}

#' @rdname summary.CDnetNPL
#' @export
"print.CDnetNPL" <- function(x, ...) {
  stopifnot(class(x) == "CDnetNPL")
  print(summary(x, ...))
}