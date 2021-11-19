#' @title Estimate Count Data Model with Social Interactions using NPL Method
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param Rbar the value of Rbar. If not provided, it is automatically set at \code{quantile(y, 0.9)}.
#' @param starting (optional) starting value of \eqn{\theta = (\lambda, \beta', \gamma')'} and \eqn{\delta = (\delta_2, ..., \delta_{\bar{R}})}{\delta = (\delta_2, ..., \delta_{Rbar})}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param yb0 (optional) expectation of y.
#' @param optimizer is either `nlm` (referring to the \link[stats]{nlm} function) or `optim` (referring to the \link[stats]{optim} function). 
#' At every step of the NPL method, the estimation is performed using \link[stats]{nlm} or \link[stats]{optim}. Other arguments 
#' of these functions such as, `control` and `method` can be defined through the argument `opt.ctr`.
#' @param npl.ctr list of controls for the NPL method (see details).
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param cov a boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `CDnetNPL` is called.
#' @return A list consisting of:
#'     \item{info}{list of general information about the model.}
#'     \item{estimate}{NPL estimator.}
#'     \item{yb}{ybar (see details), expectation of y.}
#'     \item{Gyb}{average of the expectation of y among friends.}
#'     \item{cov}{list of covariance matrices.}
#'     \item{details}{step-by-step output as returned by the optimizer.}
#' @details 
#' ## Model
#' Following Houndetoungan (2020), the count data \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i \bar{\mathbf{y}} + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*ybar + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, 1)}{\epsilon_i --> N(0, 1)}.\cr
#' Then, \eqn{y_i = r} iff \eqn{a_r \leq y_i^* \leq a_{r+1}}{a_r \le ys_i \le a_{r + 1}}, where
#' \eqn{a_0 = -\inf}{a_0 = -Inf}, \eqn{a_1 = 0}, \eqn{a_r = \sum_{k = 1}^r\delta_k}{a_r = \delta_1 + ... + \delta_r} if \eqn{1 \leq r \leq \bar{R}}{1 \le r \le Rbar}, and 
#' \eqn{a_r = (r - \bar{R})\delta_{\bar{R}} + a_{\bar{R}}}{a_r = (r - Rbar)\delta_{Rbar} + a_{Rbar}} otherwise.
#' ## \code{npl.ctr}
#' The model parameters is estimated using the Nested Partial Likelihood (NPL) method. This approach 
#' starts with a guess of \eqn{\theta} and \eqn{\bar{y}}{yb} and constructs iteratively a sequence
#' of \eqn{\theta} and \eqn{\bar{y}}{yb}. The solution converges when the \eqn{L_1}{L} distance
#' between two consecutive \eqn{\theta} and \eqn{\bar{y}}{yb} is less than a tolerance. \cr
#' The argument \code{npl.ctr} is an optional list which contain
#' \itemize{
#' \item{tol}{ the tolerance of the NPL algorithm (default 1e-4),}
#' \item{maxit}{ the maximal number of iterations allowed (default 500),}
#' \item{print}{ a boolean indicating if the estimate should be printed at each step.}
#' \item{S}{ the number of simulation performed use to compute integral in the covariance by important sampling.} 
#' }
#' @seealso \code{\link{simCDnet}}, \code{\link{SARML}} and \code{\link{SARTML}}.
#' @examples 
#' \donttest{
#' # Groups' size
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 1000))
#' n      <- sum(nvec)
#' 
#' # Parameters
#' lambda <- 0.4
#' beta   <- c(1.5, 2.2, -0.9)
#' gamma  <- c(1.5, -1.2)
#' delta  <- c(1, 0.87, 0.75, 0.6, 0.4)
#' theta  <- c(lambda, beta, gamma)
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
#' ytmp    <- simCDnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta, delta = delta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
#' 
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out   <- CDnetNPL(formula = yt ~ x1 + x2, contextual = TRUE, Glist = Glist, data = data, Rbar = 10)
#' summary(out)
#' }
#' @importFrom stats quantile
#' @export
CDnetNPL    <- function(formula,
                        contextual, 
                        Glist, 
                        Rbar      = NULL,
                        starting  = list(theta = NULL, delta = NULL), 
                        yb0       = NULL,
                        optimizer = "optim", 
                        npl.ctr   = list(), 
                        opt.ctr   = list(), 
                        cov       = TRUE,
                        data) {
  
  stopifnot(optimizer %in% c("optim", "nlm"))
  env.formula <- environment(formula)
  # controls
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  npl.S       <- npl.ctr$S
  
  thetat      <- starting$theta
  Deltat      <- starting$delta
  if (is.null(npl.print)) {
    npl.print <- TRUE
  }
  if (is.null(npl.tol)) {
    npl.tol   <- 1e-4
  }
  if (is.null(npl.maxit)) {
    npl.maxit <- 500L
  }
  if (is.null(npl.S)) {
    npl.S     <- 1e3L
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
  
  f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = thetat)
  formula  <- f.t.data$formula
  y        <- f.t.data$y
  X        <- f.t.data$X
  coln     <- c("lambda", colnames(X))
  
  maxy     <- max(y)
  K        <- ncol(X)
  
  # # simulations
  # Ksimu    <- 0
  # nsimu    <- 0
  # if(!is.null(simu1)) {
  #   Ksimu  <- 1
  #   nsimu  <- ncol(simu1)
  #   stopifnot(nrow(simu1) == n)
  #   coln   <- c(coln, "simu1")
  # } else{
  #   if(!is.null(simu2)){
  #     stop("simu1 = NULL whereas simu2 != NULL")
  #   }
  # }
  # if(!is.null(simu2)) {
  #   Ksimu  <- 2
  #   stopifnot(nrow(simu1) == n)
  #   coln   <- c(coln, "simu2")
  # }
  
  # yb 
  ybt      <- rep(0, n)
  if (!is.null(yb0)) {
    ybt    <- yb0
  }
  
  # Gybt
  Gybt     <- unlist(lapply(1:M, function(x) Glist[[x]] %*% ybt[(igr[x,1]:igr[x,2])+1]))
  
  #Rbar and delta
  if(is.null(Rbar)){
    if(is.null(Deltat)){
      Rbar   <- max(quantile(y, 0.95) + 1, 2)
      Deltat <- rep(1, Rbar - 1)
    } else {
      Rbar   <- length(Deltat) + 1
    }
  } else {
    if(is.null(Deltat)){
      Deltat <- rep(1, Rbar - 1)
    } else {
      if(Rbar != (length(Deltat) + 1)) stop("Rbar != length(delta) + 1")
    }
  }
  lDelta   <- log(Deltat)
  stopifnot(Rbar <= (max(y) + 1))
  stopifnot(Rbar >= 2)
  
  # theta
  if (!is.null(thetat)) {
    if(length(thetat) != (K + 1)) {
      stop("Length of theta is not suited.")
    }
    thetat <- c(log(thetat[1]/(1 -thetat[1])), thetat[-1])
  } else {
    Xtmp   <- cbind(f.t.data$Gy, X)
    b      <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    b      <- b/sqrt(sum((y - Xtmp%*%b)^2)/n)
    thetat <- c(log(max(b[1]/(1 - b[1]), 0.01)), b[-1])
  }
  thetat   <- c(thetat, lDelta)
  
  # other variables
  cont     <- TRUE
  t        <- 0
  Rmax     <- NULL
  theta    <- NULL
  covout   <- NULL
  REt      <- NULL
  llht     <- 0
  par0     <- NULL
  par1     <- NULL
  like     <- NULL
  steps    <- list()
  
  # arguments
  ctr      <- c(list(Gyb = Gybt, X = X, Rbar = Rbar, maxy = maxy, K = K, n = n, y = y), opt.ctr)
  # if(Ksimu == 1) {
  #   ctr    <- c(ctr, list(Simu1 = simu1, nsimu = nsimu))
  # }
  # if(Ksimu == 2) {
  #   ctr    <- c(ctr, list(Simu1 = simu1, Simu2 = simu2, nsimu = nsimu))
  # }
  
  if (optimizer == "optim") {
    # ctr    <- c(ctr, list(fn = ifelse(Ksimu == 0, foptimREM_NPL, 
    #                                   ifelse(Ksimu == 1, foptimREM_NPLncond1, 
    #                                          foptimREM_NPLncond2)), par = thetat))
    ctr    <- c(ctr, list(fn = foptimREM_NPL, par = thetat))
    
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else {
    # ctr    <- c(ctr, list(f = ifelse(Ksimu == 0, foptimREM_NPL, 
    #                                  ifelse(Ksimu == 1, foptimREM_NPLncond1, 
    #                                         foptimREM_NPLncond2)), p = thetat))
    ctr    <- c(ctr, list(f = foptimREM_NPL,  p = thetat))
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
      llht        <- -REt[[like]]
      theta       <- c(1/(1 + exp(-thetat[1])), thetat[2:(K +1)], exp(tail(thetat, Rbar-1)))
      
      # compute y
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n)
      
      # distance
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
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n)
      
      # distance
      dist        <- sum(abs(ctr[[par0]] - thetat)) + sum(abs(ybt0 - ybt))
      cont        <- (dist > npl.tol & t < (npl.maxit - 1))
      t           <- t + 1
      REt$dist    <- dist
      
      ctr[[par0]] <- thetat
      steps[[t]]  <- REt
    }
    llht          <- -REt[[like]]
    theta         <- c(1/(1 + exp(-thetat[1])), thetat[2:(K +1)], exp(tail(thetat, Rbar-1)))
  }
  
  names(theta)    <- c(coln, paste0("delta", 2:(Rbar)))
  
  environment(formula) <- env.formula
  sdata               <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist)),
    "nfriends"      = unlist(lapply(Glist, function(u) sum(u > 0))) 
  )
  if (!missing(data)) {
    sdata             <- c(sdata, list("data" = deparse(substitute(data))))
  }
  
  
  if (npl.maxit == t) {
    warning("The maximum number of iterations of the NPL algorithm has been reached.")
  }
  
  # covariance and ME
  tmp               <- fcovCDI(n, Gybt, thetat, X, Rbar, K, npl.S, Glist, igr, M, cov)
  meffects          <- c(tmp$meffects)
  covt              <- tmp$covt
  covm              <- tmp$covm
  Rmax              <- tmp$Rmax
  
  colnME            <- coln
  if("(Intercept)" %in% coln) {
    colnME          <- coln[-2]
    meffects        <- meffects[-2]
    covm            <- covm[-2, -2]
  } 
  names(meffects)   <- colnME
  
  if(!is.null(covt)) {
    colnames(covt)  <- c(coln, paste0("logdelta", 2:(Rbar)))
    rownames(covt)  <- c(coln, paste0("logdelta", 2:(Rbar)))
    colnames(covm)  <- colnME
    rownames(covm)  <- colnME
  }
  
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "nlinks"     = unlist(lapply(Glist, function(u) sum(u > 0))),
                               "formula"    = formula,
                               "Rbar"       = Rbar,
                               "Rmax"       = Rmax,
                               "log.like"   = llht, 
                               "npl.iter"   = t)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(theta = head(theta, K + 1), delta = tail(theta, Rbar - 1), marg.effects = meffects),
                               "yb"         = ybt, 
                               "Gyb"        = Gybt,
                               "cov"        = list(parms = covt, marg.effects = covm),
                               "details"    = steps)
  class(out)           <- "CDnetNPL"
  out
}


#' @title Summarize Count Data Model with Social Interactions
#' @description Summary and print methods for the class `CDnetNPL` as returned by the function \link{CDnetNPL}.
#' @param object an object of class `CDnetNPL`, output of the function \code{\link{CDnetNPL}}.
#' @param x an object of class `summary.CDnetNPL`, output of the function \code{\link{summary.CDnetNPL}},
#' class `summary.CDnetNPLs`, list of outputs of the function \code{\link{summary.CDnetNPL}} 
#' (when the model is estimated many times to control for the endogeneity) 
#' or class `CDnetNPL` of the function \code{\link{CDnetNPL}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
"summary.CDnetNPL" <- function(object,
                               ...) {
  stopifnot(class(object) == "CDnetNPL")
  out        <- c(object, list("..." = ...))
  if(is.null(object$cov$parms)){
    stop("Covariance was not computed")
  }
  class(out) <- "summary.CDnetNPL"
  out
}


#' @rdname summary.CDnetNPL
#' @export
"print.summary.CDnetNPL"  <- function(x, ...) {
  stopifnot(class(x) == "summary.CDnetNPL")
  
  M                    <- x$info$M
  n                    <- x$info$n
  iteration            <- x$info$npl.iter
  Rbar                 <- x$info$Rbar
  formula              <- x$info$formula
  coef                 <- x$estimate$theta
  K                    <- length(coef)
  meff                 <- x$estimate$marg.effects
  std                  <- sqrt(head(diag(x$cov$parms), K))
  std.meff             <- sqrt(diag(x$cov$marg.effects))
  delta                <- x$estimate$delta
  llh                  <- x$info$log.likelihood
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:6)], list(...))
  
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:6)], list(...))
  
  nfr                  <- x$info$nlinks
  cat("Count data Model with Social Interactions\n\n")
  cat("Call:\n")
  print(formula)
  cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("Rbar: ", Rbar, "\n")
  cat("delta: ", delta, "\n")
  cat("log pseudo-likelihood: ", llh, "\n")
  
  invisible(x)
}

#' @rdname summary.CDnetNPL
#' @export
"print.CDnetNPL" <- function(x, ...) {
  stopifnot(class(x) == "CDnetNPL")
  print(summary(x, ...))
}


#' @rdname summary.CDnetNPL
#' @importFrom stats cov
#' @export
"print.summary.CDnetNPLs" <- function(x, ...) {
  stopifnot(class(x) %in% c("list", "CDnetNPLs", "summary.CDnetNPLs")) 
  
  lclass        <- unique(unlist(lapply(x, class)))
  if (!all(lclass %in%c("CDnetNPL", "summary.CDnetNPL"))) {
    stop("All the components in `x` should be from `CDnetNPL` or `summary.CDnetNPL` class")
  }
  
  nsim          <- length(x)
  K             <- length(x[[1]]$estimate$theta)
  coef          <- do.call("rbind", lapply(x, function(z) t(c(z$estimate$theta, z$estimate$delta))))
  meff          <- do.call("rbind", lapply(x, function(z) t(z$estimate$marg.effects)))
  estimate      <- colSums(coef)/nsim
  meffects      <- colSums(meff)/nsim
  
  vcoef2        <- Reduce("+", lapply(x, function(z) z$cov$parms))/nsim
  vmeff2        <- Reduce("+", lapply(x, function(z) z$cov$marg.effects))/nsim
  
  vcoef1        <- cov(coef)
  vmeff1        <- cov(meff)
  
  vcoef         <- vcoef1 + vcoef2
  vmeff         <- vmeff1 + vmeff2
  
  llh           <- unlist(lapply(x, function(z) z$info$log.like))
  llh           <- c("min" = min(llh), "mean" = mean(llh), "max" = max(llh))
  
  M             <- x[[1]]$info$M
  n             <- x[[1]]$info$n
  
  coef          <- head(estimate, K)
  delta         <- estimate[-(1:K)]
  std           <- sqrt(head(diag(vcoef), K))
  std.meff      <- sqrt(diag(vmeff))
  
  tmp           <- fcoefficients(coef, std)
  out_print     <- tmp$out_print
  out           <- tmp$out
  
  tmp.meff       <- fcoefficients(meffects, std.meff)
  out_print.meff <- tmp.meff$out_print
  out.meff       <- tmp.meff$out
  
  out_print      <- c(list(out_print), x[[1]][-(1:6)], list(...))
  out_print.meff <- c(list(out_print.meff), x[[1]][-(1:6)], list(...))
  
  Rbar           <- x[[1]]$info$Rbar
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Replication of Nested pseudo-likelihood (NPL) \nReplication: ", nsim, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("Rbar: ", Rbar, "\n")
  cat("delta: ", delta, "\n")
  cat("log pseudo-likelihood: ", "\n")
  print(llh)
  
  out                  <- list("M"          = M,
                               "n"          = n,
                               "simulation" = nsim, 
                               "estimate"   = estimate, 
                               "likelihood" = llh, 
                               "cov"        = vcoef, 
                               "meffects"   = meffects,
                               "cov.me"     = vmeff,
                               ...          = ...)
  class(out)           <- "print.summary.CDnetNPLs"
  invisible(out)
} 