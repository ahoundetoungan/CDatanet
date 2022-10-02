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
#' @param print a Boolean indicating if the estimate should be printed at each step.
#' @param cov a Boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @description
#' `sar` is used to estimate peer effects continuous variables (see details). The model is presented in Lee(2004). 
#' @details 
#' ## Model
#' The variable \eqn{\mathbf{y}}{y} is given for all i as
#' \deqn{y_i = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{y_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.
#' @references  
#' Lee, L. F. (2004). Asymptotic distributions of quasi‐maximum likelihood estimators for spatial autoregressive models. \emph{Econometrica}, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}.
#' @seealso \code{\link{sart}}, \code{\link{cdnet}}, \code{\link{simsar}}.
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
#' ytmp    <- simsar(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                   theta = theta, data = data) 
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
#' 
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out     <- sar(formula = yt ~ x1 + x2, contextual = TRUE, 
#'                  Glist = Glist, optimizer = "optim", data = data)
#' summary(out)
#' }
#' @return A list consisting of:
#'     \item{info}{list of general information on the model.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{details}{outputs as returned by the optimizer.}
#' @export
sar <- function(formula,
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
  ctr        <- c(list(X = X,invXX = invXX, G = Glist, I = Ilist, n = n,
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
  
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "nlinks"     = unlist(lapply(Glist, function(u) sum(u > 0))),
                               "formula"    = formula,
                               "log.like"   = llh)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = theta, 
                               "cov"        = covout,
                               "details"    = resSAR)
  class(out)           <- "sar"
  out
}




#' @title Summarize SAR Model
#' @description Summary and print methods for the class `sar` as returned by the function \link{sar}.
#' @param object an object of class `sar`, output of the function \code{\link{sar}}.
#' @param x an object of class `summary.sar`, output of the function \code{\link{summary.sar}} or 
#' class `sar`, output of the function \code{\link{sar}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @param ... further arguments passed to or from other methods.
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
  print(formula)
  cat("\nMethod: Maximum Likelihood (ML)", "\n\n")
  
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

#' @rdname summary.sar
#' @importFrom stats cov
#' @export
"summary.sars" <- function(object, ...) {
  stopifnot(class(object) %in% c("list", "sars", "summary.sars")) 
  
  lclass        <- unique(unlist(lapply(object, class)))
  if (!all(lclass %in%c("summary.sar"))) {
    stop("All the components in `object` should be from `summary.sar` class")
  }
  
  nsim          <- length(object)
  K             <- length(object[[1]]$estimate$theta)
  coef          <- do.call("rbind", lapply(object, function(z) t(c(z$estimate$theta))))
  meff          <- do.call("rbind", lapply(object, function(z) t(z$estimate$marg.effects)))
  estimate      <- colSums(coef)/nsim
  meffects      <- colSums(meff)/nsim
  
  vcoef2        <- Reduce("+", lapply(object, function(z) z$cov$parms))/nsim
  vmeff2        <- Reduce("+", lapply(object, function(z) z$cov$marg.effects))/nsim
  
  vcoef1        <- cov(coef)
  vmeff1        <- cov(meff)
  
  vcoef         <- vcoef1 + vcoef2
  vmeff         <- vmeff1 + vmeff2
  
  llh           <- unlist(lapply(object, function(z) z$info$log.like))
  llh           <- c("min" = min(llh), "mean" = mean(llh), "max" = max(llh))
  
  M             <- object[[1]]$info$M
  n             <- object[[1]]$info$n
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "log.like"   = llh,
                               "simulation" = nsim)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(theta = estimate, marg.effects = meffects),
                               "cov"        = list(parms = vcoef, marg.effects = vmeff),
                               ...          = ...)
  class(out)           <- "summary.sars"
  out
} 


#' @rdname summary.sar
#' @importFrom stats cov
#' @export
"print.summary.sars" <- function(x, ...) {
  stopifnot(class(x) %in% c("summary.sars")) 
  
  nsim          <- x$info$simulation
  coef          <- x$estimate$theta
  meffects      <- x$estimate$marg.effects
  K             <- length(coef)
  sigma         <- tail(coef, 1)
  coef          <- head(coef, K - 1)
  RE            <- x$info$Rat.Exp
  
  vcoef         <- x$cov$parms
  vmeff         <- x$cov$marg.effects
  
  llh           <- x$info$log.like
  
  M             <- x$info$M
  n             <- x$info$n
  
  std           <- sqrt(head(diag(vcoef), K-1))
  std.meff      <- sqrt(diag(vmeff))
  
  tmp           <- fcoefficients(coef, std)
  out_print     <- tmp$out_print
  out           <- tmp$out
  
  tmp.meff       <- fcoefficients(meffects, std.meff)
  out_print.meff <- tmp.meff$out_print
  out.meff       <- tmp.meff$out
  
  out_print      <- c(list(out_print), x[-c(1:3)], list(...))
  out_print.meff <- c(list(out_print.meff), x[-c(1:3)], list(...))
  RE             <- FALSE
  cat("sar Model", ifelse(RE, "with Rational Expectation", ""), "\n\n")
  cat("Method: Replication of Maximum Likelihood (ML) \nReplication: ", nsim, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: \n")
  print(llh)
  
  invisible(x)
} 

#' @rdname summary.sar
#' @importFrom stats cov
#' @export
"print.sars" <- function(x, ...) { 
  stopifnot(class(x) %in% c("sars", "list"))
  print(summary.sars(x, ...))
} 