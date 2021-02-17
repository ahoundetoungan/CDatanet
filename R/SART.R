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
#' @param optimizer is either `nlm` (referring to the function \link[stats]{nlm}) or `optim` (referring to the function \link[stats]{optim}). 
#' Other arguments 
#' of these functions such as, the control values and the method can be defined through the argument `opt.ctr`.
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param print a boolean indicating if the estimate should be printed at each step.
#' @param cov a boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `SARTML` is called.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network, name of the object `data`,
#'      and number of zeros in `y` (see details).}
#' @details 
#' ## Model
#' The left-censored variable \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.\cr
#' The count variable \eqn{y_i} is then define that is \eqn{y_i = 0} if  
#' \eqn{y_i^* \leq 0}{ys_i \le 0} and \eqn{y_i = y_i^*}{y_i = ys_i} otherwise.
#' ## `codedata`
#' The \link[base]{class} of the output of this function is \code{SARTML}. This class has a \link[base]{summary} 
#' and \link[base]{print} \link[utils]{methods} to summarize and print the results. 
#' The adjacency matrix is needed to summarize the results. However, in order to save 
#' memory, the function does not return it. Instead, it returns `codedata` which contains the `formula` 
#' and the name of the adjacency matrix passed through the argument `Glist`.
#' `codedata` will be used to get access to the adjacency matrix. Therefore, it is
#' important to have the adjacency matrix available in \code{.GlobalEnv}. Otherwise
#' it will be necessary to provide the adjacency matrix to the \link[base]{summary} 
#' and \link[base]{print} functions.
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
  
  f.t.data   <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = theta0)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X), "sigma")
  
  K          <- ncol(X)
  
  # variables
  
  ylist      <- lapply(1:M, function(x) y[(igr[x,1]:igr[x,2]) + 1])
  
  indpos        <- (y > 1e-323)
  indzero       <- !indpos
  idpos         <- which(indpos) - 1
  idzero        <- which(indzero) - 1

  idposlis   <- lapply(ylist, function(w) which(w > 0))
  Npos       <- unlist(lapply(idposlis, length))
  
  G2list     <- lapply(1:M, function(w) Glist[[w]][idposlis[[w]], idposlis[[w]]])
  Gy         <- unlist(lapply(1:M, function(w) Glist[[w]] %*% ylist[[w]]))
  I2list     <- lapply(Npos, diag)
  Ilist      <- lapply(nvec, diag)  
  Wlist      <- lapply(1:M, function(x) (indpos[(igr[x,1]:igr[x,2]) + 1] %*% t(indpos[(igr[x,1]:igr[x,2]) + 1]))*Glist[[x]])
  
  
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
  
  if(exists("alphatde")) rm("alphatde")
  if(exists("logdetA2")) rm("logdetA2")
  alphatde   <- Inf
  logdetA2   <- 0
  
  
  # arguments
  ctr        <- c(list("X" = X, "G2" = G2list, "I2" = I2list, "K" = K, "y" = y, "Gy" = Gy,
                       "idpos" = idpos, "idzero" = idzero, "Npos" = Npos, "ngroup" = M,
                       "alphatilde" = alphatde, "logdetA2" = logdetA2, "N" = n, 
                       "I" = Ilist,  "W" = Wlist, "igroup" = igr, "indzero" = indzero,
                       "indpos" = indpos), opt.ctr)
  
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
      # fgoptimTobit <- function(theta, X, G2, I2, K, y, Gy, idpos, idzero, Npos,
      #                          ngroup, alphatilde, logdetA2, N, I,  W, igroup,
      #                          indzero, indpos) {
      #   res                   <- foptimTobit(theta = theta, X = X, logdetA2 = logdetA2,
      #                                        alphatilde = alphatilde, G2 = G2, I2 = I2,
      #                                        K = K, y = y, Gy = Gy, idpos = idpos,
      #                                        idzero = idzero, Npos = Npos, ngroup = ngroup,
      #                                        I = I, W = W, N = N, igroup = igroup,
      #                                        indzero = indzero, indpos = indpos)
      #   attr(res, "gradient") <- fgradvecTobit(theta = theta, X = X,logdetA2 = logdetA2,
      #                                          alphatilde = alphatilde, G2 = G2, I2 = I2,
      #                                          K = K, y = y, Gy = Gy, idpos = idpos,
      #                                          idzero = idzero, Npos = Npos, ngroup = ngroup,
      #                                          I = I, W = W, N = N, igroup = igroup,
      #                                          indzero = indzero, indpos = indpos)
      #   res
      # }
      ctr         <- c(ctr, list(f = foptimTobit)) 
    } else {
      # fgoptimTobit <- function(theta, X, G2, I2, K, y, Gy, idpos, idzero, Npos,
      #                          ngroup, alphatilde, logdetA2, N, I,  W, igroup) {
      #   res                   <- foptimTobit0(theta = theta, X = X, logdetA2 = logdetA2,
      #                                         alphatilde = alphatilde, G2 = G2, I2 = I2,
      #                                         K = K, y = y, Gy = Gy, idpos = idpos,
      #                                         idzero = idzero, Npos = Npos, ngroup = ngroup,
      #                                         I = I, W = W, N = N, igroup = igroup,
      #                                         indzero = indzero, indpos = indpos)
      #   attr(res, "gradient") <- c(fgradvecTobit(theta = theta, X = X, logdetA2 = logdetA2,
      #                                            alphatilde = alphatilde, G2 = G2, I2 = I2,
      #                                            K = K, y = y, Gy = Gy, idpos = idpos,
      #                                            idzero = idzero, Npos = Npos, ngroup = ngroup,
      #                                            I = I, W = W, N = N, igroup = igroup,
      #                                            indzero = indzero, indpos = indpos))
      #   res
      # }
      ctr    <- c(ctr, list(f = foptimTobit0)) 
    }
  }
  
  
  
  resTO    <- do.call(get(optimizer), ctr)
  theta    <- resTO[[par1]]
  theta    <- c(1/(1 + exp(-theta[1])), theta[2:(K + 1)], exp(theta[K + 2]))
  llh      <- -resTO[[like]]
  
  covout             <- NULL
  rm(list = c("alphatde", "logdetA2"))
  if (cov) {
    covtmp           <- fqTobit(theta, X, G2list, Ilist, Wlist, K, n, y, Gy,  indzero,
                                indpos, igr, M)
    covout           <- solve(stats::cov(covtmp))/n
    colnames(covout) <- coln
    rownames(covout) <- coln
  }
  
  names(theta)       <- coln
  
  # sdata
  environment(formula) <- env.formula
  sdata                <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist)),
    "nfriends"      = unlist(lapply(Glist, function(u) sum(u > 0)))  
  )
  if (!missing(data)) {
    sdata              <- c(sdata, list("data" = deparse(substitute(data))))
  } 
  sdata                <- c(sdata, list("pzeros" = sum(indzero)/n))
  
  out                  <- list("M"             = M,
                             "n"             = n,
                             "estimate"      = theta, 
                             "likelihood"    = llh, 
                             "cov"           = covout,
                             "optimization"  = resTO,
                             "codedata"      = sdata,
                            "y"              = y)
  class(out)            <- "SARTML"
  out
  
}



#' @title Summarize SART Model
#' @description Summary and print methods for the class `SARTML` as returned by the function \link{SARTML}.
#' @param object an object of class `SARTML`, output of the function \code{\link{SARTML}}.
#' @param x an object of class `summary.SARTML`, output of the function \code{\link{summary.SARTML}} 
#' or class `SARTML`, output of the function \code{\link{SARTML}}.
#' @param Glist the adjacency matrix or list sub-adjacency matrix. If missing make, sure that 
#' the object provided to the function \code{\link{SARTML}} is available in \code{.GlobalEnv} (see detail - codedata section of \code{\link{SARTML}}).
#' @param data dataframe containing the explanatory variables. If missing make, sure that 
#' the object provided to the function \code{\link{SARTML}} is available in \code{.GlobalEnv} (see detail - codedata section of \code{\link{SARTML}}).
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{Maximum Likelihood (ML) estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network, name of the object `data`, 
#'     and number of zeros in `y`.}
#' @export 
#' @param ... further arguments passed to or from other methods.
#' @export 
"summary.SARTML" <- function(object,
                             Glist,
                             data,
                             ...) {
  stopifnot(class(object) == "SARTML")
  
  if ((missing(Glist) & !missing(data)) | (!missing(Glist) & missing(data))) {
    stop("Glist is missing while data is provided or vice versa")
  }
  
  if(is.null(object$cov)){
    stop("Covariance was not computed")
  }
  
  # Glist and data
  codedata      <- object$codedata
  formula       <- as.formula(codedata$formula)
  
  if (missing(Glist)) {
    Glist       <- get(codedata$Glist, envir = .GlobalEnv)
  } else {
    if(!is.list(Glist)) {
      Glist     <- list(Glist)
    }
  }
  
  
  # data 
  theta         <- object$estimate
  cov           <- object$cov
  
  J             <- length(theta)
  
  lambda        <- theta[1]
  b             <- theta[2:(J - 1)]
  sigma         <- theta[J]
  
  
  M             <- length(Glist)
  nvec          <- unlist(lapply(Glist, nrow))
  n             <- sum(nvec)
  igr           <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  
  if(missing(data)) {
    if (is.null(codedata$data)) {
      data      <- environment(formula)
    } else {
      data      <- get(codedata$data, envir = .GlobalEnv)
    }
  } 
  
  f.t.data      <- formula.to.data(formula, FALSE, Glist, M, igr, data)
  
  y             <- f.t.data$y
  X             <- f.t.data$X
  Gy            <- f.t.data$Gy
  coln          <- c("lambda", colnames(X))
  
  Z             <- cbind(Gy, X)
  Zdst          <- c((Z %*% theta[-J])/sigma)
  
  phiZdst       <- dnorm(Zdst)
  PhiZdst       <- pnorm(Zdst)
  meanPhiZdst   <- mean(PhiZdst)
  meanphiZdstz  <- colSums(phiZdst*Z)/n
  
  # marginal effect
  meff              <- theta[-J]*meanPhiZdst
  tmp1              <- diag(J - 1)*meanPhiZdst + theta[-J] %*% matrix(meanphiZdstz, nrow = 1)/sigma
  tmp2              <- - mean(Zdst*phiZdst)*theta[-J]/sigma
  tmp3              <- cbind(tmp1, tmp2)
  
  covmeff           <- tmp3 %*% cov %*% t(tmp3)
  
  colnames(covmeff) <- coln
  rownames(covmeff) <- coln
  
  out           <- c(object[1:5], 
                     list("meffects" = meff,
                          "cov.me"   = covmeff),
                     object[-(1:5)],
                     list("..."       = ...)) 
  
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
  meff                 <- x$meffects
  std                  <- sqrt(diag(x$cov)[-K])
  std.meff             <- sqrt(diag(x$cov.me))
  sigma                <- estimate[K]
  llh                  <- x$likelihood
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:9)], list(...))
  
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:9)], list(...))
  

  nfr                  <- x$codedata$nfriends
  cat("SART Model\n\n")
  cat("Method: Maximum Likelihood (ML)", "\n\n")
  
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n")
  
  cat("Proportion of zeros      : ", x$codedata$pzeros, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
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


#' @rdname summary.SARTML
#' @export
"print.summary.SARTMLs" <- function(x, ...) {
  stopifnot(class(x) %in% c("list", "summary.SARTMLs", "print.summary.SARTMLs")) 
  
  type2               <- (class(x) == "print.summary.SARTMLs")
  nsim                <- NULL
  estimate            <- NULL
  meffects            <- NULL
  vcoef               <- NULL
  vmeff               <- NULL
  llh                 <- NULL
  n                   <- NULL
  M                   <- NULL
  
  if (type2) {
    nsim              <- x$simulation
    estimate          <- x$estimate
    meffects          <- x$meffects
    vcoef             <- x$cov
    vmeff             <- x$cov.me
    llh               <- x$likelihood
    n                 <- x$n
    M                 <- x$M
  } else {
    lclass            <- unique(unlist(lapply(x, class)))
    if (!all(lclass %in% "summary.SARTML")) {
      stop("All the components in `x` should be from `summary.SARTML` class")
    }
    
    nsim              <- length(x)
    coef              <- do.call("rbind", lapply(x, function(z) t(z$estimate)))
    meff              <- do.call("rbind", lapply(x, function(z) t(z$meffects)))
    estimate          <- colSums(coef)/nsim
    meffects          <- colSums(meff)/nsim
    
    vcoef2            <- Reduce("+", lapply(x, function(z) z$cov))/nsim
    vmeff2            <- Reduce("+", lapply(x, function(z) z$cov.me))/nsim
    
    vcoef1            <- cov(coef)
    vmeff1            <- cov(meff)
    
    vcoef             <- vcoef1 + vcoef2
    vmeff             <- vmeff1 + vmeff2
    
    
    llh               <- unlist(lapply(x, function(z) z$likelihood))
    llh               <- c("min" = min(llh), "mean" = mean(llh), "max" = max(llh))
    
    M                 <- x[[1]]$M
    n                 <- x[[1]]$n
  }
  
  
  
  K                   <- length(estimate)
  coef                <- estimate[-K]
  std                 <- sqrt(diag(vcoef)[-K])
  std.meff            <- sqrt(diag(vmeff))
  sigma               <- estimate[K]
  
  tmp                 <- fcoefficients(coef, std)
  out_print           <- tmp$out_print
  out                 <- tmp$out
  
  tmp.meff            <- fcoefficients(meffects, std.meff)
  out_print.meff      <- tmp.meff$out_print
  out.meff            <- tmp.meff$out
  
  
  if (type2) {
    out_print         <- c(list(out_print), x[-(1:8)], list(...))
    out_print.meff    <- c(list(out_print.meff), x[-(1:8)], list(...))
  } else {
    out_print         <- c(list(out_print), x[[1]][-(1:7)], list(...))
    out_print.meff    <- c(list(out_print.meff), x[[1]][-(1:7)], list(...))
  }
  
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
                               "estimate"   = estimate, 
                               "likelihood" = llh, 
                               "cov"        = vcoef, 
                               "meffects"   = meffects,
                               "cov.me"     = vmeff,
                               ...          = ...)
  class(out)           <- "print.summary.SARTMLs"
  invisible(out)
} 