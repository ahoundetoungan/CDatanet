#' @title Estimate SART model
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param theta0 (optional) starting value of \eqn{\theta = (\lambda, \beta, \gamma, \sigma)}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param optimizer is either `nlm` (refering to the function \link[stats]{nlm}) or `optim` (refering to the function \link[stats]{optim}). 
#' At every step of the NPL method, the estimation is performed using \link[stats]{nlm} or \link[stats]{optim}. Other arguments 
#' of these functions such as, the control values and the method can be defined through the argument `opt.ctr`.
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param print a boolean indicating if the estimate shoul be print a each step.
#' @param cov a boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `CDnetNPL` is called.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{NPL estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, formula's environment, names of the objects Glist and data (this is useful for summarizing the results, see details).}
#' @details 
#' ## Model
#' The left-censored variable \eqn{\mathbf{y}}{y} are generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
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
#' \dontrun{
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
#' ytmp    <- simCDnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                     theta = theta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
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
  
  sdata           <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist))
  )
  if (!missing(data)) {
    sdata         <- c(sdata, list("data" = deparse(substitute(data))))
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
#' @param object an object of class `SARTML`, output of the function \code{\link{SARTML}}.
#' @param x an object of class `summary.SARTML`, output of the function \code{\link{summary.SARTML}} 
#' or class `SARTML`, output of the the function \code{\link{SARTML}}.
#' @param Glist the adjacency matrix or list sub-adjacency matrix. If missing make, sure that 
#' the object provided to the function \code{\link{SARTML}} is available in \code{.GlobalEnv} (see detail - codedata section of \code{\link{SARTML}}).
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{estimate}{NPL estimator.}
#'     \item{likelihood}{likelihood value.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{optimization}{output as returned by the optimizer.}
#'     \item{codedata}{list of formula, formula's environment, names of the objects Glist and data (this is useful for summarizing the results, see details).}
#' @export 
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
"print.summary.SARTML"  <- function(x,  Glist, ...) {
  stopifnot(class(x) == "summary.SARTML")
  
  M                    <- x$M
  n                    <- x$n
  estimate             <- x$estimate
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  std                  <- sqrt(diag(x$cov))
  sigma                <- estimate[K]
  llh                  <- x$likelihood

  if (missing(Glist)) {
    Glist              <- get(x$codedata$Glist, envir = .GlobalEnv) 
  } else {
    if(!is.list(Glist)) {
      Glist            <- list(Glist)
    }
  }
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:7)], list(...))
  
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