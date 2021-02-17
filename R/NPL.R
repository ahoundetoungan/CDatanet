#' @title Estimate Count Data Model with Social Interactions using NPL Method
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
#' @param yb0 (optional) expectation of y.
#' @param optimizer is either `nlm` (referring to the \link[stats]{nlm} function) or `optim` (referring to the \link[stats]{optim} function). 
#' At every step of the NPL method, the estimation is performed using \link[stats]{nlm} or \link[stats]{optim}. Other arguments 
#' of these functions such as, `control` and `method` can be defined through the argument `opt.ctr`.
#' @param npl.ctr list of controls for the NPL method (see details).
#' @param opt.ctr list of arguments of \link[stats]{nlm} or \link[stats]{optim} (the one set in `optimizer`) such as control, method, ...
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `CDnetNPL` is called.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{iteration}{number of iterations performed by the NPL algorithm.}
#'     \item{estimate}{NPL estimator.}
#'     \item{likelihood}{pseudo-likelihood value.}
#'     \item{yb}{ybar (see details), expectation of y.}
#'     \item{Gyb}{average of the expectation of y among friends.}
#'     \item{steps}{step-by-step output as returned by the optimizer.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network and name of the object `data` (see details).}
#' @details 
#' ## Model
#' Following Houndetoungan (2020), the count data \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i \bar{\mathbf{y}} + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*ybar + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.\cr
#' The count variable \eqn{y_i} is then define by the next (greater or equal) non negative integer to 
#' \eqn{y_i^*}{ys_i}; that is \eqn{y_i = 0} if  
#' \eqn{y_i^* \leq 0}{ys_i \le 0} and \eqn{y_i = q + 1} if 
#' \eqn{q < y_i^* \leq q + 1}{q < ys_i \le q + 1}, where \eqn{q} is a non-negative integer.
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
#' }
#' ## `codedata`
#' The \link[base]{class} of the output of this function is \code{CDnetNPL}. This class has a \link[base]{summary} 
#' and \link[base]{print} \link[utils]{methods} to summarize and print the results. 
#' The adjacency matrix and the data are needed to summarize the results. However, in order to save 
#' memory, the function does not return these objects. Instead, it returns `codedata` which contains among others, the `formula` 
#' and the names of these objects passed through the argument `Glist` and `data` (if provided).
#' `codedata` will be used to get access to the adjacency matrix and the data. Therefore, it is
#' important to have the adjacency matrix and the data (or the variables) available in \code{.GlobalEnv}. Otherwise,
#' it will be necessary to provide them to the \link[base]{summary} function.
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
#' ytmp    <- simCDnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
#' 
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out   <- CDnetNPL(formula = yt ~ x1 + x2, contextual = TRUE, Glist = Glist, data = data)
#' summary(out)
#' }
#' @export
CDnetNPL    <- function(formula,
                        contextual, 
                        Glist, 
                        theta0    = NULL, 
                        yb0       = NULL,
                        optimizer = "optim", 
                        npl.ctr   = list(), 
                        opt.ctr   = list(), 
                        data) {
  
  stopifnot(optimizer %in% c("optim", "nlm"))
  env.formula <- environment(formula)
  # controls
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  
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
      
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, K, n)
      
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
      ybt0        <- ctr$yb + 0
      REt         <- do.call(get(optimizer), ctr)
      thetat      <- REt[[par1]]
      
      fL_NPL(ybt, Gybt, Glist, igr, M, X, thetat, K, n)
      
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
  out                 <- list("M"          = M,
                          "n"          = n,
                          "iteration"  = t, 
                          "estimate"   = theta, 
                          "likelihood" = llht, 
                          "yb"         = ybt, 
                          "Gyb"        = Gybt,
                          "steps"      = steps,
                          "codedata"   = sdata)
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
#' @param cov.ctr list of control values for the covariance containing two integers, `R` and `S`. The covariance summations from `0`
#' to infinity. But the summed elements decreases exponentially. The summations are approximated by summations from `0` to `R`.
#' The covariance also requires computing \eqn{\Phi(x) - \Phi(x - 1)}, where \eqn{\Phi} is 
#' the normal'  probability density function. This is done using important sampling, where `S` numbers are generated
#' form the uniform distribution.
#' @param Glist the adjacency matrix or list sub-adjacency matrix. If missing make, sure that 
#' the object provided to the function \code{\link{CDnetNPL}} is available in \code{.GlobalEnv} (see detail - codedata section of \code{\link{CDnetNPL}}).
#' @param data dataframe containing the explanatory variables. If missing make, sure that 
#' the object provided to the function \code{\link{CDnetNPL}} is available in \code{.GlobalEnv} (see detail - codedata section of \code{\link{CDnetNPL}}).
#' @param ... further arguments passed to or from other methods.
#' @return A list consisting of:
#'     \item{M}{number of sub-networks.}
#'     \item{n}{number of individuals in each network.}
#'     \item{iteration}{number of iterations performed by the NPL algorithm.}
#'     \item{estimate}{NPL estimator.}
#'     \item{likelihood}{pseudo-likelihood value.}
#'     \item{yb}{ybar (see details), expectation of y.}
#'     \item{Gyb}{average of the expectation of y among friends.}
#'     \item{steps}{step-by-step output as returned by the optimizer.}
#'     \item{cov}{covariance matrix of the estimate.}
#'     \item{meffects}{vector of marginal effects.}
#'     \item{cov.me}{covariance matrix of the marginal effects.}
#'     \item{cov.ctr}{returned value of the control values for the covariance.}
#'     \item{codedata}{list of formula, name of the object `Glist`, number of friends in the network and name of the object `data`.}
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @importFrom stats runif
#' @export 
"summary.CDnetNPL" <- function(object,
                               cov.ctr   = list(),
                               Glist,
                               data,
                               ...) {
  stopifnot(class(object) == "CDnetNPL")
  if ((missing(Glist) & !missing(data)) | (!missing(Glist) & missing(data))) {
    stop("Glist is missing while data is provided or vice versa")
  }
  codedata         <- object$codedata
  formula          <- as.formula(codedata$formula)
  
  if (missing(Glist)) {
    Glist           <- get(codedata$Glist, envir = .GlobalEnv)
  } else {
    if(!is.list(Glist)) {
      Glist         <- list(Glist)
    }
  }
  
  theta             <- object$estimate
  Gyb               <- object$Gyb
  
  J                 <- length(theta)
  
  lambda            <- theta[1]
  b                 <- theta[2:(J - 1)]
  sigma             <- theta[J]
  
  
  M                 <- length(Glist)
  nvec              <- unlist(lapply(Glist, nrow))
  n                 <- sum(nvec)
  igr               <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  
  if(missing(data)) {
    if (is.null(codedata$data)) {
      data          <- environment(formula)
    } else {
      data          <- get(codedata$data, envir = .GlobalEnv)
    }
  } 
  
  f.t.data      <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = 0)
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
  m2d           <- c(tmp[[6]])
  
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
  
  # marginal effect
  meand           <- mean(d)
  meanbz          <- colSums(b*Z)/n
  meanm2d         <- mean(m2d)
  
  meff            <- theta[-J]*meand
  tmp3            <- diag(J - 1)*meand - (theta[-J]/sigma) %*% matrix(meanbz, nrow = 1)
  tmp4            <- (meanm2d/sigma^3  - meand/sigma)*theta[-J]
  tmp5            <- cbind(tmp3, tmp4)
  tmp6            <- tmp5 %*% tmp2 %*% t(tmp5)
  
  covout            <- tmp2
  colnames(covout)  <- c(coln, "sigma")
  rownames(covout)  <- c(coln, "sigma")
  
  covmeff           <- tmp6
  colnames(covmeff) <- coln
  rownames(covmeff) <- coln
  
  cov.ctr           <- list("R0" = R0, "S0" = S0)
  
  out               <- c(object[-9], 
                         list("cov"       = covout, 
                              "meffects"  = meff,
                              "cov.me"    = covmeff,
                              "cov.ctr"   = cov.ctr, 
                              "codedata"  = codedata,
                              "..."       = ...)) 
  
  class(out)        <- "summary.CDnetNPL"
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
  meff                 <- x$meffects
  std                  <- sqrt(diag(x$cov)[-K])
  std.meff             <- sqrt(diag(x$cov.me))
  sigma                <- estimate[K]
  llh                  <- x$likelihood
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:13)], list(...))
  
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:13)], list(...))
  
  nfr                  <- x$codedata$nfriends
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log pseudo-likelihood: ", llh, "\n")
  
  out                  <- c(x[1:4], list("coefficients" = out, "coefficients.me" = out.meff), x[-(1:4)])
  class(out)           <- "print.summary.CDnetNPL"
  invisible(out)
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
  stopifnot(class(x) %in% c("list", "summary.CDnetNPLs", "print.summary.CDnetNPLs")) 
  
  type2               <- (class(x) == "print.summary.CDnetNPLs")
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
    if (!all(lclass %in% "summary.CDnetNPL")) {
      stop("All the components in `x` should be from `summary.CDnetNPL` class")
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
    out_print         <- c(list(out_print), x[[1]][-(1:13)], list(...))
    out_print.meff    <- c(list(out_print.meff), x[[1]][-(1:13)], list(...))
  }
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Replication of Nested pseudo-likelihood (NPL) \nReplication: ", nsim, "\n\n")
  
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
  class(out)           <- "print.summary.CDnetNPLs"
  invisible(out)
} 