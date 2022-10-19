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
#' @param estim.rho indicates if the parameter \eqn{\rho} should be estimated or set to zero.
#' @param starting (optional) starting value of \eqn{\theta = (\lambda, \beta', \gamma')'}, \eqn{\bar{\delta}}{deltabar}, \eqn{\delta = (\delta_2, ..., \delta_{\bar{R}})}{\delta = (\delta_2, ..., \delta_{Rbar})}, and \eqn{\rho}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param yb0 (optional) expectation of y.
#' @param optimizer is either `fastlbfgs` (L-BFGS optimization method of the package \pkg{RcppNumerical}), `nlm` (referring to the function \link[stats]{nlm}), or `optim` (referring to the function \link[stats]{optim}). 
#' Other arguments 
#' of these functions such as, `control` and `method` can be defined through the argument `opt.ctr`.
#' @param npl.ctr list of controls for the NPL method (see details).
#' @param opt.ctr list of arguments to be passed in `optim_lbfgs` of the package \pkg{RcppNumerical}, \link[stats]{nlm} or \link[stats]{optim} (the solver set in `optimizer`), such as `maxit`, `eps_f`, `eps_g`, `control`, `method`, ...
#' @param cov a Boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `cdnet` is called.
#' @return A list consisting of:
#'     \item{info}{list of general information about the model.}
#'     \item{estimate}{NPL estimator.}
#'     \item{yb}{ybar (see details), expectation of y.}
#'     \item{Gyb}{average of the expectation of y among friends.}
#'     \item{cov}{list of covariance matrices.}
#'     \item{details}{step-by-step output as returned by the optimizer.}
#' @description
#' `cdnet` is used to estimate peer effects on counting data with rational expectations (see details). The model is presented in Houndetoungan (2022). 
#' @details 
#' ## Model
#' Following Houndetoungan (2022), the count data \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i \mathbf{E}(\bar{\mathbf{y}}|\mathbf{X},\mathbf{G})  + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*E(ybar|X, G) + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, 1)}{\epsilon_i --> N(0, 1)}.\cr
#' Then, \eqn{y_i = r} iff \eqn{a_r \leq y_i^* \leq a_{r+1}}{a_r \le ys_i \le a_{r + 1}}, where
#' \eqn{a_0 = -\inf}{a_0 = -Inf}, \eqn{a_1 = 0}, \eqn{a_r = \sum_{k = 1}^r\delta_k}{a_r = \delta_1 + ... + \delta_r}. 
#' The parameter are subject to the constraints \eqn{\delta_r \geq \lambda}{\delta_r \ge \lambda} if \eqn{1 \leq r \leq \bar{R}}{1 \le r \le Rbar},  and
#' \eqn{\delta_r = (r - \bar{R})^{\rho}\bar{\delta} + \lambda}{a_r = deltabar*(r - Rbar)^{\rho} + \lambda} if \eqn{r \geq \bar{R} + 1}{r \ge Rbar + 1}. The unknown parameters to be estimated are
#' \eqn{\lambda}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta_2}, ..., \eqn{\delta_{\bar{R}}}{\delta_{Rbar}}, \eqn{\bar{\delta}}{deltabar}, and \eqn{\rho}.
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
#' @references 
#' Houndetoungan, E. A. (2022). Count Data Models with Social Interactions under Rational Expectations. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @seealso \code{\link{sart}}, \code{\link{sar}}, \code{\link{simcdnet}}.
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
#' delta  <- c(1, 0.87, 0.75, 0.6)
#' delbar <- 0.05
#' rho    <- 0.5
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
#' ytmp    <- simcdnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta,
#'                     deltabar = delbar, delta = delta, rho = rho, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))
#' 
#' data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
#' rm(list = ls()[!(ls() %in% c("Glist", "data"))])
#' 
#' out   <- cdnet(formula = yt ~ x1 + x2, contextual = TRUE, Glist = Glist, 
#'                data = data, Rbar = 6, estim.rho = TRUE)
#' summary(out)
#' }
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom utils tail
#' @export
cdnet    <- function(formula,
                     contextual, 
                     Glist, 
                     Rbar      = NULL,
                     estim.rho = FALSE,
                     starting  = list(theta = NULL, deltabar = NULL, delta = NULL, rho = NULL), 
                     yb0       = NULL,
                     optimizer = "fastlbfgs", 
                     npl.ctr   = list(), 
                     opt.ctr   = list(), 
                     cov       = TRUE,
                     data) {
  
  stopifnot(optimizer %in% c("fastlbfgs", "optim", "nlm"))
  env.formula <- environment(formula)
  
  # controls
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  npl.S       <- npl.ctr$S
  npl.incdit  <- npl.ctr$incdit
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
  if (is.null(npl.incdit)) {
    npl.incdit<- 30L
  }
  
  # starting values
  thetat      <- starting$theta
  Deltat      <- starting$delta
  Deltabart   <- starting$deltabar
  rhot        <- starting$rho
  tmp         <- c(is.null(thetat), is.null(Deltat), is.null(Deltabart))
  if(estim.rho){
    tmp       <- c(tmp, is.null(rhot))
  }
  if(all(tmp) != any(tmp)){
    stop("All the parameters in starting should be either NULL or defined. Set delta = numeric() if Rbar = 1.")
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
  
  # Rbar and delta
  if(is.null(Rbar)){
    if(is.null(Deltat)){
      Rbar   <- max(quantile(y, 0.95), 1)
    } else {
      Rbar   <- length(Deltat) + 1
    }
  } else {
    if(!is.null(Deltat)){
      if(Rbar != (length(Deltat) + 1)) stop("Rbar != length(delta) + 1")
    }
  }
  stopifnot(Rbar <= max(y))
  stopifnot(Rbar >= 1)
  
  # check starting and computed them if not defined
  lmbd0       <- NULL  
  if (!is.null(thetat)) {
    if(length(thetat) != (K + 1)) {
      stop("Length of theta is not suited.")
    }
    if(thetat[1] <= 0) stop("Starting lambda is not positive.")
    lmbd0     <- thetat[1]
    thetat    <- c(log(thetat[1]/(1 -thetat[1])), thetat[-1])
    if(Rbar > 1) if(any(Deltat < lmbd0)) stop("Multiple equilibrium issue: Starting lambda is greater than starting delta.")
    if(Deltabart <= 0) stop("Starting deltabar is not positive.")
    if(estim.rho){
      if(rhot <= 0) stop("Starting rho is not positive.")
    }
  } else {
    Xtmp      <- cbind(f.t.data$Gy, X)
    b         <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    b         <- b/sqrt(sum((y - Xtmp%*%b)^2)/n)
    b[1]      <- max(b[1], 0.01)
    lmbd0     <- b[1]
    Deltat    <- rep(b[1] + 0.01, Rbar - 1)
    thetat    <- c(log(b[1]/(1 - b[1])), b[-1])
    Deltabart <- 0.01
    rhot      <- 0.01
  }
  
  lDelta      <- log(Deltat - lmbd0)
  lDeltabart  <- log(Deltabart)
  thetat      <- c(thetat, lDelta, lDeltabart)
  if(estim.rho){
    lrho      <- log(rhot)
    thetat    <- c(thetat, lrho)
  }
  
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
  var.comp <- NULL
  steps    <- list()
  
  # arguments
  ctr      <- c(list(Gyb = Gybt, X = X, Rbar = Rbar, maxy = maxy, K = K, n = n, y = y), opt.ctr)
  # if(Ksimu == 1) {
  #   ctr    <- c(ctr, list(Simu1 = simu1, nsimu = nsimu))
  # }
  # if(Ksimu == 2) {
  #   ctr    <- c(ctr, list(Simu1 = simu1, Simu2 = simu2, nsimu = nsimu))
  # }
  
  # Arguments used in the optimizer
  if (optimizer == "fastlbfgs"){
    ctr    <- c(ctr, list(par = thetat)); 
    if(estim.rho){optimizer <- "cdnetLBFGSrho"} else{optimizer <- "cdnetLBFGS"}
    
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else if (optimizer == "optim") {
    # ctr    <- c(ctr, list(fn = ifelse(Ksimu == 0, foptimREM_NPL, 
    #                                   ifelse(Ksimu == 1, foptimREM_NPLncond1, 
    #                                          foptimREM_NPLncond2)), par = thetat))
    if(estim.rho){
      ctr  <- c(ctr, list(fn = foptimREM_NPL2, par = thetat))
    } else{
      ctr  <- c(ctr, list(fn = foptimREM_NPL, par = thetat))
    }
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else {
    # ctr    <- c(ctr, list(f = ifelse(Ksimu == 0, foptimREM_NPL, 
    #                                  ifelse(Ksimu == 1, foptimREM_NPLncond1, 
    #                                         foptimREM_NPLncond2)), p = thetat))
    if(estim.rho){
      ctr  <- c(ctr, list(f = foptimREM_NPL2,  p = thetat))
    } else{
      ctr  <- c(ctr, list(f = foptimREM_NPL,  p = thetat))
    }
    par0   <- "p"
    par1   <- "estimate"
    like   <- "minimum"
  }
  
  # Fonction used to compute new candidate for E(y)
  fL_NPLR   <- NULL
  fcovCDIR  <- NULL
  if(estim.rho){
    fL_NPLR <- fL_NPL2
    fcovCDI <- fcovCDI2
  } else{
    fL_NPLR <- fL_NPL
    fcovCDI <- fcovCDI
  }
  
  # NPL algorithm
  ninc.d   <- 0
  dist0    <- Inf
  if(npl.print) {
    while(cont) {
      # tryCatch({
        ybt0        <- ybt + 0    #copy in different memory
        # compute theta
        REt         <- do.call(get(optimizer), ctr)
        thetat      <- REt[[par1]]
        llht        <- -REt[[like]]
        theta       <- c(1/(1 + exp(-thetat[1])), thetat[2:(K +1)], exp(thetat[-(1:(K +1))]))
        if(Rbar > 1){
          theta[(K+2):(K+Rbar)] <- theta[(K+2):(K+Rbar)] + theta[1]
        }
        # compute y
        fL_NPLR(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n)
        # distance
        print(theta)
        dist        <- max(abs(c((ctr[[par0]] - thetat)/thetat, (ybt0 - ybt)/ybt)), na.rm = TRUE)
        ninc.d      <- (ninc.d + 1)*(dist > dist0) #counts the successive number of times distance increases
        dist0       <- dist
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
        if((ninc.d > npl.incdit) | (llht < -1e293)) {
          cat("** Non-convergence ** Redefining theta and computing a new yb\n")
          thetat[1] <- -4.5
          dist0     <- Inf
          if(estim.rho){
            fnewyb2(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
          } else{
            fnewyb(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
          }
          ctr[[par0]] <- thetat
        }
      # },
      # error = function(e){
      #   cat("** Non-convergence ** Redefining theta and computing a new yb\n")
      #   thetat[1]   <- -4.5
      #   dist0       <- Inf
      #   if(estim.rho){
      #     fnewyb2(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
      #   } else{
      #     fnewyb(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
      #   }
      #   cont        <- TRUE
      #   t           <- t + 1
      #   ctr[[par0]] <- thetat
      #   steps[[t]]  <- NULL
      # })
    }
  } else {
    while(cont) {
      # tryCatch({
        ybt0        <- ybt + 0    #copy in different memory
        
        # compute theta
        REt         <- do.call(get(optimizer), ctr)
        thetat      <- REt[[par1]]
        llht        <- -REt[[like]]

        # compute y
        fL_NPLR(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n)
        
        # distance
        dist        <- max(abs(c((ctr[[par0]] - thetat)/thetat, (ybt0 - ybt)/ybt)), na.rm = TRUE)
        ninc.d      <- (ninc.d + 1)*(dist > dist0) #counts the successive number of times distance increases
        dist0       <- dist
        cont        <- (dist > npl.tol & t < (npl.maxit - 1))
        t           <- t + 1
        REt$dist    <- dist
        ctr[[par0]] <- thetat
        steps[[t]]  <- REt
        
        if((ninc.d > npl.incdit) | (llht < -1e293)) {
          thetat[1] <- -4.5
          dist0     <- Inf
          if(estim.rho){
            fnewyb2(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
          } else{
            fnewyb(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
          }
          ctr[[par0]] <- thetat
        }
      # },
      # error = function(e){
      #   thetat[1]   <- -4.5
      #   dist0       <- Inf
      #   if(estim.rho){
      #     fnewyb2(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
      #   } else{
      #     fnewyb(ybt, Gybt, Glist, igr, M, X, thetat, Rbar, K, n, npl.tol, npl.maxit)
      #   }
      #   cont        <- TRUE
      #   t           <- t + 1
      #   ctr[[par0]] <- thetat
      #   steps[[t]]  <- NULL
      # })
    }
    theta           <- c(1/(1 + exp(-thetat[1])), thetat[2:(K +1)], exp(thetat[-(1:(K +1))]))
    if(Rbar > 1){
      theta[(K+2):(K+Rbar)] <- theta[(K+2):(K+Rbar)] + theta[1]
    }
  }
  
  # name theta
  namtheta          <- coln
  if(Rbar > 1){
    namtheta        <- c(namtheta, paste0("delta", 2:(Rbar)))
  }
  namtheta          <- c(namtheta, "deltabar")
  if(estim.rho){
    namtheta        <- c(namtheta, "rho")
  }
  names(theta)      <- namtheta
  
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
  var.comp          <- tmp$var.comp
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
    namecovt        <- coln
    if(Rbar > 1){
      namecovt      <- c(namecovt, paste0("log(delta", 2:(Rbar), ")"))
    }
    namecovt        <- c(namecovt, "log(deltabar)")
    if(estim.rho){
      namecovt      <- c(namecovt, "log(rho)")
    }
    colnames(covt)  <- namecovt
    rownames(covt)  <- namecovt
    colnames(covm)  <- colnME
    rownames(covm)  <- colnME
    namecovt[1]     <- "log(lambda)"
    rownames(var.comp$Sigma) <- namecovt
    rownames(var.comp$Omega) <- namecovt
    colnames(var.comp$Sigma) <- namecovt
    colnames(var.comp$Omega) <- namecovt
  }
  
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "Kz"         = K,
                               "nlinks"     = unlist(lapply(Glist, function(u) sum(u > 0))),
                               "formula"    = formula,
                               "estim.rho"  = estim.rho,
                               "Rbar"       = Rbar,
                               "Rmax"       = Rmax,
                               "log.like"   = llht, 
                               "npl.iter"   = t)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(parms = theta, marg.effects = meffects),
                               "yb"         = ybt, 
                               "Gyb"        = Gybt,
                               "cov"        = list(parms = covt, marg.effects = covm, var.comp = var.comp),
                               "details"    = steps)
  class(out)           <- "cdnet"
  out
}


#' @title Summarize Count Data Model with Social Interactions
#' @description Summary and print methods for the class `cdnet` as returned by the function \link{cdnet}.
#' @param object an object of class `cdnet`, output of the function \code{\link{cdnet}}.
#' @param x an object of class `summary.cdnet`, output of the function \code{\link{summary.cdnet}},
#' class `summary.cdnets`, list of outputs of the function \code{\link{summary.cdnet}} 
#' (when the model is estimated many times to control for the endogeneity) 
#' or class `cdnet` of the function \code{\link{cdnet}}.
#' @param Glist adjacency matrix or list sub-adjacency matrix. This is not necessary if the covariance method was computed in \link{cdnet}.
#' @param data a `dataframe` containing the explanatory variables. This is not necessary if the covariance method was computed in \link{cdnet}.
#' @param S number of simulation to be used to compute integral in the covariance by important sampling.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
"summary.cdnet" <- function(object,
                            Glist,
                            data,
                            S = 1e3L,
                            ...) {
  stopifnot(class(object) == "cdnet")
  out        <- c(object, list("..." = ...))
  if(is.null(object$cov$parms)){
    env.formula <- environment(object$info$formula)
    Rbar        <- object$info$Rbar
    Kz          <- object$info$Kz
    estim.rho   <- object$info$estim.rho
    parms       <- object$estimate$parms

    thetat      <- c(log(parms[1]/(1 -parms[1])), parms[-1])
    if(Rbar > 1){
      thetat[(2 + Kz):(Kz + Rbar)] <- thetat[(2 + Kz):(Kz + Rbar)] - parms[1]
    }
    thetat[(2 + Kz):(Kz + Rbar + 1 + estim.rho)] <- log(tail(thetat, Rbar + estim.rho))

    if(estim.rho){
      fcovCDI   <- fcovCDI2
    } else{
      fcovCDI   <- fcovCDI
    }

    Gybt        <- object$Gyb
    npl.S       <- S
    if (is.null(npl.S)) {
      npl.S     <- 1e3L
    }
    contextual  <- FALSE
    formula     <- object$info$formula
    if (!is.list(Glist)) {
      Glist     <- list(Glist)
    }
    M        <- length(Glist)
    nvec     <- unlist(lapply(Glist, nrow))
    n        <- sum(nvec)
    igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
    
    f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, theta0 = thetat)
    formula  <- f.t.data$formula
    y        <- f.t.data$y
    X        <- f.t.data$X
    K        <- ncol(X)
    coln     <- c("lambda", colnames(X))
    

    tmp      <- fcovCDI(n, Gybt, thetat, X, Rbar, K, npl.S, Glist, igr, M, TRUE)
    meffects          <- c(tmp$meffects)
    var.comp          <- tmp$var.comp
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
      namecovt        <- coln
      if(Rbar > 1){
        namecovt      <- c(namecovt, paste0("log(delta", 2:(Rbar), ")"))
      }
      namecovt        <- c(namecovt, "log(deltabar)")
      if(estim.rho){
        namecovt      <- c(namecovt, "log(rho)")
      }
      colnames(covt)  <- namecovt
      rownames(covt)  <- namecovt
      colnames(covm)  <- colnME
      rownames(covm)  <- colnME
      namecovt[1]     <- "log(lambda)"
      rownames(var.comp$Sigma) <- namecovt
      rownames(var.comp$Omega) <- namecovt
      colnames(var.comp$Sigma) <- namecovt
      colnames(var.comp$Omega) <- namecovt
    }
    
    out$cov           <- list(parms = covt, marg.effects = covm, var.comp = var.comp)
  }
  class(out) <- "summary.cdnet"
  out
}


#' @rdname summary.cdnet
#' @importFrom stats pchisq
#' @export
"print.summary.cdnet"  <- function(x, ...) {
  stopifnot(class(x) == "summary.cdnet")
  
  M                    <- x$info$M
  n                    <- x$info$n
  iteration            <- x$info$npl.iter
  Rbar                 <- x$info$Rbar
  formula              <- x$info$formula
  Kz                   <- x$info$Kz
  estim.rho            <- x$info$estim.rho
  
  parms                <- x$estimate$parms
  coef                 <- parms[1:(1 + Kz)]
  K                    <- length(coef)
  meff                 <- x$estimate$marg.effects
  std                  <- sqrt(head(diag(x$cov$parms), K))
  std.meff             <- sqrt(diag(x$cov$marg.effects))
  delta                <- tail(parms, Rbar + estim.rho)
  
  llh                  <- x$info$log.like
  
  
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
  if(Rbar > 1){
    cat("delta: ", delta[1:(Rbar - 1)], "\n")
  }
  cat("deltabar: ", delta[Rbar])
  if(estim.rho){
    cat(" -- rho: ", delta[Rbar+1], "\n")
    cat("Wald test H0: rho = 0, Prob = ", 1 - pchisq(1/x$cov$parms[Kz+Rbar+2,Kz+Rbar+2], df = 1), "\n")
  } else{
    cat("\n")
  }

  cat("log pseudo-likelihood: ", llh, "\n")
  
  invisible(x)
}

#' @rdname summary.cdnet
#' @export
"print.cdnet" <- function(x, ...) {
  stopifnot(class(x) == "cdnet")
  print(summary(x, ...))
}

#' @rdname summary.cdnet
#' @importFrom stats cov
#' @export
"summary.cdnets" <- function(object, ...) {
  stopifnot(class(object) %in% c("list", "cdnets", "summary.cdnets")) 
  
  lclass        <- unique(unlist(lapply(object, class)))
  if (!all(lclass %in%c("summary.cdnet"))) {
    stop("All the components in `object` should be from `summary.cdnet` class")
  }
  
  nsim          <- length(object)
  coef          <- do.call("rbind", lapply(object, function(z) t(z$estimate$parms)))
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
  Rbar          <- object[[1]]$info$Rbar
  Kz            <- object[[1]]$info$Kz
  estim.rho     <- object[[1]]$info$estim.rho
  
  INFO                 <- list("M"          = M,
                               "n"          = n,
                               "Kz"         = Kz,
                               "estim.rho"  = estim.rho,
                               "log.like"   = llh,
                               "Rbar"       = Rbar,
                               "simulation" = nsim)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(parms = coef, marg.effects = meffects),
                               "cov"        = list(parms = vcoef, marg.effects = vmeff),
                               ...          = ...)
  class(out)           <- "summary.cdnets"
  out
} 


#' @rdname summary.cdnet
#' @importFrom stats cov
#' @export
"print.summary.cdnets" <- function(x, ...) {
  stopifnot(class(x) %in% c("summary.cdnets")) 

  nsim          <- x$info$simulation
  
  Kz            <- x$info$Kz
  estim.rho     <- x$info$estim.rho
  
  parms         <- x$estimate$parms
  coef          <- parms[1:(1 + Kz)]
  K             <- length(coef)
  meff          <- x$estimate$marg.effects
  std           <- sqrt(head(diag(x$cov$parms), K))
  std.meff      <- sqrt(diag(x$cov$marg.effects))
  delta         <- tail(parms, Rbar + estim.rho)

  vcoef         <- x$cov$parms
  vmeff         <- x$cov$marg.effects
  
  llh           <- x$info$log.like
  
  M             <- x$info$M
  n             <- x$info$n
  Rbar          <- x$info$Rbar
  
  std           <- sqrt(head(diag(vcoef), K))
  std.meff      <- sqrt(diag(vmeff))
  
  tmp           <- fcoefficients(coef, std)
  out_print     <- tmp$out_print
  out           <- tmp$out
  
  tmp.meff       <- fcoefficients(meff, std.meff)
  out_print.meff <- tmp.meff$out_print
  out.meff       <- tmp.meff$out
  
  out_print      <- c(list(out_print), x[-c(1:3)], list(...))
  out_print.meff <- c(list(out_print.meff), x[-c(1:3)], list(...))
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Method: Replication of Nested pseudo-likelihood (NPL) \nReplication: ", nsim, "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("Rbar: ", Rbar, "\n")
  if(Rbar > 1){
    cat("delta: ", delta[1:(Rbar - 1)], "\n")
  }
  cat("deltabar: ", delta[Rbar])
  if(estim.rho){
    cat(" -- rho: ", delta[Rbar+1], "\n")
    cat("Wald test H0: rho = 0, Prob = ", 1 - pchisq(1/x$cov$parms[Kz+Rbar+2,Kz+Rbar+2], df = 1), "\n")
  } else{
    cat("\n")
  }
  
  cat("log pseudo-likelihood: ", "\n")
  print(llh)
  
  invisible(x)
} 

#' @rdname summary.cdnet
#' @importFrom stats cov
#' @export
"print.cdnets" <- function(x, ...) { 
  stopifnot(class(x) %in% c("cdnets", "list"))
  print(summary.cdnets(x, ...))
} 