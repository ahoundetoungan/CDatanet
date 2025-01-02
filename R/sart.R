#' @title Simulating Data from Tobit Models with Social Interactions
#' @description `simsart` simulates censored data with social interactions (see Xu and Lee, 2015).
#' @param formula a class object \code{\link[stats]{formula}}: a symbolic description of the model. 
#' `formula` must be, for example, \code{y ~ x1 + x2 + gx1 + gx2}, where `y` is the endogenous vector, 
#' and `x1`, `x2`, `gx1`, and `gx2` are control variables. These can include contextual variables, 
#' i.e., averages among the peers. Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist The network matrix. For networks consisting of multiple subnets, `Glist` can be a list 
#' of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes 
#' in the `m`-th subnet.
#' @param theta a vector defining the true value of \eqn{\theta = (\lambda, \Gamma, \sigma)} (see the model specification in the details).
#' @param tol the tolerance value used in the fixed-point iteration method to compute `y`. The process stops 
#' if the \eqn{\ell_1}-distance between two consecutive values of `y` is less than `tol`.
#' @param maxit the maximum number of iterations in the fixed-point iteration method.
#' @param cinfo a Boolean indicating whether information is complete (`cinfo = TRUE`) or incomplete (`cinfo = FALSE`). 
#' In the case of incomplete information, the model is defined under rational expectations.
#' @param data an optional data frame, list, or environment (or object coercible by \code{\link[base]{as.data.frame}} 
#' to a data frame) containing the variables in the model. If not found in `data`, the variables are taken 
#' from \code{environment(formula)}, typically the environment from which `simsart` is called.
#' @details 
#' For a complete information model, the outcome \eqn{y_i} is defined as:
#' \deqn{\begin{cases}
#' y_i^{\ast} = \lambda \bar{y}_i + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
#' y_i = \max(0, y_i^{\ast}),
#' \end{cases}}
#' where \eqn{\bar{y}_i} is the average of \eqn{y} among peers, 
#' \eqn{\mathbf{z}_i} is a vector of control variables, 
#' and \eqn{\epsilon_i \sim N(0, \sigma^2)}. \cr
#' 
#' In the case of incomplete information models with rational expectations, \eqn{y_i} is defined as:
#' \deqn{\begin{cases}
#' y_i^{\ast} = \lambda E(\bar{y}_i) + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
#' y_i = \max(0, y_i^{\ast}).
#' \end{cases}}
#' 
#' @return A list consisting of:
#' \describe{
#'   \item{yst}{\eqn{y^{\ast}}, the latent variable.}
#'   \item{y}{The observed censored variable.}
#'   \item{Ey}{\eqn{E(y)}, the expected value of \eqn{y}.}
#'   \item{Gy}{The average of \eqn{y} among peers.}
#'   \item{GEy}{The average of \eqn{E(y)} among peers.}
#'   \item{meff}{A list including average and individual marginal effects.}
#'   \item{iteration}{The number of iterations performed per sub-network in the fixed-point iteration method.}
#' }
#' @references 
#' Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
#' @seealso \code{\link{sart}}, \code{\link{simsar}}, \code{\link{simcdnet}}.
#' @examples 
#' \donttest{
#' # Define group sizes
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 200)) # Number of nodes per sub-group
#' n      <- sum(nvec) # Total number of nodes
#' 
#' # Define parameters
#' lambda <- 0.4
#' Gamma  <- c(2, -1.9, 0.8, 1.5, -1.2)
#' sigma  <- 1.5
#' theta  <- c(lambda, Gamma, sigma)
#' 
#' # Generate covariates (X)
#' X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))
#' 
#' # Construct network adjacency matrices
#' G      <- list()
#' for (m in 1:M) {
#'   nm           <- nvec[m] # Nodes in sub-group m
#'   Gm           <- matrix(0, nm, nm) # Initialize adjacency matrix
#'   max_d        <- 30 # Maximum degree
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1)) # Random connections
#'     Gm[i, tmp] <- 1
#'   }
#'   rs           <- rowSums(Gm) # Normalize rows
#'   rs[rs == 0]  <- 1
#'   Gm           <- Gm / rs
#'   G[[m]]       <- Gm
#' }
#' 
#' # Prepare data
#' data   <- data.frame(X, peer.avg(G, cbind(x1 = X[, 1], x2 = X[, 2])))
#' colnames(data) <- c("x1", "x2", "gx1", "gx2") # Add column names
#' 
#' # Complete information game simulation
#' ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, 
#'                    Glist = G, theta = theta, 
#'                    data = data, cinfo = TRUE)
#' data$yc <- ytmp$y # Add simulated outcome to the dataset
#' 
#' # Incomplete information game simulation
#' ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, 
#'                    Glist = G, theta = theta, 
#'                    data = data, cinfo = FALSE)
#' data$yi <- ytmp$y # Add simulated outcome to the dataset
#' }
#' @importFrom Rcpp sourceCpp
#' @export
simsart   <- function(formula,
                      Glist,
                      theta,
                      tol   = 1e-15,
                      maxit = 500,
                      cinfo = TRUE,
                      data) {
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  
  M        <- length(Glist)
  nvec     <- unlist(lapply(Glist, nrow))
  n        <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  
  f.t.data <- formula.to.data(formula, FALSE, Glist, M, igr, data, "sim", 0)
  X        <- f.t.data$X
  
  K        <- length(theta)
  if(K != (ncol(X) + 2)) {
    stop("Length of theta is not suited.")
  }
  lambda   <- theta[1]
  b        <- theta[2:(K - 1)]
  sigma    <- theta[K ]
  
  
  
  xb       <- c(X %*% b)
  eps      <- rnorm(n, 0, sigma)
  
  yst      <- numeric(n)
  y        <- NULL
  Gy       <- NULL
  Ey       <- NULL
  GEy      <- NULL
  t        <- NULL
  Ztl      <- rep(0, n)
  if(cinfo){
    y      <- rep(0, n)
    Gy     <- rep(0, n)
    t      <- fyTobit(yst, y, Gy, Ztl, Glist, eps, igr, M, xb, n, lambda, tol, maxit)
  } else {
    Ey     <- rep(0, n)
    GEy    <- rep(0, n)
    t      <- fEytbit(Ey, GEy, Glist, igr, M, xb, lambda, sigma, n, tol, maxit)
    Ztl    <- lambda*GEy + xb
    yst    <- Ztl + eps
    y      <- yst*(yst > 0)
  }
  
  # marginal effects
  coln      <- c("lambda", colnames(X))
  thetaWI   <- head(theta, K - 1)
  if("(Intercept)" %in% coln) {
    thetaWI <- thetaWI[-2]
    coln    <- coln[-2]
  }
  imeff     <- pnorm(Ztl/sigma) %*% t(thetaWI); colnames(imeff) <- coln
  
  
  list("yst"       = yst,
       "y"         = y,
       "Ey"        = Ey,
       "Gy"        = Gy,
       "GEy"       = GEy,
       "meff"      = list(ameff = apply(imeff, 2, mean), imeff = imeff),
       "iteration" = t)
}


#' @title Estimating Tobit Models with Social Interactions
#' @description
#' `sart` estimates Tobit models with social interactions based on the framework of Xu and Lee (2015). 
#' The method allows for modeling both complete and incomplete information scenarios in networks, incorporating rational expectations in the latter case.
#' 
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. The formula must follow the structure, 
#'   e.g., \code{y ~ x1 + x2 + gx1 + gx2}, where `y` is the endogenous variable, and `x1`, `x2`, `gx1`, and `gx2` are control variables. 
#'   Control variables may include contextual variables, such as peer averages, which can be computed using \code{\link{peer.avg}}.
#' @param Glist The network matrix. For networks consisting of multiple subnets, `Glist` can be a list, where the `m`-th element is 
#'   an `ns*ns` adjacency matrix representing the `m`-th subnet, with `ns` being the number of nodes in that subnet.
#' @param starting (Optional) A vector of starting values for \eqn{\theta = (\lambda, \Gamma, \sigma)}, where:
#'   \itemize{
#'     \item \eqn{\lambda} is the peer effect coefficient,
#'     \item \eqn{\Gamma} is the vector of control variable coefficients,
#'     \item \eqn{\sigma} is the standard deviation of the error term.
#'   }
#' @param Ey0 (Optional) A starting value for \eqn{E(y)}.
#' @param optimizer The optimization method to be used. Choices are:
#'   \itemize{
#'     \item `"fastlbfgs"`: L-BFGS optimization method from the \pkg{RcppNumerical} package,
#'     \item `"nlm"`: Refers to the \link[stats]{nlm} function,
#'     \item `"optim"`: Refers to the \link[stats]{optim} function.
#'   }
#'   Additional arguments for these functions, such as `control` and `method`, can be specified through the `opt.ctr` argument.
#' @param npl.ctr A list of controls for the NPL (Nested Pseudo-Likelihood) method (refer to the details in \code{\link{cdnet}}).
#' @param opt.ctr A list of arguments to be passed to the chosen solver (`fastlbfgs`, \link[stats]{nlm}, or \link[stats]{optim}), 
#'   such as `maxit`, `eps_f`, `eps_g`, `control`, `method`, etc.
#' @param cov A Boolean indicating whether to compute the covariance matrix (\code{TRUE} or \code{FALSE}).
#' @param cinfo A Boolean indicating whether the information structure is complete (\code{TRUE}) or incomplete (\code{FALSE}). 
#'   Under incomplete information, the model is defined with rational expectations.
#' @param data An optional data frame, list, or environment (or object coercible by \link[base]{as.data.frame}) containing the variables
#'   in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `sart` is called.
#' 
#' @return A list containing:
#' \describe{
#'   \item{\code{info}}{General information about the model.}
#'   \item{\code{estimate}}{The Maximum Likelihood (ML) estimates of the parameters.}
#'   \item{\code{Ey}}{\eqn{E(y)}, the expected values of the endogenous variable.}
#'   \item{\code{GEy}}{The average of \eqn{E(y)} among peers.}
#'   \item{\code{cov}}{A list including covariance matrices (if \code{cov = TRUE}).}
#'   \item{\code{details}}{Additional outputs returned by the optimizer.}
#' }
#' @details 
#' For a complete information model, the outcome \eqn{y_i} is defined as:
#' \deqn{\begin{cases}
#' y_i^{\ast} = \lambda \bar{y}_i + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
#' y_i = \max(0, y_i^{\ast}),
#' \end{cases}}
#' where \eqn{\bar{y}_i} is the average of \eqn{y} among peers, 
#' \eqn{\mathbf{z}_i} is a vector of control variables, 
#' and \eqn{\epsilon_i \sim N(0, \sigma^2)}. \cr
#' 
#' In the case of incomplete information models with rational expectations, \eqn{y_i} is defined as:
#' \deqn{\begin{cases}
#' y_i^{\ast} = \lambda E(\bar{y}_i) + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
#' y_i = \max(0, y_i^{\ast}).
#' \end{cases}}
#' @seealso \code{\link{sar}}, \code{\link{cdnet}}, \code{\link{simsart}}.
#' @references 
#' Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
#' @examples 
#' \donttest{
#' # Group sizes
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 200))
#' n      <- sum(nvec)
#' 
#' # Parameters
#' lambda <- 0.4
#' Gamma  <- c(2, -1.9, 0.8, 1.5, -1.2)
#' sigma  <- 1.5
#' theta  <- c(lambda, Gamma, sigma)
#' 
#' # Covariates (X)
#' X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))
#' 
#' # Network creation
#' G      <- list()
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
#'   Gm           <- Gm / rs
#'   G[[m]]       <- Gm
#' }
#' 
#' # Data creation
#' data   <- data.frame(X, peer.avg(G, cbind(x1 = X[, 1], x2 = X[, 2])))
#' colnames(data) <- c("x1", "x2", "gx1", "gx2")
#' 
#' ## Complete information game
#' ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, theta = theta, 
#'                    data = data, cinfo = TRUE)
#' data$yc <- ytmp$y
#' 
#' ## Incomplete information game
#' ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, theta = theta, 
#'                    data = data, cinfo = FALSE)
#' data$yi <- ytmp$y
#' 
#' # Complete information estimation for yc
#' outc1   <- sart(formula = yc ~ x1 + x2 + gx1 + gx2, optimizer = "nlm",
#'                 Glist = G, data = data, cinfo = TRUE)
#' summary(outc1)
#' 
#' # Complete information estimation for yi
#' outc1   <- sart(formula = yi ~ x1 + x2 + gx1 + gx2, optimizer = "nlm",
#'                 Glist = G, data = data, cinfo = TRUE)
#' summary(outc1)
#' 
#' # Incomplete information estimation for yc
#' outi1   <- sart(formula = yc ~ x1 + x2 + gx1 + gx2, optimizer = "nlm",
#'                 Glist = G, data = data, cinfo = FALSE)
#' summary(outi1)
#' 
#' # Incomplete information estimation for yi
#' outi1   <- sart(formula = yi ~ x1 + x2 + gx1 + gx2, optimizer = "nlm",
#'                 Glist = G, data = data, cinfo = FALSE)
#' summary(outi1)
#' }
#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' @export
sart <- function(formula,
                 Glist,
                 starting = NULL,
                 Ey0  = NULL,
                 optimizer = "fastlbfgs",
                 npl.ctr  = list(), 
                 opt.ctr = list(),
                 cov = TRUE,
                 cinfo = TRUE,
                 data) {
  stopifnot(optimizer %in% c("fastlbfgs", "optim", "nlm"))
  if(cinfo & optimizer == "fastlbfgs"){
    stop("fastlbfgs is only implemented for the rational expectation model in this version. Use another solver.")
  }
  env.formula <- environment(formula)
  # controls
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  
  #size
  if (!is.list(Glist)) {
    Glist    <- list(Glist)
  }
  M          <- length(Glist)
  nvec       <- unlist(lapply(Glist, nrow))
  n          <- sum(nvec)
  igr        <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  f.t.data   <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = starting)
  formula    <- f.t.data$formula
  y          <- f.t.data$y
  X          <- f.t.data$X
  coln       <- c("lambda", colnames(X))
  
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
  if (!is.null(starting)) {
    if(length(starting) != (K + 2)) {
      stop("Length of starting is not suited.")
    }
    thetat   <- c(log(starting[1]/(1 -starting[1])), starting[2:(K+1)], log(starting[K+2]))
  } else {
    Xtmp     <- cbind(f.t.data$Gy, X)
    b        <- solve(t(Xtmp)%*%Xtmp, t(Xtmp)%*%y)
    s        <- sqrt(sum((y - Xtmp%*%b)^2)/n)
    thetat   <- c(log(max(b[1]/(1 - b[1]), 0.01)), b[-1], log(s))
  }
  
  llh        <- NULL
  resTO      <- list()
  tmp        <- NULL
  covm       <- NULL
  var.comp   <- NULL
  t          <- NULL
  Eyt        <- NULL
  GEyt       <- NULL
  theta      <- NULL
  Ztlambda   <- NULL
  
  
  if(cinfo){
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
                         "alphatilde" = alphatde, "logdetA2" = logdetA2, "n" = n, 
                         "I" = Ilist,  "W" = Wlist, "igroup" = igr), opt.ctr)
    if (optimizer == "optim") {
      ctr    <- c(ctr, list(par = thetat))
      par1   <- "par"
      like   <- "value"
      ctr    <- c(ctr, list(fn = foptimTobit)) 
    } else {
      ctr    <- c(ctr, list(p = thetat))
      par1   <- "estimate"
      like   <- "minimum"
      ctr    <- c(ctr, list(f = foptimTobit)) 
    }
    
    resTO    <- do.call(get(optimizer), ctr)
    theta    <- resTO[[par1]]
    
    tmp      <- fcovSTC(theta = theta, X = X, G2 = G2list, I = Ilist, W = Wlist, K =  K, n = n, 
                        y = y, Gy = Gy, indzero = indzero, indpos = indpos, igroup = igr, 
                        ngroup = M, ccov = cov)
    
    theta    <- c(1/(1 + exp(-theta[1])), theta[2:(K + 1)], exp(theta[K + 2]))
    llh      <- -resTO[[like]]
    rm(list = c("alphatde", "logdetA2"))
  } else{
    # Ey 
    Eyt      <- rep(0, n)
    if (!is.null(Ey0)) {
      Eyt    <- Ey0
    }
    # GEyt
    GEyt     <- unlist(lapply(1:M, function(x) Glist[[x]] %*% Eyt[(igr[x,1]:igr[x,2])+1]))
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
    ctr      <- c(list("yidpos" = yidpos, "GEy" = GEyt, "X" = X, "npos" = sum(npos), 
                       "idpos" = idpos, "idzero" = idzero, "K" = K), opt.ctr)
    
    if (optimizer == "fastlbfgs"){
      ctr    <- c(ctr, list(par = thetat)); optimizer = "sartLBFGS"
      
      par0   <- "par"
      par1   <- "par"
      like   <- "value"
    } else if (optimizer == "optim") {
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
    
    while(cont) {
      # tryCatch({
      Eyt0        <- Eyt + 0    #copy in different memory
      
      # compute theta
      # print(optimizer)
      REt         <- do.call(get(optimizer), ctr)
      thetat      <- REt[[par1]]
      llh         <- -REt[[like]]
      
      theta       <- c(1/(1 + exp(-thetat[1])), thetat[2:(K + 1)], exp(thetat[K + 2]))
      
      # compute y
      fLTBT_NPL(Eyt, GEyt, Glist, X, thetat, igr, M, n, K)
      
      # distance
      # dist        <- max(abs(c(ctr[[par0]]/thetat, Eyt0/(Eyt + 1e-50)) - 1), na.rm = TRUE)
      dist        <- max(abs(c((ctr[[par0]] - thetat)/thetat, (Eyt0 - Eyt)/Eyt)), na.rm = TRUE)
      cont        <- (dist > npl.tol & t < (npl.maxit - 1))
      t           <- t + 1
      REt$dist    <- dist
      ctr[[par0]] <- thetat
      resTO[[t]]  <- REt
      
      if(npl.print){
        cat("---------------\n")
        cat(paste0("Step          : ", t), "\n")
        cat(paste0("Distance      : ", round(dist,3)), "\n")
        cat(paste0("Likelihood    : ", round(llh,3)), "\n")
        cat("Estimate:", "\n")
        print(theta)
      }
    }
    
    if (npl.maxit == t) {
      warning("The maximum number of iterations of the NPL algorithm has been reached.")
    }
    tmp        <- fcovSTI(n = n, GEy = GEyt, theta = thetat, X = X, K = K, G = Glist,
                          igroup = igr, ngroup = M, ccov = cov)
  }
  
  names(theta) <- c(coln, "sigma")
  
  # Marginal effects
  imeff        <- tmp$meff; colnames(imeff) <- coln
  indexWI      <- 1:(K + 1)
  if("(Intercept)" %in% coln){
    indexWI    <- indexWI[coln != "(Intercept)"]
  }
  imeff        <- imeff[, indexWI, drop = FALSE]
  meff         <- list(ameff = apply(imeff, 2, mean), imeff = imeff)
  
  # Covariances
  covm         <- tmp$covm
  covt         <- tmp$covt

  if(cov) {
    covm            <- covm[indexWI, indexWI]
    colnames(covt)  <- c(coln, "sigma")
    rownames(covt)  <- c(coln, "sigma")
    colnames(covm)  <- coln[indexWI]
    rownames(covm)  <- coln[indexWI]

  }
  
  INFO              <- list("M"          = M,
                            "n"          = n,
                            "formula"    = formula,
                            "nlinks"     = unlist(lapply(Glist, function(u) sum(u > 0))),
                            "censured"   = sum(indzero),
                            "uncensured" = n - sum(indzero),
                            "log.like"   = llh, 
                            "npl.iter"   = t)
  
  out               <- list("info"       = INFO,
                            "estimate"   = theta,
                            "meff"       = meff,
                            "Ey"         = Eyt, 
                            "GEy"        = GEyt,
                            "cov"        = list(parms = covt, ameff = covm),
                            "details"    = resTO)
  class(out)        <- "sart"
  out
}


#' @title Summary for the Estimation of Tobit Models with Social Interactions
#' @description Summary and print methods for the class `sart` as returned by the function \link{sart}.
#' @param object an object of class `sart`, output of the function \code{\link{sart}}.
#' @param x an object of class `summary.sart`, output of the function \code{\link{summary.sart}} 
#' or class `sart`, output of the function \code{\link{sart}}.
#' @param Glist adjacency matrix or list sub-adjacency matrix. This is not necessary if the covariance method was computed in \link{cdnet}.
#' @param data dataframe containing the explanatory variables. This is not necessary if the covariance method was computed in \link{cdnet}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
"summary.sart" <- function(object,
                           Glist,
                           data,
                           ...) {
  stopifnot(class(object) == "sart")
  out           <- c(object, list("..." = ...)) 
  if(is.null(object$cov$parms)){
    env.formula <- environment(object$info$formula)
    thetat      <- object$estimate
    thetat      <- c(log(thetat[1]/(1 - thetat[1])), thetat[-c(1, length(thetat))], log(thetat[length(thetat)]))
    GEyt        <- object$GEy
    formula     <- object$info$formula
    if (!is.list(Glist)) {
      Glist     <- list(Glist)
    }
    M        <- length(Glist)
    nvec     <- unlist(lapply(Glist, nrow))
    n        <- sum(nvec)
    igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
    
    f.t.data <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = thetat)
    formula  <- f.t.data$formula
    y        <- f.t.data$y
    X        <- f.t.data$X
    K        <- ncol(X)
    coln     <- c("lambda", colnames(X))
    tmp      <- NULL
    if(is.null(GEyt)){
      indpos     <- (y > 1e-323)
      indzero    <- !indpos
      idpos      <- which(indpos) - 1
      idzero     <- which(indzero) - 1
      ylist      <- lapply(1:M, function(x) y[(igr[x,1]:igr[x,2]) + 1])
      idposlis   <- lapply(ylist, function(w) which(w > 0))
      npos       <- unlist(lapply(idposlis, length))  
      G2list     <- lapply(1:M, function(w) Glist[[w]][idposlis[[w]], idposlis[[w]]])
      Gy         <- unlist(lapply(1:M, function(w) Glist[[w]] %*% ylist[[w]]))
      I2list     <- lapply(npos, diag)
      Ilist      <- lapply(nvec, diag)  
      Wlist      <- lapply(1:M, function(x) (indpos[(igr[x,1]:igr[x,2]) + 1] %*% t(indpos[(igr[x,1]:igr[x,2]) + 1]))*Glist[[x]])
      tmp        <- fcovSTC(theta = thetat, X = X, G2 = G2list, I = Ilist, W = Wlist, K =  K, n = n, 
                            y = y, Gy = Gy, indzero = indzero, indpos = indpos, igroup = igr, 
                            ngroup = M, ccov = TRUE)
    } else {
      tmp        <- fcovSTI(n = n, GEy = GEyt, theta = thetat, X = X, K = K, G = Glist,
                        igroup = igr, ngroup = M, ccov = TRUE)
    }
    
    # Marginal effects
    imeff        <- tmp$meff; colnames(imeff) <- coln
    indexWI      <- 1:(K + 1)
    if("(Intercept)" %in% coln){
      indexWI    <- indexWI[coln != "(Intercept)"]
    }
    imeff        <- imeff[, indexWI, drop = FALSE]
    meff         <- list(ameff = apply(imeff, 2, mean), imeff = imeff)
    
    # Covariances
    covm         <- tmp$covm
    covt         <- tmp$covt
    
    covm            <- covm[indexWI, indexWI]
    colnames(covt)  <- c(coln, "sigma")
    rownames(covt)  <- c(coln, "sigma")
    colnames(covm)  <- coln[indexWI]
    rownames(covm)  <- coln[indexWI]
    out$cov         <- list(parms = covt, ameff = covm)
  }
  class(out)    <- "summary.sart"
  out
}


#' @rdname summary.sart
#' @export
"print.summary.sart"  <- function(x, ...) {
  stopifnot(class(x) == "summary.sart")
  
  M                    <- x$info$M
  n                    <- x$info$n
  estimate             <- x$estimate
  iteration            <- x$info$npl.iter
  RE                   <- !is.null(iteration)
  formula              <- x$info$formula
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  meff                 <- x$meff$ameff
  std                  <- sqrt(diag(x$cov$parms)[-K])
  std.meff             <- sqrt(diag(x$cov$ameff))
  sigma                <- estimate[K]
  llh                  <- x$info$log.like
  censored             <- x$info$censured
  uncensored           <- x$info$uncensured
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:7)], list(...))
  
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:7)], list(...))
  
  
  nfr                  <- x$info$nlinks
  cat("sart Model", ifelse(RE, "with Rational Expectation", ""), "\n\n")
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

#' @rdname summary.sart
#' @export
"print.sart" <- function(x, ...) {
  stopifnot(class(x) == "sart")
  print(summary(x, ...))
}