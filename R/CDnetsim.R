#' @title Simulate data from Count Data Model with Social Interactions
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param theta the true value of the vector \eqn{\theta = (\lambda, \beta', \gamma')'}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param delta the true value of the vector \eqn{\delta = (\delta_2, ..., \delta_{\bar{R}})}{\delta = (\delta_2, ..., \delta_{Rbar})}
#' @param tol the tolerance value used in the Fixed Point Iteration Method to compute the expectancy of `y`. The process stops if the \eqn{L_1}{L} distance 
#' between two consecutive values of the expectancy of `y` is less than `tol`.
#' @param maxit the maximal number of iterations in the Fixed Point Iteration Method.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @details 
#' Following Houndetoungan (2020), the count data \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i \bar{\mathbf{y}} + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*ybar + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, 1)}{\epsilon_i --> N(0, 1)}.\cr
#' Then, \eqn{y_i = r} iff \eqn{a_r \leq y_i^* \leq a_{r+1}}{a_r \le ys_i \le a_{r + 1}}, where
#' \eqn{a_0 = -\inf}{a_0 = -Inf}, \eqn{a_1 = 0}, \eqn{a_r = \sum_{k = 1}^r\delta_k}{a_r = \delta_1 + ... + \delta_r} if \eqn{1 \leq r \leq \bar{R}}{1 \le r \le Rbar}, and 
#' \eqn{a_r = (r - \bar{R})\delta_{\bar{R}} + a_{\bar{R}}}{a_r = (r - Rbar)\delta_{Rbar} + a_{Rbar}} otherwise.
#' @seealso \code{\link{cdnet}}, \code{\link{simsart}}, \code{\link{simsar}}.
#' @return A list consisting of:
#'     \item{yst}{ys (see details), the latent variable.}
#'     \item{y}{the observed count data.}
#'     \item{yb}{ybar (see details), the expectation of y.}
#'     \item{Gyb}{the average of the expectation of y among friends.}
#'     \item{marg.effects}{the marginal effects.}
#'     \item{iteration}{number of iterations performed by sub-network in the Fixed Point Iteration Method.}
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
#' rm(list = ls()[!(ls() %in% c("Glist", "data", "theta", "delta"))])
#' 
#' ytmp    <- simcdnet(formula = ~ x1 + x2 | x1 + x2, Glist = Glist, theta = theta, 
#'                     delta = delta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y, breaks = max(y))}
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm
#' @export
simcdnet   <- function(formula,
                       contextual,
                       Glist,
                       theta,
                       delta,
                       tol   = 1e-15,
                       maxit = 500,
                       data) {
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
  
  Rbar     <- length(delta) + 1
  f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, "sim", 0)
  X        <- f.t.data$X
  
  K        <- length(theta)
  if(K != (ncol(X) + 1)) {
    stop("Length of theta is not suited.")
  }
  lambda   <- theta[1]
  b        <- theta[2:K]
  
  xb       <- c(X %*% b)
  eps      <- rnorm(n, 0, 1)
  
  yb       <- rep(0, n)
  Gyb      <- rep(0, n)

  t        <- fyb(yb, Gyb, Glist, igr, M, xb, lambda, delta, n, Rbar, tol, maxit)
  
  Ztlamda  <- lambda*Gyb + xb
  yst      <- Ztlamda + eps
  y        <- c(fy(yst, max(yst), delta, n, Rbar))
  
  coln      <- c("lambda", colnames(X))
  if("(Intercept)" %in% coln) {
    thetaWI <- theta[-2]
    coln    <- coln[-2]
  }
  meffects  <- fmeffects(Ztlamda, thetaWI, delta)
  names(meffects) <- coln
  
  
  list("yst"          = yst,
       "y"            = y,
       "yb"           = yb,
       "Gyb"          = Gyb,
       "marg.effects" = meffects,
       "iteration"    = c(t))
}

# Marginal effet
fmeffects <- function(Ztlamda, theta, delta) {
  # marginal effect
  maxZTl       <- max(Ztlamda) + 10
  avec         <- c(0, cumsum(delta))
  deltaRB      <- tail(delta, 1)
  cont         <- TRUE
  Rmax         <- length(avec)
  while (cont) {
    Rmax       <- Rmax + 1
    avec[Rmax] <- tail(avec, 1) + deltaRB
    cont       <- tail(avec, 1) < maxZTl
  }
  fir          <- sum(apply(dnorm(kronecker(Ztlamda, t(avec), "-")), 2, mean))
  theta*fir
}