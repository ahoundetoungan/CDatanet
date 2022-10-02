#' @title Simulate data from the Tobit Model with Social Interactions
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param theta the parameter value as \eqn{\theta = (\lambda, \beta, \gamma, \sigma)}. The parameter \eqn{\gamma} should be removed if the model
#' does not contain contextual effects (see details).
#' @param tol the tolerance value used in the Fixed Point Iteration Method to compute `y`. The process stops if the \eqn{L_1}{L} distance 
#' between two consecutive values of `y` is less than `tol`.
#' @param maxit the maximal number of iterations in the Fixed Point Iteration Method.
#' @param RE a Boolean which indicates if the model if under rational expectation of not.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @description
#' `simsart` is used to simulate censored data with social interactions (see details). The model is presented in Xu and Lee(2015). 
#' @details 
#' The left-censored variable \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}. 
#' The latent variable is given for all i as
#' \deqn{y_i^* = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.\cr
#' The censored variable \eqn{y_i} is then define that is \eqn{y_i = 0} if  
#' \eqn{y_i^* \leq 0}{ys_i \le 0} and \eqn{y_i = y_i^*}{y_i = ys_i} otherwise.
#' @seealso \code{\link{sart}}, \code{\link{simsar}}, \code{\link{simcdnet}}.
#' @return A list consisting of:
#'     \item{yst}{ys (see details), the latent variable.}
#'     \item{y}{the censored variable.}
#'     \item{yb}{expectation of y under rational expectation.}
#'     \item{Gy}{the average of y among friends.}
#'     \item{Gyb}{Average of expectation of y among friends under rational expectation.}
#'     \item{marg.effects}{the marginal effects.}
#'     \item{iteration}{number of iterations performed by sub-network in the Fixed Point Iteration Method.}
#' @references 
#' Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
#' @examples 
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
#' ytmp    <- simsart(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                    theta = theta, data = data)
#' 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y)
#' 
#' @importFrom Rcpp sourceCpp
#' @export
simsart   <- function(formula,
                      contextual,
                      Glist,
                      theta,
                      tol   = 1e-15,
                      maxit = 500,
                      RE = FALSE,
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
  
  
  f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, "sim", 0)
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
  yb       <- NULL
  Gyb      <- NULL
  t        <- NULL
  Ztl      <- rep(0, n)
  if(RE){
    yb     <- rep(0, n)
    Gyb    <- rep(0, n)
    t      <- fybtbit(yb, Gyb, Glist, igr, M, xb, lambda, sigma, n, tol, maxit)
    Ztl    <- lambda*Gyb + xb
    yst    <- Ztl + eps
    y      <- yst*(yst > 0)
  } else {
    y      <- rep(0, n)
    Gy     <- rep(0, n)
    t      <- fyTobit(yst, y, Gy, Ztl, Glist, eps, igr, M, xb, n, lambda, tol, maxit)
  }
  
  # marginal effects
  coln      <- c("lambda", colnames(X))
  thetaWI   <- head(theta, K - 1)
  if("(Intercept)" %in% coln) {
    thetaWI <- thetaWI[-2]
    coln    <- coln[-2]
  }
  meffects  <- thetaWI*mean(pnorm(Ztl/sigma))
  names(meffects) <- coln
  
  
  list("yst"          = yst,
       "y"            = y,
       "yb"           = yb,
       "Gy"           = Gy,
       "Gyb"          = Gyb,
       "marg.effects" = meffects,
       "iteration"    = t)
}