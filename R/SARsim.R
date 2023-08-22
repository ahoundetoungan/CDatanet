#' @title Simulate data from the linear-in-mean Model with Social Interactions
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
#' @param RE a Boolean which indicates if the model if under rational expectation of not.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @description
#' `simsar` is used to simulate continuous variables with social interactions (see details). The model is presented in Lee(2004). 
#' @references  
#' Lee, L. F. (2004). Asymptotic distributions of quasi-maximum likelihood estimators for spatial autoregressive models. \emph{Econometrica}, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}.
#' @details 
#' The variable \eqn{\mathbf{y}}{y} is given for all i as
#' \deqn{y_i = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{y_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
#' where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.
#' @seealso \code{\link{sar}}, \code{\link{simsart}}, \code{\link{simcdnet}}.
#' @return A list consisting of:
#'     \item{y}{the observed count data.}
#'     \item{Gy}{the average of y among friends.}
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
#' ytmp    <- simsar(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
#'                      theta = theta, data = data) 
#' y       <- ytmp$y
#' 
#' # plot histogram
#' hist(y)
#' 
#' @importFrom Rcpp sourceCpp
#' @export
simsar   <- function(formula,
                     contextual,
                     Glist,
                     theta,
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
  # print(colnames(X))
  K        <- length(theta)
  if(K != (ncol(X) + 2)) {
    stop("Length of theta is not suited.")
  }
  lambda   <- theta[1]
  b        <- theta[2:(K - 1)]
  sigma    <- theta[K ]
  
  
  
  xb       <- c(X %*% b)
  eps      <- rnorm(n, 0, sigma)
  
  out      <- NULL
  if(RE){
    y      <- rep(0, n)
    ye     <- rep(0, n)
    Gye    <- numeric(n)
    fySarRE(y, Gye, ye, Glist, eps, igr, M, xb, lambda)
    out    <- list("y"  = y,
                   "Gy" = Gy)
  } else {
    y      <- rep(0, n)
    Gy     <- numeric(n)
    fySar(y, Gy, Glist, eps, igr, M, xb, lambda)
    out    <- list("y"  = y,
                   "Gy" = Gy)
  }
  out
}