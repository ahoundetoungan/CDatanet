#' @title Simulate from Count Data Model With Social Interactions
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{y ~ x1 + x2 | x1 + x2}
#' where `y` is the endogenous vector, the listed variables before the pipe, `x1`, `x2` are the individual exogenous variables and
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' \code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
#' without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.
#' @param  contextual (optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
#' the `formula` as `y ~ x1 + x2` and `contextual` as `TRUE` is equivalent to set the formula as `y ~ x1 + x2 | x1 + x2`.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `mcmcARD` is called.
#' @importFrom Rcpp sourceCpp
#' @export
simCDnet   <- function(formula,
                       contextual,
                       Glist,
                       theta,
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
  
  
  f.t.data <- formula.to.data(formula, contextual, Glist, M, igr, data, "sim", 0)
  X        <- f.t.data$X
  
  K        <- length(theta)
  if(K != (ncol(X) + 2)) {
    stop("Length of theta0 is not suited.")
  }
  lambda   <- theta[1]
  b        <- theta[2:(K - 1)]
  sigma    <- theta[K ]
  

  
  xb       <- c(X %*% b)
  eps      <- rnorm(n, 0, sigma)
  
  yb       <- rep(0, n)
  Gyb      <- numeric(n)
  t        <- fyb(yb, Gyb, Glist, igr, M, xb, lambda, sigma, n, tol, maxit)
  
  
  yst      <- lambda*Gyb + xb + eps
  y        <- ceiling(yst)
  y[y < 0] <- 0
  
  list("yst"       = yst,
       "y"         = y,
       "yb"        = yb,
       "Gyb"       = Gyb,
       "iteration" = c(t))
}