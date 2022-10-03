#' @title Estimate Network Formation Model with Degree Heterogeneity as Random Effects
#' @param network matrix or list of sub-matrix of social interactions containing 0 and 1, where links are represented by 1
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{~ x1 + x2}
#' where `x1`, `x2` are explanatory variable of links formation
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `homophily` is called.
#' @param fixed.effects boolean indicating if sub-network heterogeneity as fixed effects should be included.
#' @param init (optional) list of starting values containing `beta`, an K-dimensional vector of the explanatory variables parameter, 
#' `mu` an n-dimensional vector, and `nu` an n-dimensional vector, `smu2` the variance of `mu`, 
#' and `snu2` the variance of `nu`, 
#' where K is the number of explanatory variables and n is the number of individuals.  
#' @param iteration the number of iterations to be performed. 
#' @param print boolean indicating if the estimation progression should be printed.
#' @return A list consisting of:
#'     \item{n}{number of individuals in each network.}
#'     \item{n.obs}{number of observations.}
#'     \item{n.links}{number of links.}
#'     \item{K}{number of explanatory variables.}
#'     \item{posterior}{list of simulations from the posterior distribution.}
#'     \item{iteration}{number of performed iterations.}
#'     \item{init}{returned list of starting values.}
#' @description 
#' `homophily` implements a Bayesian estimator for network formation model with homophily. The model includes degree heterogeneity as random effects (see details).
#' @details
#' Let \eqn{p_{ij}}{Pij} be a probability for a link to go from the individual \eqn{i} to the individual \eqn{j}.
#' This probability is specified as
#' \deqn{p_{ij} = F(\mathbf{x}_{ij}'\beta + \mu_j + \nu_j)}{Pij = F(Xij'*\beta + \mu_i + \nu_j),}
#' where \eqn{F} is the cumulative of the standard normal distribution. Unobserved degree heterogeneity is captured by
#' \eqn{\mu_i} and \eqn{\nu_j}. The latter are treated as random effects.  
#' @seealso \code{\link{homophily.FE}}.
#' @importFrom ddpcr quiet
#' @importFrom stats lm
#' @importFrom stats var
#' @importFrom stats cov
#' @examples 
#' \donttest{
#' set.seed(1234)
#' library(MASS)
#' M            <- 4 # Number of sub-groups
#' nvec         <- round(runif(M, 100, 500))
#' beta         <- c(.1, -.1)
#' Glist        <- list()
#' dX           <- matrix(0, 0, 2)
#' mu           <- list()
#' nu           <- list()
#' cst          <- runif(M, -1.5, 0)
#' smu2         <- 0.2
#' snu2         <- 0.2
#' rho          <- 0.8
#' Smunu        <- matrix(c(smu2, rho*sqrt(smu2*snu2), rho*sqrt(smu2*snu2), snu2), 2)
#' for (m in 1:M) {
#'   n          <- nvec[m]
#'   tmp        <- mvrnorm(n, c(0, 0), Smunu)
#'   mum        <- tmp[,1] - mean(tmp[,1])
#'   num        <- tmp[,2] - mean(tmp[,2])
#'   X1         <- rnorm(n, 0, 1)
#'   X2         <- rbinom(n, 1, 0.2)
#'   Z1         <- matrix(0, n, n)  
#'   Z2         <- matrix(0, n, n)
#'   
#'   for (i in 1:n) {
#'     for (j in 1:n) {
#'       Z1[i, j] <- abs(X1[i] - X1[j])
#'       Z2[i, j] <- 1*(X2[i] == X2[j])
#'     }
#'   }
#'   
#'   Gm           <- 1*((cst[m] + Z1*beta[1] + Z2*beta[2] +
#'                        kronecker(mum, t(num), "+") + rnorm(n^2)) > 0)
#'   diag(Gm)     <- 0
#'   diag(Z1)     <- NA
#'   diag(Z2)     <- NA
#'   Z1           <- Z1[!is.na(Z1)]
#'   Z2           <- Z2[!is.na(Z2)]
#'   
#'   dX           <- rbind(dX, cbind(Z1, Z2))
#'   Glist[[m]]   <- Gm
#'   mu[[m]]      <- mum
#'   nu[[m]]      <- num
#' }
#' 
#' mu  <- unlist(mu)
#' nu  <- unlist(nu)
#' 
#' out   <- homophily(network =  Glist, formula = ~ dX, fixed.effects = TRUE, 
#'                    iteration = 1e3)
#' 
#' # plot simulations
#' plot(out$posterior$beta[,1], type = "l")
#' abline(h = cst[1], col = "red")
#' plot(out$posterior$beta[,2], type = "l")
#' abline(h = cst[2], col = "red")
#' plot(out$posterior$beta[,3], type = "l")
#' abline(h = cst[3], col = "red")
#' plot(out$posterior$beta[,4], type = "l")
#' abline(h = cst[4], col = "red")
#' 
#' plot(out$posterior$beta[,5], type = "l")
#' abline(h = beta[1], col = "red")
#' plot(out$posterior$beta[,6], type = "l")
#' abline(h = beta[2], col = "red")
#' 
#' plot(out$posterior$sigma2_mu, type = "l")
#' abline(h = smu2, col = "red")
#' plot(out$posterior$sigma2_nu, type = "l")
#' abline(h = snu2, col = "red")
#' plot(out$posterior$rho, type = "l")
#' abline(h = rho, col = "red")
#' 
#' i <- 10
#' plot(out$posterior$mu[,i], type = "l")
#' abline(h = mu[i], col = "red")
#' plot(out$posterior$nu[,i], type = "l")
#' abline(h = nu[i], col = "red")
#' }
#' @export
homophily <- function(network,
                      formula,
                      data,
                      fixed.effects = FALSE,
                      init          = list(),
                      iteration     = 1e3,
                      print         = TRUE) {
  t1              <- Sys.time()
  # Data and dimensions
  if (!is.list(network)) {
    network       <- list(network)
  }
  
  M               <- length(network)
  nvec            <- unlist(lapply(network, nrow))
  n               <- sum(nvec)
  Nvec            <- nvec*(nvec- 1)
  N               <- sum(Nvec)
  
  network         <- unlist(lapply(network, function(x){diag(x) = NA; x}))
  network         <- network[!is.na(network)]
  quiet(gc())
  if (length(network) != N) {
    stop("network should contain only 0 and 1 (as numeric)")
  }
  if (sum(!((network == 0) | (network == 1))) != 0) {
    stop("network should contain only 0 and 1")
  } 
  tmp1            <- cumsum(unlist(lapply(nvec, function(x) rep(x - 1, x)))) - 1
  tmp2            <- c(0, tmp1[-n] + 1)
  index           <- cbind(tmp2, tmp1)
  rm(list = c("tmp1", "tmp2"))
  quiet(gc())
  indexgr         <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  INDEXgr         <- matrix(c(cumsum(c(0, Nvec[-M])), cumsum(Nvec) - 1), ncol = 2)
  # Formula to data
  f.t.data        <- formula.to.data(formula, FALSE, NULL, NULL, NULL, data,
                                     type = "network", theta0 =  NA)
  if(!missing(data)) {
    rm("data")
    quiet(gc())
  }
  
  formula         <- f.t.data$formula
  dX              <- f.t.data$X
  rm("f.t.data")
  quiet(gc())
  coln            <- colnames(dX)
  nfix            <- ifelse("(Intercept)" %in% coln, 1, 0)
  K               <- ncol(dX)
  if (fixed.effects) {
    if(M < 2){
      stop("fixed effects can be added for only one subnetwork")
    }
    K             <- K + M - nfix
    nfix          <- M
    dX            <- dX[,coln != "(Intercept)"]
    coln          <- c(paste0("(Intercept-", 1:M, ")"), coln[coln != "(Intercept)"])
  }
  
  Kx              <- ncol(dX)
  dXdX            <- crossprod(dX)
  sumnetwork      <- NULL
  invdXdX         <- NULL
  if (nfix >= 2){
    sumdX         <- do.call(cbind, lapply(1:M, function(m){colSums(dX[(INDEXgr[m,1] + 1):(INDEXgr[m,2] + 1),])}))
    sumnetwork    <- sapply(1:M, function(m){sum(network[(INDEXgr[m,1] + 1):(INDEXgr[m,2] + 1)])})
    dXdX          <- rbind(cbind(diag(Nvec), t(sumdX)), cbind(sumdX, dXdX))
    invdXdX       <- solve(as.matrix(dXdX))
    rm("sumdX")
    quiet(gc())
  } else{
    invdXdX       <- solve(as.matrix(dXdX))
  }
  rm("dXdX")
  quiet(gc())
  
  #starting value
  beta            <- init$beta
  mu              <- init$mu
  nu              <- init$nu
  smu2            <- init$smu2
  snu2            <- init$snu2
  rho             <- init$rho

  quiet(gc())
  if (is.null(beta)) {
    # print(dim(invdXdX))
    # print(length(sumnetwork))
    # print(dim(dX))
    # print(length(network))
    beta          <- c(invdXdX %*% c(sumnetwork, crossprod(dX, network)))
  } else{
    stopifnot(length(beta) == K)
  }
  
  if (is.null(mu)) {
    mu            <- rep(0, n)
  } else{
    stopifnot(length(mu) == n)
  }
  
  if (is.null(nu)) {
    nu            <- rep(0, n)
  } else{
    stopifnot(length(nu) == n)
  }
  
  if (is.null(smu2)) {
    tmp           <- var(mu)
    smu2          <- ifelse(tmp > 0, tmp, 1)
  } 
  
  if (is.null(snu2)) {
    tmp           <- var(nu)
    snu2          <- ifelse(tmp > 0, tmp, 1)
  }
  
  if (is.null(rho)) {
    rho           <- cov(mu, nu)/sqrt(smu2*snu2)
    rho           <- (rho >= 1) - (rho <= -1) + rho*((rho >= -1) & (rho <= 1))
  } 
  
  init            <- list(beta    = beta,
                          mu      = mu,
                          nu      = nu,
                          smu2    = smu2,
                          snu2    = snu2,
                          rho     = rho)
  estim           <- NULL
  quiet(gc())
  if(print) {
    estim         <- updategparms1(network, dX, invdXdX, beta, mu, nu, smu2, snu2, rho, index, indexgr,
                                   INDEXgr, nfix, N, M, K, Kx, nvec, n, iteration)
  } else {
    estim         <- updategparms2(network, dX, invdXdX, beta, mu, nu, smu2, snu2, rho, index, indexgr,
                                   INDEXgr, nfix, N, M, K, Kx, nvec, n, iteration)
  }
  
  colnames(estim$beta)     <- coln
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  nlinks      <- sum(network)
  out         <- list("n"               = nvec,
                      "n.obs"           = N,
                      "n.links"         = nlinks,
                      "K"               = K,
                      "posterior"       = estim,
                      "init"            = init)
  
  class(out)  <- "homophily"
  if(print) {
    cat("\n\n")
    cat("The program successfully executed \n")
    cat("\n")
    cat("********SUMMARY******** \n")
    cat("n.obs          : ", N, "\n")
    cat("n.links        : ", nlinks, "\n")
    cat("K              : ", K, "\n")
    cat("Fixed effects  : ", ifelse(fixed.effects, "Yes", "No"), "\n")
    cat("Iteration      : ", iteration, "\n\n")
    
    
    # Print the processing time
    nhours     <- floor(timer/3600)
    nminutes   <- floor((timer-3600*nhours)/60)%%60
    nseconds   <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  }
  
  out
}