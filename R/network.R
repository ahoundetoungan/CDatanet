#' @title Estimate Network Formation Model
#' @param network matrix or list of sub-matrix of social interactions containing 0 and 1, where links are represented by 1
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{~ x1 + x2}
#' where `x1`, `x2` are explanatory variable of links formation
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `netformation` is called.
#' @param fixed.effects boolean indicating if sub-network heterogeneity as fixed effects should be included.
#' @param init (optional) list of starting values containing `beta`, an K-dimensional vector of the explanatory variables parameter, `mu` an n-dimensional vector of unobserved parameters, `sigmau2`
#' the vector of the variances of `mu` in each sub-network (or single variance if `fixed.effects = FALSE`) and `uu` the vector of the of the means of `mu` in each sub-network (or single mean if `fixed.effects = FALSE`), 
#' where K is the number of explanatory variables and n is the number of individuals.  
#' @param mcmc.ctr (optional) list of MCMC control (see detail). 
#' @param print boolean indicating if the estimation progression should be printed.
#' @details 
#' The network formation model can be used to control the network endogeneity in the count data, Tobit, and SAR models (see the vignette for an example with the count data model). \cr
#' The MCMC control `mcmc.ctr` should be a list containing
#' \itemize{
#'    \item{jscalebeta}{ a K-dimensional vector of `beta` jumping scales.}
#'    \item{jscalemu}{ an n-dimensional vector of `mu` jumping scales.}
#'    \item{burnin}{ the number of iterations in the burn-in.}
#'    \item{iteration}{ the number of simulations.}
#'    \item{tbeta}{ target of `beta`.}
#'    \item{tmu}{ target of `mu`.}
#' }
#' The burn-in is replicated three times. The estimation is performed for each component of `beta` during the two firsts burn-in. The simulation from
#' the second burn-in are used to compute covariance of `beta`. From the third burn-in and the remaining steps of the MCMC, all the components in `beta` are jointly simulated. \cr
#' As `mu` dimension is large, the simulation is performed for each component.\cr
#' The jumping scale are also updated during the MCMC following Atchadé and Rosenthal (2005).
#' 
#' @return A list consisting of:
#'     \item{n}{number of individuals in each network.}
#'     \item{n.obs}{number of observations.}
#'     \item{n.links}{number of links.}
#'     \item{K}{number of explanatory variables.}
#'     \item{posterior}{list of simulations from the posterior distribution and the posterior density.}
#'     \item{acceptance.rate}{acceptance rate of beta and mu.}
#'     \item{mcmc.ctr}{returned list of MCMC control.}
#'     \item{init}{returned list of starting values.}
#' @importFrom ddpcr quiet
#' @examples 
#' \donttest{
#' M            <- 5 # Number of sub-groups
#' nvec         <- round(runif(M, 100, 500))
#' beta         <- c(1, -1)
#' Glist        <- list()
#' dX           <- matrix(0, 0, 2)
#' mu           <- list()
#' uu           <- runif(M, -5, 5)
#' sigma2u      <- runif(M, 0.5, 16)
#' for (m in 1:M) {
#'   n          <- nvec[m]
#'   mum        <- rnorm(n, uu[m], sqrt(sigma2u[m]))
#'   X1         <- rnorm(n)
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
#'   Gm           <- 1*((Z1*beta[1] + Z2*beta[2] +
#'                         kronecker(mum, t(mum), "+") + rlogis(n^2)) > 0)
#'   diag(Gm)     <- 0
#'   
#'   diag(Z1)     <- NA
#'   diag(Z2)     <- NA
#'   Z1           <- Z1[!is.na(Z1)]
#'   Z2           <- Z2[!is.na(Z2)]
#'   
#'   dX           <- rbind(dX, cbind(Z1, Z2))
#'   Glist[[m]]   <- Gm
#'   mu[[m]]      <- mum
#' }
#' 
#' mu  <- unlist(mu)
#' out <- netformation(network =  Glist, formula = ~ dX, fixed.effects = T,
#'                     mcmc.ctr = list(burin = 1000, iteration = 5000))
#' 
#' 
#' # plot simulations
#' plot(out$posterior$beta[,1], type = "l")
#' abline(h = beta[1], col = "red")
#' plot(out$posterior$beta[,2], type = "l")
#' abline(h = beta[2], col = "red")
#' 
#' k <- 2
#' plot(out$posterior$sigmamu2[,2], type = "l")
#' abline(h = sigma2u[2], col = "red")
#' 
#' 
#' i <- 10
#' plot(out$posterior$mu[,i], type = "l", 
#'      ylim = c(min(out$posterior$mu[,i], mu[i]), max(out$posterior$mu[,i], mu[i])))
#' abline(h = mu[i], col = "red")
#' 
#' plot(out$posterior$uu[,k] , type = "l", 
#'      ylim = c(min(out$posterior$uu[,k], uu[k]), max(out$posterior$uu[,k], uu[k])))
#' abline(h = uu[k], col = "red")
#' }
#' @export

netformation <- function(network,
                         formula,
                         data,
                         fixed.effects = TRUE,
                         init          = list(),
                         mcmc.ctr      = list(),
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
  
  if (length(network) != N) {
    stop("network should contain only 0 and 1 (as numeric)")
  }
  if (sum(!((network == 0) | (network == 1))) != 0) {
    stop("network should contain only 0 and 1")
  } 
  tmp1            <- cumsum(unlist(lapply(nvec, function(x) rep(x - 1, x)))) - 1
  tmp2            <- c(0, tmp1[-n] + 1)
  index           <- cbind(tmp2, tmp1)
  
  indexgr         <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  Msigma          <- 1
  indexsigma      <- matrix(c(0, n - 1), nrow = 1)
  possigma        <- rep(0, n)
  if(fixed.effects) {
    Msigma        <- M
    indexsigma    <- indexgr
    possigma      <- unlist(lapply(1:M, function(x) rep(0, nvec[x])))
  }
  
  # Formula to data
  f.t.data        <- formula.to.data(formula, FALSE, NULL, NULL, NULL, data,
                                     type = "network", theta0 =  NA)
  if(!missing(data)) {
    rm("data")
  }
  
  formula         <- f.t.data$formula
  dX              <- f.t.data$X
  coln            <- colnames(dX)
  rm("f.t.data")
  quiet(gc())
  
  if("(Intercept)" %in% coln) {
    dX          <- dX[,-which(coln == "(Intercept)")]
    coln        <- coln[-which(coln == "(Intercept)")]
  }
  
  K               <- ncol(dX)
  
  # MCMC control
  jsbeta          <- mcmc.ctr$jscalebeta
  jsmu            <- mcmc.ctr$jscalemu
  burnin          <- mcmc.ctr$burnin
  iteration       <- mcmc.ctr$iteration
  tbeta           <- mcmc.ctr$tbeta
  tmu             <- mcmc.ctr$tmu
  
  if (is.null(jsbeta)) {
    jsbeta        <- rep(1, K)
  }
  if (is.null(jsmu)) {
    jsmu          <- rep(1, n)
  }
  if (is.null(burnin)) {
    burnin        <- 500
  }
  if (is.null(iteration)) {
    iteration     <- 1000
  }
  if (is.null(tbeta)) {
    tbeta         <- 0.27
  }
  if (is.null(tmu)) {
    tmu           <- 0.27
  }
  
  mcmc.ctr        <- list(jscalebeta = jsbeta,
                          jscalemu   = jsmu,
                          burnin     = burnin,
                          iteration  = iteration,
                          tbeta      = tbeta,
                          tmu        = tmu)
  
  iteration1      <- 2*burnin
  iteration2      <- burnin + iteration
  iteration       <- iteration1 + iteration2

  #starting value
  beta            <- init$beta
  mu              <- init$mu
  sigmau2         <- init$sigmau2
  uu              <- init$uu
  
  if (is.null(beta)) {
    beta          <- rep(0, K)
  }
  if (is.null(mu)) {
    mu            <- rep(0, n)
  }
  if (is.null(sigmau2)) {
    sigmau2       <- rep(1, Msigma)
  }
  if (is.null(uu)) {
    uu            <- rep(0, Msigma)
  }
  
  init            <- list(beta    = beta,
                          mu      = mu,
                          uu      = uu,
                          sigmau2 = sigmau2)
  estim           <- NULL
  if(print) {
    estim         <- updategparms1(network, dX, beta, mu, sigmau2, uu, jsbeta, 
                                   jsmu, index, indexgr,  indexsigma,
                                   possigma, N, M, K,
                                   Msigma, nvec, iteration1, iteration2,
                                   tbeta, tmu)
  } else {
    estim         <- updategparms2(network, dX, beta, mu, sigmau2, uu, jsbeta, 
                                   jsmu, index, indexgr,  indexsigma,
                                   possigma, N, M, K,
                                   Msigma, nvec, iteration1, iteration2,
                                   tbeta, tmu)
  }
  
  posterior                    <- estim$posterior
  accept                       <- estim$acceptance.rate
  accept$mu                    <- c(accept$mu)
  
  colnames(posterior$beta)     <- coln
  if (Msigma == 1) {
    posterior$sigmamu2         <- c(posterior$sigmamu2)
    posterior$uu               <- c(posterior$uu)
  }
  
  t2          <- Sys.time()
  timer       <- as.numeric(difftime(t2, t1, units = "secs"))
  
  nlinks      <- sum(network)
  out         <- list("n"               = nvec,
                      "n.obs"           = N,
                      "n.links"         = nlinks,
                      "K"               = K,
                      "posterior"       = posterior,
                      "acceptance.rate" = accept,
                      "mcmc.ctr"        = mcmc.ctr,
                      "init"            = init)
  
  class(out)  <- "netformation"
  if(print) {
    cat("\n\n")
    cat("The program successfully executed \n")
    cat("\n")
    cat("********SUMMARY******** \n")
    cat("n.obs          : ", N, "\n")
    cat("n.links        : ", nlinks, "\n")
    cat("K              : ", K, "\n")
    cat("Fixed effects  : ", ifelse(fixed.effects, "Yes", "No"), "\n")
    cat("Burnin         : ", burnin, "\n")
    cat("Iteration      : ", iteration2, "\n\n")
    
    
    # Print the processing time
    nhours     <- floor(timer/3600)
    nminutes   <- floor((timer-3600*nhours)/60)%%60
    nseconds   <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
    cat("Average acceptance rate \n")
    cat("                   beta: ", accept$beta, "\n")
    cat("                     mu: ", mean(c(accept$mu)), "\n")
  }
  
  out
}