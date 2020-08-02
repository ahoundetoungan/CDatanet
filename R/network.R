#' @title Estimate Network Formation Model
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{~ x1 + x2}
#' where `x1`, `x2` are explanatory variable of links formation
#' the listed variables after the pipe, `x1`, `x2` are the contextual observable variables. Other formulas may be
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `netformation` is called.
#' @return A list consisting of:
#' @importFrom ddpcr quiet
#' @export

netformation <- function(network,
                         formula,
                         data,
                         init          = list(),
                         prior         = list(),
                         mcmc.ctr      = list(),
                         fixed.effects = FALSE,
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
  possigma        <- rep(0, N)
  if(fixed.effects) {
    Msigma        <- M
    indexsigma    <- indexgr
    possigma      <- unlist(lapply(1:M, function(x) rep(0, Nvec[x])))
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
  # prior
  a               <- prior$a
  b               <- prior$b
  invSb           <- prior$invSb
  mub             <- prior$mub
  
  if (is.null(invSb)) {
    invSb         <- diag(K)
  }
  if (is.null(mub)) {
    mub           <- rep(0, K)
  }
  
  prior           <- list(invSb = invSb,
                          mub   = mub)
  
  invSbmub        <- invSb %*% mub
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
                                   possigma, invSb, invSbmub, N, M, K,
                                   Msigma, nvec, iteration1, iteration2,
                                   tbeta, tmu)
  } else {
    estim         <- updategparms2(network, dX, beta, mu, sigmau2, uu, jsbeta, 
                                   jsmu, index, indexgr,  indexsigma,
                                   possigma, invSb, invSbmub, N, M, K,
                                   Msigma, nvec, iteration1, iteration2,
                                   tbeta, tmu)
  }
  
  posterior                    <- estim$posterior
  accept                       <- estim$acceptance.rate
  
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
                      "prior"           = prior,
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