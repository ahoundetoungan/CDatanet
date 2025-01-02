#' @title Estimating Network Formation Models with Degree Heterogeneity: the Fixed Effect Approach
#' @param network A matrix or list of sub-matrices of social interactions containing 0 and 1, where links are represented by 1.
#' @param formula An object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be, for example, \code{~ x1 + x2}, 
#' where `x1` and `x2` are explanatory variables for link formation. If missing, the model is estimated with fixed effects only.
#' @param data An optional data frame, list, or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `homophily` is called.
#' @param symmetry Indicates whether the network model is symmetric (see details).
#' @param fe.way Indicates whether it is a one-way or two-way fixed effect model. The expected value is 1 or 2 (see details).
#' @param init (optional) Either a list of starting values containing `beta`, a K-dimensional vector of the explanatory variables' parameters, 
#' `mu`, an n-dimensional vector, and `nu`, an n-dimensional vector, where K is the number of explanatory variables and n is the number of individuals; 
#' or a vector of starting values for `c(beta, mu, nu)`.  
#' @param method A character string specifying the optimization method. Expected values are `"L-BFGS"`, `"Block-NRaphson"`, or `"Mix"`. 
#' `"Block-NRaphson"` refers to the `Newton-Raphson` method applied to each subnetwork, and `"Mix"` combines the `Newton-Raphson` method for `beta` with the `L-BFGS` method for the fixed effects.
#' @param ctr (optional) A list containing control parameters for the solver. For the `optim_lbfgs` method from the \pkg{RcppNumerical} package, 
#' the list should include `maxit.opt` (corresponding to `maxit` for the `L-BFGS` method), `eps_f`, and `eps_g`. For the `Block-NRaphson` method, 
#' the list should include `maxit.nr` (corresponding to `maxit` for the `Newton-Raphson` method) and `tol`.
#' @param print A boolean indicating if the estimation progression should be printed.
#' @description 
#' `homophily.fe` implements a Logit estimator for a network formation model with homophily. The model includes degree heterogeneity using fixed effects (see details).
#' @details
#' Let \eqn{p_{ij}} be the probability for a link to go from individual \eqn{i} to individual \eqn{j}.
#' This probability is specified for two-way effect models (`fe.way = 2`) as
#' \deqn{p_{ij} = F(\mathbf{x}_{ij}'\beta + \mu_i + \nu_j),}
#' where \eqn{F} is the cumulative distribution function of the standard logistic distribution. Unobserved degree heterogeneity is captured by
#' \eqn{\mu_i} and \eqn{\nu_j}. These are treated as fixed effects (see \code{\link{homophily.re}} for random effect models). 
#' As shown by Yan et al. (2019), the estimator of the parameter \eqn{\beta} is biased. A bias correction is necessary but not implemented in this version. However, 
#' the estimators of \eqn{\mu_i} and \eqn{\nu_j} are consistent.\cr
#' 
#' For one-way fixed effect models (`fe.way = 1`), \eqn{\nu_j = \mu_j}. For symmetric models, the network is not directed, and the fixed effects need to be one-way.
#' @seealso \code{\link{homophily.re}}.
#' @references 
#' Yan, T., Jiang, B., Fienberg, S. E., & Leng, C. (2019). Statistical inference in a directed network model with covariates. \emph{Journal of the American Statistical Association}, 114(526), 857-868, \doi{https://doi.org/10.1080/01621459.2018.1448829}.
#' @return A list consisting of:
#'     \item{model.info}{A list of model information, such as the type of fixed effects, whether the model is symmetric,
#'      the number of observations, etc.}
#'     \item{estimate}{The maximizer of the log-likelihood.}
#'     \item{loglike}{The maximized log-likelihood.}
#'     \item{optim}{The returned value from the optimization solver, which contains details of the optimization. The solver used is `optim_lbfgs` from the 
#'     \pkg{RcppNumerical} package.}
#'     \item{init}{The returned list of starting values.}
#'     \item{loglike.init}{The log-likelihood at the starting values.}
#' @importFrom stats binomial
#' @importFrom matrixcalc is.symmetric.matrix
#' @examples 
#' \donttest{
#' set.seed(1234)
#' M            <- 2 # Number of sub-groups
#' nvec         <- round(runif(M, 20, 50))
#' beta         <- c(.1, -.1)
#' Glist        <- list()
#' dX           <- matrix(0, 0, 2)
#' mu           <- list()
#' nu           <- list()
#' Emunu        <- runif(M, -1.5, 0) # Expectation of mu + nu
#' smu2         <- 0.2
#' snu2         <- 0.2
#' for (m in 1:M) {
#'   n          <- nvec[m]
#'   mum        <- rnorm(n, 0.7*Emunu[m], smu2)
#'   num        <- rnorm(n, 0.3*Emunu[m], snu2)
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
#'   Gm           <- 1*((Z1*beta[1] + Z2*beta[2] +
#'                        kronecker(mum, t(num), "+") + rlogis(n^2)) > 0)
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
#' out   <- homophily.fe(network =  Glist, formula = ~ -1 + dX, fe.way = 2)
#' muhat <- out$estimate$mu
#' nuhat <- out$estimate$nu
#' plot(mu, muhat)
#' plot(nu, nuhat)
#' }
#' @export
homophily.fe <- function(network,
                         formula,
                         data,
                         symmetry   = FALSE,
                         fe.way     = 1,
                         init       = NULL,
                         method     = c("L-BFGS", "Block-NRaphson", "Mix"),
                         ctr        = list(maxit.opt = 1e4, maxit.nr = 50, eps_f = 1e-9, eps_g = 1e-9, tol = 1e-4),
                         print      = TRUE){
  t1      <- Sys.time()
  method  <- method[1]
  meth    <- tolower(method)
  if (!(meth %in% c("l-bfgs", "block-nraphson", "mix"))) stop("`method` must be eiter 'L-BFGS', 'Block-NRaphson', or 'mix'.")
  meth    <- (0:2)[c("l-bfgs", "block-nraphson", "mix") == meth]
  fe.way  <- as.numeric(fe.way[1])
  if(symmetry & fe.way == 2) stop("Two side fixed effects are not allowed for symmetric network models.")
  stopifnot(fe.way %in% (1:2))
  stopifnot(is.null(init) || is.vector(init) || is.list(init))
  
  # Data and dimensions
  if (!is.list(network)) {
    network       <- list(network)
  }
  
  M               <- length(network)
  nvec            <- unlist(lapply(network, nrow))
  n               <- sum(nvec)
  Nvec            <- NULL
  if(symmetry){
    Nvec          <- nvec*(nvec- 1)/2
    stopifnot(sapply(network, is.symmetric.matrix))
    network       <- frMtoVbyCOLsym(network, nvec, M)
  } else {
    Nvec          <- nvec*(nvec- 1)
    # network         <- unlist(lapply(network, function(x){diag(x) = NA; x}))
    # network         <- network[!is.na(network)]
    network       <- frMtoVbyCOL(network, nvec, M)
  }
  N               <- sum(Nvec)
  
  quiet(gc())
  if (sum(!((network == 0) | (network == 1))) != 0) {
    stop("Network should contain only 0 and 1.")
  }
  tmp1    <- NULL
  if(symmetry){
    tmp1  <- cumsum(unlist(lapply(nvec, function(x) (x - 1):0))) - 1
  } else {
    tmp1  <- cumsum(unlist(lapply(nvec, function(x) rep(x - 1, x)))) - 1
  }
  tmp2    <- c(0, tmp1[-n] + 1)
  index   <- cbind(tmp2, tmp1) 
  rm(list = c("tmp1", "tmp2"))
  
  quiet(gc())
  indexgr <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2) #start group, end group
  # INDEXgr         <- matrix(c(cumsum(c(0, Nvec[-M])), cumsum(Nvec) - 1), ncol = 2)
  # Formula to data
  dX              <- matrix(0, 0, 0)
  hasX            <- FALSE
  if(!missing(formula)){
    f.t.data      <- formula.to.data(formula, FALSE, NULL, NULL, NULL, data,
                                     type = "network", theta0 =  NA)
    if(!missing(data)) {
      rm("data")
      quiet(gc())
    }
    formula       <- f.t.data$formula
    dX            <- f.t.data$X
    if(nrow(dX) != N && nrow(dX) != 0) stop("The number of observations in X does not match the network.")
    rm("f.t.data")
    quiet(gc())
    hasX          <- TRUE
  }
  coln            <- colnames(dX)
  if("(Intercept)" %in% coln){stop("Fixed effect model cannot include intercept.")}
  K               <- length(coln)
  nlinks          <- sum(network)
  out             <- list()
  if(symmetry){
    out           <- homophily.LogitFESym(network = network, M = M, nvec = nvec, n = n, N = N, Nvec = Nvec, index = index, 
                                          indexgr = indexgr, formula = formula, dX = dX, coln = coln, K = K, init = init, 
                                          nlinks = nlinks, ctr = ctr, hasX = hasX, print = print, meth = meth)
  } else {
    out           <- homophily.LogitFE(network = network, fe.way = fe.way, M = M, nvec = nvec, n = n, N = N, Nvec = Nvec, 
                                       index = index, indexgr = indexgr, formula = formula, dX = dX, coln = coln, K = K, 
                                       init = init, nlinks = nlinks, ctr = ctr, hasX = hasX, print = print, meth = meth) 
  }
  
  t2              <- Sys.time()
  timer           <- as.numeric(difftime(t2, t1, units = "secs"))
  if(print) {
    cat("\n\n")
    cat("The program successfully executed \n")
    cat("\n")
    cat("********SUMMARY******** \n")
    cat("n.obs          : ", N, "\n")
    cat("n.links        : ", nlinks, "\n")
    cat("K              : ", K, "\n")
    
    
    # Print the processing time
    nhours       <- floor(timer/3600)
    nminutes     <- floor((timer-3600*nhours)/60)%%60
    nseconds     <- timer-3600*nhours-60*nminutes
    cat("Elapsed time   : ", nhours, " HH ", nminutes, " mm ", round(nseconds), " ss \n \n")
  }
  
  out
}



homophily.LogitFE <- function(network, fe.way, M, nvec, n, N, Nvec, index, indexgr, formula, 
                              dX, coln, K, init, nlinks, ctr, hasX, print, meth){
  maxit.opt       <- ctr$maxit.opt
  maxit.nr        <- ctr$maxit.nr
  tol             <- ctr$tol
  eps_f           <- ctr$eps_f
  eps_g           <- ctr$eps_g
  
  if(is.null(maxit.opt)){
    maxit.opt     <- 500
  }
  if(is.null(maxit.nr)){
    maxit.nr      <- 50
  }
  if(is.null(tol)){
    tol           <- 1e-4
  }
  if(is.null(eps_f)){
    eps_f         <- 1e-6
  }
  if(is.null(eps_g)){
    eps_g         <- 1e-5
  }
  
  #starting value
  initllh         <- NULL
  quiet(gc())
  if(is.null(init)){
    if(print) cat("starting point searching\n")
    beta          <- NULL
    mu            <- NULL
    mylogit       <- NULL
    if(hasX){
      mylogit     <- finithasX(dx = dX, theta = rep(0, K), a = network, updateb = TRUE, tol = tol, maxit = maxit.nr) 
    } else {
      mylogit     <- finit(a = network, tol = tol, maxit = maxit.nr) 
    }
    
    beta          <- c(mylogit$theta[-1])
    mu            <- rep(mylogit$theta[1], n)
    names(mu)     <- NULL
    nu            <- NULL
    if(fe.way == 2){
      nu          <- rep(0, n - M)
    } else {
      mu          <- mu/2
    }
    init          <- c(beta, mu, nu)
    initllh       <- mylogit$llh
  } else {
    if(is.list(init)){
      beta        <- c(init$beta)
      mu          <- c(init$mu)
      nu          <- c(init$nu)
      if(is.null(beta) || is.null(mu)){
        if(print) cat("starting point searching\n")
        mylogit    <- NULL
        if(hasX){
          if (is.null(beta)) {
            mylogit  <- finithasX(dx = dX, theta = rep(0, K), a = network, updateb = TRUE, tol = tol, maxit = maxit.nr) 
          } else {
            mylogit  <- finithasX(dx = dX, theta = beta, a = network, updateb = FALSE, tol = tol, maxit = maxit.nr) 
          }
        } else {
          mylogit     <- finit(a = network, tol = tol, maxit = maxit.nr) 
        }
        beta     <- c(mylogit$theta[-1])
        mu       <- rep(mylogit$theta[1], n)
        initllh  <- mylogit$llh
        if(fe.way == 2){
          nu     <- rep(0, n - M)
        }  else {
          mu     <- mu/2
        }
      }
      if(is.null(nu) & (fe.way == 2)){
        nu        <- rep(0, n - M)
      }
      stopifnot(length(beta) == K)
      stopifnot(length(mu) == n)
      if(fe.way == 2){
        stopifnot(length(nu) == (n - M))
      }
      init        <- c(beta, mu, nu)
    } else if(is.vector(init)){
      if(fe.way == 2){
        stopifnot(length(init) == (K + 2*n - M))
      } else {
        stopifnot(length(init) == (K + n))
      }
    } 
  }
  quiet(gc())
  
  theta           <- 1*init
  
  estim           <- NULL
  quiet(gc())
  
  if(print) {
    cat("maximizer searching\n")
  } 
  estim           <- NULL
  if (meth == 0) {
    if(fe.way == 2){
      estim       <- fhomobeta2f(theta = theta, a = c(network), dx = dX, nvec = nvec, index = index, indexgr = indexgr, 
                                 M = M, maxit = maxit.opt, eps_f = eps_f, eps_g = eps_g, hasX = hasX, Print = print)
    } else {
      estim       <- fhomobeta1f(theta = theta, a = c(network), dx = dX, nvec = nvec, index = index, indexgr = indexgr, 
                                 M = M, maxit = maxit.opt, eps_f = eps_f, eps_g = eps_g, hasX = hasX, Print = print)
    }
  } else if (meth == 1) {
    if(fe.way == 2){
      estim       <- NewRaph2f(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                               M = M, N = N, hasX = hasX, Print = print, tol = tol, maxit = maxit.nr)
    } else {
      estim       <- NewRaph1f(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                               M = M, N = N, hasX = hasX, Print = print, tol = tol, maxit = maxit.nr)
    }
  } else {
    if(fe.way == 2){
      estim       <- NewRaphLBFGS2f(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                                    M = M, N = N, hasX = hasX, Print = print, tol = tol, maxitNR = maxit.nr, maxitopt = maxit.opt, 
                                    eps_f = eps_f, eps_g = eps_g)
    } else {
      estim       <- NewRaphLBFGS1f(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                                    M = M, N = N, hasX = hasX, Print = print, tol = tol, maxitNR = maxit.nr, maxitopt = maxit.opt,
                                    eps_f = eps_f, eps_g = eps_g)
    }
  }
  
  # export degree
  theta           <- c(estim$estimate)
  names(theta)    <- names(init)
  beta            <- head(theta, K)
  if(hasX){
    names(beta)   <- coln
  }
  mu              <- theta[(K + 1):(K + n)]
  nu              <- NULL
  if(fe.way == 2){
    nu            <- tail(theta, n - M)
    nu            <- unlist(lapply(1:M, function(x) c(nu[(indexgr[x, 1] + 2 - x):(indexgr[x, 2] + 1 - x)], 0))) 
  }
  
  estim$estimate  <- c(estim$estimate)
  estim$gradient  <- c(estim$gradient)
  out             <- list("model.info"     = list("model"       = "logit", 
                                                  "sym.network" = FALSE,
                                                  "fe.way"      = fe.way, 
                                                  "n"           = nvec,
                                                  "n.obs"       = N,
                                                  "n.links"     = nlinks,
                                                  "K"           = K),
                          "estimate"        = list(beta = beta, mu = mu, nu = nu),
                          "loglike"         = -estim$value,
                          "optim"           = estim,
                          "init"            = init,
                          "loglike.init"    = initllh, 
                          "distance"        = estim$distance,
                          "iteration"       = estim$iteration)
  
  class(out)      <- "homophily.fe"
  out
}


homophily.LogitFESym <- function(network, M, nvec, n, N, Nvec, index, indexgr, formula, 
                                 dX, coln, K, init, nlinks, ctr, hasX, print, meth){
  maxit.opt       <- ctr$maxit.opt
  maxit.nr        <- ctr$maxit.nr
  tol             <- ctr$tol
  eps_f           <- ctr$eps_f
  eps_g           <- ctr$eps_g
  
  if(is.null(maxit.opt)){
    maxit.opt     <- 500
  }
  if(is.null(maxit.nr)){
    maxit.nr      <- 50
  }
  if(is.null(tol)){
    tol           <- 1e-4
  }
  if(is.null(eps_f)){
    eps_f         <- 1e-6
  }
  if(is.null(eps_g)){
    eps_g         <- 1e-5
  }
  
  #starting value
  initllh         <- NULL
  quiet(gc())
  if(is.null(init)){
    if(print) cat("starting point searching\n")
    beta          <- NULL
    mu            <- NULL
    mylogit       <- NULL
    if(hasX){
      mylogit     <- finithasX(dx = dX, theta = rep(0, K), a = network, updateb = TRUE, tol = tol, maxit = maxit.nr) 
    } else {
      mylogit     <- finit(a = network, tol = tol, maxit = maxit.nr) 
    }
    
    beta          <- c(mylogit$theta[-1])
    mu            <- rep(mylogit$theta[1], n)/2
    initllh       <- mylogit$llh
    names(mu)     <- NULL
    init          <- c(beta, mu)
  } else {
    if(is.list(init)){
      beta        <- c(init$beta)
      mu          <- c(init$mu)
      if(is.null(beta) || is.null(mu)){
        if(print) cat("starting point searching\n")
        mylogit       <- NULL
        if(hasX){
          if (is.null(beta)) {
            mylogit  <- finithasX(dx = dX, theta = rep(0, K), a = network, updateb = TRUE, tol = tol, maxit = maxit.nr) 
          } else {
            mylogit  <- finithasX(dx = dX, theta = beta, a = network, updateb = FALSE, tol = tol, maxit = maxit.nr) 
          }
        } else {
          mylogit     <- finit(a = network, tol = tol, maxit = maxit.nr) 
        }
        beta     <- c(mylogit$theta[-1])
        mu       <- rep(mylogit$theta[1], n)/2
        initllh  <- mylogit$llh
      }
      stopifnot(length(beta) == K)
      stopifnot(length(mu) == n)
      init        <- c(beta, mu)
    } else if(is.vector(init)){
      stopifnot(length(init) == (K + n))
    } 
  }
  
  quiet(gc())
  
  theta           <- 1*init
  
  estim           <- NULL
  quiet(gc())
  
  if(print) {
    cat("maximizer searching\n")
  } 
  
  estim           <- NULL
  if (meth == 0) {
    estim         <- fhomobetasym(theta = theta, a = c(network), dx = dX, nvec = nvec, index = index, indexgr = indexgr, 
                                  M = M, maxit = maxit.nr, eps_f = eps_f, eps_g = eps_g, hasX = hasX, Print = print)
  } else if (meth == 1) {
    estim         <- NewRaphsym(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                                M = M, N = N, hasX = hasX, Print = print, tol = tol, maxit = maxit.opt)
  } else {
    estim         <- NewRaphLBFGSsym(theta = theta, a = c(network), dx = dX, nvec = nvec, Nvec = Nvec, index = index, indexgr = indexgr, 
                                    M = M, N = N, hasX = hasX, Print = print, tol = tol, maxitNR = maxit.nr, maxitopt = maxit.opt, 
                                    eps_f = eps_f, eps_g = eps_g)
  }
  
  
  # export degree
  theta           <- c(estim$estimate)
  names(theta)    <- names(init)
  beta            <- head(theta, K)
  if(hasX){
    names(beta)   <- coln
  }
  mu              <- tail(theta, n)
  
  estim$estimate  <- c(estim$estimate)
  estim$gradient  <- c(estim$gradient)
  
  out             <- list("model.info"     = list("model"           = "logit", 
                                                  "sym.network"     = TRUE,
                                                  "n"               = nvec,
                                                  "n.obs"           = N,
                                                  "n.links"         = nlinks,
                                                  "K"               = K),
                          "estimate"        = list(beta = beta, mu = mu),
                          "loglike"         = -estim$value,
                          "optim"           = estim,
                          "init"            = init,
                          "loglike.init"    = initllh, 
                          "distance"        = estim$distance,
                          "iteration"       = estim$iteration)
  
  class(out)      <- "homophily.fe"
  out
}
