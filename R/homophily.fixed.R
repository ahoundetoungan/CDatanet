#' @title Estimate Network Formation Model with Degree Heterogeneity as Fixed Effects
#' @param network matrix or list of sub-matrix of social interactions containing 0 and 1, where links are represented by 1
#' @param formula an object of class \link[stats]{formula}: a symbolic description of the model. The `formula` should be as for example \code{~ x1 + x2}
#' where `x1`, `x2` are explanatory variable of links formation
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `homophily` is called.
#' @param init (optional) either a list of starting values containing `beta`, an K-dimensional vector of the explanatory variables parameter, 
#' `mu` an n-dimensional vector, and `nu` an n-dimensional vector, 
#' where K is the number of explanatory variables and n is the number of individuals; or a vector of starting value for `c(beta, mu, nu)`.  
#' @param opt.ctr (optimal) is a list of `maxit`, `eps_f`, and `eps_g`, which are control parameters used by the solver `optim_lbfgs`, of the package \pkg{RcppNumerical}.
#' @param print Boolean indicating if the estimation progression should be printed.
#' @description 
#' `homophily.FE` is used to estimate a network formation model with homophily. The model includes degree heterogeneity as fixed effects (see details).
#' @details
#' Let \eqn{p_{ij}}{Pij} be a probability for a link to go from the individual \eqn{i} to the individual \eqn{j}.
#' This probability is specified as
#' \deqn{p_{ij} = F(\mathbf{x}_{ij}'\beta + \mu_j + \nu_j)}{Pij = F(Xij'*\beta + \mu_i + \nu_j),}
#' where \eqn{F} is the cumulative of the standard normal distribution. Unobserved degree heterogeneity is captured by
#' \eqn{\mu_i} and \eqn{\nu_j}. The latter are treated as fixed effects. As shown by Yan et al. (2019), the estimator of 
#' the parameter \eqn{\beta} is biased. A bias correction is then necessary and is not implemented in this version. However
#' the estimator of \eqn{\mu_i} and \eqn{\nu_j} are consistent.
#' @seealso \code{\link{homophily}}.
#' @references 
#' Yan, T., Jiang, B., Fienberg, S. E., & Leng, C. (2019). Statistical inference in a directed network model with covariates. \emph{Journal of the American Statistical Association}, 114(526), 857-868, \doi{https://doi.org/10.1080/01621459.2018.1448829}.
#' @return A list consisting of:
#'     \item{n}{number of individuals in each network.}
#'     \item{n.obs}{number of observations.}
#'     \item{n.links}{number of links.}
#'     \item{K}{number of explanatory variables.}
#'     \item{estimate}{maximizer of the log-likelihood.}
#'     \item{loglike}{maximized log-likelihood.}
#'     \item{optim}{returned value of the optimization solver, which contains details of the optimization. The solver used is `optim_lbfgs` of the 
#'     package \pkg{RcppNumerical}.}
#'     \item{init}{returned list of starting value.}
#'     \item{loglike(init)}{log-likelihood at the starting value.}
#' @importFrom stats glm
#' @importFrom stats binomial
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
#' Emunu        <- runif(M, -1.5, 0) #expectation of mu + nu
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
#' out   <- homophily.FE(network =  Glist, formula = ~ -1 + dX)
#' muhat <- out$estimate$mu
#' nuhat <- out$estimate$nu
#' }
#' @export
homophily.FE <- function(network,
                         formula,
                         data,
                         init      = NULL,
                         opt.ctr   = list(maxit = 300, eps_f = 1e-6, eps_g = 1e-5),
                         print     = TRUE) {
  
  t1              <- Sys.time()
  stopifnot(is.null(init) || is.vector(init) || is.list(init))
  maxit           <- opt.ctr$maxit;
  eps_f           <- opt.ctr$eps_f
  eps_g           <- opt.ctr$eps_g
  if(is.null(maxit)){
    maxit         <- 500
  }
  if(is.null(eps_f)){
    eps_f         <- 1e-6
  }
  if(is.null(eps_g)){
    eps_g         <- 1e-5
  }
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
  # INDEXgr         <- matrix(c(cumsum(c(0, Nvec[-M])), cumsum(Nvec) - 1), ncol = 2)
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
  if("(Intercept)" %in% coln){stop("fixed effect model cannot include intercept")}
  K               <- length(coln)
  
  # dXdX            <- crossprod(dX)
  # invdXdX         <- solve(as.matrix(dXdX))
  # sumnetwork      <- sapply(1:M, function(m){sum(network[(INDEXgr[m,1] + 1):(INDEXgr[m,2] + 1)])})
  
  #starting value
  initllh         <- NULL
  quiet(gc())
  if(is.null(init)){
    if(print) cat("starting point searching\n")
    mylogit       <- glm(network ~ 1 + dX, family = binomial(link = "logit"))
    beta          <- mylogit$coefficients[-1]
    mu            <- rep(mylogit$coefficients[1], n)
    names(mu)     <- NULL
    nu            <- rep(0, n - M)
    init          <- c(beta, mu, nu)
    initllh       <- -0.5*mylogit$deviance
  } else {
    if(is.list(init)){
      beta        <- c(init$beta)
      mu          <- c(init$mu)
      nu          <- c(init$nu)
      if(is.null(beta) || is.null(mu)){
        if(print) cat("starting point searching\n")
        mylogit   <- glm(network ~ 1 + dX, family = binomial(link = "logit"))
        initllh   <- -0.5*mylogit$deviance
        if(is.null(mu)){
          mu      <- rep(mylogit$coefficients[1], n); names(mu) <- NULL
        }
        if(is.null(beta)){
          beta    <- mylogit$coefficients[-1]
        }
      }
      if(is.null(nu)){
        nu        <- rep(0, n - M)
      }
      stopifnot(length(beta) == K)
      stopifnot(length(mu) == n)
      stopifnot(length(nu) == (n - M))
      init        <- c(beta, mu, nu)
    } else if(is.vector(init)){
      stopifnot(length(init) == (K + 2*n - M))
    } 
  }
  
  quiet(gc())
  
  theta           <- init
  
  estim           <- NULL
  quiet(gc())
  
  if(print) {
    if(print){
      cat("maximizer searching\n")
    }
    estim         <- fhomobetap(theta, c(network), dX, nvec, index, indexgr, M, maxit, eps_f, eps_g)
  } else {
    estim         <- fhomobeta(theta, c(network), dX, nvec, index, indexgr, M, maxit, eps_f, eps_g)
  }
  
  # export degree
  theta           <- c(estim$estimate)
  names(theta)    <- names(init)
  beta            <- head(theta, K)
  mu              <- theta[(K + 1):(K + n)]
  nu              <- tail(theta, n - M)
  names(beta)     <- coln
  nu              <- unlist(lapply(1:M, function(x) c(nu[(indexgr[x, 1] + 2 - x):(indexgr[x, 2] + 1 - x)], 0))) 
  
  estim$estimate  <- c(estim$estimate)
  estim$gradient  <- c(estim$gradient)
  
  t2              <- Sys.time()
  timer           <- as.numeric(difftime(t2, t1, units = "secs"))
  
  nlinks          <- sum(network)
  out             <- list("n"               = nvec,
                          "n.obs"           = N,
                          "n.links"         = nlinks,
                          "K"               = K,
                          "estimate"        = list(beta = beta, mu = mu, nu = nu),
                          "loglike"         = -estim$value,
                          "optim"           = estim,
                          "init"            = init,
                          "loglike(init)"   = initllh)
  
  class(out)      <- "homophily.FE"
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