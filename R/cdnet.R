#' @title Simulating count data models with social interactions under rational expectations
#' @param formula a class object \link[stats]{formula}: a symbolic description of the model. `formula` must be as, for example, \code{y ~ x1 + x2 + gx1 + gx2}
#' where `y` is the endogenous vector and `x1`, `x2`, `gx1` and `gx2` are control variables, which can include contextual variables, i.e. averages among the peers.
#' Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogenous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param parms a vector defining the true value of \eqn{\theta = (\lambda', \Gamma', \delta')'} (see the model specification in details). 
#' Each parameter \eqn{\lambda}, \eqn{\Gamma}, or \eqn{\delta} can also be given separately to the arguments `lambda`, `Gamma`, or `delta`.
#' @param lambda the true value of the vector \eqn{\lambda}.
#' @param Gamma the true value of the vector \eqn{\Gamma}.
#' @param delta the true value of the vector \eqn{\delta}.
#' @param cost.group the vector indicating the groups in which individuals have the same cost function (same cut-points in the ordered model). The default assumes a common group.
#' @param Rmax an integer indicating the theoretical upper bound of `y`. (see the model specification in details).
#' @param Rbar an \eqn{L}-vector, where  \eqn{L} is the number of cost groups. For large `Rmax` the cost function is assumed to be semi-parametric (i.e., nonparametric from 0 to \eqn{\bar{R}} and quadratic beyond \eqn{\bar{R}}). 
#' The `l`-th element of `Rbar` indicates \eqn{\bar{R}} for the `l`-th value of `sort(unique(cost.group))` (see the model specification in details).
#' @param tol the tolerance value used in the Fixed Point Iteration Method to compute the expectancy of `y`. The process stops if the \eqn{\ell_1}-distance 
#' between two consecutive \eqn{E(y)} is less than `tol`.
#' @param maxit the maximal number of iterations in the Fixed Point Iteration Method.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `simcdnet` is called.
#' @description
#' `simcdnet` simulate the count data model with social interactions under rational expectations developed by Houndetoungan (2024).
#' @details 
#' Houndetoungan (2024) shows that the count variable \eqn{y_i}{yi} take the value \eqn{r} with probability. 
#' \deqn{P_{ir} = F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r}) - F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r + 1}).}
#' In this equation, \eqn{\mathbf{z}_i} is a vector of control variables; \eqn{F} is the distribution function of the standard normal distribution;
#' \eqn{\bar{y}_i^{e,s}} is the average of \eqn{E(y)} among peers using the `s`-th network definition;
#' \eqn{a_{h(i),r}} is the `r`-th cut-point in the cost group \eqn{h(i)}. \cr\cr
#' For identification, \eqn{a_{h(i),0} = -\infty} and \eqn{a_{h(i),1} = 0}. Moreover, \eqn{a_{h(i),r} = \infty} for any \eqn{r \geq R_{\text{max}} + 1}, 
#' which implies that \eqn{P_{ir} = 0} for any \eqn{r \geq R_{\text{max}} + 1}.
#' For any \eqn{r \geq 1}, the distance between two cut-points is \eqn{a_{h(i),r+1} - a_{h(i),r} =  \delta_{h(i),r} + \sum_{s = 1}^S \lambda_s}
#' As the number of cut-point can be large, a quadratic cost function is considered for \eqn{r \geq \bar{R}_{h(i)}}, where \eqn{\bar{R} = (\bar{R}_{1}, ..., \bar{R}_{L})}.
#' With the semi-parametric cost-function,
#' \eqn{a_{h(i),r + 1} - a_{h(i),r}= \bar{\delta}_{h(i)} + \sum_{s = 1}^S \lambda_s}.  \cr\cr
#' The model parameters are: \eqn{\lambda = (\lambda_1, ..., \lambda_S)'}, \eqn{\Gamma}, and \eqn{\delta = (\delta_1', ..., \delta_L')'}, 
#' where \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l}, \bar{\delta}_l)'} for \eqn{l = 1, ..., L}. 
#' The number of single parameters in \eqn{\delta_l} depends on  \eqn{R_{\text{max}}} and \eqn{\bar{R}_{l}}. The components \eqn{\delta_{l,2}, ..., \delta_{l,\bar{R}_l}} or/and 
#' \eqn{\bar{\delta}_l} must be removed in certain cases.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_{l} \geq 2}, then \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l})'}.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_{l} = 1} (binary models), then \eqn{\delta_l} must be empty.\cr
#' If \eqn{R_{\text{max}} > \bar{R}_{l} = 1}, then \eqn{\delta_l = \bar{\delta}_l}.
#' @seealso \code{\link{cdnet}}, \code{\link{simsart}}, \code{\link{simsar}}.
#' @return A list consisting of:
#'     \item{yst}{\eqn{y^{\ast}}, the latent variable.}
#'     \item{y}{the observed count variable.}
#'     \item{Ey}{\eqn{E(y)}, the expectation of y.}
#'     \item{GEy}{the average of \eqn{E(y)} friends.}
#'     \item{marg.effects}{the marginal effects for each individual.}
#'     \item{Rmax}{infinite sums in the marginal effects are approximated by sums up to Rmax.}
#'     \item{iteration}{number of iterations performed by sub-network in the Fixed Point Iteration Method.}
#' @references 
#' Houndetoungan, E. A. (2024). Count Data Models with Social Interactions under Rational Expectations. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @examples 
#' \donttest{
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 50, 100))
#' n      <- sum(nvec)
#' 
#' # Adjacency matrix
#' A      <- list() 
#' for (m in 1:M) {
#'   nm           <- nvec[m]
#'   Am           <- matrix(0, nm, nm)
#'   max_d        <- 30 #maximum number of friends
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
#'     Am[i, tmp] <- 1
#'   }
#'   A[[m]]       <- Am
#' }
#' Anorm  <- norm.network(A) #Row-normalization
#' 
#' # X
#' X      <- cbind(rnorm(n, 1, 3), rexp(n, 0.4))
#' 
#' # cost group defining as X1 > 1
#' crg    <- 1*(X[,1] > 0.95)
#' 
#' # Networks
#' # Let us capture heterogeneity in peer effects using crg
#' # As we have 2 crg, we can estimate 4 peer effects. 
#' # Influence of friends in 1st crg on individual in 1st crg
#' # Incluence of friends in 2nd crg on individual in 1st crg
#' # Influence of friends in 1st crg on individual in 2st crg
#' # Incluence of friends in 2nd crg on individual in 2st crg
#' G        <- list()
#' cums     <- c(0, cumsum(nvec))
#' for (m in 1:M) {
#'   tp     <- crg[(cums[m] + 1):(cums[m + 1])]
#'   Am     <- A[[m]]
#'   G[[m]] <- norm.network(list(Am * (tp %*% t(tp)), 
#'                               Am * ((1 - tp) %*% t(tp)), 
#'                               Am * (tp %*% t(1 - tp)), 
#'                               Am * ((1 - tp) %*% t(1 - tp))))
#' }
#' 
#' # Parameters
#' lambda <- c(0.2, 0.3, 0.15, 0.25) #Four definitions of networks
#' Gamma  <- c(2.5, 2.2, -0.9, 1.5, -1.2)
#' delta  <- rep(c(0.6, 0.47, 0.35, 0.2, 0.05), 2) # two cost groups.
#' 
#' # Data
#' data   <- data.frame(X, peer.avg(Anorm, cbind(x1 = X[,1], x2 =  X[,2])))
#' colnames(data) = c("x1", "x2", "gx1", "gx2")
#' 
#' ytmp   <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2), 
#'                    lambda = lambda, Gamma = Gamma, delta = delta, cost.group = crg,  
#'                    data = data)
#' y      <- ytmp$y
#' hist(y, breaks = max(y) + 1)
#' table(y)
#' }
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm
#' @export
simcdnet   <- function(formula,
                       Glist,
                       parms,
                       lambda,
                       Gamma,
                       delta,
                       cost.group,
                       Rmax,
                       Rbar,
                       tol   = 1e-10,
                       maxit = 500,
                       data) {
  if(missing(Rmax)) Rmax <- Inf
  stopifnot((Rmax >= 1) & (Rmax <= Inf))
  if(!missing(Rbar)){
    stopifnot(all(Rbar <= Rmax))
    stopifnot(all(Rbar >= 1))
  }
  stopifnot(inherits(Glist, c("list", "matrix", "array")))
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  
  if(inherits(Glist[[1]], "list")){
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, "list"))))
  } else if(inherits(Glist[[1]], c("matrix", "array"))) {
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, c("matrix", "array")))))
    Glist  <- lapply(Glist, function(x_) list(x_))
  }
  
  # Sizes
  M        <- length(Glist)
  nvec     <- sapply(Glist, function(x_) nrow(x_[[1]]))
  sumn     <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  nCl      <- length(Glist[[1]])
  
  # cost.group
  if(missing(cost.group)){
    cost.group <- rep(0, sumn)
  }
  if(length(cost.group) != sumn) stop("length(cost.group) != n")
  uCa      <- sort(unique(cost.group))
  nCa      <- length(uCa)
  lCa      <- lapply(uCa, function(x_) which(cost.group == x_) - 1)
  na       <- sapply(lCa, length)
  if(!missing(Rbar)){
    if(length(Rbar) == 1) Rbar = rep(Rbar, nCa)
    if(nCa != length(Rbar)) stop("length(Rbar) is not equal to the number of cost groups.")
  }
  
  
  # Parameters
  if(!missing(parms)){
    if(missing(lambda) & missing(Gamma) & missing(delta)) warnings("parms is defined: lambda, Gamma, or delta is ignored.")
    lambda  <- parms[1:nCl]
    delta   <- tail(parms, sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))
    Gamma   <- parms[(nCl + 1):(length(parms) - sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))]
  } else{
    if(missing(lambda) | missing(Gamma)) stop("lambda or Gamma is missing.")
    if(missing(delta)){
      if(Rmax > 1){
        stop("delta is missing.")
      } else {
        delta <- c()
      }
    }
    if(nCa == 1){
      if(!missing(Rbar)){
        if(ifelse(Rbar == Rmax, Rbar - 1, Rbar) != length(delta)) stop("length(delta) does not match Rbar and Rmax.")
      }
    }
    if(length(delta) != sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar))) stop("length(delta) does not match Rbar and Rmax.")
    if(nCl != length(lambda)) stop("length(lambda) is not equal to the number of specifications of networks.")
  }
  
  ndelta  <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  stopifnot(all(Rbar <= Rmax))
  stopifnot(all(delta > 0))
  delta   <- delta + sum(lambda)
  idelta  <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  
  # data
  f.t.data  <- formula.to.data(formula = formula, contextual = FALSE, Glist = Glist, M = M, igr = igr, 
                               data = data, type = "sim", theta0  = 0)
  X         <- f.t.data$X
  if(ncol(X) != length(Gamma)) stop("ncol(X) == length(Gamma) is not true.")
  coln      <- c("lambda", colnames(X))
  if(nCl > 1) {
    coln    <- c(paste0(coln[1], ":", 1:nCl), coln[-1])
  }
  
  # xb
  xb        <- c(X %*% Gamma)
  eps       <- rnorm(sumn, 0, 1)
  
  # E(y) and G*E(y)
  Ey        <- rep(0, sumn)
  GEy       <- matrix(0, sumn, nCl)
  t         <- fye(ye = Ey, Gye = GEy, G = Glist, lCa = lCa, nCa = nCa, nCl = nCl, igroup = igr, 
                   ngroup = M, psi = xb, lambda = lambda, delta = delta, idelta = idelta, n = na, sumn = sumn,
                   Rbar = Rbar, R = Rmax, tol = tol, maxit = maxit)
  
  # y
  Ztlamda  <- c(GEy %*% lambda + xb)
  yst      <- Ztlamda + eps
  y        <- c(fy(yst = yst, maxyst = sapply(lCa, function(x_) max(yst[x_ + 1])), lCa = lCa, nCa = nCa, delta = delta, 
                   idelta = idelta, n = na, sumn = sumn, Rbar = Rbar, R = Rmax))
  if(nCl == 1){
    GEy    <- c(GEy)
  }
  
  # marginal effects
  GamWI    <- c(lambda, Gamma)
  colWI    <- coln
  if("(Intercept)" %in% colWI){
    GamWI  <- GamWI[colWI != "(Intercept)"]
    colWI  <- colWI[colWI != "(Intercept)"]
  }
  meffects <- fmeffects(ZtLambda = Ztlamda, lbeta = GamWI, lCa = lCa, nCa = nCa, delta = delta, idelta = idelta,
                        n = na, sumn = sumn, Rbar = Rbar, R = Rmax)
  
  # rho != 0
  #   t        <- fyb2(yb, Gyb, Glist, igr, M, xb, lambda, delta, deltabar, rho, n, Rbar, tol, maxit)
  #   Ztlamda  <- lambda*Gyb + xb
  #   yst      <- Ztlamda + eps
  #   y        <- c(fy2(yst, max(yst), lambda, delta, deltabar, rho, n, Rbar))
  #   meffects <- fmeffects2(n, lambda, delta, deltabar, rho, Rbar, Ztlamda, thetaWI)
  
  
  Rmax            <- meffects$Rmax
  meffects        <- meffects$meffects; colnames(meffects) <- colWI
  
  
  list("yst"            = yst,
       "y"              = y,
       "Ey"             = Ey,
       "GEy"            = GEy,
       "marg.effects"   = meffects,
       "Rmax"           = Rmax,
       "iteration"      = c(t))
}


#' @title Estimating count data models with social interactions under rational expectations using the NPL method
#' @param formula a class object \link[stats]{formula}: a symbolic description of the model. `formula` must be as, for example, \code{y ~ x1 + x2 + gx1 + gx2}
#' where `y` is the endogenous vector and `x1`, `x2`, `gx1` and `gx2` are control variables, which can include contextual variables, i.e. averages among the peers.
#' Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogenous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param cost.group the vector indicating the groups in which individuals have the same cost function (same cut-points in the ordered model). The default assumes a common group.
#' @param Rmax an integer indicating the theoretical upper bound of `y`. (see the model specification in details).
#' @param Rbar an \eqn{L}-vector, where  \eqn{L} is the number of cost groups. For large `Rmax` the cost function is assumed to be semi-parametric (i.e., nonparametric from 0 to \eqn{\bar{R}} and quadratic beyond \eqn{\bar{R}}). 
#' @param starting (optional) a starting value for \eqn{\theta = (\lambda, \Gamma', \delta')'}, where \eqn{\lambda}, \eqn{\Gamma}, and \eqn{\delta} are the parameters to be estimated (see details).
#' @param Ey0 (optional) a starting value for \eqn{E(y)}.
#' @param optimizer is either `fastlbfgs` (L-BFGS optimization method of the package \pkg{RcppNumerical}), `nlm` (referring to the function \link[stats]{nlm}), or `optim` (referring to the function \link[stats]{optim}). 
#' Arguments for these functions such as, `control` and `method` can be set via the argument `opt.ctr`.
#' @param npl.ctr a list of controls for the NPL method (see details).
#' @param opt.ctr a list of arguments to be passed in `optim_lbfgs` of the package \pkg{RcppNumerical}, \link[stats]{nlm} or \link[stats]{optim} (the solver set in `optimizer`), such as `maxit`, `eps_f`, `eps_g`, `control`, `method`, etc.
#' @param cov a Boolean indicating if the covariance should be computed.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `cdnet` is called.
#' @return A list consisting of:
#'     \item{info}{a list of general information about the model.}
#'     \item{estimate}{the NPL estimator.}
#'     \item{Ey}{\eqn{E(y)}, the expectation of y.}
#'     \item{GEy}{the average of \eqn{E(y)} friends.}
#'     \item{cov}{a list including (if `cov == TRUE`) `parms` the covariance matrix and another list `var.comp`, which includes `Sigma`, as \eqn{\Sigma}, and `Omega`, as \eqn{\Omega}, matrices used for
#'     compute the covariance matrix.}
#'     \item{details}{step-by-step output as returned by the optimizer.}
#' @description
#' `cdnet` estimates count data models with social interactions under rational expectations using the NPL algorithm (see Houndetoungan, 2024). 
#' @details 
#' ## Model
#' The count variable \eqn{y_i}{yi} take the value \eqn{r} with probability. 
#' \deqn{P_{ir} = F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r}) - F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r + 1}).}
#' In this equation, \eqn{\mathbf{z}_i} is a vector of control variables; \eqn{F} is the distribution function of the standard normal distribution;
#' \eqn{\bar{y}_i^{e,s}} is the average of \eqn{E(y)} among peers using the `s`-th network definition;
#' \eqn{a_{h(i),r}} is the `r`-th cut-point in the cost group \eqn{h(i)}. \cr\cr
#' For identification, \eqn{a_{h(i),0} = -\infty} and \eqn{a_{h(i),1} = 0}. Moreover, \eqn{a_{h(i),r} = \infty} for any \eqn{r \geq R_{\text{max}} + 1}, 
#' which implies that \eqn{P_{ir} = 0} for any \eqn{r \geq R_{\text{max}} + 1}.
#' For any \eqn{r \geq 1}, the distance between two cut-points is \eqn{a_{h(i),r+1} - a_{h(i),r} =  \delta_{h(i),r} + \sum_{s = 1}^S \lambda_s}
#' As the number of cut-point can be large, a quadratic cost function is considered for \eqn{r \geq \bar{R}_{h(i)}}, where \eqn{\bar{R} = (\bar{R}_{1}, ..., \bar{R}_{L})}.
#' With the semi-parametric cost-function,
#' \eqn{a_{h(i),r + 1} - a_{h(i),r}= \bar{\delta}_{h(i)} + \sum_{s = 1}^S \lambda_s}.  \cr\cr
#' The model parameters are: \eqn{\lambda = (\lambda_1, ..., \lambda_S)'}, \eqn{\Gamma}, and \eqn{\delta = (\delta_1', ..., \delta_L')'}, 
#' where \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l}, \bar{\delta}_l)'} for \eqn{l = 1, ..., L}. 
#' The number of single parameters in \eqn{\delta_l} depends on  \eqn{R_{\text{max}}} and \eqn{\bar{R}_{l}}. The components \eqn{\delta_{l,2}, ..., \delta_{l,\bar{R}_l}} or/and 
#' \eqn{\bar{\delta}_l} must be removed in certain cases.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_{l} \geq 2}, then \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l})'}.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_{l} = 1} (binary models), then \eqn{\delta_l} must be empty.\cr
#' If \eqn{R_{\text{max}} > \bar{R}_{l} = 1}, then \eqn{\delta_l = \bar{\delta}_l}.
#' ## \code{npl.ctr}
#' The model parameters are estimated using the Nested Partial Likelihood (NPL) method. This approach 
#' starts with a guess of \eqn{\theta} and \eqn{E(y)} and constructs iteratively a sequence
#' of \eqn{\theta} and \eqn{E(y)}. The solution converges when the \eqn{\ell_1}-distance
#' between two consecutive \eqn{\theta} and \eqn{E(y)} is less than a tolerance. \cr
#' The argument \code{npl.ctr} must include
#' \describe{
#' \item{tol}{the tolerance of the NPL algorithm (default 1e-4),}
#' \item{maxit}{the maximal number of iterations allowed (default 500),}
#' \item{print}{a boolean indicating if the estimate should be printed at each step.}
#' \item{S}{the number of simulations performed use to compute integral in the covariance by important sampling.} 
#' }
#' @references 
#' Houndetoungan, E. A. (2024). Count Data Models with Social Interactions under Rational Expectations. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @seealso \code{\link{sart}}, \code{\link{sar}}, \code{\link{simcdnet}}.
#' @examples 
#' \donttest{
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 50, 100))
#' n      <- sum(nvec)
#' 
#' # Adjacency matrix
#' A      <- list()
#' for (m in 1:M) {
#'   nm           <- nvec[m]
#'   Am           <- matrix(0, nm, nm)
#'   max_d        <- 30 #maximum number of friends
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
#'     Am[i, tmp] <- 1
#'   }
#'   A[[m]]       <- Am
#' }
#' Anorm  <- norm.network(A)
#' 
#' # X
#' X      <- cbind(rnorm(n, 1, 3), rexp(n, 0.4))
#' 
#' # cost group defining as X1 > 1
#' crg    <- 1*(X[,1] > 0.95)
#' 
#' # Networks
#' # Let us capture heterogeneity in peer effects using crg
#' # As we have 2 crg, we can estimate 4 peer effects.
#' # Influence of friends in 1st crg on individual in 1st crg
#' # Incluence of friends in 2nd crg on individual in 1st crg
#' # Influence of friends in 1st crg on individual in 2st crg
#' # Incluence of friends in 2nd crg on individual in 2st crg
#' G        <- list()
#' cums     <- c(0, cumsum(nvec))
#' for (m in 1:M) {
#'   tp     <- crg[(cums[m] + 1):(cums[m + 1])]
#'   Am     <- A[[m]]
#'   G[[m]] <- norm.network(list(Am * (tp %*% t(tp)),
#'                               Am * ((1 - tp) %*% t(tp)),
#'                               Am * (tp %*% t(1 - tp)),
#'                               Am * ((1 - tp) %*% t(1 - tp))))
#' }
#' 
#' # Parameters
#' lambda <- c(0.2, 0.3, 0.15, 0.25) #Four definitions of networks
#' Gamma  <- c(2.5, 2.2, -0.9, 1.5, -1.2)
#' delta  <- rep(c(0.6, 0.47, 0.35, 0.2, 0.05), 2) # two cost groups.
#' 
#' # Data
#' data   <- data.frame(X, peer.avg(Anorm, cbind(x1 = X[,1], x2 =  X[,2])))
#' colnames(data) = c("x1", "x2", "gx1", "gx2")
#' 
#' ytmp   <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2),
#'                    lambda = lambda, Gamma = Gamma, delta = delta, cost.group = crg,
#'                    data = data)
#' y      <- ytmp$y
#' hist(y, breaks = max(y) + 1)
#' table(y)
#' 
#' # Estimation
#' est    <- cdnet(formula = y ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2), cost.group = crg,
#'                 optimizer = "nlm", data = data, opt.ctr = list(iterlim = 5000, steptol = 1e-15))
#' 
#' # Print summary
#' summary(est)
#' print(est)
#' 
#' # Marginal effects
#' meffects.cdnet(est, Glist = G, data = data)
#' 
#' # Standard errors of marginal effects using simulations
#' meffects.cdnet(est, Glist = G, data = data, S = 100)
#' }
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom utils tail
#' @export
cdnet    <- function(formula,
                     Glist, 
                     cost.group,
                     Rmax,
                     Rbar,
                     starting  = list(lambda = NULL, Gamma = NULL, delta = NULL), 
                     Ey0       = NULL,
                     optimizer = "fastlbfgs", 
                     npl.ctr   = list(), 
                     opt.ctr   = list(), 
                     cov       = TRUE,
                     data) {
  if(missing(Rmax)) Rmax <- Inf
  stopifnot(optimizer %in% c("fastlbfgs", "optim", "nlm"))
  stopifnot((Rmax >= 1) & (Rmax <= Inf))
  env.formula <- environment(formula)
  
  # controls
  npl.print   <- npl.ctr$print
  npl.tol     <- npl.ctr$tol
  npl.maxit   <- npl.ctr$maxit
  npl.S       <- npl.ctr$S
  npl.incdit  <- npl.ctr$incdit
  if (is.null(npl.print)) {
    npl.print <- TRUE
  }
  if (is.null(npl.tol)) {
    npl.tol   <- 1e-4
  }
  if (is.null(npl.maxit)) {
    npl.maxit <- 500L
  }
  if (is.null(npl.S)) {
    npl.S     <- 1e3L
  }
  if (is.null(npl.incdit)) {
    npl.incdit<- 30L
  }
  
  # Network
  stopifnot(inherits(Glist, c("list", "matrix", "array")))
  if (!is.list(Glist)) {
    Glist  <- list(Glist)
  }
  
  if(inherits(Glist[[1]], "list")){
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, "list"))))
  } else if(inherits(Glist[[1]], c("matrix", "array"))) {
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, c("matrix", "array")))))
    Glist  <- lapply(Glist, function(x_) list(x_))
  }
  
  # Sizes
  M        <- length(Glist)
  nvec     <- sapply(Glist, function(x_) nrow(x_[[1]]))
  sumn     <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  nCl      <- length(Glist[[1]])
  
  # cost.group
  if(missing(cost.group)){
    cost.group <- rep(0, sumn)
  }
  if(is.null(cost.group)){
    cost.group <- rep(0, sumn)
  }
  
  if(length(cost.group) != sumn) stop("length(cost.group) != n")
  uCa      <- sort(unique(cost.group))
  nCa      <- length(uCa)
  lCa      <- lapply(uCa, function(x_) which(cost.group == x_) - 1)
  na       <- sapply(lCa, length)
  if(!missing(Rbar)){
    if(length(Rbar) == 1) Rbar = rep(Rbar, nCa)
    if(nCa != length(Rbar)) stop("length(Rbar) is not equal to the number of cost groups.")
  }
  
  # starting values
  thetat   <- c(starting$lambda, starting$Gamma)
  Deltat   <- starting$delta
  if(!inherits(starting, "list") & !is.null(starting)) {
    Deltat <- tail(starting, sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))
    Lmdt   <- starting[1:nCl]
    Gamt   <- starting[(nCl + 1):(length(starting) - sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))]
    thetat <- c(Lmdt, Gamt)
  }
  if(!is.null(Deltat)){
    if(length(Deltat) != sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar))) stop("`length(delta)` does not match Rbar and Rmax.")
  }
  tmp      <- c(is.null(thetat), is.null(Deltat))
  if(all(tmp) != any(tmp)){
    stop("Starting parameters are defined, but not all of them.")
  }
  
  # data
  f.t.data <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = 0) #because G are directly included in X
  formula  <- f.t.data$formula
  y        <- f.t.data$y
  X        <- f.t.data$X
  
  if(!is.null(thetat)){
    if((ncol(X) + nCl) != length(thetat)) stop("ncol(X) == length(Gamma) is not true.")
  }
  coln     <- c("lambda", colnames(X))
  if(nCl > 1) {
    coln   <- c(paste0(coln[1], ":", 1:nCl), coln[-1])
  }
  maxy     <- sapply(lCa, function(x_) max(y[x_ + 1]))
  if(any(maxy > Rmax)) stop("`Rmax` is lower than y.")
  K        <- ncol(X)
  
  # Ey 
  Eyt      <- rep(0, sumn)
  if (!is.null(Ey0)) {
    Eyt    <- Ey0
  }
  
  # GEyt
  GEyt     <- lapply(1:nCl, function(x_2) unlist(lapply(1:M, function(x_1) Glist[[x_1]][[x_2]] %*% Eyt[(igr[x_1, 1]:igr[x_1, 2]) + 1])))
  GEyt     <- as.matrix(as.data.frame(GEyt)); colnames(GEyt) <- head(coln, nCl)
  
  # Rbar
  if(is.null(Deltat) & missing(Rbar)) {
    Rbar   <- rep(1, nCa)
  }
  stopifnot(all(Rbar >= 1))
  stopifnot(all(Rbar <= Rmax))
  if(any(sapply(1:nCa, function(x_) Rbar[x_] > max(y[lCa[[x_]]])))) {
    stop("Rbar > max(y) in at least one cost group.")
  }
  
  # check starting and computed them if not defined
  lmbd0       <- NULL  
  ndelta      <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  if (!is.null(thetat)) {
    if(length(thetat) != (K + nCl)) {
      stop("theta length is inappropriate.")
    }
    if(any(head(thetat, nCl) < 0)) stop("Negative peer effects are not supported in this version.")
    lmbd0     <- head(thetat, nCl)
    thetat    <- c(log(head(thetat, nCl)/(1 - head(thetat, nCl))), tail(thetat, K))
    if(any(Deltat <= 0)) stop("all(delta) > 0 is not true for starting values.")
  } else {
    Gy        <- lapply(1:nCl, function(x_2) unlist(lapply(1:M, function(x_1) Glist[[x_1]][[x_2]] %*% y[(igr[x_1, 1]:igr[x_1,2]) + 1])))
    Gy        <- as.matrix(as.data.frame(Gy)); colnames(Gy) <- head(coln, nCl)
    Xtmp      <- cbind(Gy, X)
    b         <- solve(crossprod(Xtmp), crossprod(Xtmp, y))
    b         <- b/sqrt(sum((y - Xtmp%*%b)^2)/sumn)
    b[1:nCl]  <- sapply(b[1:nCl], function(x_) min(max(x_, 0.01), 0.99)) 
    lmbd0     <- b[1:nCl]
    Deltat    <- rep(0.01, sum(ndelta))
    thetat    <- c(log(head(b, nCl)/(1 - head(b, nCl))), tail(b, K))
  }
  
  idelta      <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  lDelta      <- log(Deltat)
  thetat      <- c(thetat, lDelta)
  
  # other variables
  cont     <- TRUE
  t        <- 0
  theta    <- NULL
  covout   <- NULL
  REt      <- NULL
  llht     <- 0
  par0     <- NULL
  par1     <- NULL
  like     <- NULL
  var.comp <- NULL
  steps    <- list()
  
  # arguments
  ctr      <- c(list(Gye = GEyt, X = X, lCa = lCa, nCa = nCa, nCl = nCl, K = K, n = na, sumn = sumn, 
                     idelta = idelta, ndelta = ndelta, Rbar = Rbar, R = Rmax, y = y, maxy = maxy) , opt.ctr)
  
  # Arguments used in the optimizer
  if (optimizer == "fastlbfgs"){
    ctr    <- c(ctr, list(par = thetat)); 
    optimizer <- "cdnetLBFGS"
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else if (optimizer == "optim") {
    ctr    <- c(ctr, list(fn = foptimREM_NPL, par = thetat))
    par0   <- "par"
    par1   <- "par"
    like   <- "value"
  } else {
    ctr    <- c(ctr, list(f = foptimREM_NPL,  p = thetat))
    par0   <- "p"
    par1   <- "estimate"
    like   <- "minimum"
  }
  
  # NPL algorithm
  ninc.d   <- 0
  dist0    <- Inf
  
  while(cont) {
    # tryCatch({
    Eyt0        <- Eyt + 0    #copy in different memory
    
    # compute theta
    REt         <- do.call(get(optimizer), ctr)
    thetat      <- REt[[par1]]
    llht        <- -REt[[like]]
    theta       <- c(1/(1 + exp(-head(thetat, nCl))), thetat[(nCl + 1):(K + nCl)], exp(tail(thetat, sum(ndelta))) + 1e-323)
    
    # compute y
    fL_NPL(ye = Eyt, Gye = GEyt, theta = thetat, X = X, G = Glist, lCa = lCa, nCa = nCa, nCl = nCl, igroup = igr, 
           ngroup = M, K = K, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta, Rbar = Rbar, R = Rmax)
    
    # distance
    var.eps     <- abs(c((ctr[[par0]] - thetat)/thetat, (Eyt0 - Eyt)/Eyt))
    if(all(!is.finite(var.eps))){
      var.eps   <- abs(c(ctr[[par0]] - thetat, Eyt0 - Eyt))
    }
    dist        <- max(var.eps, na.rm = TRUE)
    ninc.d      <- (ninc.d + 1)*(dist > dist0) #counts the successive number of times distance increases
    dist0       <- dist
    cont        <- (dist > npl.tol & t < (npl.maxit - 1))
    t           <- t + 1
    REt$dist    <- dist
    ctr[[par0]] <- thetat
    steps[[t]]  <- REt
    
    if(npl.print){
      cat("---------------\n")
      cat(paste0("Step          : ", t), "\n")
      cat(paste0("Distance      : ", round(dist,3)), "\n")
      cat(paste0("Likelihood    : ", round(llht,3)), "\n")
      cat("Estimate:", "\n")
      print(theta)
    }
    
    if((ninc.d > npl.incdit) | (llht < -1e293)) {
      cat("** Non-convergence ** Redefining theta and computing a new E(y)\n")
      thetat[1:nCl] <- -4.5
      dist0         <- Inf
      fnewye(ye = Eyt, Gye = GEyt, theta = thetat, X = X, G = Glist, lCa = lCa, nCa = nCa, nCl = nCl, igroup = igr, 
             ngroup = M, K = K, n = na, sumn = sumn, idelta = idelta, Rbar = Rbar, tol = npl.tol, maxit = npl.maxit)
      ctr[[par0]] <- thetat
    }
  }
  
  # name theta
  namtheta          <- c(coln, unlist(lapply(1:nCa, function(x_){
    if((Rbar[x_] == 1) & (Rmax == 1)) return(c())
    if(Rbar[x_] == Rmax) return(paste0("delta:", x_, ".", 2:(Rbar[x_])))
    if(Rbar[x_] > 1) return(c(paste0("delta:", x_, ".", 2:(Rbar[x_])), paste0("deltabar:", x_)))
    paste0("deltabar:", x_)
  })))
  names(theta)      <- namtheta
  lbda              <- head(theta, nCl)
  Gam               <- theta[(nCl + 1):(nCl + K)]
  Del               <- tail(theta, sum(ndelta))
  
  environment(formula) <- env.formula
  sdata                <- list(
    "formula"       = formula,
    "Glist"         = deparse(substitute(Glist))
  )
  if (!missing(data)) {
    sdata             <- c(sdata, list("data" = deparse(substitute(data))))
  }
  
  if (npl.maxit == t) {
    warning("NPL: maximum number of iterations has been reached.")
  }
  # covariance and ME
  var.comp          <- NULL
  covt              <- NULL
  if(cov){
    tmp             <- fcovCDI(theta = thetat, Gye = GEyt, X = X, G = Glist, lCa = lCa, nCa = nCa, nCl = nCl,
                               igroup = igr, ngroup = M, K = K, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta,
                               Rbar = Rbar, R = Rmax, S = npl.S)
    var.comp        <- tmp$var.comp
    covt            <- tmp$covt
    # Rmax            <- tmp$Rmax
    
    namecovt        <- c(coln, unlist(lapply(1:nCa, function(x_){
      if((Rbar[x_] == 1) & (Rmax == 1)) return(c())
      if(Rbar[x_] == Rmax) return(paste0("log(delta:", x_, ".", 2:(Rbar[x_]), ")"))
      if(Rbar[x_] > 1) return(c(paste0("log(delta:", x_, ".", 2:(Rbar[x_]), ")"), paste0("log(deltabar:", x_, ")")))
      paste0("log(deltabar:", x_, ")")
    })))
    colnames(covt)  <- namecovt
    rownames(covt)  <- namecovt
    namecovt[1:nCl] <- "log(lambda)"; if(nCl > 1) namecovt[1:nCl] <- paste0("log(lambda:", 1:nCl, ")")
    rownames(var.comp$Sigma) <- namecovt
    rownames(var.comp$Omega) <- namecovt
    colnames(var.comp$Sigma) <- namecovt
    colnames(var.comp$Omega) <- namecovt
  }
  
  AIC                  <- 2*length(theta) - 2*llht
  BIC                  <- length(theta)*log(sumn) - 2*llht
  
  INFO                 <- list("M"          = M,
                               "n"          = nvec,
                               "Kz"         = K,
                               "formula"    = formula,
                               "n.lambda"   = nCl,
                               "Rbar"       = Rbar,
                               "Rmax"       = Rmax,
                               "cost.group" = cost.group,
                               "log.like"   = llht, 
                               "npl.iter"   = t,
                               "AIC"        = AIC,
                               "BIC"        = BIC)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(parms = theta, lambda =  lbda, Gamma = Gam, delta = Del),
                               "Ey"         = Eyt, 
                               "GEy"        = GEyt,
                               "cov"        = list(parms = covt, var.comp = var.comp),
                               "details"    = steps)
  class(out)           <- "cdnet"
  out
}


#' @title Summary for the estimation of count data models with social interactions under rational expectations
#' @description Summary and print methods for the class `cdnet` as returned by the function \link{cdnet}.
#' @param object an object of class `cdnet`, output of the function \code{\link{cdnet}}.
#' @param x an object of class `summary.cdnet`, output of the function \code{\link{summary.cdnet}}
#' or class `cdnet` of the function \code{\link{cdnet}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogenous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `summary.cdnet` is called.
#' @param S number of simulations to be used to compute integral in the covariance by important sampling.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export 
"summary.cdnet" <- function(object,
                            Glist,
                            data,
                            S = 1e3L,
                            ...) {
  stopifnot(class(object) == "cdnet")
  out           <- c(object, list("..." = ...))
  if(is.null(object$cov$parms)){
    env.formula <- environment(object$info$formula)
    Rbar        <- object$info$Rbar
    Rmax        <- object$info$Rmax
    Kz          <- object$info$Kz
    nCl         <- object$info$n.lambda
    cost.group  <- object$info$cost.group
    M           <- object$info$M
    nvec        <- object$info$n
    sumn        <- sum(nvec)
    nCa         <- length(Rbar)
    parms       <- object$estimate$parms
    lambda      <- object$estimate$lambda
    Gamma       <- object$estimate$Gamma
    delta       <- object$estimate$delta
    thetat      <- c(log(lambda/(1 - lambda)), Gamma, log(delta))
    ndelta      <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
    idelta      <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
    
    # Cost group
    uCa         <- sort(unique(cost.group))
    lCa         <- lapply(uCa, function(x_) which(cost.group == x_) - 1)
    na          <- sapply(lCa, length)
    
    # Network
    stopifnot(inherits(Glist, c("list", "matrix", "array")))
    if (!is.list(Glist)) {
      Glist     <- list(Glist)
    }
    if(inherits(Glist[[1]], "list")){
      stopifnot(all(sapply(Glist, function(x_) inherits(x_, "list"))))
    } else if(inherits(Glist[[1]], c("matrix", "array"))) {
      stopifnot(all(sapply(Glist, function(x_) inherits(x_, c("matrix", "array")))))
      Glist     <- lapply(Glist, function(x_) list(x_))
    }
    if(M != length(Glist)) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
    if(any(nvec != sapply(Glist, function(x_) nrow(x_[[1]])))) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
    if(nCl != length(Glist[[1]])) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
    igr         <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
    
    # Data
    GEyt        <- object$GEy
    npl.S       <- S
    if (is.null(npl.S)) {
      npl.S     <- 1e3L
    }
    formula     <- object$info$formula
    f.t.data    <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = 0) #because G are directly included in X
    X           <- f.t.data$X
    if(Kz != ncol(X)) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
    if(sumn != nrow(X)) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
    if(any(colnames(X) != names(Gamma))) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
    coln        <- c("lambda", colnames(X))
    if(nCl > 1) {
      coln      <- c(paste0(coln[1], ":", 1:nCl), coln[-1])
    }
    
    tmp         <- fcovCDI(theta = thetat, Gye = GEyt, X = X, G = Glist, lCa = lCa, nCa = nCa, nCl = nCl,
                           igroup = igr, ngroup = M, K = Kz, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta,
                           Rbar = Rbar, R = Rmax, S = npl.S)
    
    var.comp        <- tmp$var.comp
    covt            <- tmp$covt
    # out$info$Rmax   <- tmp$Rmax
    
    namecovt        <- c(coln, unlist(lapply(1:nCa, function(x_){
      if((Rbar[x_] == 1) & (Rmax == 1)) return(c())
      if(Rbar[x_] == Rmax) return(paste0("delta:", x_, ".", 2:(Rbar[x_])))
      if(Rbar[x_] > 1) return(c(paste0("log(delta:", x_, ".", 2:(Rbar[x_]), ")"), paste0("log(deltabar:", x_, ")")))
      paste0("log(deltabar:", x_, ")")
    })))
    colnames(covt)  <- namecovt
    rownames(covt)  <- namecovt
    namecovt[1:nCl] <- "log(lambda)"; if(nCl > 1) namecovt[1:nCl] <- paste0("log(lambda:", 1:nCl, ")")
    rownames(var.comp$Sigma) <- namecovt
    rownames(var.comp$Omega) <- namecovt
    colnames(var.comp$Sigma) <- namecovt
    colnames(var.comp$Omega) <- namecovt
    
    out$cov           <- list(parms = covt, var.comp = var.comp)
  }
  class(out) <- "summary.cdnet"
  out
}


#' @rdname summary.cdnet
#' @importFrom stats pchisq
#' @export
"print.summary.cdnet"  <- function(x, ...) {
  stopifnot(class(x) == "summary.cdnet")
  
  M                    <- x$info$M
  n                    <- x$info$n
  iteration            <- x$info$npl.iter
  Rbar                 <- x$info$Rbar
  csRbar               <- c(0, cumsum(Rbar))
  formula              <- x$info$formula
  Kz                   <- x$info$Kz
  AIC                  <- x$info$AIC
  BIC                  <- x$info$BIC
  nCa                  <- length(Rbar)
  nCl                  <- x$info$n.lambda
  
  coef                 <- c(x$estimate$lambda, x$estimate$Gamma)
  K                    <- length(coef)
  std                  <- sqrt(head(diag(x$cov$parms), K))
  delta                <- x$estimate$delta
  
  llh                  <- x$info$log.like
  
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:6)], list(...))
  
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Call:\n")
  print(formula)
  cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, sep = "", "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, sep = "", "\n")
  cat("Sample size              : ", sum(n), sep = "", "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  cat("Cost function -- Number of groups: ", nCa, "\n", sep = "")
  
  for(k in 1:nCa){
    cat("Group ", k, ", Rbar = ", Rbar[k], "\n", sep = "")
    if(Rbar[k] > 1) cat("delta:", delta[(csRbar[k] + 1):(csRbar[k + 1] - 1)], "\n", sep = " ")
    cat("deltabar: ", delta[csRbar[k + 1]], "\n\n", sep = "")
  }
  
  
  cat("log pseudo-likelihood: ", llh, sep = "", "\n")
  cat("AIC: ", AIC, " -- BIC: ", BIC, sep = "", "\n")
  
  invisible(x)
}

#' @rdname summary.cdnet
#' @export
"print.cdnet" <- function(x, ...) {
  stopifnot(class(x) == "cdnet")
  print(summary(x, ...))
}

#' @title Computing marginal effects for count data models with social interactions under rational expectations
#' @description \code{meffects.cdnet} computes and prints marginal effects for the class `cdnet` as returned by the function \link{cdnet}.
#' Standard errors for the marginal effects are computed by simulation following Krinsky and Robb (1989).
#' @param object an object of class `cdnet` or `summary.cdnet`, as returned by the function \code{\link{cdnet}} or \code{\link{summary.cdnet}}.
#' @param list.index list of observation indices to use to calculate the average marginal effect for each variable (see the example of the function \code{\link{cdnet}}).
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogenous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param data an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
#' in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `summary.cdnet` is called.
#' @param S he number of simulations to use to calculate standard errors. Only marginal effects are calculated when `S = NULL`.
#' @param ncores ncores the number of cores to use for simulations. A large number of cores should make the process faster.
#' @references 
#' Krinsky, I., & Robb, A. L. (1986). On approximating the statistical properties of elasticities. The review of economics and statistics, 715-719, \doi{10.2307/1924536}.
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doRNG "%dorng%"
#' @importFrom doParallel registerDoParallel
#' @export
"meffects.cdnet" <- function(object, Glist, data, list.index = list(), S = NULL, ncores = 1L){
  if(!is.null(S)) stopifnot(S >= 2)
  stopifnot(inherits(object, c("summary.cdnet", "cdnet")))
  env.formula  <- environment(object$info$formula)
  Rbar         <- object$info$Rbar
  Rmax         <- object$info$Rmax
  Kz           <- object$info$Kz
  nCl          <- object$info$n.lambda
  cost.group   <- object$info$cost.group
  M            <- object$info$M
  nvec         <- object$info$n
  sumn         <- sum(nvec)
  nCa          <- length(Rbar)
  parms        <- object$estimate$parms
  lambda       <- object$estimate$lambda
  Gamma        <- object$estimate$Gamma
  delta        <- object$estimate$delta
  ndelta       <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  idelta       <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  
  
  # Cost group
  uCa          <- sort(unique(cost.group))
  lCa          <- lapply(uCa, function(x_) which(cost.group == x_) - 1)
  na           <- sapply(lCa, length)
  
  # Network
  stopifnot(inherits(Glist, c("list", "matrix", "array")))
  if (!is.list(Glist)) {
    Glist      <- list(Glist)
  }
  if(inherits(Glist[[1]], "list")){
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, "list"))))
  } else if(inherits(Glist[[1]], c("matrix", "array"))) {
    stopifnot(all(sapply(Glist, function(x_) inherits(x_, c("matrix", "array")))))
    Glist      <- lapply(Glist, function(x_) list(x_))
  }
  if(M != length(Glist)) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
  if(any(nvec != sapply(Glist, function(x_) nrow(x_[[1]])))) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
  if(nCl != length(Glist[[1]])) stop("`Glist` structure does not match. Perhaps another `Glist` is used in `cdnet`.")
  igr          <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  # Data
  GEyt         <- object$GEy
  npl.S        <- S
  if (is.null(npl.S)) {
    npl.S      <- 1e3L
  }
  formula      <- object$info$formula
  f.t.data     <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = 0) #because G are directly included in X
  X            <- f.t.data$X
  if(Kz != ncol(X)) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
  if(sumn != nrow(X)) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
  if(any(colnames(X) != names(Gamma))) stop("`data` structure does not match. Perhaps another `data` is used in `cdnet` or `formula` has been changed.")
  coln         <- c("lambda", colnames(X))
  if(nCl > 1) {
    coln       <- c(paste0(coln[1], ":", 1:nCl), coln[-1])
  }
  
  # list.index
  if(length(list.index) == 0){
    if("(Intercept)" %in% coln){
      coln       <- coln[coln != "(Intercept)"]
      list.index <- rep(list(1:sumn), length(coln)); names(list.index) <- coln
    }
  } else{
    if(is.null(names(list.index))) stop("index in `list.index` must be names elements in c(lambda, Gamma)")
    if(any(!(names(list.index) %in% coln))) stop("index in `list.index` must be names elements in c(lambda, Gamma)")
  }
  
  ZtLambda <- c(GEyt%*%lambda + X%*%Gamma)
  nlindex  <- names(list.index)
  lGam     <- c(lambda, Gamma)
  meff.i   <- fmeffects(ZtLambda = ZtLambda, lbeta = lGam[nlindex], lCa = lCa, nCa = nCa, delta = delta + sum(lambda), 
                           idelta = idelta, n = na, sumn = sumn, Rbar = Rbar, R = Rmax)
  meff     <- sapply(1:length(nlindex), function(x_){mean(meff.i$meffects[list.index[[x_]], x_])})
  names(meff)  <- nlindex
  
  covm         <- NULL
  if(!is.null(S)){
    covt       <- object$cov$parms
    if(is.null(covt)) stop("`cov` was set FALSE in `cdnet`.")
    for(s in 1:nCl){
      covt[s,] <- covt[s,]/lambda[s]
      covt[,s] <- covt[,s]/lambda[s]
    }
    thetat <- c(log(lambda/(1 - lambda)), Gamma, log(delta))
    
    cl     <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
    
    covm   <- stats::cov(t(foreach(s = 1:S, .combine = cbind, .packages  = "CDatanet") %dorng% {
      fsime(thetat, covt, GEyt, X, list.index, nlindex, lCa, nCa, nCl, Kz, idelta, na, sumn, Rbar, Rmax)}))
    colnames(covm) <- nlindex
    rownames(covm) <- nlindex
  }
  
  out         <- list("info"           = object$info,
                      "marg.effects.i" = meff.i,
                      "marg.effects"   = meff,
                      "cov"            = covm)
  class(out)  <- "meffects.cdnet"
  out
}

fsime     <- function(thetat, covt, GEyt, X, list.index, nlindex, lCa, nCa, nCl, Kz, idelta, na, sumn, Rbar, Rmax){
  thet1   <- mvrnorm(n = 1, mu = thetat, Sigma = covt)
  lam1    <- 1/(1 + exp(-head(thet1, nCl)))
  Gam1    <- thet1[(nCl + 1):(nCl + Kz)]
  lGam1   <- c(lam1, Gam1)
  del1    <- exp(tail(thet1, sum(Rbar)))    
  ZtL1    <- c(GEyt%*%lam1 + X%*%Gam1)
  meff1   <- fmeffects(ZtLambda = ZtL1, lbeta = lGam1[nlindex], lCa = lCa, nCa = nCa, delta = del1 + sum(lam1), 
                       idelta = idelta, n = na, sumn = sumn, Rbar = Rbar, R = Rmax)
  sapply(1:length(nlindex), function(x_){mean(meff1$meffects[list.index[[x_]], x_])})
}


#' @title Printing marginal effects for count data models with social interactions under rational expectations
#' @description Summary and print methods for the class `meffects.cdnet` as returned by the function \link{meffects.cdnet}.
#' @param object an object of class `meffects.cdnet` as return by the function \code{\link{meffects.cdnet}}.
#' @param x an object of class `meffects.cdnet` as return by the function \code{\link{meffects.cdnet}} or \code{\link{print.meffects.cdnet}}.
#' @param ... further arguments passed to or from other methods.
#' @export
"print.meffects.cdnet" <- function(x, ...){
  stopifnot(class(x) == "meffects.cdnet")
  M                    <- x$info$M
  n                    <- x$info$n
  iteration            <- x$info$npl.iter
  Rbar                 <- x$info$Rbar
  csRbar               <- c(0, cumsum(Rbar))
  formula              <- x$info$formula
  Kz                   <- x$info$Kz
  AIC                  <- x$info$AIC
  BIC                  <- x$info$BIC
  nCa                  <- length(Rbar)
  nCl                  <- x$info$n.lambda
  
  coef                 <- x$marg.effects
  std                  <- sqrt(diag(x$cov))
  
  llh                  <- x$info$log.like
  
  if(is.null(x$cov)){
    cat("Marginal effects:\n")
    print(coef)
    return(invisible(x))
  }
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:4)], list(...))
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Call:\n")
  print(formula)
  cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, sep = "", "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, sep = "", "\n")
  cat("Sample size              : ", sum(n), sep = "", "\n\n")
  
  cat("Marginal effects:\n")
  do.call("print", out_print)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  
  cat("log pseudo-likelihood: ", llh, sep = "", "\n")
  cat("AIC: ", AIC, " -- BIC: ", BIC, sep = "", "\n")
  
  invisible(x)
}

#' @rdname print.meffects.cdnet
#' @export
"summary.meffects.cdnet" <- function(object, ...) {
  stopifnot(class(object) == "meffects.cdnet")
  object
}