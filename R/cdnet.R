#' @title Simulating Count Data Models with Social Interactions Under Rational Expectations
#' @param formula A class object of class \link[stats]{formula}: a symbolic description of the model. `formula` should be specified, for example, as \code{y ~ x1 + x2 + gx1 + gx2}, where `y` is the endogenous vector and `x1`, `x2`, `gx1`, and `gx2` are control variables. These control variables can include contextual variables, such as averages among the peers. Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist An adjacency matrix or list of adjacency matrices. For networks consisting of multiple subnets (e.g., schools), `Glist` can be a list of subnet matrices, where the \eqn{m}-th element is an \eqn{n_m \times n_m} adjacency matrix, with \eqn{n_m} representing the number of nodes in the \eqn{m}-th subnet. 
#' For heterogeneous peer effects (`length(unique(group)) = h > 1`), the \eqn{m}-th element should be a list of \eqn{h^2} \eqn{n_m \times n_m} adjacency matrices corresponding to different network specifications (see Houndetoungan, 2024). 
#' For heterogeneous peer effects in a single large network, `Glist` should be a one-item list, where the item is a list of \eqn{h^2} network specifications. The order of these networks is important and must match `sort(unique(group))` (see examples).
#' @param group A vector indicating the individual groups. By default, this assumes a common group. If there are 2 groups (i.e., `length(unique(group)) = 2`, such as `A` and `B`), four types of peer effects are defined: 
#' peer effects of `A` on `A`, `A` on `B`, `B` on `A`, and `B` on `B`.
#' @param parms A vector defining the true values of \eqn{\theta = (\lambda', \Gamma', \delta')'} (see model specification in the details section). Each parameter \eqn{\lambda}, \eqn{\Gamma}, or \eqn{\delta} can also be provided separately to the arguments `lambda`, `Gamma`, or `delta`.
#' @param lambda The true value of the vector \eqn{\lambda}.
#' @param Gamma The true value of the vector \eqn{\Gamma}.
#' @param delta The true value of the vector \eqn{\delta}.
#' @param Rmax An integer indicating the theoretical upper bound of `y` (see model specification in detail).
#' @param Rbar An \eqn{L}-vector, where \eqn{L} is the number of groups. For large `Rmax`, the cost function is assumed to be semi-parametric (i.e., nonparametric from 0 to \eqn{\bar{R}} and quadratic beyond \eqn{\bar{R}}). The \eqn{l}-th element of `Rbar` indicates \eqn{\bar{R}} for the \eqn{l}-th value of `sort(unique(group))` (see model specification in detail).
#' @param tol The tolerance value used in the Fixed Point Iteration Method to compute the expectancy of `y`. The process stops if the \eqn{\ell_1}-distance between two consecutive \eqn{E(y)} is less than `tol`.
#' @param maxit The maximum number of iterations in the Fixed Point Iteration Method.
#' @param data An optional data frame, list, or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `simcdnet` is called.
#' @description
#' `simcdnet` simulates the count data model with social interactions under rational expectations developed by Houndetoungan (2024).
#' @details 
#' The count variable \eqn{y_i} takes the value \eqn{r} with probability.
#' \deqn{P_{ir} = F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r}) - F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r + 1}).}
#' In this equation, \eqn{\mathbf{z}_i} is a vector of control variables; \eqn{F} is the distribution function of the standard normal distribution;
#' \eqn{\bar{y}_i^{e,s}} is the average of \eqn{E(y)} among peers using the `s`-th network definition;
#' \eqn{a_{h(i),r}} is the `r`-th cut-point in the cost group \eqn{h(i)}. \cr\cr
#' The following identification conditions have been introduced: \eqn{\sum_{s = 1}^S \lambda_s > 0}, \eqn{a_{h(i),0} = -\infty}, \eqn{a_{h(i),1} = 0}, and 
#' \eqn{a_{h(i),r} = \infty} for any \eqn{r \geq R_{\text{max}} + 1}. The last condition implies that \eqn{P_{ir} = 0} for any \eqn{r \geq R_{\text{max}} + 1}.
#' For any \eqn{r \geq 1}, the distance between two cut-points is \eqn{a_{h(i),r+1} - a_{h(i),r} =  \delta_{h(i),r} + \sum_{s = 1}^S \lambda_s}.
#' As the number of cut-points can be large, a quadratic cost function is considered for \eqn{r \geq \bar{R}_{h(i)}}, where \eqn{\bar{R} = (\bar{R}_{1}, ..., \bar{R}_{L})}.
#' With the semi-parametric cost function,
#' \eqn{a_{h(i),r + 1} - a_{h(i),r} = \bar{\delta}_{h(i)} + \sum_{s = 1}^S \lambda_s}.  \cr\cr
#' The model parameters are: \eqn{\lambda = (\lambda_1, ..., \lambda_S)'}, \eqn{\Gamma}, and \eqn{\delta = (\delta_1', ..., \delta_L')'}, 
#' where \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l}, \bar{\delta}_l)'} for \eqn{l = 1, ..., L}. 
#' The number of single parameters in \eqn{\delta_l} depends on  \eqn{R_{\text{max}}} and \eqn{\bar{R}_l}. The components \eqn{\delta_{l,2}, ..., \delta_{l,\bar{R}_l}} or/and 
#' \eqn{\bar{\delta}_l} must be removed in certain cases.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_l \geq 2}, then \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l})'}.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_l = 1} (binary models), then \eqn{\delta_l} must be empty.\cr
#' If \eqn{R_{\text{max}} > \bar{R}_l = 1}, then \eqn{\delta_l = \bar{\delta}_l}.
#' @seealso \code{\link{cdnet}}, \code{\link{simsart}}, \code{\link{simsar}}.
#' @return A list consisting of:
#'     \item{yst}{\eqn{y^{\ast}}, the latent variable.}
#'     \item{y}{the observed count variable.}
#'     \item{Ey}{\eqn{E(y)}, the expectation of y.}
#'     \item{GEy}{the average of \eqn{E(y)} among peers.}
#'     \item{meff}{a list including average and individual marginal effects.}
#'     \item{Rmax}{infinite sums in the marginal effects are approximated by sums up to Rmax.}
#'     \item{iteration}{number of iterations performed by sub-network in the Fixed Point Iteration Method.}
#' @references 
#' Houndetoungan, A. (2024). Count Data Models with Heterogeneous Peer Effects. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @examples 
#' \donttest{
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 200)) # Random group sizes
#' n      <- sum(nvec) # Total number of individuals
#' 
#' # Adjacency matrix for each group
#' A      <- list()
#' for (m in 1:M) {
#'   nm           <- nvec[m] # Size of group m
#'   Am           <- matrix(0, nm, nm) # Empty adjacency matrix
#'   max_d        <- 30 # Maximum number of friends
#'   for (i in 1:nm) {
#'     tmp        <- sample((1:nm)[-i], sample(0:max_d, 1)) # Sample friends
#'     Am[i, tmp] <- 1 # Set friendship links
#'   }
#'   A[[m]]       <- Am # Add to the list
#' }
#' Anorm  <- norm.network(A) # Row-normalization of the adjacency matrices
#' 
#' # Covariates (X)
#' X      <- cbind(rnorm(n, 1, 3), rexp(n, 0.4)) # Random covariates
#' 
#' # Two groups based on first covariate
#' group  <- 1 * (X[,1] > 0.95) # Assign to groups based on x1
#' 
#' # Networks: Define peer effects based on group membership
#' # The networks should capture:
#' # - Peer effects of `0` on `0`
#' # - Peer effects of `1` on `0`
#' # - Peer effects of `0` on `1`
#' # - Peer effects of `1` on `1`
#' G        <- list()
#' cums     <- c(0, cumsum(nvec)) # Cumulative indices for groups
#' for (m in 1:M) {
#'   tp     <- group[(cums[m] + 1):(cums[m + 1])] # Group membership for group m
#'   Am     <- A[[m]] # Adjacency matrix for group m
#'   # Define networks based on peer effects
#'   G[[m]] <- norm.network(list(Am * ((1 - tp) %*% t(1 - tp)),
#'                               Am * ((1 - tp) %*% t(tp)),
#'                               Am * (tp %*% t(1 - tp)),
#'                               Am * (tp %*% t(tp))))
#' }
#' 
#' # Parameters for the model
#' lambda <- c(0.2, 0.3, -0.15, 0.25) 
#' Gamma  <- c(4.5, 2.2, -0.9, 1.5, -1.2)
#' delta  <- rep(c(2.6, 1.47, 0.85, 0.7, 0.5), 2) # Repeated values for delta
#' 
#' # Prepare data for the model
#' data   <- data.frame(X, peer.avg(Anorm, cbind(x1 = X[,1], x2 = X[,2]))) 
#' colnames(data) = c("x1", "x2", "gx1", "gx2") # Set column names
#' 
#' # Simulate outcomes using the `simcdnet` function
#' ytmp   <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2),
#'                    lambda = lambda, Gamma = Gamma, delta = delta, group = group,
#'                    data = data)
#' y      <- ytmp$y
#' 
#' # Plot histogram of the simulated outcomes
#' hist(y, breaks = max(y) + 1)
#' 
#' # Display frequency table of the simulated outcomes
#' table(y)
#' }
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rnorm
#' @export
simcdnet   <- function(formula,
                       group,
                       Glist,
                       parms,
                       lambda,
                       Gamma,
                       delta,
                       Rmax,
                       Rbar,
                       tol        = 1e-10,
                       maxit      = 500,
                       data) {
  if(missing(Rmax)) Rmax <- Inf
  if(is.null(Rmax)) Rmax <- Inf
  stopifnot((Rmax >= 1) & (Rmax <= Inf))
  if(!missing(Rbar)){
    stopifnot(all(Rbar <= Rmax))
    stopifnot(all(Rbar >= 1))
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
  M        <- length(Glist)
  nvec     <- sapply(Glist, function(x_) nrow(x_[[1]]))
  sumn     <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  nCl      <- length(Glist[[1]])
  
  # group
  if(missing(group)){
    group  <- rep(0, sumn)
  }
  if(is.null(group)){
    group  <- rep(0, sumn)
  }
  if(length(group) != sumn) stop("length(group) != n")
  uCa      <- sort(unique(group))
  nCa      <- length(uCa)
  lCa      <- lapply(uCa, function(x_) which(group == x_) - 1)
  na       <- sapply(lCa, length)
  if(!missing(Rbar)){
    if(length(Rbar) == 1) Rbar = rep(Rbar, nCa)
    if(nCa != length(Rbar)) stop("length(Rbar) is not equal to the number of groups.")
  }
  if((nCa^2) != nCl) stop("The number of network specifications does not match the number of groups.")
  
  
  # Parameters
  if(!missing(parms)){
    if(missing(lambda) & missing(Gamma) & missing(delta)) warning("parms is defined: lambda, Gamma, or delta is ignored.")
    lambda  <- parms[1:nCl]
    delta   <- tail(parms, sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))
    Gamma   <- parms[(nCl + 1):(length(parms) - sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))]
  } else{
    if(missing(lambda) | missing(Gamma)) stop("lambda or Gamma is missing.")
    if(missing(delta)){
      if(Rmax > 1){
        stop("delta is missing.")
      } else {
        delta <- numeric(0)
      }
    }
    if(nCa == 1){
      if(!missing(Rbar)){
        if(ifelse(Rbar == Rmax, Rbar - 1, Rbar) != length(delta)) stop("length(delta) does not match Rbar and Rmax.")
      } else {
        Rbar <- length(delta)
      }
    }
    if(length(delta) != sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar))) stop("length(delta) does not match Rbar and Rmax.")
    if(nCl != length(lambda)) stop("length(lambda) is not equal to the number of network specifications.")
  }
  
  ndelta  <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  stopifnot(all(Rbar <= Rmax))
  if(any(Rmax > 1)) stopifnot(all(delta > 0))
  idelta  <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  delta   <- c(fdelta(deltat = delta, lambda = lambda, idelta = idelta, ndelta = ndelta, nCa = nCa))
  
  # data
  f.t.data  <- formula.to.data(formula = formula, contextual = FALSE, Glist = Glist, M = M, igr = igr, 
                               data = data, type = "sim", theta0  = 0)
  X         <- f.t.data$X
  K         <- length(Gamma)
  if(ncol(X) != K) stop("ncol(X) = ", ncol(X), " is different from length(Gamma) = ", length(Gamma))
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
  t         <- fye(ye = Ey, Gye = GEy, G = Glist, lCa = lCa, nCa = nCa, igroup = igr, ngroup = M, 
                   psi = xb, lambda = lambda, delta = delta, idelta = idelta, n = na, sumn = sumn,
                   Rbar = Rbar, R = Rmax, tol = tol, maxit = maxit)
  if(t >= maxit) warning("Point-fixed: The maximum number of iterations has been reached.")
  # y
  Ztlamda  <- c(GEy %*% lambda + xb)
  yst      <- Ztlamda + eps
  y        <- c(fy(yst = yst, maxyst = sapply(lCa, function(x_) max(yst[x_ + 1])), lCa = lCa, nCa = nCa, delta = delta, 
                   idelta = idelta, n = na, sumn = sumn, Rbar = Rbar, R = Rmax))
  if(nCl == 1){
    GEy    <- c(GEy)
  }
  
  # marginal effects
  Gamma2   <- Gamma
  colnme   <- coln
  if("(Intercept)" %in% coln){
    Gamma2 <- Gamma[tail(coln, K) != "(Intercept)"]
    colnme <- colnme[coln != "(Intercept)"]
  }
  meff     <- fmeffects(ZtLambda = Ztlamda, lambda = lambda, Gamma2 = Gamma2, lCa = lCa, nCa = nCa,
                        delta = delta, idelta = idelta, sumn = sumn, Rbar = Rbar, R = Rmax)
  # Rmax     <- meff$Rmax
  imeff    <- meff$imeff; colnames(imeff) <- colnme
  
  
  list("yst"       = yst,
       "y"         = y,
       "Ey"        = Ey,
       "GEy"       = GEy,
       "meff"      = list(ameff = apply(imeff, 2, mean), imeff = imeff),
       "Rmax"      = Rmax,
       "iteration" = c(t))
}


#' @title Estimating Count Data Models with Social Interactions under Rational Expectations Using the NPL Method
#' @param formula a class object \link[stats]{formula}: a symbolic description of the model. The `formula` must be, for example, \code{y ~ x1 + x2 + gx1 + gx2}, where `y` is the endogenous vector, and `x1`, `x2`, `gx1`, and `gx2` are control variables, which may include contextual variables (i.e., averages among the peers). Peer averages can be computed using the function \code{\link{peer.avg}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets (e.g., schools), `Glist` can be a list of subnets, with the `m`-th element being an \eqn{n_m \times n_m} adjacency matrix, where \eqn{n_m} is the number of nodes in the `m`-th subnet. For heterogeneous peer effects (i.e., when `length(unique(group)) = h > 1`), the `m`-th element must be a list of \eqn{h^2} \eqn{n_m \times n_m} adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024, Section 2.1). For heterogeneous peer effects in the case of a single large network (a single school), `Glist` must be a one-item list (since there is one school). This item must be a list of \eqn{h^2} network specifications. The order in which the networks are specified is important and must match the order of the groups in `sort(unique(group))` (see argument `group` and examples).
#' @param group a vector indicating the individual groups. The default assumes a common group. For two groups, i.e., `length(unique(group)) = 2` (e.g., `A` and `B`), four types of peer effects are defined: peer effects of `A` on `A`, of `A` on `B`, of `B` on `A`, and of `B` on `B`. In this case, in the argument `Glist`, the networks must be defined in this order: `AA`, `AB`, `BA`, `BB`.
#' @param Rmax an integer indicating the theoretical upper bound of `y` (see model specification in detail).
#' @param Rbar an \eqn{L}-vector, where \eqn{L} is the number of groups. For large `Rmax`, the cost function is assumed to be semi-parametric (i.e., nonparametric from 0 to \eqn{\bar{R}} and quadratic beyond \eqn{\bar{R}}).
#' @param starting (optional) a starting value for \eqn{\theta = (\lambda, \Gamma', \delta')}, where \eqn{\lambda}, \eqn{\Gamma}, and \eqn{\delta} are the parameters to be estimated (see details).
#' @param Ey0 (optional) a starting value for \eqn{E(y)}.
#' @param optimizer specifies the optimization method, which can be one of: `fastlbfgs` (L-BFGS optimization method from the \pkg{RcppNumerical} package), `nlm` (from the function \link[stats]{nlm}), or `optim` (from the function \link[stats]{optim}). Arguments for these functions, such as `control` and `method`, can be set via the argument `opt.ctr`.
#' @param npl.ctr a list of controls for the NPL method (see details).
#' @param opt.ctr a list of arguments to be passed to `optim_lbfgs` from the \pkg{RcppNumerical} package, or to \link[stats]{nlm} or \link[stats]{optim} (the solver specified in `optimizer`), such as `maxit`, `eps_f`, `eps_g`, `control`, `method`, etc.
#' @param cov a Boolean indicating whether the covariance should be computed.
#' @param ubslambda a positive value indicating the upper bound of \eqn{\sum_{s = 1}^S \lambda_s > 0}.
#' @param data an optional data frame, list, or environment (or an object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model. If not found in `data`, the variables are taken from \code{environment(formula)}, typically the environment from which `cdnet` is called.
#' @return A list consisting of:
#'     \item{info}{a list containing general information about the model.}
#'     \item{estimate}{the NPL estimator.}
#'     \item{Ey}{\eqn{E(y)}, the expectation of \eqn{y}.}
#'     \item{GEy}{the average of \eqn{E(y)} across peers.}
#'     \item{cov}{a list that includes (if `cov == TRUE`): `parms`, the covariance matrix, and another list, `var.comp`, which contains `Sigma` (\eqn{\Sigma}) and `Omega` (\eqn{\Omega}), the matrices used to compute the covariance matrix.}
#'     \item{details}{step-by-step output returned by the optimizer.}
#' @description
#' `cdnet` estimates count data models with social interactions under rational expectations using the NPL algorithm (see Houndetoungan, 2024).
#' @details 
#' ## Model
#' The count variable \eqn{y_i} takes the value \eqn{r} with probability.
#' \deqn{P_{ir} = F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r}) - F(\sum_{s = 1}^S \lambda_s \bar{y}_i^{e,s}  + \mathbf{z}_i'\Gamma - a_{h(i),r + 1}).}
#' In this equation, \eqn{\mathbf{z}_i} is a vector of control variables; \eqn{F} is the distribution function of the standard normal distribution;
#' \eqn{\bar{y}_i^{e,s}} is the average of \eqn{E(y)} among peers using the `s`-th network definition;
#' \eqn{a_{h(i),r}} is the `r`-th cut-point in the cost group \eqn{h(i)}. \cr\cr
#' The following identification conditions have been introduced: \eqn{\sum_{s = 1}^S \lambda_s > 0}, \eqn{a_{h(i),0} = -\infty}, \eqn{a_{h(i),1} = 0}, and 
#' \eqn{a_{h(i),r} = \infty} for any \eqn{r \geq R_{\text{max}} + 1}. The last condition implies that \eqn{P_{ir} = 0} for any \eqn{r \geq R_{\text{max}} + 1}.
#' For any \eqn{r \geq 1}, the distance between two cut-points is \eqn{a_{h(i),r+1} - a_{h(i),r} =  \delta_{h(i),r} + \sum_{s = 1}^S \lambda_s}.
#' As the number of cut-points can be large, a quadratic cost function is considered for \eqn{r \geq \bar{R}_{h(i)}}, where \eqn{\bar{R} = (\bar{R}_{1}, ..., \bar{R}_{L})}.
#' With the semi-parametric cost function,
#' \eqn{a_{h(i),r + 1} - a_{h(i),r} = \bar{\delta}_{h(i)} + \sum_{s = 1}^S \lambda_s}.  \cr\cr
#' The model parameters are: \eqn{\lambda = (\lambda_1, ..., \lambda_S)'}, \eqn{\Gamma}, and \eqn{\delta = (\delta_1', ..., \delta_L')'}, 
#' where \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l}, \bar{\delta}_l)'} for \eqn{l = 1, ..., L}. 
#' The number of single parameters in \eqn{\delta_l} depends on  \eqn{R_{\text{max}}} and \eqn{\bar{R}_l}. The components \eqn{\delta_{l,2}, ..., \delta_{l,\bar{R}_l}} or/and 
#' \eqn{\bar{\delta}_l} must be removed in certain cases.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_l \geq 2}, then \eqn{\delta_l = (\delta_{l,2}, ..., \delta_{l,\bar{R}_l})'}.\cr
#' If \eqn{R_{\text{max}} = \bar{R}_l = 1} (binary models), then \eqn{\delta_l} must be empty.\cr
#' If \eqn{R_{\text{max}} > \bar{R}_l = 1}, then \eqn{\delta_l = \bar{\delta}_l}.

#' ## \code{npl.ctr}
#' The model parameters are estimated using the Nested Partial Likelihood (NPL) method. This approach 
#' begins with an initial guess for \eqn{\theta} and \eqn{E(y)} and iteratively refines them. 
#' The solution converges when the \eqn{\ell_1}-distance between two consecutive estimates of 
#' \eqn{\theta} and \eqn{E(y)} is smaller than a specified tolerance.
#' 
#' The argument \code{npl.ctr} must include the following parameters:
#' \describe{
#' \item{tol}{the tolerance level for the NPL algorithm (default is 1e-4).}
#' \item{maxit}{the maximum number of iterations allowed (default is 500).}
#' \item{print}{a boolean value indicating whether the estimates should be printed at each step.}
#' \item{S}{the number of simulations performed to compute the integral in the covariance using importance sampling.}
#' }
#' @references 
#' Houndetoungan, A. (2024). Count Data Models with Heterogeneous Peer Effects. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @seealso \code{\link{sart}}, \code{\link{sar}}, \code{\link{simcdnet}}.
#' @examples 
#' \donttest{
#' set.seed(123)
#' M      <- 5 # Number of sub-groups
#' nvec   <- round(runif(M, 100, 200))
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
#' # Two group:
#' group  <- 1*(X[,1] > 0.95)
#' 
#' # Networks
#' # length(group) = 2 and unique(sort(group)) = c(0, 1)
#' # The networks must be defined as to capture:
#' # peer effects of `0` on `0`, peer effects of `1` on `0`
#' # peer effects of `0` on `1`, and peer effects of `1` on `1`
#' G        <- list()
#' cums     <- c(0, cumsum(nvec))
#' for (m in 1:M) {
#'   tp     <- group[(cums[m] + 1):(cums[m + 1])]
#'   Am     <- A[[m]]
#'   G[[m]] <- norm.network(list(Am * ((1 - tp) %*% t(1 - tp)),
#'                               Am * ((1 - tp) %*% t(tp)),
#'                               Am * (tp %*% t(1 - tp)),
#'                               Am * (tp %*% t(tp))))
#' }
#' 
#' # Parameters
#' lambda <- c(0.2, 0.3, -0.15, 0.25) 
#' Gamma  <- c(4.5, 2.2, -0.9, 1.5, -1.2)
#' delta  <- rep(c(2.6, 1.47, 0.85, 0.7, 0.5), 2) 
#' 
#' # Data
#' data   <- data.frame(X, peer.avg(Anorm, cbind(x1 = X[,1], x2 =  X[,2])))
#' colnames(data) = c("x1", "x2", "gx1", "gx2")
#' 
#' ytmp   <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2),
#'                    lambda = lambda, Gamma = Gamma, delta = delta, group = group,
#'                    data = data)
#' y      <- ytmp$y
#' hist(y, breaks = max(y) + 1)
#' table(y)
#' 
#' # Estimation
#' est    <- cdnet(formula = y ~ x1 + x2 + gx1 + gx2, Glist = G, Rbar = rep(5, 2), group = group,
#'                 optimizer = "fastlbfgs", data = data,
#'                 opt.ctr = list(maxit = 5e3, eps_f = 1e-11, eps_g = 1e-11))
#' summary(est)
#' }
#' @importFrom stats quantile
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom stats runif
#' @export
cdnet    <- function(formula,
                     Glist, 
                     group,
                     Rmax,
                     Rbar,
                     starting  = list(lambda = NULL, Gamma = NULL, delta = NULL), 
                     Ey0       = NULL,
                     ubslambda = 1L,
                     optimizer = "fastlbfgs", 
                     npl.ctr   = list(), 
                     opt.ctr   = list(), 
                     cov       = TRUE,
                     data) {
  if(missing(Rmax)) Rmax <- Inf
  stopifnot(optimizer %in% c("fastlbfgs", "optim", "nlm"))
  stopifnot((Rmax >= 1) & (Rmax <= Inf))
  stopifnot(length(ubslambda) == 1)
  stopifnot(ubslambda > 0)
  if(Rmax == 1) Rbar <- 1
  lb_sl       <- 0
  ub_sl       <- ubslambda
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
  M        <- length(Glist)
  nvec     <- sapply(Glist, function(x_) nrow(x_[[1]]))
  sumn     <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  nCl      <- length(Glist[[1]])
  
  # group
  if(missing(group)){
    group <- rep(0, sumn)
  }
  if(is.null(group)){
    group <- rep(0, sumn)
  }
  
  if(length(group) != sumn) stop("length(group) != n")
  uCa      <- sort(unique(group))
  nCa      <- length(uCa)
  lCa      <- lapply(uCa, function(x_) which(group == x_) - 1)
  na       <- sapply(lCa, length)
  if(!missing(Rbar)){
    if(length(Rbar) == 1) Rbar = rep(Rbar, nCa)
    if(nCa != length(Rbar)) stop("The length of Rbar must equal the number of groups.")
  }
  if((nCa^2) != nCl) stop("The number of network specifications must match the number of groups.")
  
  # starting values
  thetat   <- NULL
  Deltat   <- NULL
  if(!inherits(starting, "list") & !is.null(starting)) {
    Deltat <- tail(starting, sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))
    Lmdt   <- starting[1:nCl]
    Gamt   <- starting[(nCl + 1):(length(starting) - sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar)))]
    thetat <- c(Lmdt, Gamt)
  }
  if(inherits(starting, "list")){
    if(!is.null(starting$lambda)){
      thetat <- c(starting$lambda, starting$Gamma)
      Deltat <- starting$delta
      if(Rmax == 1) Deltat <- numeric(0)
    }
  }
  if(!is.null(Deltat)){
    if(length(Deltat) != sum(ifelse(Rbar == Rmax, Rbar - 1, Rbar))) stop("The length of `delta` must match `Rbar` and `Rmax`.")
  }
  tmp      <- c(is.null(thetat), is.null(Deltat))
  if((all(tmp) != any(tmp))){
    stop("Starting parameters are defined, but not all required parameters are provided.")
  }
  
  # data
  f.t.data <- formula.to.data(formula, FALSE, Glist, M, igr, data, theta0 = 0) #because G are directly included in X
  formula  <- f.t.data$formula
  y        <- f.t.data$y
  X        <- f.t.data$X
  
  if(!is.null(thetat)){
    if((ncol(X) + nCl) != length(thetat)) stop("ncol(X) != length(Gamma).")
  }
  coln     <- c("lambda", colnames(X))
  if(nCl > 1) {
    coln   <- c(paste0(coln[1], ":", 1:nCl), coln[-1])
  }
  maxy     <- sapply(lCa, function(x_) max(y[x_ + 1]))
  if(any(maxy > Rmax)) stop("`Rmax` is lower than the maximum value of y.")
  K        <- ncol(X)
  
  # Ey 
  Eyt      <- runif(sumn, y*0.98, y*1.02)
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
    stop("`Rbar` is greater than the maximum value of y in at least one cost group.")
  }
  
  # check starting and computed them if not defined
  lmbd0       <- NULL  
  ndelta      <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  if (!is.null(thetat)) {
    if(length(thetat) != (K + nCl)) {
      stop("The length of `theta` is inappropriate.")
    }
    t0il      <- c(0, cumsum(rep(nCa, nCa)))
    if(any(sapply(1:nCa, function(x_) sum(thetat[(t0il[x_] + 1):t0il[x_ + 1]])) < 0)) stop("Negative peer effects are not supported in this version.")
    lmbd0     <- head(thetat, nCl)
    thetat    <- c(fcdlambdat(head(thetat, nCl), nCa, lb_sl, ub_sl), tail(thetat, K))
    if(any(Deltat <= 0)) stop("all(delta) > 0 is not true for starting values.")
  } else {
    Gy        <- lapply(1:nCl, function(x_2) unlist(lapply(1:M, function(x_1) Glist[[x_1]][[x_2]] %*% y[(igr[x_1, 1]:igr[x_1,2]) + 1])))
    Gy        <- as.matrix(as.data.frame(Gy)); colnames(Gy) <- head(coln, nCl)
    Xtmp      <- cbind(Gy, X)
    b         <- solve(crossprod(Xtmp), crossprod(Xtmp, y))
    b         <- b/sqrt(sum((y - Xtmp%*%b)^2)/sumn)
    b[1:nCl]  <- sapply(b[1:nCl], function(x_) min(max(x_, 0.01), 0.99)) 
    lmbd0     <- b[1:nCl]
    Deltat    <- rep(1, sum(ndelta))
    thetat    <- c(fcdlambdat(head(b, nCl), nCa, lb_sl, ub_sl), tail(b, K))
  }
  
  idelta      <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  lDelta      <- log(Deltat)
  thetat      <- c(thetat, lDelta)
  thetat0     <- thetat
  
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
  ctr      <- c(list(lb_sl = lb_sl, ub_sl = ub_sl, Gye = GEyt, X = X, lCa = lCa, nCa = nCa, 
                     K = K, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta, Rbar = Rbar, 
                     R = Rmax, y = y, maxy = maxy) , opt.ctr)
  
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
    # print(ctr$Gye)
    # print(ctr[[par]])
    REt         <- do.call(get(optimizer), ctr)
    thetat      <- REt[[par1]]
    llht        <- -REt[[like]]
    theta       <- c(fcdlambda(head(thetat, nCl), nCa, lb_sl, ub_sl), thetat[(nCl + 1):(K + nCl)], exp(tail(thetat, sum(ndelta))) + 1e-323)
    
    # compute y
    fL_NPL(ye = Eyt, Gye = GEyt, theta = theta, X = X, G = Glist, lCa = lCa, nCa = nCa, igroup = igr, 
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
    
    if((ninc.d > npl.incdit) | (!is.finite(llht)) | (llht == 1e250)) {
      cont          <- (t < (npl.maxit - 1))
      cat("** Non-convergence ** Redefining theta and computing a new E(y)\n")
      if(t > 1){
        thetat      <- steps[[t - 1]][[par1]]
        if(any(!is.finite(thetat))){
          thetat        <- thetat0
          thetat[1:nCl] <- 1e-7
          theta         <- c(thetat[1:(K + nCl)], exp(tail(thetat, sum(ndelta))) + 1e-323)
          thetat        <- c(fcdlambdat(head(theta, nCl), nCa, lb_sl, ub_sl), theta[(nCl + 1):(K + nCl)], log(tail(theta, sum(ndelta))))
        }
      } else {
        thetat          <- thetat0
        thetat[1:nCl]   <- 1e-7
        theta           <- c(thetat[1:(K + nCl)], exp(tail(thetat, sum(ndelta))) + 1e-323)
        thetat          <- c(fcdlambdat(head(theta, nCl), nCa, lb_sl, ub_sl), theta[(nCl + 1):(K + nCl)], log(tail(theta, sum(ndelta))))
      }
      Eyt           <- runif(sumn, y*0.98, y*1.02)
      GEyt          <- lapply(1:nCl, function(x_2) unlist(lapply(1:M, function(x_1) Glist[[x_1]][[x_2]] %*% Eyt[(igr[x_1, 1]:igr[x_1, 2]) + 1])))
      GEyt          <- as.matrix(as.data.frame(GEyt)); colnames(GEyt) <- head(coln, nCl)
      dist0         <- Inf
      # fnewye(ye = Eyt, Gye = GEyt, theta = theta, X = X, G = Glist, lCa = lCa, nCa = nCa, igroup = igr, 
      #        ngroup = M, K = K, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta, Rbar = Rbar, tol = npl.tol, 
      #        maxit = npl.maxit, R = Rmax)
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
    warning("NPL: The maximum number of iterations has been reached.")
  }
  # covariance and ME
  Gamma2   <- Gam
  colnme   <- coln
  ixWi     <- 0:(K - 1)
  if("(Intercept)" %in% coln){
    Gamma2 <- Gam[tail(coln, K) != "(Intercept)"]
    colnme <- colnme[coln != "(Intercept)"]
    ixWi   <- ixWi[tail(coln, K) != "(Intercept)"]
  }
  
  tmp      <- fcovCDI(theta = theta, Gamma2 = Gamma2, Gye = GEyt, X = X, ixWi = ixWi, G = Glist, lCa = lCa, nCa = nCa, 
                      igroup = igr, ngroup = M, K = K, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta,
                      Rbar = Rbar, R = Rmax, S = npl.S, ccov = cov)
  # Marginal effects
  imeff    <- tmp$imeff; colnames(imeff) <- colnme
  meff     <- list(ameff = apply(imeff, 2, mean), imeff = imeff)
  
  # Covariances
  covt     <- tmp$covt
  covm     <- tmp$covm
  if(cov){
    colnames(covt)  <- namtheta
    rownames(covt)  <- namtheta
    colnames(covm)  <- colnme
    rownames(covm)  <- colnme
  }
  
  AIC                  <- 2*length(theta) - 2*llht
  BIC                  <- length(theta)*log(sumn) - 2*llht
  
  INFO                 <- list("M"          = M,
                               "n"          = nvec,
                               "Kz"         = K,
                               "formula"    = formula,
                               "n.lambda"   = nCl,
                               "bslambda"   = c(lb_sl, ub_sl),
                               "Rbar"       = Rbar,
                               "Rmax"       = Rmax,
                               "group"      = group,
                               "log.like"   = llht, 
                               "npl.iter"   = t,
                               "AIC"        = AIC,
                               "BIC"        = BIC)
  
  out                  <- list("info"       = INFO,
                               "estimate"   = list(parms = theta, lambda =  lbda, Gamma = Gam, delta = Del),
                               "meff"       = meff,
                               "Ey"         = Eyt, 
                               "GEy"        = GEyt,
                               "cov"        = list(parms = covt, ameff = covm),
                               "details"    = steps)
  class(out)           <- "cdnet"
  out
}


#' @title Summary for the Estimation of Count Data Models with Social Interactions under Rational Expectations
#' @description Summary and print methods for the class `cdnet` as returned by the function \link{cdnet}.
#' @param object an object of class `cdnet`, output of the function \code{\link{cdnet}}.
#' @param x an object of class `summary.cdnet`, output of the function \code{\link{summary.cdnet}} or class `cdnet`, output of the function \code{\link{cdnet}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogeneous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param data an optional data frame, list, or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `summary.cdnet` is called.
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
    group       <- object$info$group
    M           <- object$info$M
    nvec        <- object$info$n
    sumn        <- sum(nvec)
    nCa         <- length(Rbar)
    theta       <- object$estimate$parms
    lambda      <- object$estimate$lambda
    Gamma       <- object$estimate$Gamma
    delta       <- object$estimate$delta
    ndelta      <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
    idelta      <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
    
    # Cost group
    uCa         <- sort(unique(group))
    lCa         <- lapply(uCa, function(x_) which(group == x_) - 1)
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
    if(M != length(Glist)) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
    if(any(nvec != sapply(Glist, function(x_) nrow(x_[[1]])))) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
    if(nCl != length(Glist[[1]])) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
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
    
    namtheta    <- c(coln, unlist(lapply(1:nCa, function(x_){
      if((Rbar[x_] == 1) & (Rmax == 1)) return(c())
      if(Rbar[x_] == Rmax) return(paste0("delta:", x_, ".", 2:(Rbar[x_])))
      if(Rbar[x_] > 1) return(c(paste0("delta:", x_, ".", 2:(Rbar[x_])), paste0("deltabar:", x_)))
      paste0("deltabar:", x_)
    })))
    
    Gamma2      <- Gamma
    ixWi        <- 0:(Kz - 1)
    if("(Intercept)" %in% coln){
      Gamma2    <- Gamma[tail(coln, Kz) != "(Intercept)"]
      ixWi      <- ixWi[tail(coln, Kz) != "(Intercept)"]
    }
    tmp         <- fcovCDI(theta = theta, Gamma2 = Gamma2, Gye = GEyt, X = X, ixWi = ixWi, G = Glist, lCa = lCa, nCa = nCa, 
                           igroup = igr, ngroup = M, K = Kz, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta,
                           Rbar = Rbar, R = Rmax, S = npl.S, ccov = TRUE)
    
    # Covariances
    covt        <- tmp$covt
    covm        <- tmp$covm
    colnames(covt)  <- namtheta
    rownames(covt)  <- namtheta
    colnames(covm)  <- names(out$meff$ameff)
    rownames(covm)  <- names(out$meff$ameff)
    
    out$cov         <- list(parms = covt, ameff = covm)
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
  Rmax                 <- x$info$Rmax
  csRbar               <- c(0, cumsum(Rbar - (Rbar == Rmax)))
  formula              <- x$info$formula
  Kz                   <- x$info$Kz
  AIC                  <- x$info$AIC
  BIC                  <- x$info$BIC
  nCa                  <- length(Rbar)
  nCl                  <- x$info$n.lambda
  
  coef                 <- c(x$estimate$lambda, x$estimate$Gamma)
  K                    <- length(coef)
  meff                 <- x$meff$ameff
  std                  <- sqrt(head(diag(x$cov$parms), K))
  std.meff             <- sqrt(diag(x$cov$ameff))
  delta                <- x$estimate$delta
  
  llh                  <- x$info$log.like
  
  tmp                  <- fcoefficients(coef, std)
  out_print            <- tmp$out_print
  out                  <- tmp$out
  out_print            <- c(list(out_print), x[-(1:7)], list(...))
  
  tmp.meff             <- fcoefficients(meff, std.meff)
  out_print.meff       <- tmp.meff$out_print
  out.meff             <- tmp.meff$out
  out_print.meff       <- c(list(out_print.meff), x[-(1:7)], list(...))
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Call:\n")
  print(formula)
  cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, sep = "", "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, sep = "", "\n")
  cat("Sample size              : ", sum(n), sep = "", "\n\n")
  
  cat("Coefficients:\n")
  do.call("print", out_print)
  
  cat("\nMarginal Effects:\n")
  do.call("print", out_print.meff)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  cat("Cost function -- Number of groups: ", nCa, "\n", sep = "")
  
  for(k in 1:nCa){
    cat("Group ", k, ", Rbar = ", Rbar[k], "\n", sep = "")
    if(Rbar[k] > 1){
      cat("delta:", delta[(csRbar[k] + 1):(csRbar[k + 1] - 1)], "\n", sep = " ")
      if(Rmax > Rbar[k]){
        cat("deltabar: ", delta[csRbar[k + 1]], "\n\n", sep = "")
      } else {
        cat("deltabar not defined: Rbar = Rmax\n\n", sep = "")
      }
    } else {
      if(Rmax > Rbar[k]){
        cat("delta not defined: Rbar = 1\n", sep = " ")
        cat("deltabar: ", delta[k], "\n\n", sep = "")
      } else{
        cat("delta and deltabar not defined: Rbar = Rmax = 1\n\n", sep = " ")
      }
    }
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

#' @title Counterfactual Analyses with Count Data Models and Social Interactions
#' @param object an object of class `summary.cdnet`, output of the function \code{\link{summary.cdnet}} or class `cdnet`, output of the function \code{\link{cdnet}}.
#' @param Glist adjacency matrix. For networks consisting of multiple subnets, `Glist` can be a list of subnets with the `m`-th element being an `ns*ns` adjacency matrix, where `ns` is the number of nodes in the `m`-th subnet.
#' For heterogeneous peer effects (e.g., boy-boy, boy-girl friendship effects), the `m`-th element can be a list of many `ns*ns` adjacency matrices corresponding to the different network specifications (see Houndetoungan, 2024).
#' For heterogeneous peer effects in the case of a single large network, `Glist` must be a one-item list. This item must be a list of many specifications of large networks.
#' @param data an optional data frame, list, or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which `summary.cdnet` is called.
#' @param S number of simulations to be used to compute integral in the covariance by important sampling.
#' @param tol the tolerance value used in the Fixed Point Iteration Method to compute the expectancy of `y`. The process stops if the \eqn{\ell_1}-distance between two consecutive \eqn{E(y)} is less than `tol`.
#' @param maxit the maximal number of iterations in the Fixed Point Iteration Method.
#' @param group the vector indicating the individual groups (see function \code{\link{cdnet}}). If missing, the former group saved in `object` will be used.
#' @description
#' `simcdpar` computes the average expected outcomes for count data models with social interactions and standard errors using the Delta method. 
#' This function can be used to examine the effects of changes in the network or in the control variables.
#' @seealso \code{\link{simcdnet}}
#' @return A list consisting of:
#'     \item{Ey}{\eqn{E(y)}, the expectation of y.}
#'     \item{GEy}{the average of \eqn{E(y)} friends.}
#'     \item{aEy}{the sampling mean of \eqn{E(y)}.}
#'     \item{se.aEy}{the standard error of the sampling mean of \eqn{E(y)}.}
#' @export
simcdEy <- function(object,
                    Glist,
                    data,
                    group,
                    tol          = 1e-10,
                    maxit        = 500,
                    S            = 1e3){
  stopifnot(inherits(object, c("cdnet", "summary.cdnet")))
  covparms    <- object$cov$parms
  if(is.null(covparms)) stop("`object` does not include a covariance matrix")
  env.formula <- environment(object$info$formula)
  Rbar        <- object$info$Rbar
  Rmax        <- object$info$Rmax
  Kz          <- object$info$Kz
  nCl         <- object$info$n.lambda
  if(missing(group)){
    group     <- object$info$group
  }
  uCa         <- sort(unique(group))
  nCa         <- length(uCa)
  lCa         <- lapply(uCa, function(x_) which(group == x_) - 1)
  na          <- sapply(lCa, length)
  if(nCa != length(Rbar)) stop("length(Rbar) is not equal to the number of groups.")
  if((nCa^2) != nCl) stop("The number of network specifications does not match the number of groups.")
  M           <- object$info$M
  nvec        <- object$info$n
  sumn        <- sum(nvec)
  theta       <- object$estimate$parms
  lambda      <- object$estimate$lambda
  Gamma       <- object$estimate$Gamma
  delta       <- object$estimate$delta
  ndelta      <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  idelta      <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  delta       <- c(fdelta(deltat = delta, lambda = lambda, idelta = idelta, ndelta = ndelta, nCa = nCa))
  

  
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
  if(M != length(Glist)) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
  if(any(nvec != sapply(Glist, function(x_) nrow(x_[[1]])))) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
  if(nCl != length(Glist[[1]])) stop("`Glist` structure does not match. Perhaps another `Glist` structure is used in `cdnet`.")
  igr         <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  # Data
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
  
  # E(y) and G*E(y)
  Ey          <- rep(0, sumn)
  GEy         <- matrix(0, sumn, nCl)
  xb          <- c(X %*% Gamma)
  t           <- fye(ye = Ey, Gye = GEy, G = Glist, lCa = lCa, nCa = nCa, igroup = igr, ngroup = M, 
                   psi = xb, lambda = lambda, delta = delta, idelta = idelta, n = na, sumn = sumn,
                   Rbar = Rbar, R = Rmax, tol = tol, maxit = maxit)
  if(t >= maxit) warning("Point-fixed: The maximum number of iterations has been reached.")
  
  dEy         <- fcddEy(theta = theta, Gye = GEy, X = X, psi = xb, G = Glist, lCa = lCa, nCa = nCa, igroup = igr, 
                        ngroup = M, K = Kz, n = na, sumn = sumn, idelta = idelta, ndelta = ndelta, Rbar = Rbar,
                        R = Rmax, S = npl.S)
  
  if(nCl == 1){
    GEy       <- c(GEy)
  }
  Ey          <- c(Ey)
  aEy         <- mean(Ey)
  
  # if(fixed.lambda){
  #   dEy       <- dEy[, -(1:nCl), drop = FALSE]
  #   covparms  <- covparms[-(1:nCl), -(1:nCl), drop = FALSE]
  # }
  adEy        <- apply(dEy, 2, mean)
  se.aEy      <- sqrt(c(t(adEy) %*% covparms %*% adEy))
  
  out         <- list(Ey     = Ey,
                      GEy    = GEy,
                      aEy    = aEy,
                      se.aEy = se.aEy)
  class(out)  <- "simcdEy"
  out
}

#' @title Printing the Average Expected Outcomes for Count Data Models with Social Interactions
#' @description Summary and print methods for the class `simcdEy` as returned by the function \link{simcdEy}.
#' @param object an object of class `simcdEy`, output of the function \code{\link{simcdEy}}.
#' @param x an object of class `summary.simcdEy`, output of the function \code{\link{summary.simcdEy}} or class `simcdEy`, output of the function \code{\link{simcdEy}}.
#' @param ... further arguments passed to or from other methods.
#' @return A list of the same objects in `object`.
#' @export
print.simcdEy <- function(x, ...){
  stopifnot(class(x) == "simcdEy")
  cat("Average expected outcome\n\n")
  tp           <- cbind("Estimate" = x$aEy, "Std. Error" = x$se.aEy)
  rownames(tp) <- ""
  print(tp, ...)
  invisible(x)
}

#' @rdname print.simcdEy
#' @export
"summary.simcdEy"  <- function(object, ...) return(object)

#' @rdname print.simcdEy
#' @export
"print.summary.simcdEy"  <- function(x, ...) print.simcdEy(x, ...)
