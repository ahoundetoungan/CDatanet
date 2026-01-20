#' @title Marginal Effects for Count Data Models and Tobit Models with Social Interactions
#' @description
#' `meffects` computes marginal effects for count data and Tobit models with social interactions.
#' It is a generic function which means that new printing methods can be easily added for new classes.
#' @param model an object of class `cdnet` (`summary.cdnet`) or `sart` (`summary.sart`), output of the function 
#' \code{\link{cdnet}} or \code{\link{sart}}, respectively.
#' @param ... Additional arguments passed to methods.
#' @param cont.var A character vector of continuous variable names for which the marginal effects should be computed.
#' @param bin.var A character vector of binary variable names for which the marginal effects should be computed.
#' @param Glist The network matrix used to obtain `model`. Typically, this is the `Glist` argument supplied to 
#'   the function \code{\link{cdnet}} or \code{\link{sart}}.
#' @param Glist.contextual The network matrix used to compute contextual variables, if any are specified in the `type.var` argument. 
#'   For networks consisting of multiple subnets, `Glist` can be a list of subnets, where the `m`-th element is an 
#'   `ns*ns` adjacency matrix, with `ns` denoting the number of nodes in the `m`-th subnet.
#' @param type.var A list indicating "own" and contextual variables that appear in the `cont.var` and `bin.var` arguments. 
#'   The list contains pairs of variable names, with the first element being the "own" variable and the second being the 
#'   contextual variable. When a variable has no associated contextual variable, only the variable name is included. 
#'   For example, `type.var = list(c("x1", "gx1"), c("x2", "gx2"), "x3")` means that `gx1` is the contextual variable for `x1`, 
#'   `gx2` is the contextual variable for `x2`, and `x3` has no contextual variable. This information is used to compute the 
#'   indirect and total marginal effects for `x1`, `x2`, and `x3`.
#' @param ncores Number of CPU cores (threads) used to run the bootstrap process in parallel.
#' @param data An optional data frame, list, or environment (or object coercible by \code{\link[base]{as.data.frame}} 
#'   to a data frame) containing the variables in the model. If not found in `data`, the variables are taken from 
#'   \code{environment(model)}, typically the environment from which `meffects` is called.
#' @param tol The tolerance value used in the fixed-point iteration method to compute `y`. The process stops if the 
#'   \eqn{\ell_1}-distance between two consecutive values of `y` is less than `tol`.
#' @param maxit The maximum number of iterations in the fixed-point iteration method.
#' @param boot The number of bootstrap simulations to compute standard errors and confidence intervals.
#' @param progress A logical value indicating whether the progress of the bootstrap simulations should be printed to the console.
#' @return A list containing:
#' \describe{
#'   \item{\code{info}}{General information about the model.}
#'   \item{\code{estimate}}{The Maximum Likelihood (ML) estimates of the parameters.}
#'   \item{\code{Ey}}{\eqn{E(y)}, the expected values of the endogenous variable.}
#'   \item{\code{GEy}}{The average of \eqn{E(y)} among peers.}
#'   \item{\code{cov}}{A list containing covariance matrices (if \code{cov = TRUE}).}
#'   \item{\code{details}}{Additional outputs returned by the optimizer.}
#'   \item{\code{meffects}}{A list containing the marginal effects.}
#' }
#' 
#' @examples
#' \donttest{
#' #' set.seed(123)
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
#' 
#' meffects(est, Glist = G, data = data, cont.var = c("x1", "x2", "gx1", "gx2"),
#'          type.var = list(c("x1", "gx1"), c("x2", "gx2")), Glist.contextual = Anorm,
#'          boot = 100, ncores = 2)
#' }
#' @export
meffects <- function(model, ...) {
  UseMethod("meffects")
}

#' @rdname meffects
#' @export
meffects.cdnet <- function(model, 
                           Glist, 
                           cont.var, 
                           bin.var,
                           type.var,
                           Glist.contextual,
                           data,
                           tol        = 1e-10,
                           maxit      = 500,
                           boot       = 1000,
                           progress   = TRUE, 
                           ncores     = 1,
                           ...) {
  stopifnot(class(model) == "cdnet")
  stopifnot(ncores >= 1)
  
  if (is.null(model$cov) & (boot > 0)) {
    model    <- summary.cdnet(model, Glist = Glist, data = data)
  } else if (is.null(model$cov)) {
    model$cov <- matrix(0)
  }
  
  # cont.var and bin.var
  if (missing(cont.var)) {
    cont.var <- NULL
  } 
  
  if (missing(bin.var)) {
    bin.var <- NULL
  } 
  
  meff.var   <- c(cont.var, bin.var)
  if (any(duplicated(meff.var))) {
    stop("At least one variable is declared as both continuous and binary.")
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
  nvec     <- model$info$n
  M        <- model$info$M
  Rbar     <- model$info$Rbar
  Rmax     <- model$info$Rmax
  sumn     <- sum(nvec)
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  group    <- model$info$group
  uCa      <- sort(unique(group))
  nCa      <- length(uCa)
  nCl      <- nCa^2
  lCa      <- lapply(uCa, function(x_) which(group == x_) - 1)
  na       <- sapply(lCa, length)
  ndelta   <- ifelse(Rbar == Rmax, Rbar - 1, Rbar)
  idelta   <- matrix(c(0, cumsum(ndelta)[-length(ndelta)], cumsum(ndelta) - 1), ncol = 2); idelta[ndelta == 0,] <- NA
  
  # data
  formula   <- model$info$formula
  f.t.data  <- formula.to.data(formula = formula, contextual = FALSE, Glist = Glist, M = M, igr = igr, 
                               data = data, type = "model", theta0  = 0)
  X         <- f.t.data$X
  K         <- ncol(X)
  
  # X columns
  cnames   <- colnames(X)
  if (!all(meff.var %in% cnames)) {
    stop("Variables indicated for the marginal effects are not included in the model.")
  }
  
  # indices
  tp        <- fmeffect_aux(type.var, meff.var, cont.var, bin.var, cnames)
  conti     <- tp$conti
  dis0      <- tp$dis0
  dis1      <- tp$dis1
  idmarg    <- tp$idmarg
  indexX    <- tp$indexX
  hasCont   <- tp$hasCont
  indexGX   <- tp$indexGX
  idinmarg  <- tp$idinmarg
  
  # parameter
  lb_sl      <- model$info$bslambda[1]
  ub_sl      <- model$info$bslambda[2]
  if (is.null(lb_sl)) {
    lb_sl    <- 0
  }
  if (is.null(ub_sl)) {
    ub_sl    <- 1
  }
  thetat     <- c(fcdlambdat(model$estimate$lambda, nCa, lb_sl, ub_sl),
                  model$estimate$Gamma,
                  log(model$estimate$delta))
  # Other variables
  Gye        <- model$GEy
  if (any(hasCont == 1)) {
    if (missing(Glist.contextual)) {
      stop("Glist.contextual is missing without a default value.")
    }
    if (!is.list(Glist)) {
      Glist.contextual  <- list(Glist.contextual)
    }
  } else {
    Glist.contextual    <- rep(list(matrix(0)), M)
  }
  # print(conti)
  # print(dis0)
  # print(dis1)
  # print(idmarg)
  # print(indexX)
  # print(hasCont)
  # print(indexGX)
  # print(idinmarg)
  nparms <- length(thetat)
  ME     <- cdmeffects(thetat, Gye = as.matrix(Gye), X = X, conti = conti, dis0 = dis0, dis1 = dis1,
                       indexmarg = idmarg, indexX = indexX, hasCont = hasCont, indexGX = indexGX,
                       indexinmarg = idinmarg, G = Glist, Gcont = Glist.contextual, lCa = lCa, nCa = nCa,
                       igroup = igr, ngroup = M, idelta = idelta, ndelta = ndelta, sumn = sumn, Rbar = Rbar,
                       lb_sl = lb_sl, ub_sl = ub_sl, n = na, R = Rmax, tol = tol, maxit = maxit,
                       covparm = model$cov, simNorm = matrix(rnorm(nparms*boot), nparms), boot = boot,
                       print = progress, nthreads = ncores)
  # Direct effects
  meff     <- list(direct = ME$direct[1,])
  cov.meff <- list(direct = cov(ME$direct[ifelse(boot == 0, 1, -1),, drop = FALSE]))
  qt.meff  <- list(direct = apply(ME$direct[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
    quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
  }))
  names(meff$direct) <- colnames(cov.meff$direct) <-
    rownames(cov.meff$direct) <- colnames(qt.meff$direct) <- 
    c(names(model$estimate$lambda), cnames[idmarg + 1])
  if (length(conti) > 0) {
    # Indirect effects and total
    meff$indirect     <- ME$indirect[1,]
    meff$total        <- ME$total[1,]
    cov.meff$indirect <- cov(ME$indirect[ifelse(boot == 0, 1, -1),, drop = FALSE])
    cov.meff$total    <- cov(ME$total[ifelse(boot == 0, 1, -1),, drop = FALSE])
    qt.meff$indirect  <- apply(ME$indirect[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
      quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
    })
    qt.meff$total     <- apply(ME$total[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
      quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
    })
    
    names(meff$indirect) <- names(meff$total) <- 
      colnames(cov.meff$indirect) <- colnames(cov.meff$total) <-
      rownames(cov.meff$indirect) <- rownames(cov.meff$total) <-
      colnames(qt.meff$indirect) <- colnames(qt.meff$total) <- 
      cnames[indexX + 1]
  }
  model$meffects <- list(estimate = meff, cov = cov.meff, pctl = qt.meff)
  model$info$boot<- boot
  class(model)   <- "meffects.cdnet"
  model
}


#' @rdname meffects
#' @export
meffects.summary.cdnet <- function(model, 
                                   Glist, 
                                   cont.var, 
                                   bin.var,
                                   type.var,
                                   Glist.contextual, 
                                   data,
                                   tol        = 1e-10,
                                   maxit      = 500,
                                   boot       = 1000,
                                   progress   = TRUE, 
                                   ncores     = 1,
                                   ...) {
  stopifnot(class(model) == "summary.cdnet")
  class(model) <- "cdnet"
  meffects.cdnet(model = model, Glist = Glist, cont.var = cont.var, bin.var = bin.var,
                 type.var = type.var, Glist.contextual = Glist.contextual, data = data,
                 tol = tol, maxit = maxit, boot = boot, progress = progress, 
                 ncores = ncores, ...)
}

#' @export
"print.meffects.cdnet"  <- function(x, ...) {
  stopifnot(class(x) == "meffects.cdnet")
  
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
  
  llh                  <- x$info$log.like
  
  tmp                  <- fcoefficients(x$meffects$estimate$direct, sqrt(diag(x$meffects$cov$direct)))
  out_print_d          <- tmp$out_print
  out_d                <- tmp$out
  out_print_d          <- c(list(out_print_d), x[-(1:7)], list(...))
  
  out_print_i          <- NULL
  out_i                <- NULL
  out_print_t          <- NULL
  out_t                <- NULL
  if (!is.null(x$meffects$estimate$indirect)) {
    tmp                <- fcoefficients(x$meffects$estimate$indirect, sqrt(diag(x$meffects$cov$indirect)))
    out_print_i        <- tmp$out_print
    out_i              <- tmp$out
    out_print_i        <- c(list(out_print_i), x[-(1:7)], list(...))
    
    tmp                <- fcoefficients(x$meffects$estimate$total, sqrt(diag(x$meffects$cov$total)))
    out_print_t        <- tmp$out_print
    out_t              <- tmp$out
    out_print_t        <- c(list(out_print_t), x[-(1:7)], list(...))
  }
  
  cat("Count data Model with Social Interactions\n\n")
  cat("Call:\n")
  print(formula)
  cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, sep = "", "\n\n")
  cat("Network:\n")
  cat("Number of groups         : ", M, sep = "", "\n")
  cat("Sample size              : ", sum(n), sep = "", "\n\n")
  
  cat("Direct Marginal Effects:\n")
  do.call("print", out_print_d)
  
  if (!is.null(x$meffects$estimate$indirect)) {
    cat("\nIndirect Marginal Effects:\n")
    do.call("print", out_print_i)
    
    cat("\nTotal Marginal Effects:\n")
    do.call("print", out_print_t)
  }
  
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  
  cat("log pseudo-likelihood: ", llh, sep = "", "\n")
  cat("AIC: ", AIC, " -- BIC: ", BIC, sep = "", "\n")
  x$meffects$table <- list(direct = out_d, indirect = out_i, total = out_t)
  invisible(x)
}



########################################## Sart model
#' @rdname meffects
#' @export
meffects.sart <- function(model, 
                          Glist, 
                          cont.var, 
                          bin.var,
                          type.var,
                          Glist.contextual, 
                          data,
                          tol        = 1e-10,
                          maxit      = 500,
                          boot       = 1000,
                          progress   = TRUE, 
                          ncores     = 1,
                          ...) {
  stopifnot(class(model) == "sart")
  iteration <- model$info$npl.iter
  cinfo     <- is.null(iteration)
  
  # cont.var and bin.var
  if (missing(cont.var)) {
    cont.var <- NULL
  } else if (cinfo) {
    stop("Marginal effects are not provided for the complete information game in this version. 
         You can use simulations to obtain them instead. 
         Feel free to contribute to the package by adding marginal effects for the complete information model.")
  }
  
  if (missing(bin.var)) {
    bin.var <- NULL
  } else if (cinfo) {
    stop("Marginal effects are not provided for the complete information game in this version. 
         You can use simulations to obtain them instead. 
         Feel free to contribute to the package by adding marginal effects for the complete information model.")
  }
  meff.var   <- c(cont.var, bin.var)
  if (any(duplicated(meff.var))) {
    stop("At least one variable is declared as both continuous and binary.")
  }
  
  
  if (is.null(model$cov) & (boot > 0)) {
    model     <- summary.sart(model, Glist = Glist, data = data)
  } else if (is.null(model$cov)) {
    model$cov <- matrix(0)
  }
  
  # cont.var and bin.var
  if (missing(cont.var)) {
    cont.var <- NULL
  } 
  
  if (missing(bin.var)) {
    bin.var <- NULL
  } 
  
  meff.var   <- c(cont.var, bin.var)
  if (any(duplicated(meff.var))) {
    stop("At least one variable is declared as both continuous and binary.")
  }
  
  # Network
  if (!is.list(Glist)) {
    Glist    <- list(Glist)
  }
  
  
  # Sizes
  nvec     <- unlist(lapply(Glist, nrow))
  M        <- model$info$M
  sumn     <- model$info$n
  igr      <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2)
  
  # data
  formula   <- model$info$formula
  f.t.data  <- formula.to.data(formula = formula, contextual = FALSE, Glist = Glist, M = M, igr = igr, 
                               data = data, type = "model", theta0  = 0)
  X         <- f.t.data$X
  K         <- ncol(X)
  
  # X columns
  cnames   <- colnames(X)
  if (!all(meff.var %in% cnames)) {
    stop("Variables indicated for the marginal effects are not included in the model.")
  }
  
  # indices
  tp        <- fmeffect_aux(type.var, meff.var, cont.var, bin.var, cnames)
  conti     <- tp$conti
  dis0      <- tp$dis0
  dis1      <- tp$dis1
  idmarg    <- tp$idmarg
  indexX    <- tp$indexX
  hasCont   <- tp$hasCont
  indexGX   <- tp$indexGX
  idinmarg  <- tp$idinmarg
  
  # parameter
  thetat     <- c(log(model$estimate[1]) - log(1 - model$estimate[1]),
                  model$estimate[2:(K + 1)],
                  log(model$estimate[K + 2]))
  
  # Other variables
  Gye        <- model$GEy
  if (any(hasCont == 1)) {
    if (missing(Glist.contextual)) {
      Glist.contextual  <- Glist
    }
    if (!is.list(Glist)) {
      Glist.contextual  <- list(Glist.contextual)
    }
  } else {
    Glist.contextual    <- rep(list(matrix(0)), M)
  }
  
  nparms <- length(thetat)
  ME     <- SImeffects(thetat, Gye = Gye, X = X, conti = conti, dis0 = dis0, dis1 = dis1,
                       indexmarg = idmarg, indexX = indexX, hasCont = hasCont, indexGX = indexGX,
                       indexinmarg = idinmarg, G = Glist, Gcont = Glist.contextual, igroup = igr, 
                       ngroup = M, sumn = sumn, tol = tol, maxit = maxit, covparm = model$cov, 
                       simNorm = matrix(rnorm(nparms*boot), nparms), boot = boot, print = progress,
                       nthreads = ncores)
  # Direct effects
  meff     <- list(direct = ME$direct[1,])
  cov.meff <- list(direct = cov(ME$direct[ifelse(boot == 0, 1, -1),, drop = FALSE]))
  qt.meff  <- list(direct = apply(ME$direct[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
    quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
  }))
  names(meff$direct) <- colnames(cov.meff$direct) <-
    rownames(cov.meff$direct) <- colnames(qt.meff$direct) <- 
    c(names(model$estimate)[1], cnames[idmarg + 1])
  if (length(conti) > 0) {
    # Indirect effects and total
    meff$indirect     <- ME$indirect[1,]
    meff$total        <- ME$total[1,]
    cov.meff$indirect <- cov(ME$indirect[ifelse(boot == 0, 1, -1),, drop = FALSE])
    cov.meff$total    <- cov(ME$total[ifelse(boot == 0, 1, -1),, drop = FALSE])
    qt.meff$indirect  <- apply(ME$indirect[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
      quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
    })
    qt.meff$total     <- apply(ME$total[ifelse(boot == 0, 1, -1),, drop = FALSE], 2, function(x) {
      quantile(x, probs = c(0.005, 0.025, 0.5, 0.95, 0.975, 0.995))
    })
    
    names(meff$indirect) <- names(meff$total) <- 
      colnames(cov.meff$indirect) <- colnames(cov.meff$total) <-
      rownames(cov.meff$indirect) <- rownames(cov.meff$total) <-
      colnames(qt.meff$indirect) <- colnames(qt.meff$total) <- 
      cnames[indexX + 1]
  }
  model$meffects <- list(estimate = meff, cov = cov.meff, pctl = qt.meff)
  model$info$boot<- boot
  class(model)   <- "meffects.sart"
  model
}

#' @rdname meffects
#' @export
meffects.summary.sart <- function(model, 
                                  Glist, 
                                  cont.var, 
                                  bin.var,
                                  type.var,
                                  Glist.contextual, 
                                  data,
                                  tol        = 1e-10,
                                  maxit      = 500,
                                  boot       = 1000,
                                  progress   = TRUE, 
                                  ncores     = 1,
                                  ...) {
  stopifnot(class(model) == "summary.sart")
  class(model) <- "sart"
  meffects.sart(model = model, Glist = Glist, cont.var = cont.var, bin.var = bin.var,
                type.var = type.var, Glist.contextual = Glist.contextual, data = data,
                tol = tol, maxit = maxit, boot = boot, progress = progress, 
                ncores = ncores, ...)
}

#' @export
"print.meffects.sart"  <- function(x, ...) {
  stopifnot(class(x) == "meffects.sart")
  
  M                    <- x$info$M
  n                    <- x$info$n
  estimate             <- x$estimate
  iteration            <- x$info$npl.iter
  RE                   <- !is.null(iteration)
  formula              <- x$info$formula
  K                    <- length(estimate)
  coef                 <- estimate[-K]
  std                  <- sqrt(diag(x$cov)[-K])
  sigma                <- estimate[K]
  llh                  <- x$info$log.like
  censored             <- x$info$censured
  uncensored           <- x$info$uncensured
  
  tmp                  <- fcoefficients(x$meffects$estimate$direct, sqrt(diag(x$meffects$cov$direct)))
  out_print_d          <- tmp$out_print
  out_d                <- tmp$out
  out_print_d          <- c(list(out_print_d), x[-(1:7)], list(...))
  
  out_print_i          <- NULL
  out_i                <- NULL
  out_print_t          <- NULL
  out_t                <- NULL
  if (!is.null(x$meffects$estimate$indirect)) {
    tmp                <- fcoefficients(x$meffects$estimate$indirect, sqrt(diag(x$meffects$cov$indirect)))
    out_print_i        <- tmp$out_print
    out_i              <- tmp$out
    out_print_i        <- c(list(out_print_i), x[-(1:7)], list(...))
    
    tmp                <- fcoefficients(x$meffects$estimate$total, sqrt(diag(x$meffects$cov$total)))
    out_print_t        <- tmp$out_print
    out_t              <- tmp$out
    out_print_t        <- c(list(out_print_t), x[-(1:7)], list(...))
  }
  
  nfr                  <- x$info$nlinks
  cat("sart Model", ifelse(RE, "with Rational Expectation", ""), "\n\n")
  cat("Call:\n")
  print(formula)
  if(RE){
    cat("\nMethod: Nested pseudo-likelihood (NPL) \nIteration: ", iteration, "\n\n")
  } else{
    cat("\nMethod: Maximum Likelihood (ML)", "\n\n")
  }
  
  cat("Network:\n")
  cat("Number of groups         : ", M, "\n")
  cat("Sample size              : ", n, "\n")
  cat("Average number of friends: ", sum(nfr)/n, "\n\n")
  
  cat("No. censored             : ", paste0(censored, "(", round(censored/n*100, 0), "%)"), "\n")
  cat("No. uncensored           : ", paste0(uncensored, "(", round(uncensored/n*100, 0), "%)"), "\n\n")
  
  cat("Direct Marginal Effects:\n")
  do.call("print", out_print_d)
  
  if (!is.null(x$meffects$estimate$indirect)) {
    cat("\nIndirect Marginal Effects:\n")
    do.call("print", out_print_i)
    
    cat("\nTotal Marginal Effects:\n")
    do.call("print", out_print_t)
  }
  
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat("sigma: ", sigma, "\n")
  cat("log likelihood: ", llh, "\n")
  x$meffects$table <- list(direct = out_d, indirect = out_i, total = out_t)
  invisible(x)
}



##############################################
fmeffect_aux <- function(type.var, meff.var, cont.var, bin.var, cnames) {
  # type.var
  if (missing(type.var)) {
    type.var <- as.list(unlist(meff.var))
  } else {
    type.var <- type.var[sapply(type.var, function(xx) !is.null(xx))]
  }
  if (is.null(type.var)) {
    type.var <- list()
  }
  if (any(duplicated(unlist(type.var)))) {
    stop("Variables in `type.var` are duplicated.")
  }
  if (any(!(meff.var %in% unlist(type.var))) | any(!(unlist(type.var) %in% meff.var))) {
    stop("Variables in `type.var` are not found in `cont.var` or `bin.var`, or vice versa.")
  }
  if (!is.list(type.var)) {
    stop("`type.var` is not a list.")
  }
  if (!is.null(meff.var)) {
    if(any(!(sapply(type.var, length) %in% 1:2))) {
      stop("Elements in `type.var` should each have a length of one or two.")
    }
  }
  
  # indices
  idmcon     <- numeric()
  idmdis     <- numeric()
  ncont      <- length(cont.var)
  nbin       <- length(bin.var)
  if (ncont > 0){
    idmcon   <- sapply(cont.var, function(xx) which(xx == cnames)) - 1
  }
  if (nbin > 0){
    idmdis   <- sapply(bin.var, function(xx) which(xx == cnames)) - 1
  }
  idmarg     <- c(idmcon, idmdis)
  conti      <- c(rep(1, ncont), rep(0, nbin))
  
  ntypeVar   <- length(type.var)
  indexX     <- numeric()
  hasCont    <- numeric()
  indexGX    <- rep(0, ntypeVar)
  if (ntypeVar) {
    indexX   <- sapply(sapply(type.var, function(xx) xx[1]),  function(zz) which(zz == cnames)) - 1
    hasCont  <- sapply(type.var, function(xx) ifelse(length(xx) == 2, 1, 0))
    indexGX  <- sapply(1:ntypeVar, function(xx) {
      ifelse(hasCont[xx] == 1, which(type.var[[xx]][2] == cnames) - 1, 0)})
  }
  
  # Sorting
  Order      <- order(idmarg)
  idmarg     <- idmarg[Order]
  conti      <- conti[Order]
  
  Order      <- order(indexX)
  indexX     <- indexX[Order]
  hasCont    <- hasCont[Order]
  indexGX    <- indexGX[Order]
  idinmarg   <- numeric() 
  if (length(indexX) > 0) {
    idinmarg <- sapply(indexX, function(xx)  which(xx == idmarg)) - 1
  }
  
  K          <- length(cnames)
  dis0       <- rep(0, K)
  dis1       <- rep(1, K)
  list(conti = conti, dis0 = dis0, dis1 = dis1, idmarg = idmarg, 
       indexX = indexX, hasCont = hasCont, indexGX = indexGX,
       idinmarg = idinmarg)
}
