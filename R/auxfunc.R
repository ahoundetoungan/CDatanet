#' @title Simulating network data
#' @description `simnetwork` simulates adjacency matrices.
#' @param dnetwork is a list of sub-network matrices, where the (i, j)-th position of the m-th matrix is the probability that i be connected to j, with i and j individuals from the m-th network.
#' @param normalise boolean takes `TRUE` if the returned matrices should be row-normalized and `FALSE` otherwise.
#' @return list of (row-normalized) adjacency matrices.
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N         <- c(250, 370, 120)  
#' ## distribution
#' dnetwork  <- lapply(N, function(x) matrix(runif(x^2), x))
#' ## network
#' G         <- simnetwork(dnetwork)
#' @export
simnetwork <- function(dnetwork, normalise = FALSE) {
  trsf          <- FALSE
  if (!is.list(dnetwork)) {
    if (is.matrix(dnetwork)) {
      dnetwork  <- list(dnetwork)
      trsf      <- TRUE
    } else {
      stop("dnetwork is neither a matrix nor a list of matrices")
    }
  }
  
  M        <- length(dnetwork)
  N        <- unlist(lapply(dnetwork, nrow))
  out      <- NULL
  if (normalise) {
    out    <- simGnorm(dnetwork = dnetwork, N = N, M = M)
  } else {
    out    <- simG(dnetwork = dnetwork, N = N, M = M)
  }
  if (trsf) {
    out    <- out[[1]]
  }
  out
}


#' @rdname vec.to.mat
#' @export
norm.network <- function(W) {
  trsf     <- FALSE
  if (!is.list(W)) {
    if (is.matrix(W)) {
      W    <- list(W)
      trsf <- TRUE
    } else {
      stop("W is neither a matrix nor a list of matrices")
    }
  }
  stopifnot(inherits(W, "list"))
  
  M        <- length(W)
  out      <- fGnormalise(W, M)
  if (trsf) {
    out    <- out[[1]]
  }
  out
}


#' @title Creating objects for network models
#' @description  `vec.to.mat` creates a list of square matrices from a given vector.
#' The elements of the generated matrices are taken from the vector and placed column-wise (ie. the first column is filled up before filling the second column) and from the first matrix of the list to the last matrix of the list. 
#' The diagonal of the generated matrices are zeros.
#' `mat.to.vec` creates a vector from a given list of square matrices .
#' The elements of the generated vector are taken from column-wise and from the first matrix of the list to the last matrix of the list, while dropping the diagonal entry.
#' `norm.network` row-normalizes matrices in a given list.
#' @param u numeric vector to convert.
#' @param W matrix or list of matrices to convert. 
#' @param N vector of sub-network sizes  such that `length(u) == sum(N*(N - 1))`.
#' @param normalise Boolean takes `TRUE` if the returned matrices should be row-normalized and `FALSE` otherwise.
#' @param ceiled Boolean takes `TRUE` if the given matrices should be ceiled before conversion and `FALSE` otherwise.
#' @param byrow Boolean takes `TRUE` is entries in the matrices should be taken by row and `FALSE` if they should be taken by column.
#' @return a vector of size `sum(N*(N - 1))` or list of `length(N)` square matrices. The sizes of the matrices are `N[1], N[2], ...`
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N <- c(250, 370, 120)  
#' ## rate of friendship
#' p <- c(.2, .15, .18)   
#' ## network data
#' u <- unlist(lapply(1: 3, function(x) rbinom(N[x]*(N[x] - 1), 1, p[x])))
#' W <- vec.to.mat(u, N)
#' 
#' # Convert G into a list of row-normalized matrices
#' G <- norm.network(W)
#' 
#' # recover u
#' v <- mat.to.vec(G, ceiled = TRUE)
#' all.equal(u, v)
#' @seealso 
#' \code{\link{simnetwork}}, \code{\link{peer.avg}}.
#' @export
vec.to.mat <- function(u, N, normalise = FALSE, byrow = FALSE) {
  M        <- length(N)
  stopifnot(length(u) == sum(N*(N - 1)))
  out      <- NULL
  if (normalise) {
    out    <- frVtoMnorm(u, N, M)
  } else {
    out    <- frVtoM(u, N, M)
  }
  
  if(byrow) {
    out    <- lapply(out, t)
  }
  
  out
}


#' @rdname vec.to.mat
#' @export
mat.to.vec <- function(W, ceiled = FALSE, byrow = FALSE) {
  if (!is.list(W)) {
    if (is.matrix(W)) {
      W    <- list(W)
    } else {
      stop("W is neither a matrix nor a list")
    }
  }
  
  M        <- length(W)
  N        <- unlist(lapply(W, nrow))
  
  out      <- W
  if(byrow) {
    out    <- lapply(W, t)
  }
  if (ceiled) {
    out    <- frMceiltoV(out, N, M)
  } else {
    out    <- frMtoV(out, N, M)
  }
  
  out
}


#' @title Computing peer averages
#' @description `peer.avg` computes peer average value using network data (as a list) and observable characteristics.
#' @param Glist the adjacency matrix or list sub-adjacency matrix.
#' @param V vector or matrix of observable characteristics.
#' @param export.as.list (optional) boolean to indicate if the output should be a list of matrices or a single matrix.
#' @return the matrix product `diag(Glist[[1]], Glist[[2]], ...) %*% V`, where `diag()` is the block diagonal operator.
#' @examples 
#' # Generate a list of adjacency matrices
#' ## sub-network size
#' N  <- c(250, 370, 120)  
#' ## rate of friendship
#' p  <- c(.2, .15, .18)   
#' ## network data
#' u  <- unlist(lapply(1: 3, function(x) rbinom(N[x]*(N[x] - 1), 1, p[x])))
#' G  <- vec.to.mat(u, N, normalise = TRUE)
#' 
#' # Generate a vector y
#' y  <- rnorm(sum(N))
#' 
#' # Compute G%*%y
#' Gy <- peer.avg(Glist = G, V = y)
#' @seealso 
#' \code{\link{simnetwork}}
#' @export
peer.avg    <- function(Glist, V, export.as.list = FALSE) {
  if (!is.list(Glist)) {
    if (is.matrix(Glist)) {
      Glist <- list(Glist)
    } else {
      stop("Glist is neither a matrix nor a list")
    }
  }
  
  v.is.mat  <- !is.null(dim(V))
  V         <- as.matrix(V)
  
  cnames    <- colnames(V)
  if(!is.null(cnames)) {
    cnames  <-  paste0("G.", cnames)
  }
  
  M         <- length(Glist)
  N         <- unlist(lapply(Glist, ncol))
  
  if (sum(N) != nrow(V)) {
    stop("Glist and V do not match")
  }
  Ncum      <- c(0, cumsum(N))
  
  out       <- lapply(1:M, peer.avg.single, Glist = Glist, V = V, Ncum = Ncum, cnames = cnames, v.is.mat = v.is.mat)
  if (!export.as.list) {
    out     <- do.call(rbind, out)
  }
  if(!v.is.mat) out <- c(out)
  out
}


peer.avg.single <- function(m, Glist, V, Ncum, cnames, v.is.mat) {
  if (!is.matrix(Glist[[m]])) {
    stop("All components in Glist are not matrices")
  }
  out <- Glist[[m]]%*%V[(Ncum[m] + 1):Ncum[m + 1],,drop = FALSE]
}


#' @title Removing IDs with NA from Adjacency Matrices Optimally
#' @description
#' `remove.ids` optimally removes identifiers with NA from adjacency matrices. Many combinations of rows and columns can be deleted
# before getting rid of NA. This function removes the smallest number of rows and columns.   
#' removing many rows and column
#' @param network is a list of adjacency matrices
#' @param ncores is the number of cores to be used to run the program in parallel
#' @return List of adjacency matrices without missing values and a list of vectors of retained indeces
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%" 
#' @importFrom doRNG "%dorng%"
#' @examples 
#' A <- matrix(1:25, 5)
#' A[1, 1] <- NA
#' A[4, 2] <- NA
#' remove.ids(A)
#' 
#' B <- matrix(1:100, 10)
#' B[1, 1] <- NA
#' B[4, 2] <- NA
#' B[2, 4] <- NA
#' B[,8]   <-NA
#' remove.ids(B)
#' @export
remove.ids <- function(network, ncores = 1L){
  stopifnot(inherits(network, c("matrix", "data.frame", "list")))
  if(inherits(network, c("matrix", "data.frame"))) network <- list(network)
  
  # Construct cluster
  cl   <- makeCluster(ncores)
  registerDoParallel(cl)
  M    <- length(network)
  m    <- NULL
  out  <- foreach(m = 1:M, .packages  = "CDatanet") %dorng% {rem_non_fin(as.matrix(network[[m]]))}
  stopCluster(cl)
  network <- lapply(1:M, function(m) out[[m]]$net)
  id      <- lapply(1:M, function(m) c(out[[m]]$id))
  list(network = network, id = id)
}


fcoefficients          <- function(coef, std) {
  cnames               <- names(coef)
  tval                 <- coef/std
  pval                 <- 2*(1 - pnorm(abs(tval)))
  
  pval_print           <- unlist(lapply(pval, function(u){
    ifelse(u < 2e-16, "<2e-16", format(u, digit = 3))
  }))
  
  refprob              <- c(0.001, 0.01, 0.05, 0.1)
  refstr               <- c("***",  "**", "*", ".", "")
  str                  <- unlist(lapply(pval, function(u) refstr[1 + sum(u > refprob)]))
  
  
  out_print            <- data.frame("X1" = round(coef, 6),
                                     "X2" = round(std, 6),
                                     "X3" = round(tval, 2),
                                     "X4" = pval_print,
                                     "X5" = str)
  
  out                  <- data.frame("X1" = coef,
                                     "X2" = std,
                                     "X3" = tval,
                                     "X4" = pval)
  
  rownames(out_print)  <- cnames
  colnames(out_print)  <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "")
  rownames(out)        <- cnames
  colnames(out)        <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  
  list(out_print = out_print, out = out)
}



#' @importFrom Formula as.Formula
#' @importFrom formula.tools env
#' @importFrom stats model.frame
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom stats model.response
#' @importFrom stats model.matrix
#' @importFrom stats delete.response
#' @importFrom stats as.formula
#' @importFrom Matrix rankMatrix
#' @importFrom ddpcr quiet
formula.to.data <- function(formula,
                            contextual, 
                            Glist, 
                            M, 
                            igr, 
                            data, 
                            type          = "model",
                            theta0        = NULL,
                            fixed.effects = FALSE) {
  
  ## Extract data from the formula
  if (missing(data)) {
    data           <- env(formula)
  }
  formula          <- as.Formula(formula)
  
  if (type == "model") {
    stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  } else {
    stopifnot(length(formula)[1] == 0L, length(formula)[2] %in% 1:2)
  }
  
  
  # try to handle dots in formula
  has_dot          <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if(has_dot(formula)) {
    f1             <- formula(formula, rhs = 1)
    f2             <- formula(formula, lhs = 0, rhs = 2)
    if(!has_dot(f1) & has_dot(f2)) {
      formula      <- as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
    }
  }
  
  ## call model.frame()
  mf               <- model.frame(formula, data = data)
  ## extract response, terms, model matrices
  y                <- model.response(mf, "numeric")
  mtXone           <- terms(formula, data = data, rhs = 1)
  Xone             <- model.matrix(mtXone, mf) ## X before pipe
  cnames           <- colnames(Xone)
  
  mtX              <- NULL
  X                <- NULL
  if(length(formula)[2] >= 2L) { ## X after pipe
    mtX            <- delete.response(terms(formula, data = data, rhs = 2))
    X              <- model.matrix(mtX, mf)
  }
  
  if (contextual) {
    if (!is.null(X)) {
      stop("contextual cannot be TRUE while contextual variable are declared after the pipe.")
    }
    X              <- Xone
    tmpx           <- as.character.default(formula(formula, lhs = 1, rhs = 1))
    formula        <- Formula::as.Formula(paste(c(tmpx[c(2, 1, 3)], "|", tmpx[3]), collapse = " "))
  } 
  
  cnames.x         <- colnames(X) 
  intercept        <- "(Intercept)" %in% cnames.x
  
  if (intercept) {
    X              <- X[,-1, drop = FALSE]
  }  
  
  if(("(Intercept)" %in% cnames) & fixed.effects){
    Xone           <- Xone[, -1, drop = FALSE]
    cnames         <- cnames[-1]
  }
  # GX and Gy
  Gy                 <- NULL
  if (is.null(theta0)) {
    if(!is.null(X)) {
      GXlist         <- list()
      Gylist         <- list()
      for (m in 1:M) {
        n1           <- igr[m,1] + 1
        n2           <- igr[m,2] + 1
        GXlist[[m]]  <- Glist[[m]] %*% X[n1:n2,]
        Gylist[[m]]  <- Glist[[m]] %*% y[n1:n2]
        
        if(fixed.effects){
          y[n1:n2]     <- y[n1:n2] - mean(y[n1:n2])
          Gylist[[m]]  <- Gylist[[m]] - mean(Gylist[[m]])
          GXlist[[m]]  <- apply(GXlist[[m]], 2, function(x) x - mean(x))
          Xone[n1:n2,] <- apply(Xone[n1:n2, ,drop = FALSE], 2, function(x) x - mean(x))
        }
      }
      
      GX           <- do.call("rbind", GXlist)
      Gy           <- unlist(Gylist)
      Xone         <- cbind(Xone, GX)
      cnames       <- c(cnames,  paste0("G: ", colnames(X)))
    } else {
      Gylist         <- list()
      for (m in 1:M) {
        n1           <- igr[m,1] + 1
        n2           <- igr[m,2] + 1
        Gylist[[m]]  <- Glist[[m]] %*% y[n1:n2]
        
        if(fixed.effects){
          y[n1:n2]     <- y[n1:n2] - mean(y[n1:n2])
          Gylist[[m]]  <- Gylist[[m]] - mean(Gylist[[m]])
          Xone[n1:n2,] <- apply(Xone[n1:n2, ,drop = FALSE], 2, function(x) x - mean(x))
        }  
      }
      Gy             <- unlist(Gylist)
    }
  } else {
    if(!is.null(X)) {
      GXlist         <- list()
      for (m in 1:M) {
        n1           <- igr[m,1] + 1
        n2           <- igr[m,2] + 1
        GXlist[[m]]  <- Glist[[m]] %*% X[n1:n2,]
      }
      
      GX             <- do.call("rbind", GXlist)
      Xone           <- cbind(Xone, GX)
      cnames         <- c(cnames,  paste0("G: ", colnames(X)))
    }
  }
  if(type != "network") {
    if(rankMatrix(Xone)[1] != ncol(Xone))  {
      stop("X or [X, GX] is not a full rank matrix. May be there is an intercept in X and in GX.")
    }
  } else {
    if(rankMatrix(Xone)[1] != ncol(Xone))  {
      stop("X is not a full rank matrix. May be there is an intercept in X and you add intercept in the formula or fixed effects.")
    }
  }
  
  
  colnames(Xone)   <- cnames
  
  list("formula" = formula, 
       "X"       = Xone, 
       "y"       = y,
       "Gy"      = Gy)
}
