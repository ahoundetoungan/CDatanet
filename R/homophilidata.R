#' @title Converting Data between Directed Network Models and Symmetric Network Models.
#' @param data A `matrix` or `data.frame` of the explanatory variables of the network formation model. This 
#' corresponds to the \code{X} matrix in \code{\link{homophily.fe}} or \code{\link{homophily.re}}.
#' @param nvec A vector of the number of individuals in the networks.
#' @param to Indicates the direction of the conversion. For a matrix of explanatory variables `X` (`n*(n-1)` rows), one can 
#' select lower triangular entries (`to = "lower"`) or upper triangular entries (`to = "upper"`).
#' For a triangular `X` (`n*(n-1)/2` rows), one can convert to a full matrix of `n*(n-1)` rows by using symmetry (`to = "symmetric"`).
#' @description 
#' `homophili.data` converts the matrix of explanatory variables between directed network models and symmetric network models.
#' @return The transformed `data.frame`.
#' @export
homophili.data <- function(data, nvec, to = c("lower", "upper", "symmetric")){
  to      <- tolower(to[1])
  stopifnot(to %in% c("lower", "upper", "symmetric"))
  
  M       <- length(nvec)
  n       <- sum(nvec)
  tmp1    <- NULL
  if(to == "symmetric"){
    stopifnot(nrow(data) == sum(nvec*(nvec- 1)/2))
    tmp1  <- cumsum(unlist(lapply(nvec, function(x) (x - 1):0))) - 1
  } else {
    stopifnot(nrow(data) == sum(nvec*(nvec- 1)))
    tmp1  <- cumsum(unlist(lapply(nvec, function(x) rep(x - 1, x)))) - 1
  }
  tmp2    <- c(0, tmp1[-n] + 1)
  index   <- cbind(tmp2, tmp1) 
  rm(list = c("tmp1", "tmp2"))
  indexgr <- matrix(c(cumsum(c(0, nvec[-M])), cumsum(nvec) - 1), ncol = 2) 
  out     <- NULL
  if(to == "symmetric"){
    out   <- hdata2S(as.matrix(data), nvec, index, indexgr, M)
  } 
  if(to == "lower"){
    out   <- hdataF2L(as.matrix(data), nvec, index, M)
  }
  if(to == "upper"){
    out   <- hdataF2U(as.matrix(data), nvec, index, indexgr, M)
  }
  colnames(out) <- colnames(data)
  as.data.frame(out)
}