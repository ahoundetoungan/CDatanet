#' @importFrom Formula as.Formula
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
                            type   = "model",
                            theta0 = NULL) {
  
  ## Extract data from the formula
  if (missing(data)) {
    data           <- environment(formula)
  }
  cl               <- match.call()
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
  Xone             <- model.matrix(mtXone, mf)
  cnames           <- colnames(Xone)
  
  mtX              <- NULL
  X                <- NULL
  if(length(formula)[2] >= 2L) {
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
  
  # GX and Gy
  Gy               <- NULL
  if (is.null(theta0)) {
    if(!is.null(X)) {
      GXlist         <- list()
      Gylist         <- list()
      for (m in 1:M) {
        n1           <- igr[m,1] + 1
        n2           <- igr[m,2] + 1
        GXlist[[m]]  <- Glist[[m]] %*% X[n1:n2,]
        Gylist[[m]]  <- Glist[[m]] %*% y[n1:n2]
      }
      
      GX             <- do.call("rbind", GXlist)
      Gy             <- unlist(Gylist)
      Xone           <- cbind(Xone, GX)
      cnames         <- c(cnames,  paste0("G: ", colnames(X)))
    } else {
      Gylist         <- list()
      for (m in 1:M) {
        n1           <- igr[m,1] + 1
        n2           <- igr[m,2] + 1
        Gylist[[m]]  <- Glist[[m]] %*% y[n1:n2]
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