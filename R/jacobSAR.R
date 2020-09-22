jacobSAR <- function(alpha, X, invXX, XX, y, N, G, I, Gy, ngroup, igroup, cov = TRUE) {
  Ay           <- (y - alpha * Gy)
  beta         <- invXX %*% (t(X) %*% Ay)
  Xbeta        <- X %*% beta
  s2           <- sum((Ay - Xbeta) ^ 2) / N
  nbeta        <- length(beta)
  
  covout       <- NULL
  if (cov) {
  #   GinvA        <- lapply(1:ngroup, function(x)
  #     t(solve(t(I[[x]] - alpha * G[[x]]), t(G[[x]]))))
  #   
  #   
  #   GinvAXb      <- unlist(lapply(1:ngroup, function(x)
  #     GinvA[[x]] %*% Xbeta[(igroup[x, 1]:igroup[x, 2]) + 1]))
  #   
  #   
  #   trGinA2      <- sum(unlist(lapply(1:ngroup, function(x)
  #     sum(diag(GinvA[[x]] %*% GinvA[[x]])))))
  #   trGinvAGinvA <- sum(unlist(lapply(1:ngroup, function(x)
  #     sum(diag(t(GinvA[[x]]) %*% GinvA[[x]])))))
  #   
  #   JAC          <- matrix(NA, 1 + nbeta, 1 + nbeta)
  #   JAC[1, 1]    <- trGinA2 + trGinvAGinvA + sum(GinvAXb ^ 2) / s2
  #   JAC[-1, 1]   <- t(X) %*% GinvAXb / s2
  #   JAC[1,-1]    <- t(JAC[-1, 1])
  #   JAC[-1,-1]   <- XX / s2
    covout   <- fSARjac(alpha, s2, X, XX, Xbeta, G,
                        I, igroup, ngroup, N, nbeta)
  }

  
  list("alpha"   = alpha,
       "beta"    = beta,
       "sigma2"  = s2,
       "cov"     = covout)
}