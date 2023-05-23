jacobSAR <- function(alpha, X, invXX, XX, y, N, G, I, Gy, ngroup, igroup, 
                     cov = TRUE, fixed.effects = FALSE) {
  Ay           <- (y - alpha * Gy)
  beta         <- invXX %*% (t(X) %*% Ay)
  Xbeta        <- X %*% beta
  s2           <- sum((Ay - Xbeta) ^ 2) / ifelse(fixed.effects, N - ngroup, N)
  nbeta        <- length(beta)
  
  covout       <- NULL
  if (cov) {
    covout   <- fSARjac(alpha, s2, X, XX, Xbeta, G, I, igroup, ngroup, N, 
                        nbeta, fixed.effects)
  }

  
  list("alpha"   = alpha,
       "beta"    = beta,
       "sigma2"  = s2,
       "cov"     = covout)
}