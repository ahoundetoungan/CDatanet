############### THIS CODE REPLICATES THE MONTE CARLO RESULTS FOR DGP B ############### 
rm(list = ls())
library(CDatanet)
library(doParallel)
pdir <- c("~/Dropbox/Papers - In progress/CountDNtw/Code/Monte Carlo/_output",
          "~/CountDNtw/Code/Monte Carlo/_output")
setwd(pdir[sapply(pdir, dir.exists)])

# Parameters
lambda <- 0.25
Gamma  <- c(0.5, 1.5, -1.2, 0.5, -0.9)
delta  <- c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005)
Rbar   <- length(delta)
lGamma <- c(lambda, Gamma); 
names(lGamma) <- c("lambda", "(Intercept)", "x1", "x2", "gx1", "gx2")

# Functions
# This function performs the iteration in the Monte Carlo
f.estim  <- function(n, optimizer = "fastlbfgs"){
  G      <- matrix(0, n, n)
  for (i in 1:n) {
    max_d        <- 10
    tmp          <- sample((1:n)[-i], sample(0:max_d, 1))
    G[i, tmp]    <- 1
  }
  G       <- norm.network(G)
  
  X       <- cbind(runif(n, 0, 5), rpois(n, 2))
  data    <- data.frame(X, peer.avg(G, X)); colnames(data) <- c("x1", "x2", "gx1", "gx2")
  
  ytmp    <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, lambda = lambda, 
                      Gamma = Gamma, delta = delta, Rbar = Rbar, data = data)
  data$y  <- ytmp$y
  # hist(data$y, breaks = max(data$y) + 1)
  ameff   <- ytmp$meff$ameff; names(ameff) <- paste("ameff", names(ameff))
  
  cont    <- TRUE
  Rbh     <- 0
  BIC     <- Inf
  ecd     <- list()
  while(cont){
    Rbh         <- Rbh + 1
    ecd[[Rbh]]  <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = optimizer, 
                     npl.ctr   = list(print = FALSE, maxit = 1e3, tol = 0.01), 
                     opt.ctr = list(maxit = 1e9, eps_f = 1e-10, eps_g = 1e-10),
                     Rbar = Rbh, data = data, cov = FALSE)
    cont        <- (BIC > ecd[[Rbh]]$info$BIC)
    BIC         <- ecd[[Rbh]]$info$BIC
  }
  ccd0    <- c(ecd[[1]]$estimate$lambda, ecd[[1]]$estimate$Gamma) 
  names(ccd0) <- paste0("cd.coef.R0", names(ccd0))
  mcd0    <- ecd[[1]]$meff$ameff; names(mcd0) <- paste0("cd.meff.R0.", names(mcd0))
  ccd     <- c(ecd[[Rbh]]$estimate$lambda, ecd[[Rbh]]$estimate$Gamma) 
  names(ccd) <- paste0("cd.coef.Rh.", names(ccd))
  mcd     <- ecd[[Rbh]]$meff$ameff; names(mcd) <- paste0("cd.meff.Rh.", names(mcd))
  
  esart   <- sart(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = optimizer, 
                  npl.ctr   = list(print = FALSE, maxit = 1e3, tol = 0.001), 
                  opt.ctr = list(maxit = 1e9, eps_f = 1e-12, eps_g = 1e-12),
                  cinfo = FALSE, data = data, cov = FALSE)
  
  cTo     <- esart$estimate; names(cTo) <- paste0("To.coef.", names(cTo))
  mTo     <- esart$meff$ameff; names(mTo) <- paste0("To.meff.", names(mTo))
  
  c(lGamma, ameff, ccd0, mcd0, "Rh" = Rbh, ccd, mcd, cTo, mTo)
}

# The summary functions
sum.func     <- function(x) {
  out        <- c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.")
  return(out)
}


# Simulations
nsimu     <- 1000
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
n         <- 500
out500    <- mclapply(1:nsimu, function(x) f.estim(n = n), mc.cores = 2L)

n         <- 2000
out2000   <- mclapply(1:nsimu, function(x) f.estim(n = n), mc.cores = 2L)
write.csv(cbind(t(apply(do.call(cbind, out500), 1, sum.func)), 
                t(apply(do.call(cbind, out2000), 1, sum.func))), file = "dgpB.csv")
