############### THIS CODE REPLICATES THE MONTE CARLO RESULTS FOR DGP A ############### 
rm(list = ls())
library(CDatanet)
library(doParallel)
pdir <- c("~/Dropbox/Papers - In progress/CountDNtw/Code/Monte Carlo/_output",
          "~/CountDNtw/Code/Monte Carlo/_output")
setwd(pdir[sapply(pdir, dir.exists)])

# Parameters
lambda <- 0.25
Gamma  <- c(2, 1.5, -1.2, 0.5, -0.9)
delta  <- 0.3
parms  <- c(lambda, Gamma, delta); 
lGamma <- c(lambda, Gamma)
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
  
  ytmp    <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = G, 
                      parms = parms, Rbar = 1, data = data)
  data$y  <- ytmp$y
  # hist(data$y, breaks = max(data$y) + 1)
  ameff   <- ytmp$meff$ameff; names(ameff) <- paste("ameff", names(ameff))
  
  ecdnet  <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = optimizer, 
                   npl.ctr   = list(print = FALSE, maxit = 5e2, tol = 0.001), 
                   opt.ctr = list(maxit = 1e9, eps_f = 1e-20, eps_g = 1e-20),
                   Rbar = 1, data = data, cov = FALSE)
  ccd     <- c(ecdnet$estimate$lambda, ecdnet$estimate$Gamma); names(ccd) <- paste0("cd.coef.", names(ccd))
  mcd     <- ecdnet$meff$ameff; names(mcd) <- paste0("cd.meff.", names(mcd))
  
  esart   <- sart(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = optimizer, 
                npl.ctr   = list(print = FALSE, maxit = 5e2, tol = 0.001), 
                opt.ctr = list(maxit = 1e9, eps_f = 1e-20, eps_g = 1e-20),
                cinfo = FALSE, data = data, cov = FALSE)
  
  cTo     <- esart$estimate; cTo <- cTo[-c(2, length(cTo))]; names(cTo) <- paste0("To.coef.", names(cTo))
  mTo     <- esart$meff$ameff; names(mTo) <- paste0("To.meff.", names(mTo))
  
  c(lGamma, ameff, ccd, mcd, cTo, mTo)
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
               t(apply(do.call(cbind, out2000), 1, sum.func))), file = "dgpA.csv")
