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
lGamma <- c(lambda, Gamma)
names(lGamma) <- c("lambda", "(Intercept)", "x1", "x2", "gx1", "gx2")

opt.ctr0 <- list(maxit = 5e3, eps_f = 1e-13, eps_g = 1e-13)
opt.ctr1 <- list(control = list(abstol = 1e-11, reltol = 1e-11, maxit = 5e3),
                 method  = "Nelder-Mead")
npl.ctr  <- list(maxit   = 1e3, tol = 5e-4, print = FALSE)

# Functions
# This function performs the iteration in the Monte Carlo
f.estim  <- function(n){
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
    tp          <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = G, npl.ctr = npl.ctr, 
                         opt.ctr = opt.ctr0, Rbar = Rbh, data = data, cov = FALSE)
    starting    <- NULL
    Ey0         <- NULL
    if(!is.null(tp$estimate$parms)){
      starting  <- tp$estimate
    }
    if(!is.null(tp$Ey)){
      Ey0       <- tp$Ey
    }
    ecd[[Rbh]]  <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = "optim", 
                         npl.ctr = npl.ctr, opt.ctr = opt.ctr1, Rbar = Rbh, data = data, 
                         cov = FALSE, starting = starting, Ey0 = Ey0)
    cont      <- (BIC > ecd[[Rbh]]$info$BIC)
    BIC       <- ecd[[Rbh]]$info$BIC
  }
  ccd0        <- c(ecd[[1]]$estimate$lambda, ecd[[1]]$estimate$Gamma) 
  names(ccd0) <- paste0("cd.coef.R0", names(ccd0))
  mcd0        <- ecd[[1]]$meff$ameff; names(mcd0) <- paste0("cd.meff.R0.", names(mcd0))
  ccd         <- c(ecd[[Rbh]]$estimate$lambda, ecd[[Rbh]]$estimate$Gamma) 
  names(ccd)  <- paste0("cd.coef.Rh.", names(ccd))
  mcd         <- ecd[[Rbh]]$meff$ameff; names(mcd) <- paste0("cd.meff.Rh.", names(mcd))
  
  tp          <- sart(y ~ x1 + x2 + gx1 + gx2, Glist = G, npl.ctr = npl.ctr, 
                  opt.ctr = opt.ctr0, cinfo = FALSE, data = data, cov = FALSE)
  starting    <- NULL
  Ey0         <- NULL
  if(!is.null(tp$estimate)){
    starting  <- tp$estimate
  }
  if(!is.null(tp$Ey)){
    Ey0       <- tp$Ey
  }
  esart       <- sart(y ~ x1 + x2 + gx1 + gx2, Glist = G, optimizer = "optim", npl.ctr = npl.ctr, 
                  opt.ctr = opt.ctr1, cinfo = FALSE, data = data, cov = FALSE, 
                  starting = starting, Ey0 = Ey0)
  
  cTo         <- esart$estimate; names(cTo) <- paste0("To.coef.", names(cTo))
  mTo         <- esart$meff$ameff; names(mTo) <- paste0("To.meff.", names(mTo))
  
  c(lGamma, ameff, ccd0, mcd0, "Rh" = Rbh, ccd, mcd, cTo, mTo)
}

# The summary functions
sum.func     <- function(x) {
  out        <- c("Mean" = mean(x, na.rm = TRUE), 
                  "Sd."  = sd(x, na.rm = TRUE),
                  quantile(x, na.rm = TRUE, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))
  return(out)
}


# Simulations
nsimu     <- 1000
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
n         <- 500
out500    <- mclapply(1:nsimu, function(x){cat(x, "\n"); f.estim(n = n)}, mc.cores = 2L)

n         <- 2000
out2000   <- mclapply(1:nsimu, function(x){cat(x, "\n"); f.estim(n = n)}, mc.cores = 2L)

save(out500, out2000, file = "dgpB.rda")
write.csv(cbind(t(apply(do.call(cbind, out500), 1, sum.func)), 
                t(apply(do.call(cbind, out2000), 1, sum.func))), file = "dgpB.csv")
