############### THIS CODE REPLICATES THE MONTE CARLO RESULTS FOR DGP C ############### 
rm(list = ls())
library(CDatanet)
library(doParallel)
pdir <- c("~/Dropbox/Academy/1.Papers/CountDNtw/Code/Monte Carlo/_output",
          "~/CountDNtw/Code/Monte Carlo/_output")
setwd(pdir[sapply(pdir, dir.exists)])

# Parameters
lambda <- c(0.3, 0.15, 0.1, 0.15)
Gamma  <- c(2, 1.5, -1.2, 0.5, -0.9)
delta  <- c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005)
Rbar   <- length(delta)
lGamma <- c(lambda, Gamma); 
names(lGamma) <- c(paste0("lambda:", 1:4), "(Intercept)", "x1", "x2", "gx1", "gx2")

opt.ctr0 <- list(maxit = 5e3, eps_f = 1e-13, eps_g = 1e-13)
opt.ctr1 <- list(control = list(abstol = 1e-11, reltol = 1e-11, maxit = 5e3),
                 method  = "Nelder-Mead")
npl.ctr  <- list(maxit   = 2e3, tol = 5e-4, print = FALSE)

# Functions
# This function performs the iteration in the Monte Carlo
f.estim  <- function(nvec){
  n      <- sum(nvec)
  S      <- length(nvec)
  
  G      <- list()
  for(s in 1:S){
    ns           <- nvec[s]
    Gs           <- matrix(0, ns, ns)
    for (i in 1:ns) {
      max_d      <- 10
      tmp        <- sample((1:ns)[-i], sample(0:max_d, 1))
      Gs[i, tmp] <- 1
    }
    G[[s]]       <- Gs
  }
  
  X       <- cbind(runif(n, 0, 5), rpois(n, 2))
  nc      <- c(0, cumsum(nvec))
  grp     <- lapply(1:S, function(s) 1*(X[(nc[s] + 1):nc[s + 1], 1] > 2.5))
  Gmu     <- lapply(1:S, function(s) norm.network(list(G[[s]] * ((1 - grp[[s]]) %*% t(1 - grp[[s]])), 
                                                   G[[s]] * ((1 - grp[[s]]) %*% t(grp[[s]])), 
                                                   G[[s]] * (grp[[s]] %*% t(1 - grp[[s]])), 
                                                   G[[s]] * (grp[[s]] %*% t(grp[[s]])))))
  grp     <- unlist(grp)
 
  data    <- data.frame(X, peer.avg(norm.network(G), X)); colnames(data) <- c("x1", "x2", "gx1", "gx2")
  
  ytmp    <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = Gmu, group = grp, lambda = lambda, 
                      Gamma = Gamma, delta = rep(delta, 2), Rbar = rep(Rbar, 2), data = data, Rmax = 100,
                      cont.var = c("x1", "x2", "gx1", "gx2"))
  data$y  <- ytmp$y
  # hist(data$y, breaks = max(data$y) + 1)
  ameff   <- ytmp$meff$ameff; names(ameff) <- paste("ameff", names(ameff))
  
  cont    <- TRUE
  Rbh     <- 0
  BIC     <- Inf
  ecd     <- list()
  while(cont){
    Rbh         <- Rbh + 1
    tp          <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = Gmu, group = grp, npl.ctr = npl.ctr, 
                         opt.ctr = opt.ctr0, Rbar = rep(Rbh, 2), data = data, cov = FALSE, Rmax = 100)
    starting    <- NULL
    Ey0         <- NULL
    if(!is.null(tp$estimate$parms)){
      starting  <- tp$estimate
    }
    if(!is.null(tp$Ey)){
      Ey0       <- tp$Ey
    }
    ecd[[Rbh]]  <- cdnet(y ~ x1 + x2 + gx1 + gx2, Glist = Gmu, group = grp, npl.ctr = npl.ctr, 
                         opt.ctr = opt.ctr1, Rbar = rep(Rbh, 2), data = data, cov = FALSE, 
                         optimizer = "optim", starting = starting, Ey0 = Ey0, Rmax = 100)
    cont        <- (BIC > ecd[[Rbh]]$info$BIC)
    BIC         <- ecd[[Rbh]]$info$BIC
  }
  
  ccd0          <- c(ecd[[1]]$estimate$lambda, ecd[[1]]$estimate$Gamma) 
  names(ccd0)   <- paste0("cd.coef.R0", names(ccd0))
  ccd           <- c(ecd[[Rbh]]$estimate$lambda, ecd[[Rbh]]$estimate$Gamma) 
  names(ccd)    <- paste0("cd.coef.Rh.", names(ccd))
  
  # Marginal effects
  mcd0          <- meffects(ecd[[1]], Glist = Gmu, cont.var = c("x1", "x2", "gx1", "gx2"), data = data, 
                            boot = 0, progress = FALSE)$meffects$estimate$direct
  names(mcd0)   <- paste0("cd.meff.R0.", names(mcd0))
  mcd           <- meffects(ecd[[Rbh]], Glist = Gmu, cont.var = c("x1", "x2", "gx1", "gx2"), data = data, 
                            boot = 0, progress = FALSE)$meffects$estimate$direct
  names(mcd)    <- paste0("cd.meff.Rh.", names(mcd))
  
  c(lGamma, ameff, ccd0, mcd0, "Rh" = Rbh, ccd, mcd)
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
nvec      <- rep(250, 2)
out500    <- mclapply(1:nsimu, function(x){cat(x, "\n"); f.estim(nvec)}, mc.cores = 4L)

nvec      <- rep(250, 8)
out2000   <- mclapply(1:nsimu, function(x){cat(x, "\n"); f.estim(nvec)}, mc.cores = 4L)

save(out500, out2000, file = "dgpC.rda")
write.csv(cbind(t(apply(do.call(cbind, out500), 1, sum.func)), 
                t(apply(do.call(cbind, out2000), 1, sum.func))), file = "dgpC.csv")
