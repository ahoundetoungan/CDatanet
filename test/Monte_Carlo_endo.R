############### THIS CODE REPLICATES THE MONTE CARLO RESULTS WITH ENDOGENEITY ############### 
rm(list = ls())
library(PartialNetwork)
library(CDatanet)
library(ggplot2)

# Parameters peer effect model
lambda <- 0.3
beta   <- c(2.5, 1.5, -1.2)
betast <- c(0.8, 0.5)
gamma  <- c(0.5, -0.9)
theta  <- c(lambda, beta, betast, gamma)

# Parameters network model
dbeta0 <- c(-2.3, -1.8, -1.2, -2)
dbeta  <- c(.1, .1)
smu2   <- 0.2
snu2   <- 0.3
rho    <- 0.2
Smunu  <- matrix(c(smu2, rho*sqrt(smu2*snu2), rho*sqrt(smu2*snu2), snu2), 2)

# Functions
# This function performs the one iteration of the Monte Carlo
# 'nvec' is the vector of sample sizes. 
# 'delta' is the parameter delta
# 'Rbar' specifies Rbar
f.estim <- function(nvec, delta, Rbar) {
  M      <- length(nvec)
  mu     <- c()
  nu     <- c()
  X1     <- c()
  X2     <- c()
  distr  <- list()
  dX1    <- list()
  dX2    <- list()
  for (m in 1:M) {
    n          <- nvec[m]
    tmp        <- MASS::mvrnorm(n, c(0, 0), Smunu)
    mum        <- tmp[,1]
    num        <- tmp[,2]
    X1m        <- rnorm(n, 1, 1)
    X2m        <- rpois(n, 2)
    dX1m       <- abs(kronecker(X1m, t(X1m), "-"))
    dX2m       <- abs(kronecker(X2m, t(X2m), "-"))
    distr[[m]] <- pnorm(dbeta0[m] + dbeta[1]*dX1m + dbeta[1]*dX2m + kronecker(mum, t(num), "+"))
    diag(distr[[m]]) <- 0
    dX1[[m]]   <- dX1m
    dX2[[m]]   <- dX2m
    X1         <- c(X1, X1m)
    X2         <- c(X2, X2m)
    mu         <- c(mu, mum)
    nu         <- c(nu, num)
  }
  
  dX      <- cbind(mat.to.vec(dX1), mat.to.vec(dX2))
  network <- sim.network(dnetwork = distr)
  Glist   <- norm.network(network)
  
  ytmp    <- simCDnet(formula = ~ X1 + X2 + mu + nu| X1 + X2, Glist = Glist, 
                      theta = theta, delta = delta)
  y       <- ytmp$y
  yb      <- ytmp$yb
  
  tmp         <- c()
  tmp[1:length(delta)]                <- delta
  tmp[(length(delta) + 1):(Rbar - 1)] <- tail(delta,1)
  
  start1 <- list(theta = c(lambda, beta, gamma), delta = head(tmp, Rbar - 1))
  est1   <- CDnetNPL(y ~ X1 + X2| X1 + X2, Glist =  Glist, starting = start1, optimizer = "optim", 
                     npl.ctr   = list(print = FALSE, maxit = 1e3), Rbar = Rbar, yb = yb, cov = FALSE)
  
  init   <- list(beta = c(dbeta0, dbeta), mu = mu, nu = nu, smu2 = smu2, snu2 = snu2, rho = rho)
  out.n  <- homophily(network =  network, formula = ~ dX, iteration = 1e3, fixed.effects = TRUE,
                      print = FALSE)
  
  mus    <- colMeans(tail(out.n$posterior$mu, 5e2)) 
  nus    <- colMeans(tail(out.n$posterior$nu, 5e2)) 
  
  start2 <- list(theta = c(lambda, beta, betast,  gamma), delta = head(tmp, Rbar - 1))
  
  est2   <- CDnetNPL(y ~ X1 + X2 + mus + nus| X1 + X2, Glist =  Glist, starting = start2, optimizer = "optim", 
                     npl.ctr   = list(print = FALSE, maxit = 1e3), Rbar = Rbar, yb = yb, cov = FALSE)
  
  c(est1$estimate$theta, est1$estimate$delta, est2$estimate$theta, est2$estimate$delta,
    ytmp$marg.effects, est1$estimate$marg.effects, est2$estimate$marg.effects, 
    quantile(y, probs = c(0.90, 0.95, 0.975, 0.99)))
}

# This function loops on f.estim
floop.estim <- function(iter, nvec, Rbar, delta) {
  out       <- matrix(NA, 0, 35 + 2*Rbar)
  for (i in 1:iter) {
    out     <- rbind(out, t(f.estim(nvec, delta, Rbar)))
    cat("Iter: ", i, "/", iter, "\n")
    print(apply(out, 2, mean, na.rm = TRUE))
  }
  out
}

## PERFORM THE MONTE CARLO
nvec     <- rep(250, 4)
delta    <- c(1, 0.87, 0.75, 0.55, 0.35)
Rbar     <- 8 #correspond to quantile 90
ncores   <- 2
set.seed(123)
out      <- floop.estim(1000, nvec, Rbar, delta)

# The summary functions
sum.func     <- function(x) {
  out        <- c(mean(x, na.rm = TRUE),
                  sd(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.")
  return(out)
}

apply(out, 2, sum.func)