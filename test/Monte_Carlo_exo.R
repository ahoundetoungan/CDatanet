############### THIS CODE REPLICATES THE MONTE CARLO RESULTS FOR EXOGENEITY ############### 
rm(list = ls())
library(PartialNetwork)
library(CDatanet)
library(ggplot2)

# Parameters
lambda <- 0.3
beta   <- c(2.5, 1.5, -1.2)
gamma  <- c(0.5, -0.9)
theta  <- c(lambda, beta, gamma)

# Functions
# This function performs the one iteration of the Monte Carlo
# 'nvec' is the vector of sample sizes. 
# 'delta' is the parameter delta
# 'Rbar' specifies Rbar
f.estim  <- function(nvec, delta, Rbar){
  n      <- sum(nvec)
  M      <- length(nvec)
  
  # X
  X1     <- rnorm(n, 1, 1)
  X2     <- rpois(n, 2)
  
  # Network
  Glist  <- list()
  for (m in 1:M) {
    nm           <- nvec[m]
    Gm           <- matrix(0, nm, nm)
    max_d        <- 30
    for (i in 1:nm) {
      tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
      Gm[i, tmp] <- 1
    }
    rs           <- rowSums(Gm); rs[rs == 0] <- 1
    Gm           <- Gm/rs
    Glist[[m]]   <- Gm
  }
  ytmp    <- simcdnet(formula = ~ X1 + X2 | X1 + X2, Glist = Glist, 
                      theta = theta, delta = delta)
  y       <- ytmp$y
  yb      <- ytmp$yb
  print(ggplot(data = data.frame(y = y), aes(x = y)) +
          geom_bar(color = "black", fill = "cyan") +
          theme_bw() + xlab("") + ylab("Frequency"))
  
  tmp     <- c()
  tmp[1:length(delta)]                <- delta
  tmp[(length(delta) + 1):(Rbar - 1)] <- tail(delta,1)
  start   <- list(theta = theta, delta = head(tmp, Rbar - 1))
  out1    <- cdnet(y ~ X1 + X2 | X1 + X2, Glist =  Glist, starting = start,  optimizer = "optim", 
                   npl.ctr   = list(print = FALSE),  opt.ctr   = list(), Rbar = Rbar, yb0 = yb,
                   cov = FALSE)
  start   <- list(theta = theta, delta = mean(start$delta))
  out2    <- cdnet(y ~ X1 + X2 | X1 + X2, Glist =  Glist, starting = start,  optimizer = "optim", 
                   npl.ctr   = list(print = FALSE),  opt.ctr   = list(), Rbar = 2, yb0 = yb,
                   cov = FALSE)
  
  outsart <- sart(y ~ X1 + X2 | X1 + X2, Glist =  Glist, optimizer = "optim", 
                  print = FALSE,  opt.ctr   = list(), RE = TRUE, theta0 = c(theta, 1),
                  cov = FALSE)
  
  c(out1$estimate$theta, out1$estimate$delta, 
    out2$estimate$theta, out2$estimate$delta, 
    outsart$estimate$theta, 
    ytmp$marg.effects, 
    out1$estimate$marg.effects, out2$estimate$marg.effects, outsart$estimate$marg.effects, 
    quantile(y, probs = c(0.90, 0.95, 0.975, 0.99)))
}

# This function loops on f.estim
floop.estim <- function(iter, nvec, Rbar, delta) {
  out       <- matrix(NA, 0, 43 + Rbar)
  for (i in 1:iter) {
    out     <- rbind(out, t(f.estim(nvec, delta, Rbar)))
    cat("Iter: ", i, "/", iter, "\n")
    print(apply(out, 2, mean, na.rm = TRUE))
  }
  out
}

## PERFORM THE MONTE CARLO
# Model A
delta    <- c(1, 0.87, 0.75, 0.55, 0.35)
Rbar     <- 8 #correspond to quantile 90
####### n = 500
nvec     <- 500
set.seed(123)
outA500  <- floop.estim(1000, nvec, Rbar, delta)

####### n = 1500
nvec     <- 1500
set.seed(123)
outA1500 <- floop.estim(1000, nvec, Rbar, delta)

# Model B
delta    <- c(1.2, 0.7, 0.55, rep(0.5, 2), rep(0.4, 2), rep(0.3, 2), rep(0.25, 2), 0.2)
Rbar     <- 9 #correspond to quantile 90
####### n = 500
nvec     <- 500
set.seed(123)
outB500  <- floop.estim(1000, nvec, Rbar, delta)

####### n = 1500
nvec     <- 1500
set.seed(123)
outB1500 <- floop.estim(1000, nvec, Rbar, delta)

# Model C (True Rbar = 1)
delta    <- 1
Rbar     <- 5 #correspond to quantile 90
####### n = 500
nvec     <- 500
set.seed(123)
outC500  <- floop.estim(1000, nvec, Rbar, delta)

####### n = 1500
nvec     <- 1500
set.seed(123)
outC1500 <- floop.estim(1000, nvec, Rbar, delta)

# The summary functions
sum.func     <- function(x) {
  out        <- c(mean(x, na.rm = TRUE),
                  sd(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.")
  return(out)
}


# Print monte Carlo results for each sample size, dispersion and type
apply(outA500, 2, sum.func)
apply(outA1500, 2, sum.func)
apply(outB500, 2, sum.func)
apply(outB1500, 2, sum.func)
apply(outC500, 2, sum.func)
apply(outC1500, 2, sum.func)
