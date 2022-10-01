############### THIS CODE REPLICATES THE MONTE CARLO RESULTS FOR EXOGENEITY ############### 
rm(list = ls())
library(PartialNetwork)
library(CDatanet)
library(ggplot2)

# Parameters
lambda <- 0.25
beta   <- c(2.5, 1.5, -1.2)
gamma  <- c(0.5, -0.9)
theta  <- c(lambda, beta, gamma)

# Functions
# This function performs the one iteration of the Monte Carlo
# 'nvec' is the vector of sample sizes. 
# 'delta' is the parameter delta[2:Rbar]
# 'deltabar' is the parameter bar{delta}
# 'rho' is the parameter rho
# 'Rbarvec' is the vector of Rbar
# optimizer is the optimizer used
f.estim  <- function(nvec, delta, deltabar, rho, Rbarvec, optimizer){
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
                      theta = theta, delta = delta, deltabar = deltabar, rho = rho)
  y       <- ytmp$y
  yb      <- ytmp$yb
  tmp0    <- ytmp$marg.effects; names(tmp0) <- paste("Meff", names(tmp0))
  # print(ggplot(data = data.frame(y = y), aes(x = y)) +
  #         geom_bar(color = "black", fill = "cyan") +
  #         theme_bw() + xlab("") + ylab("Frequency"))
  
  out1    <- list()
  out2    <- list()
  rr      <- 0
  for (Rbar in Rbarvec) {
    rr    <- rr + 1
    # Use true values as starting
    dlt   <- numeric(0)
    dlt[1:length(delta)]                <- delta
    dltb                                <- deltabar
    ro                                  <- rho
    if((Rbar - 1) <= length(dlt)){
      dlt                               <- head(dlt, Rbar - 1)
      dltb                              <- 0.5*(deltabar + ifelse(is.na(delta[Rbar]), 0, delta[Rbar]))
    } else{
      dlt[(length(dlt) + 1):(Rbar - 1)] <- deltabar*(1:(Rbar - length(dlt) - 1))^rho + theta[1]
    }
    start      <- list(theta = theta, delta = dlt, deltabar = dltb, rho = max(c(ro, 1e-5)))
    # flexible rho
    # this is fast but the solution may not be good. I use it as starting in optim
    out1[[rr]] <- cdnet(y ~ X1 + X2 | X1 + X2, starting = start, Glist = Glist, optimizer = "optim", 
                        npl.ctr   = list(print = FALSE, maxit = 1e3, tol = 0.01), 
                        opt.ctr = list(control = list(maxit = 1e5, abstol = 1e-11, reltol = 1e-11)), 
                        Rbar = Rbar, yb0 = yb, cov = FALSE, estim.rho = TRUE)
    # start      <- tmp$estimate$parms
    # start      <- list(theta = start[1:6], delta = start[-c(1:6, Rbar + 6, Rbar + 7)], deltabar = start[Rbar + 6], rho = start[Rbar + 7])
    # out1[[rr]] <- cdnet(y ~ X1 + X2 | X1 + X2, Glist =  Glist, starting = start, optimizer = "optim", 
    #                     npl.ctr = list(print = TRUE, maxit = 1e9), 
    #                     opt.ctr = list(control = list(maxit = 1e5, abstol = 1e-13, reltol = 1e-13)), 
    #                     Rbar = Rbar, yb0 = tmp$yb, cov = FALSE, estim.rho = TRUE)
    
    # now set rho to 0
    # this is fast but the solution may not be good. I use it as starting in optim
    out2[[rr]] <- cdnet(y ~ X1 + X2 | X1 + X2, starting = start, Glist = Glist, optimizer = "optim", 
                        npl.ctr = list(print = FALSE, maxit = 1e3, tol = 0.01), 
                        opt.ctr = list(control = list(maxit = 1e5, abstol = 1e-11, reltol = 1e-11)), 
                        Rbar = Rbar, yb0 = yb, cov = FALSE, estim.rho = FALSE)
    # start      <- tmp$estimate$parms
    # start      <- list(theta = start[names(start) != "deltabar"], deltabar = start["deltabar"])
    # out2[[rr]] <- cdnet(y ~ X1 + X2 | X1 + X2, Glist =  Glist, starting = start, optimizer = "optim", 
    #                     npl.ctr = list(print = TRUE, maxit = 1e9),  
    #                     opt.ctr = list(control = list(maxit = 1e5, abstol = 1e-13, reltol = 1e-13)), 
    #                     Rbar = Rbar, yb0 = tmp$yb, cov = FALSE, estim.rho = FALSE)
  }
  outsart <- sart(y ~ X1 + X2 | X1 + X2, Glist =  Glist, optimizer = optimizer, 
                  print = FALSE,  opt.ctr   = list(), RE = TRUE, theta0 = c(theta, 1),
                  cov = FALSE)
  tmp1    <- outsart$estimate$theta; names(tmp1) <- paste("Coef.Tobit", names(tmp1))
  tmp2    <- outsart$estimate$marg.effects; names(tmp2) <- paste("Meff.Tobit", names(tmp2))
  
  c(unlist(c(lapply(1:length(Rbarvec), function(rr){
    tmp3 <- out1[[rr]]$estimate$parms; names(tmp3) <- paste0("Coef.CD.Rbar=", Rbarvec[rr], ".rho ", names(tmp3))
    tmp4 <- out2[[rr]]$estimate$parms; names(tmp4) <- paste0("Coef.CD.Rbar=", Rbarvec[rr], " ", names(tmp4))
    c(tmp3, tmp4)
  }))), tmp1, tmp0, 
  unlist(c(lapply(1:length(Rbarvec), function(rr){
    tmp3 <- out1[[rr]]$estimate$marg.effects; names(tmp3) <- paste0("Meff.CD.Rbar=", Rbarvec[rr], ".rho ", names(tmp3))
    tmp4 <- out2[[rr]]$estimate$marg.effects; names(tmp4) <- paste0("Meff.CD.Rbar=", Rbarvec[rr], " ", names(tmp4))
    c(tmp3, tmp4)
  }))), tmp2, quantile(y, probs = c(0.90, 0.95, 0.975, 0.99)))
}

# This function loops on f.estim
floop.estim <- function(iter, nvec, delta, deltabar, rho, Rbarvec, optimizer = "optim") {
  tmp       <- f.estim(nvec, delta, deltabar, rho, Rbarvec, optimizer)
  out       <- matrix(tmp, nrow = 1)
  for (i in 2:iter) {
    out     <- rbind(out, t(f.estim(nvec, delta, deltabar, rho, Rbarvec, optimizer)))
    cat("Iter: ", i, "/", iter, "\n")
    print(apply(out, 2, mean, na.rm = TRUE))
  }
  out
}

nvec     <- 1500
iter     <- 1e3
# PERFORM THE MONTE CARLO
# Model A
delta    <- c(1, 0.87, 0.75, 0.55)
deltabar <- 0.05
rho      <- 0.3
Rbarvec  <- c(1, 8, 9) #true Rbar is 5 and 90th percentile is 8
set.seed(123)
outA     <- floop.estim(iter, nvec, delta, deltabar, rho, Rbarvec, "nlm")
saveRDS(outA, file = "outA.RDS")

# Model B
# delta    <- c(1, 0.6, 0.52, 0.48, 0.4, rep(0.38, 4), rep(0.35, 2), rep(0.33, 2), rep(0.26,2))
delta    <- c(1.2, 0.7, 0.55, 0.5, 0.5, 0.4, 0.4, 0.3, 0.3, 0.27, 0.27, 0.25)
deltabar <- 0.005
rho      <- 0
Rbarvec  <- c(1, 8, 9) #true Rbar is 13 and 90th percentile is 8
set.seed(123)
outB     <- floop.estim(iter, nvec, delta, deltabar, rho, Rbarvec, "nlm")
saveRDS(outB, file = "outB.RDS")

# Model C
delta    <- numeric(0)
deltabar <- 0.4
rho      <- 0
Rbarvec  <- c(1, 7, 8) #true Rbar is 4 and 90th percentile is 7
set.seed(123)
outC     <- floop.estim(iter, nvec, delta, deltabar, rho, Rbarvec, "nlm")
saveRDS(outC, file = "outC.RDS")

# The summary functions
sum.func     <- function(x) {
  out        <- c(mean(x), sd(x))
  names(out) <- c("Mean", "Sd.")
  return(out)
}

# Print monte Carlo results for each sample size, dispersion and type
(resA <- t(apply(outA, 2, sum.func)))
(resB <- t(apply(outB, 2, sum.func)))
(resC <- t(apply(outC, 2, sum.func)))

# export 
write.csv(resA, file = "resA.csv")
write.csv(resB, file = "resB.csv")
write.csv(resC, file = "resC.csv")