rm(list = ls())
library(CDatanet)
library(doParallel)

set.seed(123)

# summary function
my.sum <- function(x) {
  out <- c(mean(x, na.rm = TRUE),
           sd(x, na.rm = TRUE),
           min(x, na.rm = TRUE),
           quantile(x, 0.25, na.rm = TRUE),
           median(x, na.rm = TRUE),
           quantile(x, 0.75, na.rm = TRUE),
           max(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.", "Min", "1st Qu.", "Median", "3rd Qu.", "Max")
  return(out)
}

# estimation on one data 
f      <- function(N, disp, type) {
  lambda         <- 0.4
  
  blist          <- rbind(c(1, -1.9, 0.5),
                          c(2, -3.8, 1.3),
                          c(3, -5.5, 1.9),
                          c(2, -0.5, 0.2),
                          c(2, 1.8, 0.9),
                          c(3, 3.9, 1.7))
  
  b              <- NULL
  if (type == 1) {
    b            <- blist[disp,]
  } else {
    b            <- blist[disp + 3,]
  }
  
  
  sigma          <- 1.5
  
  theta          <- c(lambda, b, sigma)
  
  # X
  X              <- cbind(rep(1, N), rnorm(N, 1, 1), rexp(N, 0.4))

  # Network
  G              <- matrix(0, N, N)
  nmax_f         <- 30
  for (i in 1:N) {
    tmp          <- sample((1:N)[-i], sample(0:nmax_f, 1))
    G[i, tmp]    <- 1
  }
  rs             <- rowSums(G); rs[rs == 0] <- 1
  G              <- G/rs
  Glist          <- list(G)
  
  # data
  y              <- simCDnet(X, Glist, theta)$y
  
  # solver
  opt.ctr1       <- list(gradtol = 1e-12, steptol = 1e-12, iterlim = 2e3)
  opt.ctr2       <- list(iterlim = 2e3)
  opt.ctr3       <- list(control = list(method = "BFGS"))
  
  # init 
  theta          <- c(alpha, b, sigma^2)
  
  
  # RE
  resRE          <- CDnet(X = X, Glist = Glist, y = y, yb = yb, opt.ctr = opt.ctr1,
                          theta0 = theta, npl.ctr = list(print = FALSE), cov = FALSE,
                          optimizer = "nlm")
  
  
  # TOBIT
  resTO          <- SARTML(y, X, Glist, opt.ctr = opt.ctr2, theta0 = theta,
                           print =  FALSE, cov = FALSE, optimizer = "nlm")
  
  
  #SAR
  resSAR         <- SARML(y, X, Glist, opt.ctr = opt.ctr3, alpha0 = alpha,
                          print =  FALSE, cov = FALSE, optimizer = "optim")
  
  cat("RE", "\n")
  print(resRE$estimate)
  
  cat("TO", "\n")
  print(resTO$estimate)
  
  cat("SAR", "\n")
  print(resSAR$estimate)
  
  c(resRE$estimate, resTO$estimate, resSAR$estimate)
}


fMC       <- function(iteration, N, disp, type) {
  out.mc  <- mclapply(1:iteration, function(m) 
    f(m, iteration, N, disp, type), mc.cores = 4L)
}
tmp      <- fMC(iteration = 1000, N = 750, disp = 1, type = 1)



out                <-  readRDS("Output/MC_N=1500_disp=1_type=1.RDS")
out                <- do.call(rbind, out)
out[,c(5, 10, 15)] <- sqrt(out[,c(5, 10, 15)])
r.out              <- t(apply(out, 2, my.sum))
sum(is.na(out))
print(r.out)
# write.csv(r.out, file = "MC.yst.out.csv")
