# This code replicates the Monte Carlo Results

rm(list = ls())
library(CDatanet)
library(doParallel)

# You need de define the number of cores
# I set the number of available cores - 1
n.cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)-1

# This function does one iteration in the Monte Carlo
# 'm' stand for the m-th iteration to perform out of 'iteration'
# 'N' is the sample size. 'disp' is the dispersion: low or high
#  type specifies the type of data: A or B
f      <- function(m, iteration, N, type, dispersion) {
  # parameters
  thetal         <- rbind(c(0.4, 2, -2.5, 0.6, 1.3, -1.2, 1.5),
                          c(0.4, 3, -5.5, 1.9, -0.5, -0.9, 1.5),
                          c(0.4, 1.1, 0.5, 0.5, 0.4, -0.6, 1.5),
                          c(0.4, 2, 1.8, 0.9, 0.5, 1.6, 1.5))
  
  theta          <- NULL
  
  disp           <- (1:2)[c("LOW", "HIGH") == toupper(dispersion)]
  if (type == "A") {
    theta        <- thetal[disp,]
  } else {
    theta        <- thetal[disp + 2,]
  }

  # X
  X              <- cbind(rnorm(N, 1, 1), rexp(N, 0.4))
  
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
  ytmp           <- simCDnet(formula = ~ X|X, Glist = Glist, theta = theta)
  y              <- ytmp$y
  
  # CD
  CDest          <- CDnetNPL(formula = y ~ X|X, Glist = Glist, npl.ctr = list(print = F))
  
  
  # TOBIT
  TOest          <- SARTML(formula = y ~ X|X, Glist = Glist, print =  F, cov = F, optimizer = "nlm")
  
  
  #SAR
  SARest         <- SARML(formula = y ~ X|X, Glist = Glist, print =  F, cov = F, optimizer = "nlm")
  
  
  cat("Iteration : ", m, "/", iteration, "\n")
  cat("CDSI", "\n")
  print(CDest$estimate)
  
  cat("SART", "\n")
  print(TOest$estimate)
  
  cat("SAR", "\n")
  print(SARest$estimate)
  
  c(CDest$estimate, TOest$estimate, SARest$estimate)
}


# This function performs 'iteration' iteration one one model
# defined by N, dispersion and type
# and saves the in the working directory 
fMC       <- function(iteration, N, dispersion, type) {
  set.seed(123)
  out.mc  <- mclapply(1:iteration, function(m) 
    f(m, iteration, N, type, dispersion), mc.cores = n.cores)
  saveRDS(out.mc, file = paste0("_output/MC_N=", N, "_type=", type, "_disp=", dispersion, ".RDS"))
}

## Create a folder _output to save the results
if (dir.exists("_output")) {
  unlink("_output", recursive = TRUE)
}
dir.create("_output")

## RUN THE MONTE CARLO
iteration <- 10

fMC(iteration, N = 250, disp = "low", type = "A")
fMC(iteration, N = 250, disp = "low", type = "B")
fMC(iteration, N = 250, disp = "high", type = "A")
fMC(iteration, N = 250, disp = "high", type = "B")

fMC(iteration, N = 750, disp = "low", type = "A")
fMC(iteration, N = 750, disp = "low", type = "B")
fMC(iteration, N = 750, disp = "high", type = "A")
fMC(iteration, N = 750, disp = "high", type = "B")

fMC(iteration, N = 1500, disp = "low", type = "A")
fMC(iteration, N = 1500, disp = "low", type = "B")
fMC(iteration, N = 1500, disp = "high", type = "A")
fMC(iteration, N = 1500, disp = "high", type = "B")



# The sumary functions
my.sum <- function(x) {
  out        <- c(mean(x, na.rm = TRUE),
                  sd(x, na.rm = TRUE))
  names(out) <- c("Mean", "Sd.")
  return(out)
}
MC.summary <- function(N, type, dispersion) {
  cat("Sample size: ", N, "\n")
  cat("Data type  : ", type, "\n")
  cat("Dispersion : ", dispersion, "\n")
  out                <- readRDS(paste0("_output/MC_N=", N, "_type=", type, "_disp=", dispersion, ".RDS"))
  out                <- do.call(rbind, out)
  s.out              <- t(apply(out, 2, my.sum))
  s.out              <- cbind(s.out[1:7,], s.out[8:14,], s.out[15:21,])
  colnames(s.out)    <- paste(c("Mean", "Sd."), rep(c("CDSI", "SART", "SAR"), c(2, 2, 2)))
  s.out
}


# Print monte Carlo results for each sample size, dispersion and type
MC.summary(N = 250, type = "A", dispersion = "low")
MC.summary(N = 250, type = "A", dispersion = "high")
MC.summary(N = 250, type = "B", dispersion = "low")
MC.summary(N = 250, type = "B", dispersion = "high")


MC.summary(N = 750, type = "A", dispersion = "low")
MC.summary(N = 750, type = "A", dispersion = "high")
MC.summary(N = 750, type = "B", dispersion = "low")
MC.summary(N = 750, type = "B", dispersion = "high")


MC.summary(N = 1500, type = "A", dispersion = "low")
MC.summary(N = 1500, type = "A", dispersion = "high")
MC.summary(N = 1500, type = "B", dispersion = "low")
MC.summary(N = 1500, type = "B", dispersion = "high")