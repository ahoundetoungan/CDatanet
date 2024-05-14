# This code plots examples of simulated data
rm(list = ls())
set.seed(1)
library(CDatanet)
library(ggplot2)
library(latex2exp)

# Parameters
Gamma   <- c(0.5, 1.5, -1.2, 0.5, -0.9)
llambda <- list(0.25, 0.25,
                c(0.3, 0.15, 0.1, 0.15),
                c(0.4, -0.1, 0.2, 0.1))
ldelta  <- list(0.3,
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005),
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005),
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005))

n       <- 2000

fsim    <- function(j) {
  delta <- ldelta[[j]]
  lamb  <- llambda[[j]]
  Rbar  <- length(delta)
  G     <- matrix(0, n, n)
  for (i in 1:n) {
    max_d        <- 10
    tmp          <- sample((1:n)[-i], sample(0:max_d, 1))
    G[i, tmp]    <- 1
  }
  
  X     <- cbind(runif(n, 0, 5), rpois(n, 2))
  Gmu   <- G
  grp   <- rep(1, n)
  if(j >= 3){
    grp <- 1*(X[,1] > 2.5)
    Gmu <- list(norm.network(list(G * ((1 - grp) %*% t(1 - grp)), 
                                  G * ((1 - grp) %*% t(grp)), 
                                  G * (grp %*% t(1 - grp)), 
                                  G * (grp %*% t(grp)))))
    Rbar  <- rep(Rbar, 2)
    delta <- rep(delta, 2)
  } else {
    Gmu <- norm.network(G)
  }
  
  data  <- data.frame(X, peer.avg(norm.network(G), X)); colnames(data) <- c("x1", "x2", "gx1", "gx2")
  ytmp  <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = Gmu, group = grp, lambda = lamb, 
                    Gamma = Gamma, delta = delta, Rbar = Rbar, data = data)
  data.frame(y = ytmp$y, DGP = paste0("DGP ", LETTERS[j]))
}

out     <- do.call(rbind, lapply(1:4, fsim))

# Figure 2: Simulated data using the count data model with social interactions
(graph  <- ggplot(out, aes(x = y)) +
  geom_bar(color = "black", fill = "#eeeeee") + 
  theme_bw() + xlab("y") + ylab("Frequency") + 
  facet_wrap(~ DGP, ncol = , scales = "free") + 
    theme(strip.text = element_text(face = "italic"), 
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 12, family = "Palatino")))

ggsave("mc_plot.pdf", path = "~/Dropbox/Papers - In progress/CountDNtw/Code/Monte Carlo/_output", 
       plot = graph, device = "pdf", width = 7, height = 4)
