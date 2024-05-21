# This code plots examples of simulated data
rm(list = ls())
set.seed(1234)
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
nvec    <- rep(250, 8)
n       <- sum(nvec)
S       <- length(nvec)

fsim    <- function(j) {
  delta <- ldelta[[j]]
  lamb  <- llambda[[j]]
  Rbar  <- length(delta)
  G     <- list()
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
  
  X     <- cbind(runif(n, 0, 5), rpois(n, 2))
  Gmu   <- G
  grp   <- rep(1, n)
  if(j >= 3){
    nc  <- c(0, cumsum(nvec))
    grp <- lapply(1:S, function(s) 1*(X[(nc[s] + 1):nc[s + 1], 1] > 2.5))
    Gmu <- lapply(1:S, function(s) norm.network(list(G[[s]] * ((1 - grp[[s]]) %*% t(1 - grp[[s]])), 
                                                     G[[s]] * ((1 - grp[[s]]) %*% t(grp[[s]])), 
                                                     G[[s]] * (grp[[s]] %*% t(1 - grp[[s]])), 
                                                     G[[s]] * (grp[[s]] %*% t(grp[[s]])))))
    Rbar  <- rep(Rbar, 2)
    delta <- rep(delta, 2)
    grp   <- unlist(grp)
  } else {
    Gmu <- norm.network(G)
  }
  
  data  <- data.frame(X, peer.avg(norm.network(G), X)); colnames(data) <- c("x1", "x2", "gx1", "gx2")
  ytmp  <- simcdnet(formula = ~ x1 + x2 + gx1 + gx2, Glist = Gmu, group = grp, lambda = lamb, 
                    Gamma = Gamma, delta = delta, Rbar = Rbar, data = data, Rmax = 100)
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
