# This code plots examples of simulated data
rm(list = ls())
set.seed(1234)
library(CDatanet)
library(ggplot2)
library(latex2exp)

# Parameters
Gamma   <- c(0.5, 1.5, -1.2, 0.5, -0.9)
Xlab    <- c("", "", "y", "y")
Ylab    <- c("Frequency", "", "Frequency", "")
llambda <- list(0.25, 0.25,
                c(0.3, 0.15, 0.1, 0.15),
                c(0.4, -0.1, 0.2, 0.1))
ldelta  <- list(0.3,
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005),
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005),
                c(1.8, 1, 0.6, 0.45, 0.25, 0.15, 0.08, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005))

PAT     <- LETTERS[1:4]

n       <- 2000

fgraph  <- function(j) {
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
  
  y     <- ytmp$y
  
  ggplot(data = data.frame(y = y), aes(x = y)) +
    geom_bar(color = "black", fill = "#eeeeee") + 
    theme_bw() + xlab(Xlab[j]) + ylab(Ylab[j]) + 
    ggtitle(TeX(sprintf(paste0("DGP ", PAT[j])))) + 
    ylab("") + theme(plot.title = element_text(size = 8, vjust = -12, hjust = 0.96))
}

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Graph <- lapply(1:4, fgraph)

multiplot(Graph[[1]], Graph[[3]], Graph[[2]], Graph[[4]], cols = 2)
# save with dimension 5 x 6
