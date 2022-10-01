# This code plots an example of simulated data following the count data model
# with social interactions
rm(list = ls())
set.seed(12)
library(CDatanet)
library(ggplot2)
library(latex2exp)

# Parameters
lambda <- 0.25
beta   <- c(2.5, 1.5, -1.2)
gamma  <- c(0.5, -0.9)
theta  <- c(lambda, beta, gamma)

Xlab    <- c("y", "y", "y")
Ylab    <- c("Frequency", "", "")

deltal  <- list(c(1, 0.87, 0.75, 0.55, 0.05, 0.3),
                c(1.2, 0.7, 0.55, 0.5, 0.5, 0.4, 0.4, 0.3, 0.3, 0.27, 0.27, 0.25, 0.005, 0),
                c(0.4, 0))

PAT     <- LETTERS[1:3]

n       <- 1500

fgraph  <- function(j) {
  delta <- head(deltal[[j]], length(deltal[[j]]) - 2)
  bdlta <- tail(deltal[[j]], 2)[1]
  rho   <- tail(deltal[[j]], 2)[2]
  
  # X
  X1    <- rnorm(n, 1, 1)
  X2    <- rpois(n, 2)
  
  # Network
  G              <- matrix(0, n, n)
  nmax_f         <- 30
  for (i in 1:n) {
    tmp          <- sample((1:n)[-i], sample(0:nmax_f, 1))
    G[i, tmp]    <- 1
  }
  rs             <- rowSums(G); rs[rs == 0] <- 1
  G              <- G/rs
  Glist          <- list(G)
  
  # data
  ytmp           <- simcdnet(formula = ~ X1 + X2 | X1 + X2, Glist = Glist, 
                             theta = theta, deltabar = bdlta, delta = delta, rho = rho)
  y              <- ytmp$y
  cat("Rbar = ", length(delta), "\n")
  print(quantile(y, prob = c(0.8, 0.9, 0.95, 0.975, 0.99, 1)))
  
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

Graph <- lapply(1:3, fgraph)

multiplot(Graph[[1]], Graph[[2]], Graph[[3]], cols = 3)
# save with dimension 3 x 7.5
