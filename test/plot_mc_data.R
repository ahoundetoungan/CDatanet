# This code plots an example of simulated data following the cout data model 
# with social interactions

rm(list = ls())
set.seed(123)

library(CDatanet)
library(ggplot2)
library(latex2exp)


thetal  <- rbind(c(0.4, 2, -2.5, 0.6, 1.3, -1.2, 1.5),
                        c(0.4, 3, -5.5, 1.9, -0.5, -0.9, 1.5),
                        c(0.4, 1.1, 0.5, 0.5, 0.4, -0.6, 1.5),
                        c(0.4, 2, 1.8, 0.9, 0.5, 1.6, 1.5))



TYPE    <- c("A", "B")
DIS     <- c("Low", "High")
N       <- 1500

# j should be in 1:4
# 1: type A dispersion low
# 2: type A dispersion high
# 3: type B dispersion low
# 4: type B dispersion high

fgraph  <- function(j) {
  # parameters
  theta          <- thetal[,j]
  
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
  ytmp           <- simCDnet(formula = ~ X | X, Glist = Glist, theta = theta)
  y              <- ytmp$y
  
  ggplot(data = data.frame(y = y), aes(x = y)) +
    geom_bar(color = "black", fill = "#eeeeee") + 
    theme_bw() + xlab("") + ylab("Frequency")  + 
    
    ggtitle(TeX(sprintf(paste0("Type ", TYPE[ceiling(j/2)], "\n Dispersion: ",  DIS[((j - 1)%%2) + 1])))) + 
    ylab("") + theme(plot.title = element_text(size = 8, vjust = -12, hjust = 0.96))
}

# This function puts the graphics on a same loyout
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

# save with de dimension 7.42 × 4 