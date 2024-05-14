# Counterfactual Analysis 
rm(list = ls())
library(CDatanet)
library(dplyr)
library(splines)
library(ggplot2)
library(doParallel)

# Data
type  <- "FE" 
load("MydataCount.rda")
load("_output/AH_fixed.rda")
rm(list = c("datanet", "GX")); gc()

SCD_EndoFE <- readRDS(file = paste0("_output/AH_Endo.FE.RDS"))
load(paste0("_output/munu.", type, ".rda"))  # use munu.FE for fixed effects

mu     <- mu/max(abs(mu))
nu     <- nu/max(abs(nu))
Zmu    <- bs(mu, degree = 3L, knots = quantile(mu, probs = seq(0.05, 0.95, 0.05))); colnames(Zmu) <- paste0("mu", 1:ncol(Zmu))
Znu    <- bs(nu, degree = 3L, knots = quantile(nu, probs = seq(0.05, 0.95, 0.05))); colnames(Znu) <- paste0("nu", 1:ncol(Znu))
KZmu   <- ncol(Zmu)
KZnu   <- ncol(Znu)

mydata <- mydata %>% bind_cols(Zmu) %>% bind_cols(Znu)

# Settings
expmunu  <- c(expvc, paste0("mu", 1:KZmu), paste0("nu", 1:KZnu))
n        <- sapply(Gnet, nrow)
ncum     <- c(0, cumsum(n))
sumn     <- sum(n)

# This function simulate the dependent variable considering two situations.
# The first one is when nobody interacts and the second one consider the observed network
# pmale is the proportion of boys
fsimu    <- function(pmale){
  mydata$male <- rbinom(sumn, 1, pmale)
  if(length(unique(mydata$male)) == 1){
    mydata$male[1] <- 1 - mydata$male[1]
  }
  group            <- 1 - mydata$male
  
  # Constructing network matrices
  # A network in which nobody interacts
  G0          <- lapply(Gnet, function(x) x*0)
  
  # A single network
  Gnetnorm    <- norm.network(Gnet)
  
  # Heterogeneity
  GHet        <- vector("list", n.school)
  for(m in 1:n.school){
    grp       <- group[(ncum[m] + 1):ncum[m + 1]]
    GHet[[m]] <- norm.network(list(Gnet[[m]] * ((grp == 0) %*% t(grp == 0)), 
                                    Gnet[[m]] * ((grp == 0) %*% t(grp == 1)), 
                                    Gnet[[m]] * ((grp == 1) %*% t(grp == 0)), 
                                    Gnet[[m]] * ((grp == 1) %*% t(grp == 1))))
  }
  

  # Simulations
  Ey0 <- simcdEy(object = SCD_f, Glist = G0, data = mydata)
  Ey1 <- simcdEy(object = SCD_f0, Glist = Gnetnorm, data = mydata)
  Ey2 <- simcdEy(object = SCD_f, Glist = Gnetnorm, data = mydata)
  Ey3 <- simcdEy(object = SCD_EndoFE, Glist = GHet, data = mydata, group = group)
  data.frame(pmale = pmale, 
             model = rep(0:3, each = sumn),
             Ey    = c(Ey0$Ey, Ey1$Ey, Ey2$Ey, Ey3$Ey),
             seaEy = rep(c(Ey0$se.aEy, Ey1$se.aEy, Ey2$se.aEy, Ey3$se.aEy), each = sumn))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
# Run on 7 cores and it took 2 days
simu   <- mclapply(seq(0, 1, 0.01), function(x){cat(x, "\n"); fsimu(x)}, mc.cores = 7L)

saveRDS(simu, "_output/AH_CFA.RDS")
simu   <- readRDS("_output/AH_CFA.RDS")

out    <- do.call(rbind, simu) %>% group_by(pmale, model) %>% 
  summarise(m = mean(Ey), s = mean(seaEy)) %>%
  mutate(model = factor(model, labels =  c("(A) Model 4 without interations (Semparametric cost)", 
                                           "(B) Model 5 (Quadratic cost)",
                                           "(C) Model 4 (Semiparametric cost)",
                                           "(D) Model 8 (Semparametric cost, Heter., and Endog.)")))

# Figure 4: Counterfactual analysis results
(graph <- ggplot(out, aes(x = pmale, y = m)) +
    # geom_errorbar(width = .02, aes(ymin = m - qnorm(0.975)*s, ymax = m + qnorm(0.975)*s)) + 
    geom_ribbon(aes(ymin =  m - qnorm(0.975)*s, ymax = m + qnorm(0.975)*s), 
                color = "#c90000", fill = "#d0ffff", linewidth = 0.3, linetype = 2) +
    geom_line(linewidth = 0.5, linetype = 1, color = "#240423") + theme_bw() + 
    xlab("Proportion of male students") + ylab("Number of activities") + 
    facet_wrap(~ model, ncol = 2) + 
    theme(strip.text = element_text(face = "italic"), 
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 12, family = "Palatino")))

ggsave("plot:AH:CFA.pdf", path = "_output", plot = graph, device = "pdf", width = 7, height = 5)

write.csv(out, row.names = FALSE, file = "_output/plot:AH:CFA.csv")
