# Estimations with heterogeneity in peer effects, while addressing network endogeneity (Model 7)
rm(list = ls())
library(CDatanet)
library(dplyr)
library(splines)
setwd("/home/aristide/Dropbox/Academy/1.Papers/CountDNtw/Code/Application")

# Data
type  <- "FE" # use RE for fixed effects
load("MydataCount.rda")
load("_output/AH_hete.rda")
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
form.en  <- as.formula(paste0(c("nclubs ~ -1 + ", paste0(c(eff, expmunu), collapse = " + ")), collapse = ""))
opt.ctr0 <- list(maxit = 5e3, eps_f = 1e-13, eps_g = 1e-13)
opt.ctr1 <- list(control = list(abstol = 1e-16, reltol = 1e-16, maxit = 5e3),
                     method  = "Nelder-Mead")
npl.ctr  <- list(maxit   = 1e4, tol = 5e-4, print = TRUE)

# list of continuous variables to compute marginal effects
cont.var     <- c("age", "yearinschl", "gage", "gmale", "ghispanic", "graceblack", 
                  "graceasian", "graceother", "gyearinschl", "gwithbothpar", "gmelhigh", "gmemhigh",
                  "gmemiss", "gmjprof", "gmjother", "gmjmiss") 

# list of binary variables to compute marginal effects
bin.var      <- c("male", "hispanic", "raceblack", "raceasian", "raceother",  "withbothpar",  
                  "melhigh", "memhigh", "memiss", "mjprof", "mjother", "mjmiss") 

# Variable type
type.var     <- list(c("age", "gage"), c("male", "gmale"), c("hispanic", "ghispanic"),
                     c("raceblack", "graceblack"), c("raceasian", "graceasian"), c("raceother", "graceother"),
                     c("yearinschl", "gyearinschl"), c("withbothpar", "gwithbothpar"),
                     c("melhigh", "gmelhigh"), c("memhigh", "gmemhigh"), c("memiss", "gmemiss"),
                     c("mjprof", "gmjprof"), c("mjother", "gmjother"), c("mjmiss", "gmjmiss"))

# Constructing network matrices
Gnetnorm     <- norm.network(Gnet)
group        <- mydata$female
GHet         <- vector("list", n.school)
n            <- sapply(Gnet, nrow)
ncum         <- c(0, cumsum(n))
for(m in 1:n.school){
  grp        <- mydata$female[(ncum[m] + 1):ncum[m + 1]]
  GHet[[m]]  <- norm.network(list(Gnet[[m]] * ((grp == 0) %*% t(grp == 0)), 
                                  Gnet[[m]] * ((grp == 0) %*% t(grp == 1)), 
                                  Gnet[[m]] * ((grp == 1) %*% t(grp == 0)), 
                                  Gnet[[m]] * ((grp == 1) %*% t(grp == 1))))
}


Rbh      <- 10
starting <- SCD_hf$estimate; starting$Gamma <- c(starting$Gamma, rep(0, KZmu + KZnu))
Ey0      <- SCD_hf$Ey

SCD_en   <- cdnet(formula = form.en, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr0, cov = FALSE, group = group, 
                  starting = starting, Ey0 = Ey0)
SCD_en   <- cdnet(formula = form.en, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = TRUE, group = group,
                  optimizer = "optim", starting = SCD_en$estimate, Ey0 = SCD_en$Ey)

(SCD_en  <- meffects(SCD_en, Glist = GHet, data = mydata, cont.var = cont.var, 
                     bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                     Glist.contextual = Gnetnorm))

saveRDS(SCD_en, file = paste0("_output/AH_Endo.", type, ".RDS"))

################ Figure B.1: Histograms of the observed and simulated dependent variables
SCD_EndoFE  <- readRDS(file = paste0("_output/AH_Endo.FE.RDS"))
set.seed(123)
sform.fix   <- as.formula(paste0(c(" ~ -1 +", paste0(c(eff, expmunu), collapse = " + ")), collapse = ""))
Rbar        <- 10
y           <- simcdnet(formula = sform.fix, Glist = GHet, parms = SCD_EndoFE$estimate$parms, 
                        data = mydata, Rbar = Rbar, Rmax = 33, group = group)$y

simuAH_hete <- data.frame(model = "D", y = y, Data  = 7) %>%
  mutate(Datacode = factor(Data, labels = c(expression(paste("(F) Het with Endo and FE, ", bar(R), " = 10")))))
saveRDS(simuAH_hete, file = "_output/simuD.RDS")


dataplot    <- rbind(readRDS(file = "_output/simuA.RDS"),
                     readRDS(file = "_output/simuB.RDS"),
                     readRDS(file = "_output/simuC.RDS"),
                     readRDS(file = "_output/simuD.RDS")) %>% filter(Data %in% c(1, 4, 5, 7)) %>%
  mutate(Datacode = factor(Data, labels = c("(A) Observed Data",
                                            "(B) Model 5 (Quadratic cost)", 
                                            "(C) Model 4 (Semiparametric cost)",
                                            "(D) Model 8 (Semparametric cost, Heter., and Endog.)")))

library(ggplot2)
(graph       <- ggplot(dataplot, aes(x = y)) + geom_bar(color = "black", fill = "#eeeeee") + 
  theme_bw() + facet_wrap(~ Datacode, ncol = 2) + 
    theme(strip.text = element_text(face = "italic"), 
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 12, family = "Palatino")) +  
  ylab("") + xlab("Number of activities"))
ggsave("plot:AH:simu.pdf", path = "_output", plot = graph, device = "pdf", width = 7, height = 4)
