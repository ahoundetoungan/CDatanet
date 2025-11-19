# Estimations without fixed effects (Models 1-3)
rm(list = ls())
library(CDatanet)
library(ggplot2)
library(dplyr)
setwd("")

# Data
load("MydataCount.rda")

# Settings
form.nofix   <- as.formula(paste0(c("nclubs ~ ", paste0(expvc, collapse = " + ")), collapse = ""))
opt.ctr0     <- list(maxit = 5e3, eps_f = 1e-10, eps_g = 1e-10)
opt.ctr1     <- list(control = list(abstol = 1e-13, reltol = 1e-13, maxit = 5e3),
                     method  = "Nelder-Mead")
npl.ctr      <- list(maxit   = 1e4, tol = 5e-4, print = TRUE)
Gnetnorm     <- norm.network(Gnet)

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

## Tobit model (Model 3)
SART_nf      <- sart(formula = form.nofix, Glist = Gnetnorm, data = mydata, cov = FALSE,
                     opt.ctr = opt.ctr0, npl.ctr = npl.ctr, cinfo = FALSE)
SART_nf      <- sart(formula = form.nofix, Glist = Gnetnorm, data = mydata, cov = FALSE,
                     opt.ctr = opt.ctr1, npl.ctr = npl.ctr, cinfo = TRUE, optimizer = "optim",
                     starting = SART_nf$estimate, Ey0 = SART_nf$Ey)
(SSART_nf    <- meffects(SART_nf, Glist = Gnetnorm, data = mydata, cont.var = cont.var, 
                         bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                         ncores = 5))

## Count data model (Models 1-2)
cont  <- TRUE
BIC   <- Inf
Rbh   <- 0
CDtmp <- vector("list", max(mydata$nclubs))
while (cont) {
  Rbh          <- Rbh + 1
  starting     <- NULL
  Ey0          <- NULL
  if(Rbh > 1){
    starting   <- CDtmp[[Rbh - 1]]$estimate; starting$delta <- c(starting$delta, starting$delta[length(starting$delta)])
    Ey0        <- CDtmp[[Rbh - 1]]$Ey
  }
  CDtmp[[Rbh]] <- cdnet(formula = form.nofix, Glist =  Gnetnorm, Rbar = Rbh, data = mydata, 
                        npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = FALSE, optimizer = "optim", 
                        starting = starting, Ey0 = Ey0)
  cont         <- (BIC > CDtmp[[Rbh]]$info$BIC)
  BIC          <- CDtmp[[Rbh]]$info$BIC
}
save(SART_nf, SSART_nf, CDtmp, Rbh, file = "_output/AH_nofixed.rda")

# summary with Rbar = 0 and optimal Rbar
(SCD0_nf       <- meffects(CDtmp[[1]], Glist = Gnetnorm, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm, ncores = 5))
(SCD1_nf       <- meffects(CDtmp[[Rbh - 1]], Glist = Gnetnorm, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm, ncores = 5))
(SCD2_nf       <- meffects(CDtmp[[Rbh]], Glist = Gnetnorm, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm, ncores = 5))

save(SART_nf, SSART_nf, CDtmp, SCD0_nf, SCD1_nf, SCD2_nf, Rbh, file = "_output/AH_nofixed.rda")

################ simulate data using these results 
# (the simulations are used in file `D3_Add_Health_Endo.R` to construct Figure Figure B.1)
load(file = "_output/AH_nofixed.rda")
set.seed(123)
sform.nofix <- as.formula(paste0(c("~ ", paste0(expvc, collapse = " + ")), collapse = ""))
Rbar        <- 1
y0          <- simcdnet(formula = sform.nofix, Glist = Gnetnorm, parms = CDtmp[[Rbar]]$estimate$parms, 
                        data = mydata, Rbar = Rbar, Rmax = 33)$y
Rbar        <- 13
y1          <- simcdnet(formula = sform.nofix, Glist = Gnetnorm, parms = CDtmp[[Rbar]]$estimate$parms, 
                        data = mydata, Rbar = Rbar, Rmax = 33)$y

simuAH_nf   <- data.frame(model = "A",
                          y     = c(mydata$nclubs, y0, y1),
                          Data  = rep(1:3, each = length(y0))) %>%
  mutate(Datacode = factor(Data, labels = c(expression(paste("(A) Real data")),
                                            expression(paste("(B) Without FE, ", bar(R), " =  1")), 
                                            expression(paste("(C) Without FE, ", bar(R), " = 13")))))
saveRDS(simuAH_nf, file = "_output/simuA.RDS")