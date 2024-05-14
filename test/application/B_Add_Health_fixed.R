# Estimations with fixed effects (Models 4-6)
rm(list = ls())
library(CDatanet)
library(ggplot2)
library(dplyr)

# Data
load("MydataCount.rda")
load("_output/AH_nofixed.rda")

# Settings
form.fix     <- as.formula(paste0(c("nclubs ~ -1 + ", paste0(c(eff, expvc), collapse = " + ")), collapse = ""))
opt.ctr0     <- list(maxit = 5e3, eps_f = 1e-10, eps_g = 1e-10)
opt.ctr1     <- list(control = list(abstol = 1e-13, reltol = 1e-13, maxit = 5e3),
                     method  = "Nelder-Mead")
npl.ctr      <- list(maxit   = 1e4, tol = 5e-4, print = TRUE)
Gnetnorm     <- norm.network(Gnet)

## Tobit model (Model 6)
starting    <- SART_nf$estimate
starting    <- c(starting[1], rep(starting[2], n.school), starting[-(1:2)]); names(starting) <- NULL
SART_f      <- sart(formula = form.fix, Glist = Gnetnorm, data = mydata, starting = starting,
                    cov = FALSE, npl.ctr = npl.ctr, cinfo = FALSE, Ey0 = SART_nf$Ey, 
                    opt.ctr = opt.ctr0)
SART_f      <- sart(formula = form.fix, Glist = Gnetnorm, data = mydata, starting = SART_f$estimate,
                    cov = FALSE, opt.ctr = opt.ctr1, npl.ctr = npl.ctr, cinfo = FALSE, optimizer = "optim",
                    Ey0 = SART_nf$Ey)

(SSART_f    <- summary(SART_f, Glist = Gnetnorm, data = mydata))
save(SART_f, SSART_f, file = "_output/AH_fixed.rda")

## Count data
# Model 4
Rbh      <- 1
starting <- CDtmp[[Rbh]]$estimate; starting$Gamma <- c(rep(starting$Gamma[1], n.school), starting$Gamma[-1])
Ey0      <- CDtmp[[Rbh]]$Ey
SCD_f0   <- cdnet(formula = form.fix, Glist =  Gnetnorm, Rbar = Rbh, data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr0, cov = FALSE, 
                  starting = starting, Ey0 = Ey0)
SCD_f0   <- cdnet(formula = form.fix, Glist =  Gnetnorm, Rbar = Rbh, data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = TRUE, optimizer = "optim", 
                  starting = SCD_f0$estimate, Ey0 = SCD_f0$Ey)

# Model 5
Rbh      <- 13
starting <- CDtmp[[Rbh]]$estimate; starting$Gamma <- c(rep(starting$Gamma[1], n.school), starting$Gamma[-1])
Ey0      <- CDtmp[[Rbh]]$Ey
SCD_f    <- cdnet(formula = form.fix, Glist =  Gnetnorm, Rbar = Rbh, data = mydata, 
                npl.ctr = npl.ctr, opt.ctr = opt.ctr0, cov = FALSE, 
                starting = starting, Ey0 = Ey0)
SCD_f    <- cdnet(formula = form.fix, Glist =  Gnetnorm, Rbar = Rbh, data = mydata, 
                npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = TRUE, optimizer = "optim", 
                starting = SCD_f$estimate, Ey0 = SCD_f$Ey)

save(SSART_f, SCD_f0, SCD_f, file = "_output/AH_fixed.rda")

################ simulate data using these results
# (the simulations are used in file `D3_Add_Health_Endo.R` to construct Figure Figure B.1)
load(file = "_output/AH_fixed.rda")
set.seed(123)
sform.fix   <- as.formula(paste0(c("~ -1 + ", paste0(c(eff, expvc), collapse = " + ")), collapse = ""))
Rbar        <- 1
y0          <-  simcdnet(formula = sform.fix, Glist = Gnetnorm, parms = SCD_f0$estimate$parms, 
                         data = mydata, Rbar = Rbar, Rmax = 33)$y 

Rbar        <- 13
y1          <- simcdnet(formula = sform.fix, Glist = Gnetnorm, parms = SCD_f$estimate$parms, 
                        data = mydata, Rbar = Rbar, Rmax = 33)$y

simuAH_f    <- data.frame(model = "B",
                          y     = c(y0, y1),
                          Data  = rep(4:5, each = length(y0))) %>%
  mutate(Datacode = factor(Data, labels = c(expression(paste("(B) Without FE, ", bar(R), " =  1")), 
                                            expression(paste("(C) Without FE, ", bar(R), " = 13")))))
saveRDS(simuAH_f, file = "_output/simuB.RDS")