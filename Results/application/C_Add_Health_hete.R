# Estimations with heterogeneity in peer effects without addressing network endogeneity (Model 7)
rm(list = ls())
library(CDatanet)
library(ggplot2)
library(dplyr)
setwd("")

# Data
load("MydataCount.rda")

################ Without fixed effects (This is not presented in the paper)
# Settings
form.nofix   <- as.formula(paste0(c("nclubs ~ ", paste0(expvc, collapse = " + ")), collapse = ""))
opt.ctr1     <- list(control = list(abstol = 1e-13, reltol = 1e-13, maxit = 5e3),
                     method  = "Nelder-Mead")
npl.ctr      <- list(maxit   = 1e4, tol = 5e-4, print = TRUE)

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

# Estimation
cont  <- TRUE
BIC   <- Inf
Rbh   <- 0
CDtmp <- vector("list", max(mydata$nclubs))
while (cont) {
  Rbh          <- Rbh + 1
  starting     <- NULL
  Ey0          <- NULL
  if(Rbh > 1){
    starting   <- CDtmp[[Rbh - 1]]$estimate
    starting$delta <- starting$delta[c(1:(Rbh - 1), Rbh - 1, Rbh:(2*Rbh - 2), 2*Rbh - 2)]
    Ey0        <- CDtmp[[Rbh - 1]]$Ey
  }
  CDtmp[[Rbh]] <- cdnet(formula = form.nofix, Glist = GHet, Rbar = rep(Rbh, 2), data = mydata, 
                        npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = FALSE, optimizer = "optim", 
                        starting = starting, Ey0 = Ey0, group = group)
  cont         <- (BIC > CDtmp[[Rbh]]$info$BIC)
  BIC          <- CDtmp[[Rbh]]$info$BIC
}

# summary with Rbar = 0 and optimal Rbar
(SCD0_hnf      <- meffects(CDtmp[[1]], Glist = GHet, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm))
(SCD1_hnf      <- meffects(CDtmp[[Rbh - 1]], Glist = GHet, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm))
(SCD2_hnf      <- meffects(CDtmp[[Rbh]], Glist = GHet, data = mydata, cont.var = cont.var, 
                           bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                           Glist.contextual = Gnetnorm))


save(CDtmp, SCD0_hnf, SCD1_hnf, SCD2_hnf, Rbh, file = "_output/AH_hete.rda")

################ With fixed effects (Model 7)
load("_output/AH_hete.rda")
# Settings
form.fix     <- as.formula(paste0(c("nclubs ~ -1 + ", paste0(c(eff, expvc), collapse = " + ")), collapse = ""))
opt.ctr0     <- list(maxit = 5e3, eps_f = 1e-10, eps_g = 1e-10)
opt.ctr1     <- list(control = list(abstol = 1e-13, reltol = 1e-13, maxit = 5e3),
                     method  = "Nelder-Mead")
npl.ctr      <- list(maxit   = 1e4, tol = 5e-4, print = TRUE)

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

Rbh      <- 1
starting <- SCD0_hnf$estimate; starting$Gamma <- c(rep(starting$Gamma[1], n.school), starting$Gamma[-1])
Ey0      <- CDtmp[[Rbh]]$Ey

SCD0_hf  <- cdnet(formula = form.fix, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr0, cov = FALSE, group = group, 
                  starting = starting, Ey0 = Ey0)
SCD0_hf  <- cdnet(formula = form.fix, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = TRUE, group = group,
                  optimizer = "optim", starting = SCD0_hf$estimate, Ey0 = SCD0_hf$Ey)

(SCD0_hf <- meffects(SCD0_hf, Glist = GHet, data = mydata, cont.var = cont.var, 
                     bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                     Glist.contextual = Gnetnorm))

Rbh      <- 10
starting <- SCD1_hnf$estimate; starting$Gamma <- c(rep(starting$Gamma[1], n.school), starting$Gamma[-1])
Ey0      <- CDtmp[[Rbh]]$Ey

SCD_hf   <- cdnet(formula = form.fix, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr0, cov = FALSE, group = group, 
                  starting = starting, Ey0 = Ey0)
SCD_hf   <- cdnet(formula = form.fix, Glist =  GHet, Rbar = rep(Rbh, 2), data = mydata, 
                  npl.ctr = npl.ctr, opt.ctr = opt.ctr1, cov = TRUE, group = group,
                  optimizer = "optim", starting = SCD_hf$estimate, Ey0 = SCD_hf$Ey)

(SCD_hf  <- meffects(SCD_hf, Glist = GHet, data = mydata, cont.var = cont.var, 
                     bin.var = bin.var, type.var = type.var, boot = 100, progress = TRUE,
                     Glist.contextual = Gnetnorm))

save(CDtmp, SCD0_hnf, SCD1_hnf, SCD2_hnf, SCD0_hf, SCD_hf, file = "_output/AH_hete.rda")


################ simulate data using these results
# (the simulations are used in file `D3_Add_Health_Endo.R` to construct Figure Figure B.1)
load(file = "_output/AH_hete.rda")
set.seed(123)
sform.fix   <- as.formula(paste0(c("~ -1 + ", paste0(c(eff, expvc), collapse = " + ")), collapse = ""))
Rbar        <- 10
y           <- simcdnet(formula = sform.fix, Glist = GHet, parms = SCD_hf$estimate$parms, 
                        data = mydata, Rbar = Rbar, Rmax = 33, group = group)$y

simuAH_hete <- data.frame(model = "C", y = y, Data  = 6) %>%
  mutate(Datacode = factor(Data, labels = c(expression(paste("(E) Het with FE, ", bar(R), " = 10")))))
saveRDS(simuAH_hete, file = "_output/simuC.RDS")