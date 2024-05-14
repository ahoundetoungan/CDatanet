# The Logit estimation of the dyadic linking model.
rm(list = ls())
library(PartialNetwork)
library(CDatanet)
library(ggplot2)

# Data
proot <- c("~/tmpah/SchoolData",
           "~/Dropbox/Papers - In progress/CountDNtw/Code/Application",
           "~/CountDNtw/Code/Application/")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)
load("MydataCount.rda")

form.net       <- as.formula(paste0(c("~ -1 +", paste0(va.net, collapse = " + ")), collapse = ""))
init           <- list(beta = c(-0.7574687,  0.2630633,  1.6065785, -0.9774744, 3.3870634, 2.7455576, -0.5066828, 0.1571598),
                       mu = rep(0, nrow(mydata)),
                       nu = rep(0, nrow(mydata) - n.school))

rm(list = ls()[!(ls() %in% c("Gnet", "form.net", "datanet", "init"))])
gc()

out            <- homophily.FE(network =  Gnet, formula = form.net, data = datanet, init = init,
                               opt.ctr = list(maxit = 1e9, eps_f = 1e-20, eps_g = 1e-20))
 
save(out, file = "_output/Net.FE.rda")
mu             <- out$estimate$mu
nu             <- out$estimate$nu
save(mu, nu, file = "_output/munu.FE.rda")

