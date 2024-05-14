# The Bayesian estimation of the dyadic linking model (Online Appendix B.3.3)
rm(list = ls())
library(PartialNetwork)
library(CDatanet)
library(ggplot2)
library(dplyr)

# Data
load("MydataCount.rda")
obj            <- ls()
rm(list = obj[!(obj %in%c("Gnet", "form.net", "datanet"))])
gc()

form.net       <- as.formula(paste0(c("~ ", paste0(va.net, collapse = " + ")), collapse = ""))

set.seed(123)
out            <- homophily(network =  Gnet, formula = form.net, fixed.effects = TRUE,
                            iteration = 2e4, data = datanet)

save(out, file = "Net.RE.rda")

load("_output/Net.RE.rda")
set.seed(123)
index    <- sample(1:72291, 9)
tbeta    <- c("Age", "Same sex", "Hispanic", "White", "Black", "Asian", "Year in school", "Mom Job Prof.")

dataplot <- data.frame(iter = rep(1:nrow(out$posterior$beta), 20),
                       simu = c(out$posterior$beta[, 121:128],
                                out$posterior$mu[,index[1:5]],
                                out$posterior$nu[,index[6:9]],
                                out$posterior$sigma2_mu,
                                out$posterior$sigma2_nu,
                                out$posterior$rho),
                       parms = rep(1:20, each = nrow(out$posterior$beta))) %>% 
  mutate(parms = factor(parms, labels = c(expr(paste(psi, ": ", !!tbeta[1])),
                                          expr(paste(psi, ": ", !!tbeta[2])),
                                          expr(paste(psi, ": ", !!tbeta[3])),
                                          expr(paste(psi, ": ", !!tbeta[4])),
                                          expr(paste(psi, ": ", !!tbeta[5])),
                                          expr(paste(psi, ": ", !!tbeta[6])),
                                          expr(paste(psi, ": ", !!tbeta[7])),
                                          expr(paste(psi, ": ", !!tbeta[8])),
                                          expr(paste(mu[!!(index[1])])),
                                          expr(paste(mu[!!(index[2])])),
                                          expr(paste(mu[!!(index[3])])),
                                          expr(paste(mu[!!(index[4])])),
                                          expr(paste(mu[!!(index[5])])),
                                          expr(paste(nu[!!(index[6])])),
                                          expr(paste(nu[!!(index[7])])),
                                          expr(paste(nu[!!(index[8])])),
                                          expr(paste(nu[!!(index[9])])),
                                          expr(sigma[mu]^2), expr(sigma[nu]^2), expr(rho[paste(mu, ",", nu)])))) %>% 
  group_by(parms) %>% mutate(min = quantile(tail(simu, 1e3), 0.025), max = quantile(tail(simu, 1e3), 0.975)) %>%
  ungroup() %>% mutate(simu = ifelse(iter <= 1e3, NA, simu)) 


(graph <- ggplot(dataplot, aes(x = iter, y = simu)) +
    # geom_errorbar(width = .02, aes(ymin = m - qnorm(0.975)*s, ymax = m + qnorm(0.975)*s)) + 
    # geom_ribbon(aes(ymin =  m - qnorm(0.975)*s, ymax = m + qnorm(0.975)*s), 
    #             color = "#c90000", fill = "#d0ffff", linewidth = 0.3, linetype = 2) +
    geom_line(linewidth = 0.25, linetype = 1, color = "blue") + theme_bw() +
    geom_line(aes(x = iter, y = min), color = "red", linewidth = 0.5, linetype = 2) +
    geom_line(aes(x = iter, y = max), color = "red", linewidth = 0.5, linetype = 2) +
    xlab("") + ylab("") + facet_wrap(~ parms, ncol = 4, scales = "free", labeller = label_parsed) + 
    scale_x_continuous(breaks = seq(0, 2e4, 1e4)) + 
    theme(strip.text = element_text(face = "italic", margin = margin(-0.07, 0, -0.07, 0, "cm")), 
          axis.text.x = element_text(hjust = 1),
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 7, family = "Palatino")))


ggsave("plot:Bayes.pdf", path = "_output", plot = graph, device = "pdf", width = 8, height = 8)

# compute posterior means of mu and nu (they are used in file `D3_Add_Health_Endo` to address network endogeneity)
load("_output/Net.RE.rda")
gc()
nu <- colMeans(tail(out$posterior$nu, 1e4))
gc()
mu <- colMeans(tail(out$posterior$mu, 1e4))
gc()
save(mu, nu, file = "_output/munu.RE.rda")
rm(list = ls())
gc()
