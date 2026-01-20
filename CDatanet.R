rm(list = ls())
library(RcppArmadillo)
library(Rcpp)
library(roxygen2)
library(devtools)
library(pkgbuild)
library(usethis)

setwd("~/Dropbox/Academy/1.Papers/CountDNtw/Code/Package/CDatanet")

roxygenise()
compile_dll(force = TRUE)
roxygenise()
compileAttributes()
use_build_ignore(c("doc", "Meta", "test", ".gitignore", ".Rbuildignore", ".Rhistory",
                   "QuantilePeer.Rproj, README.md",  "vignettes/.gitignore", "scr/*.so", "scr/*.o",
                   ".tar.gz", "vignettes/.pdf"))
document()

## tests
# usethis::use_testthat(3)


# Build and install
system("R CMD build --compact-vignettes='both' .") 


system("R CMD check --as-cran CDatanet_2.2.3.tar.gz")

system("R CMD INSTALL CDatanet_2.2.3.tar.gz")

# sum(cranlogs::cran_downloads(packages = "CDatanet", from = format(Sys.Date() - 365, "%Y-%m-%d"))[,2])
# system("R CMD check --as-cran CDatanet_2.2.2.tar.gz")
#tools::showNonASCIIfile
#system("R CMD INSTALL Rd2pdf CDatanet --no-clean")
# check_rhub("CDatanet")
# check_win_devel(pkg = "CDatanet")

# knitr::write_bib(c('knitr', 'Rcpp',  'Formula', 
#                    'formula.tools', 'ddpcr', 'Matrix',
#                    'RcppArmadillo', 'RcppProgress',
#                    'microbenchmark', 'ggplot2', 'devtools',
#                    'doParallel'), 
#                  'Documentation/Packages.bib', width = 60)
