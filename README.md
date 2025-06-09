## An R package for econometrics of network data

<!-- badges: start -->
  [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![R-CMD-check](https://github.com/ahoundetoungan/CDatanet/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ahoundetoungan/CDatanet/actions/workflows/R-CMD-check.yml)
  
  [![R-universe](https://ahoundetoungan.r-universe.dev/badges/CDatanet)](https://ahoundetoungan.r-universe.dev/CDatanet)
  [![CRAN](https://www.r-pkg.org/badges/version/CDatanet)](https://CRAN.R-project.org/package=CDatanet)
  [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/CDatanet)](https://cran.r-project.org/package=CDatanet)
   
<!-- badges: end -->

The **CDatanet** package simulates and estimates peer effect models and network formation models. The class of peer effect models includes the linear-in-means models ([Lee, 2004](https://doi.org/10.1111/j.1468-0262.2004.00558.x); [Lee et al., 2010](https://doi.org/10.1111/j.1368-423X.2010.00310.x)),  the Tobit models ([Xu and Lee, 2015](https://doi.org/10.1016/j.jeconom.2015.05.004)), and discrete numerical data models ([Houndetoungan, 2024](https://doi.org/10.2139/ssrn.3721250)).  The network formation models include pair-wise regressions with degree heterogeneity ([Graham, 2017](https://doi.org/10.3982/ECTA12679); [Yan et al., 2019](https://doi.org/10.1080/01621459.2018.1448829)) and exponential random graph models ([Mele, 2017](https://doi.org/10.3982/ECTA10400)).

**CDatanet** package includes all functions for the replication of the results in [Houndetoungan (2024)](https://doi.org/10.2139/ssrn.3721250). The exact replication codes are located in the folder [**Results**](https://github.com/ahoundetoungan/CDatanet/tree/master/Results).

See CRAN version [here](https://CRAN.R-project.org/package=CDatanet).

## Installation
### CRAN version
**CDatanet** can be directly installed from CRAN.
```R
install.packages("CDatanet")
```

### GitHub version
It may be possible that I updated the package without submitting the new version to CRAN. The latest version (*but not necessary stable*) of **CDatanet** can be installed from this GitHub repos.
```R
library(remotes)
install_github("ahoundetoungan/CDatanet")
```
