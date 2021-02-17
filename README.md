# An R package for estimating a Count Data Model With Social Interactions
The **CDatanet** package implements a count data model with social interactions and a dyadic linking model to deal with the endogeneity of the social network. It includes all functions for the replication of the results in Houndetoungan (2020). The exact replication codes of the Monte Carlo study are located in the folder [**test**](https://github.com/ahoundetoungan/CDatanet/tree/master/test).

## Installation
### Requirements
- **CDatanet** package needs [**R**](https://cran.r-project.org/) version 3.0.0 or later which can be installed on Linux, Mac and Windows. See [**CRAN**](https://cran.r-project.org/) for installation details.
- [**devtools**](https://cran.r-project.org/package=devtools) package should be installed on [**R**](https://cran.r-project.org/). If not already done, install [**devtools**](https://cran.r-project.org/package=devtools) using the code ` install.packages("devtools") `.
- (*Only for windows users*) Windows users should install  [**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) compatible with their [**R**](https://cran.r-project.org/) version.

### How to install
**CDatanet** package can be installed from this GitHub repos using the `install_github` function of the [**devtools**](https://cran.r-project.org/package=devtools) package. All the dependencies will also be installed automatically.
```R
library(devtools)
install_github("ahoundetoungan/CDatanet")
```
The option `build_vignettes = TRUE` can be added if one desires to install the vignettes (this will take several minutes).
### Load CDatanet
Once the installation is done, **CDatanet** can be loaded as a common package in [**R**](https://cran.r-project.org/).
```R
library(CDatanet)
```
### How to use CDatanet
I provide a [vignette](https://nbviewer.jupyter.org/github/ahoundetoungan/CDatanet/blob/master/doc/CDatanet.pdf) that fully documents the package through simple and practical examples. Moreover, each function of the package has a help page accessible using `?name-of-the-function`. 