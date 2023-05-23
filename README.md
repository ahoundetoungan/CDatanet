# An R package for estimating a Count Data Model With Social Interactions
The **CDatanet** package implements a count data model with social interactions and a dyadic linking model to deal with the endogeneity of the social network. It includes all functions for the replication of the results in Houndetoungan (2022). The exact replication codes of the Monte Carlo study are located in the folder [**test**](https://github.com/ahoundetoungan/CDatanet/tree/master/test).

See CRAN version [here](https://CRAN.R-project.org/package=CDatanet).

## Installation
### CRAN version
**CDatanet** can be directly installed from CRAN.
```R
install.packages("CDatanet")
```

### GitHub version
It may be possible that I updated the package without submitting the new version to CRAN. The latest version (*but not necessary stable*) of **CDatanet** can be installed from this GitHub repos using the `install_github` function of the [**devtools**](https://cran.r-project.org/package=devtools) package. All the dependencies will also be installed automatically.
```R
library(remotes)
install_github("ahoundetoungan/CDatanet", build_vignettes = TRUE)
```
