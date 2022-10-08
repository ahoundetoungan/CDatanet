# Changes in version 1.0.1
## The new version follows the major revision of the paper in April 2022
- The count data model has changed. It is now possible to have several gamma parameters instead of gamma = 1 as assumed in the previous version.
- In the count data model identification, sigma is now set to 1, because there are several gamma's.
- The network formation model now includes two unobserved effects instead of one.

## Some function name has changed
- `SARML` has been replaced by `sar`.
- `SARTML` has been replaced by `sart`.
- `CDnetNPL` has been replaced by `cdnet`.
- `simSARnet` has been replaced by `simsar`.
- `simTobitnet` has been replaced by `simsart`.
- `simCDnet` has been replaced by `simcdnet`.
- `netformation` has been replaced by `homophily`.

## SART model under rational expectations
It is now possible to estimate the SART model under rational expectations. In the previous version, the SART model is only available under complete information.

# Changes in version 2.0.1
This version follows the major revision of the paper in September 2022. 
- The count data model includes a more flexible specification. Especially, it is possible to assume that the cut points are not equally spaced for large values of the dependent variable. 
- I also implement a network formation model with degree heterogeneity as fixed effects (see [Yan et al., 2019](https://doi.org/10.1080/01621459.2018.1448829)).
- Models under incomplete information are now estimated using LBFGS algorithm of the package RcppNumerical. Thus, the optimization is performed in C++ and is very fast compared to the version 1.0.1.