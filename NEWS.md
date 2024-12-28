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

# Changes in versions 2.0.2 and 2.0.3
Note and Warning found in the check for MACOS have been fixed

# Changes in version 2.1.0
R defaulted to C++11 in R 4.0.0, to C++14 in R 4.2.0 and to C++17.

# Changes in version 2.1.1
Fixed effect is allowed in the model SAR.

# Changes in version 2.1.2
- I address the problem of single-agent subnetwork.
- AIC and BIC are added to the output of cdnet. They can be used to choose Rbar.

# Changes in version 2.1.3
- `homophily` has been changed to `homophily.re` for the random effect models.
- `homophily.FE` has heen changed to `homophily.fe`.
- Random effects in `homophily.re` can be one side or two sides. 
- Fixed effects in `homophily.fe` can be one side or two sides.
- Symmetric network models are included in `homophily.re`.
- Symmetric network models are included in `homophily.fe`.
- The function `homophili.data` is added to convert data between directed network models and symmetric network models.

# Changes in version 2.2.0
- Following the revision of the paper Count Data Model with Social Interaction under Rational Expectations, `cdnet` now allows heterogeneity in peer effects. For example peer effects can be estimated between Blacks and Whites, Girls and Boys etc.
- `cdnet` allows heterogeneity in the cut-points. The cost function in the utility function is nonparametric and depends on some characteristics X such as gender, race, etc.

# Changes in version 2.2.1
- The Newton-Raphson method has been added for estimating the homophily model.
