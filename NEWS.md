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