#' @title The CDatanet package
#' @description The \pkg{CDatanet} package implements the count data model with social interactions and the dyadic linking model developed in Houndetoungan (2022). 
#' It also simulates data from the count data model and implements the Spatial Autoregressive Tobit model (LeSage, 2000; Xu and Lee, 2015) for left censored data and the Spatial Autoregressive Model (Lee, 2004). 
#' Network formation models, such as that studied by Yan et al. (2019), are also implemented.
#' To make the computations faster \pkg{CDatanet} uses \code{C++} through the \pkg{Rcpp} package (Eddelbuettel et al., 2011). 
#'
#' @references 
#' Eddelbuettel, D., & François, R. (2011). \pkg{Rcpp}: Seamless \R and \code{C++} integration. \emph{Journal of Statistical Software}, 40(8), 1-18, \doi{10.18637/jss.v040.i08}.
#' @references 
#' Houndetoungan, E. A. (2022). Count Data Models with Social Interactions under Rational Expectations. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}.
#' @references  
#' Lee, L. F. (2004). Asymptotic distributions of quasi‐maximum likelihood estimators for spatial autoregressive models. \emph{Econometrica}, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}.
#' @references 
#' Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
#' @references 
#' Yan, T., Jiang, B., Fienberg, S. E., & Leng, C. (2019). Statistical inference in a directed network model with covariates. \emph{Journal of the American Statistical Association}, 114(526), 857-868, \doi{https://doi.org/10.1080/01621459.2018.1448829}.
#' 
#' @useDynLib CDatanet, .registration = TRUE
"_PACKAGE"
NULL