#' @title The CDatanet package
#' @description The \pkg{CDatanet} package implements the count data model with social interactions and the dyadic linking model developed in Houndetoungan (2020). 
#' It also simulates data from the count data model and implements the Spatial Autoregressive Tobit model (LeSage, 2000; Xu and Lee, 2015) for left censored data and the Spatial Autoregressive Model (Lee, 2004). 
#' To make the computations faster \pkg{CDatanet} uses \code{C++} through the \pkg{Rcpp} package (Eddelbuettel et al., 2011). 
#'
#' @seealso \code{\link{simCDnet}}, \code{\link{CDnetNPL}}, \code{\link{SARML}}, \code{\link{SARTML}} and \code{\link{netformation}}.
#' 
#' @references 
#' Atchade, Y. F., & Rosenthal, J. S. (2005). On adaptive markov chain monte carlo algorithms, \emph{Bernoulli}, 11(5), 815-828, \doi{10.3150/bj/1130077595}
#' @references Eddelbuettel, D., François, R., Allaire, J., Ushey, K., Kou, Q., Russel, N., ... & Bates, D., 2011,
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration, \emph{Journal of Statistical Software}, 40(8), 1-18, \doi{10.18637/jss.v040.i08}
#' @references 
#' Houndetoungan, E. A., 2020, A Count Data Model with Social Interactions. Available at SSRN 3721250, \doi{10.2139/ssrn.3721250}
#' @references  
#' Lee, L. F., 2004, Asymptotic distributions of quasi-maximum likelihood estimators for spatial autoregressive models. Econometrica, 72(6), 1899-1925, \doi{10.1111/j.1468-0262.2004.00558.x}
#' @references  
#' LeSage, J. P., 2000, Bayesian estimation of limited dependent variable spatial autoregressive models, \emph{Geographical Analysis}, 32(1), 19-35, \doi{10.1111/j.1538-4632.2000.tb00413.x}
#' @references 
#' Xu, X., & Lee, L. F., 2015, Maximum likelihood estimation of a spatial autoregressive Tobit model, Journal of Econometrics, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}
#' 
#' @useDynLib CDatanet, .registration = TRUE
"_PACKAGE"
NULL