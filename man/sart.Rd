% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SARTestim.R
\name{sart}
\alias{sart}
\title{Estimate sart model}
\usage{
sart(
  formula,
  contextual,
  Glist,
  theta0 = NULL,
  yb0 = NULL,
  optimizer = "fastlbfgs",
  npl.ctr = list(),
  opt.ctr = list(),
  print = TRUE,
  cov = TRUE,
  RE = FALSE,
  data
)
}
\arguments{
\item{formula}{an object of class \link[stats]{formula}: a symbolic description of the model. The \code{formula} should be as for example \code{y ~ x1 + x2 | x1 + x2}
where \code{y} is the endogenous vector, the listed variables before the pipe, \code{x1}, \code{x2} are the individual exogenous variables and
the listed variables after the pipe, \code{x1}, \code{x2} are the contextual observable variables. Other formulas may be
\code{y ~ x1 + x2} for the model without contextual effects, \code{y ~ -1 + x1 + x2 | x1 + x2} for the model
without intercept or \code{ y ~ x1 + x2 | x2 + x3} to allow the contextual variable to be different from the individual variables.}

\item{contextual}{(optional) logical; if true, this means that all individual variables will be set as contextual variables. Set the
\code{formula} as \code{y ~ x1 + x2} and \code{contextual} as \code{TRUE} is equivalent to set the formula as \code{y ~ x1 + x2 | x1 + x2}.}

\item{Glist}{the adjacency matrix or list sub-adjacency matrix.}

\item{theta0}{(optional) starting value of \eqn{\theta = (\lambda, \beta, \gamma, \sigma)}. The parameter \eqn{\gamma} should be removed if the model
does not contain contextual effects (see details).}

\item{yb0}{(optional) expectation of y.}

\item{optimizer}{is either \code{fastlbfgs} (L-BFGS optimization method of the package \pkg{RcppNumerical}), \code{nlm} (referring to the function \link[stats]{nlm}), or \code{optim} (referring to the function \link[stats]{optim}).
Other arguments
of these functions such as, \code{control} and \code{method} can be defined through the argument \code{opt.ctr}.}

\item{npl.ctr}{list of controls for the NPL method (see \code{\link{cdnet}}).}

\item{opt.ctr}{list of arguments to be passed in \code{optim_lbfgs} of the package \pkg{RcppNumerical}, \link[stats]{nlm} or \link[stats]{optim} (the solver set in \code{optimizer}), such as \code{maxit}, \code{eps_f}, \code{eps_g}, \code{control}, \code{method}, ...}

\item{print}{a Boolean indicating if the estimate should be printed at each step.}

\item{cov}{a Boolean indicating if the covariance should be computed.}

\item{RE}{a Boolean which indicates if the model if under rational expectation of not.}

\item{data}{an optional data frame, list or environment (or object coercible by \link[base]{as.data.frame} to a data frame) containing the variables
in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{sart} is called.}
}
\value{
A list consisting of:
\item{info}{list of general information on the model.}
\item{estimate}{Maximum Likelihood (ML) estimator.}
\item{yb}{ybar (see details), expectation of y.}
\item{Gyb}{average of the expectation of y among friends.}
\item{cov}{List of covariances.}
\item{details}{outputs as returned by the optimizer.}
}
\description{
\code{sart} is used to estimate peer effects on censored data (see details). The model is presented in Xu and Lee(2015).
}
\details{
\subsection{Model}{

The left-censored variable \eqn{\mathbf{y}}{y} is generated from a latent variable \eqn{\mathbf{y}^*}{ys}.
The latent variable is given for all i as
\deqn{y_i^* = \lambda \mathbf{g}_i y + \mathbf{x}_i'\beta + \mathbf{g}_i\mathbf{X}\gamma + \epsilon_i,}{ys_i = \lambda g_i*y + x_i'\beta + g_i*X\gamma + \epsilon_i,}
where \eqn{\epsilon_i \sim N(0, \sigma^2)}{\epsilon_i --> N(0, \sigma^2)}.\cr
The count variable \eqn{y_i} is then define that is \eqn{y_i = 0} if
\eqn{y_i^* \leq 0}{ys_i \le 0} and \eqn{y_i = y_i^*}{y_i = ys_i} otherwise.
}
}
\examples{
\donttest{
# Groups' size
M      <- 5 # Number of sub-groups
nvec   <- round(runif(M, 100, 1000))
n      <- sum(nvec)

# Parameters
lambda <- 0.4
beta   <- c(2, -1.9, 0.8)
gamma  <- c(1.5, -1.2)
sigma  <- 1.5
theta  <- c(lambda, beta, gamma, sigma)

# X
X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))

# Network
Glist  <- list()

for (m in 1:M) {
  nm           <- nvec[m]
  Gm           <- matrix(0, nm, nm)
  max_d        <- 30
  for (i in 1:nm) {
    tmp        <- sample((1:nm)[-i], sample(0:max_d, 1))
    Gm[i, tmp] <- 1
  }
  rs           <- rowSums(Gm); rs[rs == 0] <- 1
  Gm           <- Gm/rs
  Glist[[m]]   <- Gm
}


# data
data    <- data.frame(x1 = X[,1], x2 =  X[,2])

rm(list = ls()[!(ls() \%in\% c("Glist", "data", "theta"))])

ytmp    <- simsart(formula = ~ x1 + x2 | x1 + x2, Glist = Glist,
                   theta = theta, data = data)

y       <- ytmp$y

# plot histogram
hist(y)

opt.ctr <- list(method  = "Nelder-Mead", 
                control = list(abstol = 1e-16, abstol = 1e-11, maxit = 5e3))
data    <- data.frame(yt = y, x1 = data$x1, x2 = data$x2)
rm(list = ls()[!(ls() \%in\% c("Glist", "data"))])

out     <- sart(formula = yt ~ x1 + x2, optimizer = "nlm",
                  contextual = TRUE, Glist = Glist, data = data)
summary(out)
}
}
\references{
Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
}
\seealso{
\code{\link{sar}}, \code{\link{cdnet}}, \code{\link{simsart}}.
}
