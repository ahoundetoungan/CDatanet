% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sart.R
\name{simsart}
\alias{simsart}
\title{Simulating Data from Tobit Models with Social Interactions}
\usage{
simsart(formula, Glist, theta, tol = 1e-15, maxit = 500, cinfo = TRUE, data)
}
\arguments{
\item{formula}{a class object \code{\link[stats]{formula}}: a symbolic description of the model.
\code{formula} must be, for example, \code{y ~ x1 + x2 + gx1 + gx2}, where \code{y} is the endogenous vector,
and \code{x1}, \code{x2}, \code{gx1}, and \code{gx2} are control variables. These can include contextual variables,
i.e., averages among the peers. Peer averages can be computed using the function \code{\link{peer.avg}}.}

\item{Glist}{The network matrix. For networks consisting of multiple subnets, \code{Glist} can be a list
of subnets with the \code{m}-th element being an \code{ns*ns} adjacency matrix, where \code{ns} is the number of nodes
in the \code{m}-th subnet.}

\item{theta}{a vector defining the true value of \eqn{\theta = (\lambda, \Gamma, \sigma)} (see the model specification in the details).}

\item{tol}{the tolerance value used in the fixed-point iteration method to compute \code{y}. The process stops
if the \eqn{\ell_1}-distance between two consecutive values of \code{y} is less than \code{tol}.}

\item{maxit}{the maximum number of iterations in the fixed-point iteration method.}

\item{cinfo}{a Boolean indicating whether information is complete (\code{cinfo = TRUE}) or incomplete (\code{cinfo = FALSE}).
In the case of incomplete information, the model is defined under rational expectations.}

\item{data}{an optional data frame, list, or environment (or object coercible by \code{\link[base]{as.data.frame}}
to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken
from \code{environment(formula)}, typically the environment from which \code{simsart} is called.}
}
\value{
A list consisting of:
\describe{
\item{yst}{\eqn{y^{\ast}}, the latent variable.}
\item{y}{The observed censored variable.}
\item{Ey}{\eqn{E(y)}, the expected value of \eqn{y}.}
\item{Gy}{The average of \eqn{y} among peers.}
\item{GEy}{The average of \eqn{E(y)} among peers.}
\item{meff}{A list including average and individual marginal effects.}
\item{iteration}{The number of iterations performed per sub-network in the fixed-point iteration method.}
}
}
\description{
\code{simsart} simulates censored data with social interactions (see Xu and Lee, 2015).
}
\details{
For a complete information model, the outcome \eqn{y_i} is defined as:
\deqn{\begin{cases}
y_i^{\ast} = \lambda \bar{y}_i + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
y_i = \max(0, y_i^{\ast}),
\end{cases}}
where \eqn{\bar{y}_i} is the average of \eqn{y} among peers,
\eqn{\mathbf{z}_i} is a vector of control variables,
and \eqn{\epsilon_i \sim N(0, \sigma^2)}. \cr

In the case of incomplete information models with rational expectations, \eqn{y_i} is defined as:
\deqn{\begin{cases}
y_i^{\ast} = \lambda E(\bar{y}_i) + \mathbf{z}_i'\Gamma + \epsilon_i, \\ 
y_i = \max(0, y_i^{\ast}).
\end{cases}}
}
\examples{
\donttest{
# Define group sizes
set.seed(123)
M      <- 5 # Number of sub-groups
nvec   <- round(runif(M, 100, 200)) # Number of nodes per sub-group
n      <- sum(nvec) # Total number of nodes

# Define parameters
lambda <- 0.4
Gamma  <- c(2, -1.9, 0.8, 1.5, -1.2)
sigma  <- 1.5
theta  <- c(lambda, Gamma, sigma)

# Generate covariates (X)
X      <- cbind(rnorm(n, 1, 1), rexp(n, 0.4))

# Construct network adjacency matrices
G      <- list()
for (m in 1:M) {
  nm           <- nvec[m] # Nodes in sub-group m
  Gm           <- matrix(0, nm, nm) # Initialize adjacency matrix
  max_d        <- 30 # Maximum degree
  for (i in 1:nm) {
    tmp        <- sample((1:nm)[-i], sample(0:max_d, 1)) # Random connections
    Gm[i, tmp] <- 1
  }
  rs           <- rowSums(Gm) # Normalize rows
  rs[rs == 0]  <- 1
  Gm           <- Gm / rs
  G[[m]]       <- Gm
}

# Prepare data
data   <- data.frame(X, peer.avg(G, cbind(x1 = X[, 1], x2 = X[, 2])))
colnames(data) <- c("x1", "x2", "gx1", "gx2") # Add column names

# Complete information game simulation
ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, 
                   Glist = G, theta = theta, 
                   data = data, cinfo = TRUE)
data$yc <- ytmp$y # Add simulated outcome to the dataset

# Incomplete information game simulation
ytmp    <- simsart(formula = ~ x1 + x2 + gx1 + gx2, 
                   Glist = G, theta = theta, 
                   data = data, cinfo = FALSE)
data$yi <- ytmp$y # Add simulated outcome to the dataset
}
}
\references{
Xu, X., & Lee, L. F. (2015). Maximum likelihood estimation of a spatial autoregressive Tobit model. \emph{Journal of Econometrics}, 188(1), 264-280, \doi{10.1016/j.jeconom.2015.05.004}.
}
\seealso{
\code{\link{sart}}, \code{\link{simsar}}, \code{\link{simcdnet}}.
}
