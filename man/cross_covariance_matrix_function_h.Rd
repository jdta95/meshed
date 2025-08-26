\name{cross_covariance_matrix_function_h}
\alias{cross_covariance_matrix_function_h}
\title{Cross-covariance or cross-correlation matrix MCMC chain for models based on MGPs}
\description{Calculates the cross-covariance or cross-correlation matrix among outcomes at each MCMC iteration using the output from spmeshed.}
\usage{
cross_covariance_matrix_function_h(spmeshed_object, distance = 0, correl = FALSE, num_threads = 1)
}
\arguments{
    \item{spmeshed_object}{Object output from \code{spmeshed}.}
    \item{distance}{Distance at which to compute the cross-covariance/cross-correlation matrix.}
    \item{correl}{Logical indicating whether to compute the cross-correlation matrix instead of the cross-covariance matrix.}
    \item{num_threads}{Integer number of OpenMP threads. This is ineffective if \code{meshed} was not compiled with OpenMP support.}
}
\value{
    An array of size \eqn{(q, q, n_thin*n_samples)} containing the cross-covariance/cross-correlation matrix at each MCMC iteration The full chain minus burn-in is returned NOT thinned since \code{q} is relatively small.
    }

\references{
  Peruzzi, M., Banerjee, S., and Finley, A.O. (2022)
  Highly Scalable Bayesian Geostatistical Modeling via Meshed Gaussian Processes on Partitioned Domains. \emph{Journal of the American Statistical Association}, 117(538):969-982. \doi{10.1080/01621459.2020.1833889}
  
  Peruzzi, M., Dunson, D.B. (2024)
  Spatial meshing for general Bayesian multivariate models. \emph{Journal of Machine Learning Research}, 25(87):1−49. \url{https://arxiv.org/abs/2201.10080}
  
  Peruzzi, M., Banerjee, S., Dunson, D.B., and Finley, A.O.  (2025)
  Gridding and Parameter Expansion for Scalable Latent Gaussian Models of Spatial Multivariate Data. \emph{Bayesian Analysis} \doi{10.1214/25-BA1515}
}

\author{
  Michele Peruzzi \email{michele.peruzzi@umich.edu}
}

\examples{
# toy example with small dataset and short MCMC
# on a bivariate outcome
library(tidyverse)
library(meshed, lib.loc = "../temp")

set.seed(2021)

SS <- 12 # number of coordinates for the jth dimension
n <- SS ^ 2 # total n locations
q <- 2 # number of outcomes

# create coordinates with SS points in each dimension
coords <- expand.grid(
  x = seq(0, 1, length.out = SS),
  y = seq(0, 1, length.out = SS)
)

# generate data
sigmasq <- c(2.3, 1.5)
phi <- c(6, 4)
tausq <- c(0.1, 0.2)
Blist <- list(c(-1,.5,1), c(1,-.5,2))

distmat = as.matrix(dist(coords))
LClist = lapply(1:q, function(i) {
    t(chol(exp(-phi[i] * distmat)))
})
wlist = 1:q \%>\% lapply(function(i) sigmasq[i]^0.5 * LClist[[i]] \%*\% rnorm(n))
p <- length(Blist[[1]])
X <- rnorm(n * p) \%>\% matrix(ncol=p)
ylist <- 1:q \%>\% lapply(function(i) X \%*\% Blist[[i]] + wlist[[i]] + tausq[i]^0.5 * rnorm(n))
y <- do.call(cbind, ylist)

# MCMC setup
mcmc_keep <- 500
mcmc_burn <- 100
mcmc_thin <- 2

# run MCMC
meshout <- spmeshed(y, X, coords,
                    axis_partition=c(4,4),
                    n_samples = mcmc_keep, 
                    n_burn = mcmc_burn, 
                    n_thin = mcmc_thin, 
                    settings = list(cache=FALSE),
                    prior=list(phi=c(1,15)),
                    verbose = 0,
                    n_threads = 1)

# calculate cross-covariance matrix MCMC chain at distance 0
cross_covariance_matrix_function_h(meshout, distance = 0, correl = FALSE, num_threads = 1)

# calculate cross-correlation matrix MCMC chain at distance 0.1
cross_covariance_matrix_function_h(meshout, distance = 0.1, correl = TRUE, num_threads = 1)
}