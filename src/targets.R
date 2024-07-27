# Title: targets.R
# Description: R file containing different target distributions and models

# NOTE: Recall that for HMC, when the target distribution is p(x),
# the potential energy function is U := -log(p(x))

################################################################################
# Imports
################################################################################
library(mvtnorm) # For drawing samples from multivariate normal
library(posteriordb) # For accessing different target distributions
library(rstan) # Has useful functions to interface with Stan (for posteriordb)

################################################################################
# Simple 1d standard Gaussian case
################################################################################

# Potential energy
# gaussian_1d_U <- function(x){
#   return(0.5 * x**2)
# }
# 
# # Gradient of potential energy
# grad_gaussian_1d_U <- function(x){
#   return(x)
# }

################################################################################
# N-d standard Gaussian case
################################################################################

# Potential energy
gaussian_nd_U <- function(x){return(-dmvnorm(x, log=T))}
 
# Gradient of potential energy
grad_gaussian_nd_U <- function(x){return(x)}

################################################################################
# Ill-conditioned Gaussian based on NUTS paper (variance sampled from Wishart)
################################################################################
# Set number of dimensions
illgaussian_d <- 20

# Generate covariance matrix from Wishart distribution
set.seed(1234) # Set seed to ensure same target produced everytime
A <- rWishart(n=1, df=illgaussian_d, Sigma=diag(illgaussian_d))
A <- matrix(A, nrow=illgaussian_d, ncol=illgaussian_d)

# Potential energy
illgaussian_U <- function(x){return(0.5 * x %*% A %*% x)}

# Gradient of potential energy
grad_illgaussian_U <- function(x){return(c(A %*% x))}

# Function to plot marginal density of first dimension
illgaussian_d1_density <- function(x){
  return(dnorm(x, mean=0, sd=sqrt(solve(A)[1,1])))
}
################################################################################
# Neal's funnel
################################################################################
# NOTE: First coordinate corresponds to y

# Potential energy
funnel_U <- function(x){
  d <- length(x)
  y <- x[1]
  x <- x[-1]
  neglogp <- 0.5*log(18*pi) + (d-1)/2 *log(2*pi) +
    (d-1)/2 * y +
    (y^2) / 18 +
    sum(x^2) / (2*exp(y))
  return(neglogp)
}

# Gradient of potential energy
funnel_grad_U <- function(x){
  d <- length(x)
  y <- x[1]
  x <- x[-1]
  grad <- rep(NA, d)
  grad[1] <- (d-1)/2 + y/9 - sum(x^2)/(2*exp(y))
  grad[-1] <- x/exp(y)
  return(grad)
}

# Sampling from funnel distribution
funnel_sampler <- function(n=1, d=10){
  samples <- matrix(NA, nrow=n, ncol=d)
  for(i in 1:n){
    samples[i,1] <- rnorm(1, mean=0, sd=3)
    samples[i,-1] <- rnorm(d-1, mean=0, sd=exp(samples[i,1]/2))
  }
  return(samples)
}

# Function to plot first dimension of Neal's funnel, y ~ N(0, 9)
funnel_y_density <- function(x){return(dnorm(x, mean=0, sd=3))}

################################################################################
# Posteriordb setup
################################################################################
# Connect to posteriordb database
my_pdb <- pdb_github()