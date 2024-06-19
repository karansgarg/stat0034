# Title: targets.R
# Description: R file containing different target distributions and models

# NOTE: Recall that for HMC, when the target distribution is p(x),
# the potential energy function is U := -log(p(x))

################################################################################
# Imports
################################################################################
library(mvtnorm)

################################################################################
# Simple 1d standard Gaussian case
################################################################################

# Potential energy
gaussian_1d_U <- function(x){
  return(0.5 * x**2)
}

# Gradient of potential energy
grad_gaussian_1d_U <- function(x){
  return(x)
}

################################################################################
# 10-d standard Gaussian case
################################################################################

# Potential energy
gaussian_nd_U <- function(x){
  return(-dmvnorm(x, log=T))
}

# Gradient of potential energy
grad_gaussian_nd_U <- function(x){
  return(c(-solve(diag(1, nrow=length(x), ncol=length(x))) %*% x))
}

################################################################################
# 1d Laplace tails case
################################################################################