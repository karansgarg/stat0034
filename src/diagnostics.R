# Title: diagnostics.R
# Description: R file containing functions to evaluate output of MCMC algorithms

################################################################################
# Imports
################################################################################
library(coda)

################################################################################
# Finite difference method to evaluate derivatives
################################################################################
finite_diff <- function(target, # Function that we are deriving
                        grad_target, # Analytical derivative of function
                        epsilon=0.00000001, # Epsilon to use in finite diff method
                        eval_value=rep(0, 10)){ # Test value for function
  differences <- rep(NA, length(eval_value))
  deriv <- grad_target(eval_value)
  for(i in 1:length(eval_value)){
    val_plus <- eval_value
    val_neg <- eval_value
    val_plus[i] <- eval_value[i] + epsilon
    val_neg[i] <- eval_value[i] - epsilon
    numerical_approx <- (target(val_plus) - target(val_neg)) / (2*epsilon)
    diff <- abs(numerical_approx - deriv[i])
    differences[i] <- diff
  }
  return(differences)
}