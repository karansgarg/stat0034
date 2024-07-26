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

################################################################################
# Collect HMC output into singular matrix
################################################################################
collect_data <- function(hmc_output){
  #XXX
}

################################################################################
# Output acceptance and divergence statistics
################################################################################
accept_diverge_stats <- function(hmc_output_data, accept_col=11, diverge_col=13){
  
  cat("Acceptance rate: ", mean(hmc_output_data[,accept_col]), "\n",
      "Divergence rate: ", mean(hmc_output_data[,diverge_col]==1), "\n",
      "Number of divergences: ", sum(hmc_output_data[,diverge_col]==1), "\n",
      "Divergence conditional acceptance rate: ",
      mean(hmc_output_data[hmc_output_data[,diverge_col]==1,accept_col]), "\n")
}

################################################################################
# Compute rolling bias of MCMC estimates
################################################################################
rolling_bias <- function(samples_vector, true_mean, true_sd){
  bias_vector <- rep(NA, length(samples_vector))
  for(i in 1:length(samples_vector)){
    bias_vector[i] <- mean(samples_vector[1:i])
  }
  bias_vector <- (bias_vector - true_mean) / true_sd
  return(bias_vector)
}