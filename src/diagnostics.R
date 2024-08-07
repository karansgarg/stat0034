# Title: diagnostics.R
# Description: R file containing functions to evaluate output of MCMC algorithms

################################################################################
# Imports
################################################################################
#library(coda)

################################################################################
# Finite difference method to evaluate derivatives
################################################################################
finite_diff <- function(target, # Function that we are deriving
                        grad_target, # Analytical derivative of function
                        epsilon=0.00000001, # Epsilon to use in finite diff method
                        eval_value=rep(0, 10)){ # Test value for function
  # Vector to store differences
  differences <- rep(NA, length(eval_value))
  
  # Gradient vector according to supplied function
  deriv <- grad_target(eval_value)
  
  # For each dimension numerically approximate gradient to compute gradient vector
  # And then compare against supplied function to get a difference
  for(i in 1:length(eval_value)){
    val_plus <- eval_value
    val_neg <- eval_value
    val_plus[i] <- eval_value[i] + epsilon
    val_neg[i] <- eval_value[i] - epsilon
    numerical_approx <- (target(val_plus) - target(val_neg)) / (2*epsilon)
    diff <- abs(numerical_approx - deriv[i])
    differences[i] <- diff
  }
  # Return vector of differences
  return(differences)
}

################################################################################
# Collect HMC output into singular matrix
################################################################################
collect_data <- function(hmc_output){
  
  # Store output into single matrix
  data <- cbind(hmc_output$chain,
                hmc_output$acceptance,
                hmc_output$divergences)
  # Add U-turn-length if GIST sampler used
  if(any(!is.na(hmc_output$U_turn_length))){
    data <- cbind(data, hmc_output$U_turn_length)
  }
  
  print("Collecting data from HMC output...")
  
  # Convert to data frame
  data <- data.frame(data)
  
  # Extract dimension of target distribution
  d <- ncol(hmc_output$chain)
  
  # Name columns
  for(i in 1:d){
    name <- paste0("dim", i)
    colnames(data)[i] <- c(name)
  }

  colnames(data)[d+1] <- "acceptance"
  colnames(data)[d+2] <- "divergence"
  # Add U-turn-length if GIST sampler used
  if(any(!is.na(hmc_output$U_turn_length))){
    colnames(data)[d+3] <- "U_turn_length"
  }
  
  print("Data collection complete!")
  return(data)
}

################################################################################
# Output acceptance and divergence statistics
################################################################################
accept_diverge_stats <- function(hmc_output_data){
  
  # Calculate and print acceptance and divergence rates
  cat("Acceptance rate: ", mean(hmc_output_data[,"acceptance"]), "\n",
      "Divergence rate: ", mean(hmc_output_data[,"divergence"]==1), "\n",
      "Number of divergences: ", sum(hmc_output_data[,"divergence"]==1), "\n",
      "Divergence conditional acceptance rate: ",
      mean(hmc_output_data[hmc_output_data[,"divergence"]==1,"acceptance"]), "\n")
}

################################################################################
# Compute rolling bias of MCMC estimates
################################################################################
# rolling_bias_old <- function(samples_vector, true_mean=0, true_sd=3){
#   bias_vector <- rep(NA, length(samples_vector))
#   for(i in 1:length(samples_vector)){
#     bias_vector[i] <- mean(samples_vector[1:i])
#   }
#   bias_vector <- abs(bias_vector - true_mean) / true_sd
#   return(bias_vector)
# }

rolling_bias <- function(samples_vector, true_mean=0){
  bias_vector <- rep(NA, length(samples_vector))
  for(i in 1:length(samples_vector)){
    bias_vector[i] <- mean(samples_vector[1:i])
  }
  bias_vector <- abs(bias_vector - true_mean)
  return(bias_vector)
}

rolling_bias_sq <- function(samples_vector, true_mean=0, true_sd=3){
  bias_vector <- rep(NA, length(samples_vector))
  samples_vector <- samples_vector^2
  for(i in 1:length(samples_vector)){
    bias_vector[i] <- mean(samples_vector[1:i])
  }
  bias_vector <- abs(bias_vector - (true_sd^2 + true_mean^2))
  return(bias_vector)
}

################################################################################
# Plot rolling bias
################################################################################