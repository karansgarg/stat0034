################################################################################
# HMC sampler (old, as of 10/07/24) # NOTE: Bug found in this (current_q not updated)
################################################################################
HMC <- function(U, grad_U, # Potential energy and its gradient
                K='gaussian', # Kinetic energy (can be 'gaussian' or 'laplace')
                current_q, # Starting point of chain
                L, # Path length
                epsilon=NULL,
                epsilon_dist='exponential',
                epsilon_params=list(rate=1),
                chain_length=5000,
                warmup=5000){
  
  # Extract dimensionality of target
  d <- length(current_q)
  
  # Assign data structure for chain
  chain <- matrix(NA, nrow=chain_length, ncol=d)
  
  # Set parameters/functions for kinetic energy
  if(K == 'gaussian'){
    sample_p <- function(x){return(rnorm(x))}
    K <- function(x){return(sum(x^2)/2)}
    grad_K <- function(x){return(x)}
  }
  else if(K == 'laplace'){
    sample_p <- function(x){return(rlaplace(x))}
    K <- function(x){return(sum(abs(x)))}
    grad_K <- function(x){return(sign(x))}
  }
  
  # Perform HMC steps in loop to generate chain
  for(t in 1:chain_length){
    
    # Draw value for epsilon if distribution specified
    if(epsilon_dist=='exponential'){
      epsilon <- rexp(1, epsilon_params)
    }
    
    # Assign variable for updated position variables
    q <- current_q
    
    # Sample fresh momentum
    p <- sample_p(d) # Sample p based on given kinetic energy
    current_p <- p # Set initial p based on what was just sampled
    
    # Make a half-step for momentum at the beginning
    p <- p - epsilon*grad_U(q)/2
    
    # Alternate full-steps for position and momentum
    for(i in 1:L){
      # Full step for position variables
      q <- q + epsilon*grad_K(p)
      
      # Full step for momentum, except at the end
      if(i!=L){p <- p - epsilon*grad_U(q)}
    }
    
    # Final half step for momentum
    p <- p - epsilon*grad_U(q)/2
    
    # Negate momentum to get symmetric proposal
    p <- -p
    
    # Evaluate kinetic energies at start and end of trajectory
    current_K <- K(current_p)
    proposed_K <- K(p)
    
    # Evaluate potential energy at end of trajectory (and start for 1st iteration)
    if(t==1){current_U <- U(current_q)}
    proposed_U <- U(q)
    
    # Create vector of energies for Inf/NaN check
    energies <- c(current_U, proposed_U, current_K, proposed_K)
    
    # Metropolis accept/reject step
    # Check if all energy values are finite and not NaN and also
    # standard M-H accept/reject step
    if(all(is.finite(energies)) & all(!is.nan(energies)) &
       runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
      chain[t,] <- q # Accept new proposal
      current_U <- proposed_U # Store potential energy for next iteration
    }
    else{
      chain[t,] <- current_q # Reject new proposal
    }
  }
  return(chain)
}