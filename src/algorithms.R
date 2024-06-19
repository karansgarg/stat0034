# Title: algorithms.R
# Description: R file containing HMC algorithms and supplementary functions

################################################################################
# Vanilla HMC function taken from Neal (2011)
################################################################################
HMC_Neal_2011 <- function(U, grad_U, epsilon, L, current_q){
  
  q <- current_q
  p <- rnorm(length(q), 0, 1)
  current_p <- p
  
  # Make a half-step for momentum at the beginning
  p <- p - epsilon*grad_U(q)/2
  
  # Alternate full-steps for position and momentum
  for(i in 1:L){
    # Full step for position variables
    q <- q + epsilon*p
    
    # Full step for momentum, except at the end
    if(i!=L){p <- p - epsilon*grad_U(q)}
  }
  # Final half step for momentum
  p <- p - epsilon*grad_U(q)/2
  
  # Negate momentum to get symmetric proposal
  p <- -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- sum(current_p^2)/2
  proposed_U <- U(q)
  proposed_K <- sum(p^2)/2
  
  # Metropolis accept/reject step
  if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
    return(q) # Accept new proposal
  }
  else{ return(current_q) } # Reject new proposal
}

################################################################################
# Modified Neal function to produce chain
################################################################################
HMC_sampler_Neal <- function(U, grad_U, epsilon, L, current_q, chain_length,...){
  chain <- matrix(NA, nrow=chain_length, ncol=length(current_q))
  for(t in 1:chain_length){
    q <- current_q
    p <- rnorm(length(q), 0, 1)
    current_p <- p
    
    # Make a half-step for momentum at the beginning
    p <- p - epsilon*grad_U(q,...)/2
    
    # Alternate full-steps for position and momentum
    for(i in 1:L){
      # Full step for position variables
      q <- q + epsilon*p
      
      # Full step for momentum, except at the end
      if(i!=L){p <- p - epsilon*grad_U(q,...)}
    }
    # Final half step for momentum
    p <- p - epsilon*grad_U(q,...)/2
    
    # Negate momentum to get symmetric proposal
    p <- -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- U(current_q)
    current_K <- sum(current_p^2)/2
    proposed_U <- U(q)
    proposed_K <- sum(p^2)/2
    
    #print(exp(current_U-proposed_U+current_K-proposed_K))
    #print(current_U)
    #print(current_K)
    #print(proposed_U)
    #print(propsed_K)
    
    # Metropolis accept/reject step
    if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
      chain[t,] <- q # Accept new proposal
    }
    else{
      chain[t,] <- current_q # Reject new proposal
    } 
  }
  return(chain)
}