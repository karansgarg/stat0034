# Title: algorithms.R
# Description: R file containing HMC algorithms and supplementary functions

################################################################################
# Imports
################################################################################
library(VGAM) # For drawing samples from Laplace distribution
library(HyperbolicDist) # For using hyperbolic distribution sampler

################################################################################
# Leapfrog integrators for HMC sampler
################################################################################

################# Standard Leapfrog integrator ("vanilla HMC") #################

leapfrog_fixed_L <- function(q, # Position variables
                             p, # Momentum variables
                             U, # Potential energy
                             K, # Kinetic energy
                             grad_U, # Gradient of potential energy
                             grad_K, # Gradient of kinetic energy
                             epsilon, # Step size
                             L, # Number of leapfrog steps
                             divergence_limit=NULL # Divergence in Hamiltonians
                             ){ 
  # Set flag for divergence
  divergence <- 0
  
  # Store starting values
  current_q <- q
  current_p <- p
  
  # Sample path length based on set maximum
  L_sampled <- sample(L, size=1)
  
  # Make a half-step for momentum at the beginning
  p <- p - epsilon*grad_U(q)/2
  
  # Alternate full-steps for position and momentum
  for(i in 1:L_sampled){
    # Full step for position variables
    q <- q + epsilon*grad_K(p)
    
    # Full step for momentum, except at the end
    if(i!=L){p <- p - epsilon*grad_U(q)}
  }
  
  # Final half step for momentum
  p <- p - epsilon*grad_U(q)/2
  
  # Negate momentum to get symmetric proposal
  p <- -p
  
  # Indicate divergence reached if Inf/Nan value seen
  if(any(is.na(c(q, p, U(q), K(p)))) | 
     any(is.infinite(c(q, p, U(q), K(p))))){divergence <- 1}
  
  # Indicate divergence if difference in Hamiltonians exceeds threshold
  if(!is.null(divergence_limit) & divergence==0){
    if(abs(U(current_q) + K(current_p) - U(q) - K(p)) > divergence_limit){
      divergence <- 1
    }
  }
  
  # Return q and p
  return(list(q=q,
              p=p,
              divergence=divergence))
}

################ GIST Leapfrog integrator (no-u-turn condition) ################

leapfrog_GIST <- function(q, # Position variables
                          p, # Momentum variables
                          U, # Potential energy
                          K, # Kinetic energy
                          grad_U, # Gradient of potential energy
                          grad_K, # Gradient of kinetic energy
                          epsilon, # Step size
                          L_max=1000, # Maxiumum path length
                          divergence_limit=NULL # Divergence in Hamiltonians
                          ){
  
  L_max_forward <- 0 # Variable for storing u-turn step number
  q_0_forward <- q # Store first q for use in checking condition
  u_turn_condition <- FALSE # Set variable for u-turn condition
  divergence <- 0 # Set flag for divergence
  
  # Data structures to store q and p
  q_store <- matrix(data=q, nrow=1, ncol=length(q))
  p_store <- matrix(data=p, nrow=1, ncol=length(p))
  
  # Compute Hamiltonian at start to assess divergence
  H_0 <- U(q) + K(p)
  
  while(u_turn_condition == FALSE){
    L_max_forward <- L_max_forward + 1 # Add 1 timestep to counter
    
    # Perform leapfrog update for one timestep
    p <- p - epsilon*grad_U(q)/2
    q <- q + epsilon*grad_K(p)
    p <- p - epsilon*grad_U(q)/2
    
    # Store updated results
    q_store <- rbind(q_store, q)
    p_store <- rbind(p_store, p)
    

    # Indicate divergence if difference in Hamiltonians exceeds threshold
    if(!is.null(divergence_limit)){
      H_t <- U(q) + K(p)
      # If Hamiltonian has diverged to inf then we stop as its a divergence anyway
      if(!is.nan(H_t) & !is.infinite(H_t)){
        if(abs(H_t - H_0) > divergence_limit){
          u_turn_condition <- TRUE
          divergence <- 1
          break
        }
      }
      else{
        u_turn_condition <- TRUE
        divergence <- 1
        break
      }
    }
    
    # Stop early if L_max reaches specified limit (default is 1000)
    if(!is.null(L_max)){if(L_max_forward==L_max){break}}
    
    # Catch Inf/NaN values and stop integration
    if(any(is.na(c(q, q_0_forward, p, grad_K(p)))) | 
       any(is.infinite(c(q, q_0_forward, p, grad_K(p))))){
      L_max_forward <- L_max_forward - 1
      u_turn_condition <- TRUE
      divergence <- 1
      break
    }
    
    # Assess u-turn condition
    # Break if Inf/NaN
    u_turn_grad <- (q - q_0_forward) %*% grad_K(p)
    if(is.infinite(u_turn_grad) | is.nan(u_turn_grad)){
      L_max_forward <- L_max_forward - 1
      u_turn_condition <- TRUE
      divergence <- 1
      break
    }
    # Stop if u-turn reached
    if(u_turn_grad < 0){u_turn_condition <- TRUE}
  }
  
  # Sample L discrete uniformly based on u-turn condition
  L <- sample(L_max_forward, size=1)
  
  # Get probability of choosing the sampled L
  prob_L <- 1/L_max_forward
  
  # Determine chosen q and p and negate momentum
  q <- q_store[L+1,]
  p <- p_store[L+1,] * -1
  
  # Only need to do backwards trajectory if we met a stopping critera
  # No need to do backwards check if we stopped because we hit L_max
  if(u_turn_condition == TRUE){
    
    # Reset u-turn condition, L_max and starting q for reverse trajectory
    u_turn_condition <- FALSE
    L_max_backward <- 0
    q_0_backward <- q
    
    # Compute Hamiltonian at start to assess divergence
    H_0 <- U(q) + K(p)
    
    while(u_turn_condition == FALSE){
      L_max_backward <- L_max_backward + 1 # Add 1 timestep to counter
      
      # Perform leapfrog update for one timestep
      p <- p - epsilon*grad_U(q)/2
      q <- q + epsilon*grad_K(p)
      p <- p - epsilon*grad_U(q)/2
      
      # Indicate divergence if difference in Hamiltonians exceeds threshold
      if(!is.null(divergence_limit)){
        H_t <- U(q) + K(p)
        # If Hamiltonian has diverged to inf then we stop as its a divergence anyway
        if(!is.nan(H_t) & !is.infinite(H_t)){
          if(abs(H_t - H_0) > divergence_limit){
            u_turn_condition <- TRUE
            divergence <- 1
            break
          }
        }
        else{
          u_turn_condition <- TRUE
          divergence <- 1
          break
        }
      }
      
      # Catch Inf/NaN values and stop integration
      if(any(is.na(c(q, q_0_backward, p, grad_K(p)))) | 
         any(is.infinite(c(q, q_0_backward, p, grad_K(p))))){
        L_max_backward <- L_max_backward - 1
        divergence <- 1
        break
      }
      
      # Assess u-turn condition
      u_turn_grad <- (q - q_0_backward) %*% grad_K(p)
      if(is.infinite(u_turn_grad) | is.nan(u_turn_grad)){
        L_max_backward <- L_max_backward - 1
        divergence <- 1
        break
      }
      if(u_turn_grad < 0){u_turn_condition <- TRUE}
    }
    
    # Get probability of choosing L within reverse trajectory (0 if start not reached by reverse u-turn)
    if(L-L_max_backward > 0){reverse_prob_L <- 0}
    else{reverse_prob_L <- 1/L_max_backward}
    
    # Recover chosen q and p (with negated momentum for symmetric proposal)
    q <- q_store[L+1,]
    p <- p_store[L+1,] * -1
  }
  else{reverse_prob_L <- prob_L}
  
  return(list(q=q, 
              p=p, 
              prob_L=prob_L, 
              reverse_prob_L=reverse_prob_L,
              N=L_max_forward, 
              divergence=divergence))
}

################################################################################
# HMC sampler
################################################################################
HMC <- function(U, # Potential energy
                grad_U, # Gradient of potential energy
                K='gaussian', # Kinetic energy (can be 'gaussian' or 'laplace')
                current_q, # Starting point of chain
                L=NULL, # Path length (defaults to GIST if null)
                divergence_auto_reject=FALSE, # Auto-reject divergent tracjectories
                epsilon=NULL, # Leapfrog step size (defaults to RHMC if null)
                epsilon_dist='exponential', # Used for RHMC (only exp allowed)
                epsilon_params=c(rate=1), # Starting value for lambda_epsilon
                chain_length=5000,
                ...){
  
  # Additional arguments in ... that can be supplied:
  # L_max=1000, # Maximum path length to limit algorithm runtime
  # divergence_limit=NULL, # Defines limit in Hamiltonian differences allowed
  
  # Extract additional supplied parameters to pass to function
  args <- list(...)
  L_max <- args$L_max
  divergence_limit <- args$divergence_limit
  
  # Extract dimensionality of target
  d <- length(current_q)
  
  # Assign data structures for chain
  chain <- matrix(NA, nrow=chain_length, ncol=d)
  acceptance <- rep(NA, chain_length) # Boolean vector of acceptances
  divergences <- rep(NA, chain_length) # Boolean vector indicating divergences
  N_store <- rep(NA, chain_length) # Store u-turn condition
  
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
  else if(K == 'hyperbolic'){
    # Convert to parameterisation required by HyperbolicDist package
    k_params <- hyperbChangePars(from=2, to=1, Theta=c(1, 0, 1, 0))
    sample_p <- function(x){return(rhyperb(x, Theta=k_params))}
    K <- function(x){return(sum(sqrt(1+x^2)))}
    grad_K <- function(x){return(x/(sqrt(1+x^2)))}
  }
  
  # Perform HMC steps in loop to generate chain
  for(t in 1:chain_length){
    
    # Draw value for epsilon if distribution specified
    if(!is.null(epsilon_dist)){
      if(epsilon_dist=='exponential'){
        epsilon <- rexp(1, epsilon_params)
      }
    }
    
    # Assign variable for updated position variables
    q <- current_q
    
    # Sample fresh momentum
    p <- sample_p(d) # Sample p based on given kinetic energy
    current_p <- p # Set initial p based on what was just sampled
    
    # Perform leapfrog integration to get proposal q and p
    # If L is NULL we use GIST leapfrog integrator
    if(is.null(L)){
      proposal_results <- leapfrog_GIST(q=q,
                                        p=p,
                                        U=U,
                                        K=K,
                                        grad_U=grad_U,
                                        grad_K=grad_K,
                                        epsilon=epsilon,
                                        L_max=L_max,
                                        divergence_limit=divergence_limit)
      q <- proposal_results$q
      p <- proposal_results$p
      prob_L <- proposal_results$prob_L
      reverse_prob_L <- proposal_results$reverse_prob_L
      N_store[t] <- proposal_results$N
      divergences[t] <- proposal_results$divergence
    }
    # If L is supplied we use that in the leapfrog integrator
    else{
      proposal_results <- leapfrog_fixed_L(q=q, 
                                           p=p,
                                           U=U,
                                           K=K,
                                           grad_U=grad_U,
                                           grad_K=grad_K,
                                           epsilon=epsilon,
                                           L=L,
                                           divergence_limit=divergence_limit)
      q <- proposal_results$q
      p <- proposal_results$p
      divergences[t] <- proposal_results$divergence
    }
    
    # Evaluate kinetic energies at start and end of trajectory
    current_K <- K(current_p)
    proposed_K <- K(p)
    
    # Evaluate potential energy at end of trajectory (and start for 1st iteration)
    if(t==1){current_U <- U(current_q)}
    proposed_U <- U(q)
    
    # Create vector of energies for Inf/NaN check
    energies <- c(current_U, proposed_U, current_K, proposed_K)
    
    # Automatically reject proposals from divergent trajectory if option selected
    if(divergence_auto_reject==TRUE){
      if(proposal_results$divergence==1){
        chain[t,] <- current_q # Reject new proposal
        acceptance[t] <- 0
        next
      }
    }
    
    # Check energies are valid quantities - if not we automatically reject
    if(all(is.finite(energies)) & all(!is.nan(energies))){
      
      # Different M-H condition based on choice of L - first is GIST sampler
      if(is.null(L)){
        
        if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K)*reverse_prob_L/prob_L){
          chain[t,] <- q # Accept new proposal
          acceptance[t] <- 1
          current_U <- proposed_U # Store potential energy for next iteration
          current_q <- q # Replace current state in chain
        }
        else{
          chain[t,] <- current_q # Reject new proposal
          acceptance[t] <- 0
        }
      }
      
      # Standard HMC M-H step if L hard-coded supplied by main function call
      else if(!is.null(L)){
        
        if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
          chain[t,] <- q # Accept new proposal
          acceptance[t] <- 1
          current_U <- proposed_U # Store potential energy for next iteration
          current_q <- q # Replace current state in chain
        }
        else{
          chain[t,] <- current_q # Reject new proposal
          acceptance[t] <- 0
        }
      }
      
    }
    else{
      chain[t,] <- current_q # Reject new proposal
      acceptance[t] <- 0
    } 
  }
  return(list(chain=chain,
              acceptance=acceptance,
              divergences=divergences,
              U_turn_length=N_store))
}

################################################################################
# HMC sampler with warmup
################################################################################

HMC_with_warmup <- function(warmup_batch_length=500, # Length of each warmup batch
                            warmup_batches=20, # Number of warmup batches
                            target_a_rate=0.675,
                            tolerance=0.025,
                            ...){
  
  # Extract parameters needed for function
  args <- list(...)
  U <- args$U
  grad_U <- args$grad_U
  K <- args$K
  current_q <- args$current_q
  L <- args$L
  L_max <- args$L_max
  epsilon <- args$epsilon
  epsilon_dist <- args$epsilon_dist
  epsilon_params <- args$epsilon_params
  chain_length <- args$chain_length
  divergence_auto_reject <- args$divergence_auto_reject
  divergence_limit <- args$divergence_limit
  
  # Print start time
  cat("Starting warmup at: ")
  print(Sys.time())
  
  # Extract parameters for adjustment
  if(!is.null(epsilon_dist)){
    if(epsilon_dist=='exponential'){
      rate <- epsilon_params[['rate']]
      cat("Starting lambda for epsilon: ", rate, "\n")
    }
  }
  else{
    cat("Starting epsilon: ", epsilon, "\n")
  }
  
  # Set starting point
  q <- current_q
  
  # Run warmup batches to tune parameters for epsilon
  for(i in 1:warmup_batches){
    results <- HMC(U=U,
                   grad_U=grad_U,
                   K=K,
                   current_q=q,
                   L=L,
                   divergence_auto_reject=divergence_auto_reject,
                   epsilon=epsilon,
                   epsilon_dist=epsilon_dist,
                   epsilon_params=epsilon_params,
                   chain_length=warmup_batch_length,
                   L_max=L_max,
                   divergence_limit=divergence_limit)
    
    a_rate <- mean(results$acceptance)
    
    # Check if epsilon distribution is supplied or if epsilon is supplied
    if(!is.null(epsilon_dist)){
      
      # If epsilon distribution is exponential, we update rate parameter (lambda)
      if(epsilon_dist=='exponential'){
        if(a_rate - target_a_rate < -tolerance){
          # Increase lambda (decrease mean epsilon) if not accepting enough
          rate <- rate * 1.2
        }
        else if(a_rate - target_a_rate > tolerance){
          # Decrease lambda (increase mean epsilon) if not accepting enough
          rate <- rate * 0.8
        }
        else{break}
      }
    }
    # If epsilon distribution is null we update epsilon directly
    else{
      if(a_rate - target_a_rate < -tolerance){
        # Decrease epsilon if not accepting enough
        epsilon <- epsilon * 0.8
      }
      else if(a_rate - target_a_rate > tolerance){
        # Increase epsilon if not accepting enough
        epsilon <- epsilon * 1.2
      }
      else{break}
    }
    
    # Update starting point for next batch
    q <- results$chain[nrow(results$chain),]
  }
  
  # Set parameters after adjustment
  if(!is.null(epsilon_dist)){
    if(epsilon_dist=='exponential'){
      epsilon_params <- c(rate=rate)
      cat("Final lambda for epsilon: ", rate, "\n")
    }
  }
  else{
    cat("Final epsilon: ", epsilon, "\n")
  }
  
  # Print time when starting final chain
  cat("Warmup complete! Starting final chain at time: ")
  print(Sys.time())
  
  # Run HMC with tuned parameters for epsilon
  results <- HMC(U=U,
                 grad_U=grad_U,
                 K=K,
                 current_q=q,
                 L=L,
                 divergence_auto_reject=divergence_auto_reject,
                 epsilon=epsilon,
                 epsilon_dist=epsilon_dist,
                 epsilon_params=epsilon_params,
                 chain_length=chain_length,
                 L_max=L_max,
                 divergence_limit=divergence_limit)
  
  # Print time when final chain complete
  cat("Final chain complete at time: ")
  print(Sys.time())
  
  return(results)
}

################################################################################
# Run sampler for multiple chains with different starting points in parallel
################################################################################
