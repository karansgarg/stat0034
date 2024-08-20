# Title: figs.R
# Description: Produce (generic, non-results) figures needed for report

# NOTE: Recall that for HMC, when the target distribution is p(x),
# the potential energy function is U := -log(p(x))

################################################################################
# Imports
################################################################################
library(tidyverse)
library(ggforce)

################################################################################
# Plots to demonstrate choice of step size
################################################################################
plot_trajectory <- function(epsilon, L){
  
  # Define the true trajectory as an ellipse/circle
  theta <- seq(0, 2*pi, length.out = 100)
  true_x <- sqrt(1.0) * cos(theta)
  true_p <- sqrt(1.0) * sin(theta)
  
  # Initialize the simulated trajectory
  x <- numeric(L)
  p <- numeric(L)
  x[1] <- 1.0 # initial position
  p[1] <- 0.0 # initial momentum
  
  # Perform the Leapfrog integration
  for (i in 1:(L-1)) {
    # Half-step momentum update
    p_half <- p[i] - (epsilon / 2) * x[i]
    
    # Full-step position update
    x[i+1] <- x[i] + epsilon * p_half
    
    # Final half-step momentum update
    p[i+1] <- p_half - (epsilon / 2) * x[i+1]
  }
    
  # Create a data frame for plotting
  simulated_trajectory <- data.frame(x = x, p = p)
  true_trajectory <- data.frame(x = true_x, p = true_p)
  
  # Plot using ggplot2
  plot <- ggplot() + geom_circle(aes(x0 = 0, y0 = 0, r = 1),
                                 color = "grey", linetype = "dotted", linewidth = 1) +
    geom_point(data = simulated_trajectory, aes(x = x, y = p), color = "black") +
    geom_line(data = simulated_trajectory, aes(x = x, y = p), color = "black") +
    labs(x = "Position (q)", y = "Momentum (p)") +
    coord_fixed(ratio = 1) + # Fix aspect ratio to ensure circle is not distorted
    theme_minimal()
  
  return(plot)
}