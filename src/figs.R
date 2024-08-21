# Title: figs.R
# Description: Produce (generic, non-results) figures needed for report

# NOTE: Recall that for HMC, when the target distribution is p(x),
# the potential energy function is U := -log(p(x))

################################################################################
# Imports
################################################################################
library(tidyverse) # Mainly for ggplot2
library(ggforce) # For plotting circles
library(VGAM) # For dlaplace function
library(HyperbolicDist) # For dhyperb function

################################################################################
# Plots to demonstrate choice of step size
################################################################################
plot_trajectory <- function(epsilon, L, text_size=24){
  
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
    theme_minimal() + theme(text = element_text(size=text_size))
  
  return(plot)
}

small_step_size_plot <- plot_trajectory(epsilon=0.1, L=5)
big_step_size_plot <- plot_trajectory(epsilon=1.1, L=5)
diverge_step_size_plot <- plot_trajectory(epsilon=2.1, L=5, text_size=12)

png(filename='./out/graphics/small_step_size_trajectory.png')
small_step_size_plot
dev.off()

png(filename='./out/graphics/big_step_size_trajectory.png')
big_step_size_plot
dev.off()

png(filename='./out/graphics/diverge_step_size_trajectory.png')
diverge_step_size_plot
dev.off()
################################################################################
# Plots to show kinetic energy tail behaviour
################################################################################

kinetic_energies_plot <- function(x_limit=6){
  
  # Get samples from each distribution
  x <- seq(-x_limit, x_limit, length.out=1000000)
  gaussian <- dnorm(x)
  laplace <- dlaplace(x)
  k_params <- hyperbChangePars(from=2, to=1, Theta=c(1, 0, 1, 0))
  hyperbolic <- dhyperb(x, k_params)

  plot(x, gaussian, type='l', xlim=c(-x_limit, x_limit), ylim=c(0, 0.51),
       xlab="", ylab="Density", col="red", cex.axis=1.3, cex.lab=1.5)
  lines(x, laplace, col="blue")
  lines(x, hyperbolic, col="green")
  legend('topright',
         legend=c("Gaussian Density", "Laplace Density", "Hyperbolic Density"),
         col=c("red", "blue", "green"), pch=16, cex=1.5)
}

png(filename='./out/graphics/kinetic-energy-tails.png', width = 12, height = 6, units = 'in', res=300)
kinetic_energies_plot()
dev.off()