# Title: 03-funnel-graphs.R
# Description: Analysing results of sampling: creating graphs

################################################################################
# Setup
################################################################################
# Load utils functions
source('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034/src/utils.R')
setup() # Import relevant scripts

################################################################################
# Import data
################################################################################
# Not needed at the moment since all data is still in memory - need to add later
# so that scripts can be run top to tail from a fresh R session with no dependency

y_rolling_bias <- `02-funnel-diagnostics_results`$y_rolling_bias
y_rolling_bias_sq <- `02-funnel-diagnostics_results`$y_rolling_bias_sq

stan_rolling_bias <- read_csv(file='./out/diagnostics/rolling_bias/stan.csv') %>% data.frame
stan_rolling_bias <- stan_rolling_bias$x

stan_rolling_bias_sq <- read_csv(file='./out/diagnostics/rolling_bias_sq/stan.csv') %>% data.frame
stan_rolling_bias_sq <- stan_rolling_bias_sq$x
################################################################################
# Graph histogram of y samples
################################################################################

for(name in names(y_data)){
  
  # Set filename
  filename <- paste0('./out/graphics/histograms/', name, '.png')
  
  # Create graphics device
  png(filename=filename)

  # Plot histogram
  hist(y_data[[name]], freq = F, breaks=100, xlim=c(-10,10), ylim=c(0, 0.15),
       xlab="y", ylab="Density", main="",
       cex.axis = 1.3, cex.lab = 1.3)
  curve(funnel_y_density, from=-10, to=10, add=T, col='red', lwd=2)

  # Turn off graphics device
  dev.off()
  rm(filename)
}

# Create histogram for Stan samples
png(filename='./out/graphics/histograms/stan.png')
hist(stan_samples, freq = F, breaks=100, xlim=c(-10,10), ylim=c(0, 0.15),
     xlab="y", ylab="Density", main="", 
     cex.axis = 1.3, cex.lab = 1.3)
curve(funnel_y_density, from=-10, to=10, add=T, col='red', lwd=2)
dev.off()
################################################################################
# Graph rolling bias
################################################################################

# NUTS/Stan
png(filename='./out/graphics/rolling_bias/stan.png')
plot(stan_rolling_bias, type='l', xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators', ylim=c(0, 4),
     cex.axis = 1.3, cex.lab = 1.3)
lines(stan_rolling_bias_sq, col='blue')
legend('topright', legend=c('f(y)=y', as.expression('f(y)='~y^2)), col=c('black', 'blue'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, fixedE, k1
png(filename='./out/graphics/rolling_bias/hmc_fixedE_k1.png')
plot(y_rolling_bias$hmc_fixedE_k1_div_mouth, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_fixedE_k1_div_neck, col='blue')
lines(y_rolling_bias$hmc_fixedE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias$hmc_fixedE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, randomE, k1
png(filename='./out/graphics/rolling_bias/hmc_randomE_k1.png')
plot(y_rolling_bias$hmc_randomE_k1_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k1_div_neck, col='blue')
lines(y_rolling_bias$hmc_randomE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias$hmc_randomE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, fixedE, k2
png(filename='./out/graphics/rolling_bias/hmc_fixedE_k2.png')
plot(y_rolling_bias$hmc_fixedE_k2_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_fixedE_k2_div_neck, col='blue')
lines(y_rolling_bias$hmc_fixedE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias$hmc_fixedE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, randomE, k2
png(filename='./out/graphics/rolling_bias/hmc_randomE_k2.png')
plot(y_rolling_bias$hmc_randomE_k2_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k2_div_neck, col='blue')
lines(y_rolling_bias$hmc_randomE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias$hmc_randomE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, fixedE, k1
png(filename='./out/graphics/rolling_bias/gist_fixedE_k1.png')
plot(y_rolling_bias$gist_fixedE_k1_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_fixedE_k1_div_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias$gist_fixedE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, randomE, k1
png(filename='./out/graphics/rolling_bias/gist_randomE_k1.png')
plot(y_rolling_bias$gist_randomE_k1_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_randomE_k1_div_neck, col='blue')
lines(y_rolling_bias$gist_randomE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias$gist_randomE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, fixedE, k2
png(filename='./out/graphics/rolling_bias/gist_fixedE_k2.png')
plot(y_rolling_bias$gist_fixedE_k2_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_fixedE_k2_div_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias$gist_fixedE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, randomE, k2
png(filename='./out/graphics/rolling_bias/gist_randomE_k2.png')
plot(y_rolling_bias$gist_randomE_k2_div_mouth, type='l', ylim=c(0,2),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_randomE_k2_div_neck, col='blue')
lines(y_rolling_bias$gist_randomE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias$gist_randomE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

################################################################################
# Graph rolling bias for square function
################################################################################

# Plots for HMC, fixedE, k1
png(filename='./out/graphics/rolling_bias_sq/hmc_fixedE_k1.png')
plot(y_rolling_bias_sq$hmc_fixedE_k1_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$hmc_fixedE_k1_div_neck, col='blue')
lines(y_rolling_bias_sq$hmc_fixedE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$hmc_fixedE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, randomE, k1
png(filename='./out/graphics/rolling_bias_sq/hmc_randomE_k1.png')
plot(y_rolling_bias_sq$hmc_randomE_k1_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$hmc_randomE_k1_div_neck, col='blue')
lines(y_rolling_bias_sq$hmc_randomE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$hmc_randomE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, fixedE, k2
png(filename='./out/graphics/rolling_bias_sq/hmc_fixedE_k2.png')
plot(y_rolling_bias_sq$hmc_fixedE_k2_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$hmc_fixedE_k2_div_neck, col='blue')
lines(y_rolling_bias_sq$hmc_fixedE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$hmc_fixedE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for HMC, randomE, k2
png(filename='./out/graphics/rolling_bias_sq/hmc_randomE_k2.png')
plot(y_rolling_bias_sq$hmc_randomE_k2_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$hmc_randomE_k2_div_neck, col='blue')
lines(y_rolling_bias_sq$hmc_randomE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$hmc_randomE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, fixedE, k1
png(filename='./out/graphics/rolling_bias_sq/gist_fixedE_k1.png')
plot(y_rolling_bias_sq$gist_fixedE_k1_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$gist_fixedE_k1_div_neck, col='blue')
lines(y_rolling_bias_sq$gist_fixedE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$gist_fixedE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, randomE, k1
png(filename='./out/graphics/rolling_bias_sq/gist_randomE_k1.png')
plot(y_rolling_bias_sq$gist_randomE_k1_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$gist_randomE_k1_div_neck, col='blue')
lines(y_rolling_bias_sq$gist_randomE_k1_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$gist_randomE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, fixedE, k2
png(filename='./out/graphics/rolling_bias_sq/gist_fixedE_k2.png')
plot(y_rolling_bias_sq$gist_fixedE_k2_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$gist_fixedE_k2_div_neck, col='blue')
lines(y_rolling_bias_sq$gist_fixedE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$gist_fixedE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plots for GIST, randomE, k2
png(filename='./out/graphics/rolling_bias_sq/gist_randomE_k2.png')
plot(y_rolling_bias_sq$gist_randomE_k2_div_mouth, type='l', ylim=c(0,5),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias_sq$gist_randomE_k2_div_neck, col='blue')
lines(y_rolling_bias_sq$gist_randomE_k2_nodiv_mouth, col='green')
lines(y_rolling_bias_sq$gist_randomE_k2_nodiv_neck, col='red')
legend('topright', legend=c('Div limit, start in mouth',
                            'Div limit, start in neck',
                            'No div limit, start in mouth',
                            'No div limit, start in neck'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

################################################################################
# Print acceptance divergence statistics
################################################################################
for(name in names(sampling_data)){
  print(name)
  accept_diverge_stats(sampling_data[[name]])
}

# Consider amount of times L_max was hit during a trajectory
for(name in names(sampling_data)[17:32]){
  print(name)
  print(sum(sampling_data[[name]][,13]==1000))
}
