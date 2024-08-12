# Title: 04-high-dimensional-funnel.R
# Description: Check if results hold when sampling from high-dimensional funnel too
# NOTE: Was taking way too long to run - over 20 hours for a single sampler, so
      # instead will include some justification for not testing higher D and 
      # moving onto the double funnel instead

################################################################################
# Setup
################################################################################
# Load utils functions
source('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034/src/utils.R')
setup() # Import relevant scripts

# For this version we will only check GIST + RHMC to see if we can still sample well
# If it can handle a 100D funnel then we can move on to the double funnel

################################################################################
# Sampling
################################################################################

sink(file='./log/04-high-dimensional-funnel-script/sampling.log', split=T)
# Create data structure to hold all data
sampling_data <- list()

# GIST, random E, Gaussian KE, divergence limit = 1000, starting in mouth
sampling_data$gist_randomE_k1_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='gaussian',
                                                           current_q=rep(5,100),
                                                           L=NULL,
                                                           epsilon=NULL,
                                                           epsilon_dist='exponential',
                                                           epsilon_params=c(rate=1),
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=1000,
                                                           divergence_auto_reject=F) %>% collect_data

# GIST, random E, Gaussian KE, divergence limit = NULL, starting in neck
sampling_data$gist_randomE_k1_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                            grad_U=funnel_grad_U,
                                                            K='gaussian',
                                                            current_q=rep(1,100),
                                                            L=NULL,
                                                            epsilon=NULL,
                                                            epsilon_dist='exponential',
                                                            epsilon_params=c(rate=1),
                                                            chain_length=500000,
                                                            L_max=1000,
                                                            divergence_limit=NULL,
                                                            divergence_auto_reject=F) %>% collect_data

# GIST, random E, Gaussian KE, divergence limit = NULL, starting in mouth
sampling_data$gist_randomE_k1_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                             grad_U=funnel_grad_U,
                                                             K='gaussian',
                                                             current_q=rep(5,100),
                                                             L=NULL,
                                                             epsilon=NULL,
                                                             epsilon_dist='exponential',
                                                             epsilon_params=c(rate=1),
                                                             chain_length=500000,
                                                             L_max=1000,
                                                             divergence_limit=NULL,
                                                             divergence_auto_reject=F) %>% collect_data

# GIST, random E, Laplace KE, divergence limit = 1000, starting in neck
sampling_data$gist_randomE_k2_div_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='laplace',
                                                          current_q=rep(1,100),
                                                          L=NULL,
                                                          epsilon=NULL,
                                                          epsilon_dist='exponential',
                                                          epsilon_params=c(rate=1),
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=1000,
                                                          divergence_auto_reject=F) %>% collect_data

# GIST, random E, Laplace KE, divergence limit = 1000, starting in mouth
sampling_data$gist_randomE_k2_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='laplace',
                                                           current_q=rep(5,100),
                                                           L=NULL,
                                                           epsilon=NULL,
                                                           epsilon_dist='exponential',
                                                           epsilon_params=c(rate=1),
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=1000,
                                                           divergence_auto_reject=F) %>% collect_data

# GIST, random E, Laplace KE, divergence limit = NULL, starting in neck
sampling_data$gist_randomE_k2_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                            grad_U=funnel_grad_U,
                                                            K='laplace',
                                                            current_q=rep(1,100),
                                                            L=NULL,
                                                            epsilon=NULL,
                                                            epsilon_dist='exponential',
                                                            epsilon_params=c(rate=1),
                                                            chain_length=500000,
                                                            L_max=1000,
                                                            divergence_limit=NULL,
                                                            divergence_auto_reject=F) %>% collect_data

# GIST, random E, Laplace KE, divergence limit = NULL, starting in mouth
sampling_data$gist_randomE_k2_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                             grad_U=funnel_grad_U,
                                                             K='laplace',
                                                             current_q=rep(5,100),
                                                             L=NULL,
                                                             epsilon=NULL,
                                                             epsilon_dist='exponential',
                                                             epsilon_params=c(rate=1),
                                                             chain_length=500000,
                                                             L_max=1000,
                                                             divergence_limit=NULL,
                                                             divergence_auto_reject=F) %>% collect_data

# Export data from batch 1
# Directory to save the CSV files
output_dir <- "./out/100d-funnel/samples/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(sampling_data)){
  
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write the dataframe to a CSV file
  write.csv(sampling_data[[name]], file = output_file, row.names = FALSE)
}

# Print acceptance divergence statistics
for(name in names(sampling_data)){
  print(name)
  accept_diverge_stats(sampling_data[[name]])
  print(sum(sampling_data[[name]][,13]==1000))
}
sink()

################################################################################
# Extract y values (first dimension of funnel samples)
################################################################################
y_data <- lapply(sampling_data, function(sampled_data_df){return(sampled_data_df[,1])})

################################################################################
# Compute rolling bias
################################################################################
print(Sys.time())
y_rolling_bias <- lapply(y_data, rolling_bias)
print(Sys.time())

# Export rolling bias data
# Directory to save the CSV files
output_dir <- "./out/100d-funnel/diagnostics/rolling_bias/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(y_rolling_bias)){
  
  # Create the full path for the output file and write dataframe to CSV
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  write.csv(y_rolling_bias[[name]], file = output_file, row.names = FALSE)
}

################################################################################
# Graph histogram of y samples
################################################################################

for(name in names(y_data)){
  
  # Set filename
  filename <- paste0('./out/100d-funnel/graphics/histograms/', name, '.png')
  
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

################################################################################
# Plots for rolling bias
################################################################################

# Plots for GIST, randomE, k1
png(filename='./out/100d-funnel/graphics/rolling_bias/gist_randomE_k1.png')
plot(y_rolling_bias$gist_randomE_k1_div_mouth, type='l',
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

# Plots for GIST, randomE, k2
png(filename='./out/100d-funnel/graphics/rolling_bias/gist_randomE_k2.png')
plot(y_rolling_bias$gist_randomE_k2_div_mouth, type='l',
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