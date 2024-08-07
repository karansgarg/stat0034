# Title: 01-funnel-sampling.R
# Description: Sampling from 10D Neal's funnel using different variants of HMC

################################################################################
# Setup
################################################################################
# Load utils functions
source('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034/src/utils.R')
setup() # Import relevant scripts

# For all we don't autoreject
# Different variants:
# GIST vs HMC
# Randomised epsilon vs fixed
# Gaussian vs Laplace kinetic energy
# Include or don't include divergence break at 1000
# Starting at neck of funnel vs mouth

# Overall 32 variants of HMC

################################################################################
# BATCH 1 (HMC, Fixed Epsilon)
################################################################################
sink(file='./log/01-funnel-sampling-script/batch1.log', split=T)
# Create data structure to hold all data
sampling_data <- list()

# HMC, fixed E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$hmc_fixedE_k1_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(1,10),
                                                        L=100,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Gaussian KE, divergence limit = 1000, starting in mouth
sampling_data$hmc_fixedE_k1_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(5,10),
                                                        L=100,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Gaussian KE, divergence limit = NULL, starting in neck
sampling_data$hmc_fixedE_k1_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(1,10),
                                                        L=100,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=NULL,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Gaussian KE, divergence limit = NULL, starting in mouth
sampling_data$hmc_fixedE_k1_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(5,10),
                                                        L=100,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=NULL,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Laplace KE, divergence limit = 1000, starting in neck
sampling_data$hmc_fixedE_k2_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='laplace',
                                                        current_q=rep(1,10),
                                                        L=100,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Laplace KE, divergence limit = 1000, starting in mouth
sampling_data$hmc_fixedE_k2_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='laplace',
                                                         current_q=rep(5,10),
                                                         L=100,
                                                         epsilon=1,
                                                         epsilon_dist=NULL,
                                                         epsilon_params=NULL,
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Laplace KE, divergence limit = NULL, starting in neck
sampling_data$hmc_fixedE_k2_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='laplace',
                                                          current_q=rep(1,10),
                                                          L=100,
                                                          epsilon=1,
                                                          epsilon_dist=NULL,
                                                          epsilon_params=NULL,
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=NULL,
                                                          divergence_auto_reject=F) %>% collect_data

# HMC, fixed E, Laplace KE, divergence limit = NULL, starting in mouth
sampling_data$hmc_fixedE_k2_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='laplace',
                                                           current_q=rep(5,10),
                                                           L=100,
                                                           epsilon=1,
                                                           epsilon_dist=NULL,
                                                           epsilon_params=NULL,
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=NULL,
                                                           divergence_auto_reject=F) %>% collect_data

# Export data from batch 1
# Directory to save the CSV files
output_dir <- "./out/samples/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(sampling_data)){
  
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write the dataframe to a CSV file
  write.csv(sampling_data[[name]], file = output_file, row.names = FALSE)
}
sink()

################################################################################
# BATCH 2 (HMC, Randomised Epsilon)
################################################################################
sink(file='./log/01-funnel-sampling-script/batch2.log', split=T)
# Create data structure to hold all data
sampling_data <- list()

# HMC, random E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$hmc_randomE_k1_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(1,10),
                                                        L=100,
                                                        epsilon=NULL,
                                                        epsilon_dist='exponential',
                                                        epsilon_params=c(rate=1),
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, random E, Gaussian KE, divergence limit = 1000, starting in mouth
sampling_data$hmc_randomE_k1_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='gaussian',
                                                         current_q=rep(5,10),
                                                         L=100,
                                                         epsilon=NULL,
                                                         epsilon_dist='exponential',
                                                         epsilon_params=c(rate=1),
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# HMC, random E, Gaussian KE, divergence limit = NULL, starting in neck
sampling_data$hmc_randomE_k1_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='gaussian',
                                                          current_q=rep(1,10),
                                                          L=100,
                                                          epsilon=NULL,
                                                          epsilon_dist='exponential',
                                                          epsilon_params=c(rate=1),
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=NULL,
                                                          divergence_auto_reject=F) %>% collect_data

# HMC, random E, Gaussian KE, divergence limit = NULL, starting in mouth
sampling_data$hmc_randomE_k1_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='gaussian',
                                                           current_q=rep(5,10),
                                                           L=100,
                                                           epsilon=NULL,
                                                           epsilon_dist='exponential',
                                                           epsilon_params=c(rate=1),
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=NULL,
                                                           divergence_auto_reject=F) %>% collect_data

# HMC, random E, Laplace KE, divergence limit = 1000, starting in neck
sampling_data$hmc_randomE_k2_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='laplace',
                                                        current_q=rep(1,10),
                                                        L=100,
                                                        epsilon=NULL,
                                                        epsilon_dist='exponential',
                                                        epsilon_params=c(rate=1),
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# HMC, random E, Laplace KE, divergence limit = 1000, starting in mouth
sampling_data$hmc_randomE_k2_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='laplace',
                                                         current_q=rep(5,10),
                                                         L=100,
                                                         epsilon=NULL,
                                                         epsilon_dist='exponential',
                                                         epsilon_params=c(rate=1),
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# HMC, random E, Laplace KE, divergence limit = NULL, starting in neck
sampling_data$hmc_randomE_k2_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='laplace',
                                                          current_q=rep(1,10),
                                                          L=100,
                                                          epsilon=NULL,
                                                          epsilon_dist='exponential',
                                                          epsilon_params=c(rate=1),
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=NULL,
                                                          divergence_auto_reject=F) %>% collect_data

# HMC, random E, Laplace KE, divergence limit = NULL, starting in mouth
sampling_data$hmc_randomE_k2_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='laplace',
                                                           current_q=rep(5,10),
                                                           L=100,
                                                           epsilon=NULL,
                                                           epsilon_dist='exponential',
                                                           epsilon_params=c(rate=1),
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=NULL,
                                                           divergence_auto_reject=F) %>% collect_data

# Export data from batch 1
# Directory to save the CSV files
output_dir <- "./out/samples/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(sampling_data)){
  
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write the dataframe to a CSV file
  write.csv(sampling_data[[name]], file = output_file, row.names = FALSE)
}
sink()

################################################################################
# BATCH 3 (GIST, Fixed Epsilon)
################################################################################
sink(file='./log/01-funnel-sampling-script/batch3.log', split=T)
# Create data structure to hold all data
sampling_data <- list()

# GIST, fixed E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$gist_fixedE_k1_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='gaussian',
                                                        current_q=rep(1,10),
                                                        L=NULL,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Gaussian KE, divergence limit = 1000, starting in mouth
sampling_data$gist_fixedE_k1_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='gaussian',
                                                         current_q=rep(5,10),
                                                         L=NULL,
                                                         epsilon=1,
                                                         epsilon_dist=NULL,
                                                         epsilon_params=NULL,
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Gaussian KE, divergence limit = NULL, starting in neck
sampling_data$gist_fixedE_k1_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='gaussian',
                                                          current_q=rep(1,10),
                                                          L=NULL,
                                                          epsilon=1,
                                                          epsilon_dist=NULL,
                                                          epsilon_params=NULL,
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=NULL,
                                                          divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Gaussian KE, divergence limit = NULL, starting in mouth
sampling_data$gist_fixedE_k1_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='gaussian',
                                                           current_q=rep(5,10),
                                                           L=NULL,
                                                           epsilon=1,
                                                           epsilon_dist=NULL,
                                                           epsilon_params=NULL,
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=NULL,
                                                           divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Laplace KE, divergence limit = 1000, starting in neck
sampling_data$gist_fixedE_k2_div_neck <- HMC_with_warmup(U=funnel_U,
                                                        grad_U=funnel_grad_U,
                                                        K='laplace',
                                                        current_q=rep(1,10),
                                                        L=NULL,
                                                        epsilon=1,
                                                        epsilon_dist=NULL,
                                                        epsilon_params=NULL,
                                                        chain_length=500000,
                                                        L_max=1000,
                                                        divergence_limit=1000,
                                                        divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Laplace KE, divergence limit = 1000, starting in mouth
sampling_data$gist_fixedE_k2_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='laplace',
                                                         current_q=rep(5,10),
                                                         L=NULL,
                                                         epsilon=1,
                                                         epsilon_dist=NULL,
                                                         epsilon_params=NULL,
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Laplace KE, divergence limit = NULL, starting in neck
sampling_data$gist_fixedE_k2_nodiv_neck <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='laplace',
                                                          current_q=rep(1,10),
                                                          L=NULL,
                                                          epsilon=1,
                                                          epsilon_dist=NULL,
                                                          epsilon_params=NULL,
                                                          chain_length=500000,
                                                          L_max=1000,
                                                          divergence_limit=NULL,
                                                          divergence_auto_reject=F) %>% collect_data

# GIST, fixed E, Laplace KE, divergence limit = NULL, starting in mouth
sampling_data$gist_fixedE_k2_nodiv_mouth <- HMC_with_warmup(U=funnel_U,
                                                           grad_U=funnel_grad_U,
                                                           K='laplace',
                                                           current_q=rep(5,10),
                                                           L=NULL,
                                                           epsilon=1,
                                                           epsilon_dist=NULL,
                                                           epsilon_params=NULL,
                                                           chain_length=500000,
                                                           L_max=1000,
                                                           divergence_limit=NULL,
                                                           divergence_auto_reject=F) %>% collect_data

# Export data from batch 1
# Directory to save the CSV files
output_dir <- "./out/samples/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(sampling_data)){
  
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write the dataframe to a CSV file
  write.csv(sampling_data[[name]], file = output_file, row.names = FALSE)
}
sink()

################################################################################
# BATCH 4 (HMC, Randomised Epsilon)
################################################################################
sink(file='./log/01-funnel-sampling-script/batch4.log', split=T)
# Create data structure to hold all data
sampling_data <- list()

# GIST, random E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$gist_randomE_k1_div_neck <- HMC_with_warmup(U=funnel_U,
                                                         grad_U=funnel_grad_U,
                                                         K='gaussian',
                                                         current_q=rep(1,10),
                                                         L=NULL,
                                                         epsilon=NULL,
                                                         epsilon_dist='exponential',
                                                         epsilon_params=c(rate=1),
                                                         chain_length=500000,
                                                         L_max=1000,
                                                         divergence_limit=1000,
                                                         divergence_auto_reject=F) %>% collect_data

# GIST, random E, Gaussian KE, divergence limit = 1000, starting in mouth
sampling_data$gist_randomE_k1_div_mouth <- HMC_with_warmup(U=funnel_U,
                                                          grad_U=funnel_grad_U,
                                                          K='gaussian',
                                                          current_q=rep(5,10),
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
                                                           current_q=rep(1,10),
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
                                                            current_q=rep(5,10),
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
                                                         current_q=rep(1,10),
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
                                                          current_q=rep(5,10),
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
                                                           current_q=rep(1,10),
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
                                                            current_q=rep(5,10),
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
output_dir <- "./out/samples/"

# Loop over the list and write each dataframe to a CSV file
for(name in names(sampling_data)){
  
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(name, ".csv"))
  
  # Write the dataframe to a CSV file
  write.csv(sampling_data[[name]], file = output_file, row.names = FALSE)
}
sink()