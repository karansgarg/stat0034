# Title: 05-time-normalisation.R
# Description: Assess time normalised quality of sampling

################################################################################
# Setup
################################################################################
# Load utils functions
source('/Users/karansgarg/Documents/MSc Statistics/STAT0034/Project/git/stat0034/src/utils.R')
setup() # Import relevant scripts

################################################################################
# Run 1 version of each sampler for timing (HMC, RHMC, GIST, RHMC-GIST)
################################################################################
sink(file='./log/05-time-normalisation-script/sampling.log', split=T)

# Create data structure to hold all data
sampling_data <- list()

# HMC, fixed E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$hmc_fixedE <- HMC_with_warmup(U=funnel_U,
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

# HMC, random E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$hmc_randomE <- HMC_with_warmup(U=funnel_U,
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

# GIST, fixed E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$gist_fixedE <- HMC_with_warmup(U=funnel_U,
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

# GIST, random E, Gaussian KE, divergence limit = 1000, starting in neck
sampling_data$gist_randomE <- HMC_with_warmup(U=funnel_U,
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

sink()

################################################################################
# Import data
################################################################################

# Import sampling results
sampling_data <- import_data()

# Extract y samples
y_data <- lapply(sampling_data, function(sampled_data_df){return(sampled_data_df[,1])})

# Import rolling bias data
y_rolling_bias <- import_data('./out/diagnostics/rolling_bias/')

# y_rolling_bias imports as dataframes so we need to convert to vector
for(name in names(y_rolling_bias)){
  y_rolling_bias[[name]] <- pull(y_rolling_bias[[name]], var=1)
}

################################################################################
# Histograms using time-normalised data
################################################################################

restricted_rhmc_gist <- y_data$gist_randomE_k1_div_neck[1:9268]

png(filename='./out/graphics/histograms/restricted-rhmc-gist.png')

# Plot histogram of restricted samples
hist(restricted_rhmc_gist, freq = F, breaks=100, xlim=c(-10,10), ylim=c(0, 0.15),
     xlab="y", ylab="Density", main="",
     cex.axis = 1.3, cex.lab = 1.3)
curve(funnel_y_density, from=-10, to=10, add=T, col='red', lwd=2)

dev.off()