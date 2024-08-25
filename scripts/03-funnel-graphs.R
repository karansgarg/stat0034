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
#png(filename='./out/graphics/histograms/stan.png')
#hist(stan_samples, freq = F, breaks=100, xlim=c(-10,10), ylim=c(0, 0.15),
#     xlab="y", ylab="Density", main="", 
#     cex.axis = 1.3, cex.lab = 1.3)
#curve(funnel_y_density, from=-10, to=10, add=T, col='red', lwd=2)
#dev.off()

################################################################################
# QQ plots comparing y sampling (GIST/RHMC comparisons)
################################################################################
# Overall settings for QQ plots
# Set colours for lines
qq_colors <- c("Fixed path length and step size" = "black",
                   "Fixed path length, random step size" = "blue",
                   "GIST, fixed step size" = "green",
                   "GIST, random step size" = "red")


######## k1-div-neck (Comparing step size and L choices, 4 lines total) ########
# Combine vectors into a dataframe
k1_div_neck_samples <- data.frame(sample = c(y_data$hmc_fixedE_k1_div_neck,
                                             y_data$hmc_randomE_k1_div_neck,
                                             y_data$gist_fixedE_k1_div_neck,
                                             y_data$gist_randomE_k1_div_neck),
                                  group = rep(c("Fixed path length and step size",
                                                "Fixed path length, random step size",
                                                "GIST, fixed step size",
                                                "GIST, random step size"),
                                              each = 500000))

# Create the overlayed QQ plot
k1_div_neck_qq_plot <- ggplot(k1_div_neck_samples, aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2))

k1_div_neck_qq_plot
ggsave(filename='k1_div_neck.png', path='./out/graphics/qq_plots/', width=9, height=9)

####### k1-nodiv-neck (Comparing step size and L choices, 4 lines total) #######
# Combine vectors into a dataframe
k1_nodiv_neck_samples <- data.frame(sample = c(y_data$hmc_fixedE_k1_nodiv_neck,
                                             y_data$hmc_randomE_k1_nodiv_neck,
                                             y_data$gist_fixedE_k1_nodiv_neck,
                                             y_data$gist_randomE_k1_nodiv_neck),
                                  group = rep(c("Fixed path length and step size",
                                                "Fixed path length, random step size",
                                                "GIST, fixed step size",
                                                "GIST, random step size"),
                                              each = 500000))

# Create the overlayed QQ plot
k1_nodiv_neck_qq_plot <- ggplot(k1_nodiv_neck_samples, aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2))

k1_nodiv_neck_qq_plot
ggsave(filename='k1_nodiv_neck.png', path='./out/graphics/qq_plots/', width=9, height=9)

######## k2-div-neck (Comparing step size and L choices, 4 lines total) ########
# Combine vectors into a dataframe
k2_div_neck_samples <- data.frame(sample = c(y_data$hmc_fixedE_k2_div_neck,
                                             y_data$hmc_randomE_k2_div_neck,
                                             y_data$gist_fixedE_k2_div_neck,
                                             y_data$gist_randomE_k2_div_neck),
                                  group = rep(c("Fixed path length and step size",
                                                "Fixed path length, random step size",
                                                "GIST, fixed step size",
                                                "GIST, random step size"),
                                              each = 500000))

# Create the overlayed QQ plot
k2_div_neck_qq_plot <- ggplot(k2_div_neck_samples, aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2))

k2_div_neck_qq_plot
ggsave(filename='k2_div_neck.png', path='./out/graphics/qq_plots/', width=9, height=9)

######## k3-div-neck (Comparing step size and L choices, 4 lines total) ########
# Combine vectors into a dataframe
k3_div_neck_samples <- data.frame(sample = c(y_data$hmc_fixedE_k3_div_neck,
                                             y_data$hmc_randomE_k3_div_neck,
                                             y_data$gist_fixedE_k3_div_neck,
                                             y_data$gist_randomE_k3_div_neck),
                                  group = rep(c("Fixed path length and step size",
                                                "Fixed path length, random step size",
                                                "GIST, fixed step size",
                                                "GIST, random step size"),
                                              each = 500000))

# Create the overlayed QQ plot
k3_div_neck_qq_plot <- ggplot(k3_div_neck_samples, aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2))

k3_div_neck_qq_plot
ggsave(filename='k3_div_neck.png', path='./out/graphics/qq_plots/', width=9, height=9)

################################################################################
# QQ plots comparing y sampling (KE comparisons)
################################################################################
# Set colours for lines
qq_colors_ke <- c("Gaussian" = "black",
               "Laplace" = "green",
               "Hyperbolic" = "red")

############## hmc-fixedE-div-neck (Comparing KEs, 3 lines total) ##############
# Combine vectors into a dataframe
hmc_fixedE_div_neck_samples <- data.frame(sample = c(y_data$hmc_fixedE_k1_div_neck,
                                             y_data$hmc_fixedE_k2_div_neck,
                                             y_data$hmc_fixedE_k3_div_neck),
                                  group = rep(c("Gaussian",
                                                "Laplace",
                                                "Hyperbolic"),
                                              each = 500000))

# Create the overlayed QQ plot
hmc_fixedE_div_neck_qq_plot <- ggplot(hmc_fixedE_div_neck_samples,
                                      aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors_ke) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))

hmc_fixedE_div_neck_qq_plot
ggsave(filename='hmc_fixedE_div_neck.png',
       path='./out/graphics/qq_plots/', width=9, height=9)

############## hmc-randomE-div-neck (Comparing KEs, 3 lines total) ##############
# Combine vectors into a dataframe
hmc_randomE_div_neck_samples <- data.frame(sample = c(y_data$hmc_randomE_k1_div_neck,
                                                     y_data$hmc_randomE_k2_div_neck,
                                                     y_data$hmc_randomE_k3_div_neck),
                                          group = rep(c("Gaussian",
                                                        "Laplace",
                                                        "Hyperbolic"),
                                                      each = 500000))

# Create the overlayed QQ plot
hmc_randomE_div_neck_qq_plot <- ggplot(hmc_randomE_div_neck_samples,
                                      aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors_ke) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))

hmc_randomE_div_neck_qq_plot
ggsave(filename='hmc_randomE_div_neck.png',
       path='./out/graphics/qq_plots/', width=9, height=9)

############## gist-fixedE-div-neck (Comparing KEs, 3 lines total) ##############
# Combine vectors into a dataframe
gist_fixedE_div_neck_samples <- data.frame(sample = c(y_data$gist_fixedE_k1_div_neck,
                                                     y_data$gist_fixedE_k2_div_neck,
                                                     y_data$gist_fixedE_k3_div_neck),
                                          group = rep(c("Gaussian",
                                                        "Laplace",
                                                        "Hyperbolic"),
                                                      each = 500000))

# Create the overlayed QQ plot
gist_fixedE_div_neck_qq_plot <- ggplot(gist_fixedE_div_neck_samples,
                                      aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors_ke) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))

gist_fixedE_div_neck_qq_plot
ggsave(filename='gist_fixedE_div_neck.png',
       path='./out/graphics/qq_plots/', width=9, height=9)

############## gist-randomE-div-neck (Comparing KEs, 3 lines total) ##############
# Combine vectors into a dataframe
gist_randomE_div_neck_samples <- data.frame(sample = c(y_data$gist_randomE_k1_div_neck,
                                                      y_data$gist_randomE_k2_div_neck,
                                                      y_data$gist_randomE_k3_div_neck),
                                           group = rep(c("Gaussian",
                                                         "Laplace",
                                                         "Hyperbolic"),
                                                       each = 500000))

# Create the overlayed QQ plot
gist_randomE_div_neck_qq_plot <- ggplot(gist_randomE_div_neck_samples,
                                       aes(sample = sample, color = group)) +
  stat_qq(distribution = qnorm, dparams = list(mean = 0, sd = 3)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  scale_color_manual(values = qq_colors_ke) +
  labs(title = "", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(text = element_text(size=20), legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))

gist_randomE_div_neck_qq_plot
ggsave(filename='gist_randomE_div_neck.png',
       path='./out/graphics/qq_plots/', width=9, height=9)

################################################################################
# Graph rolling bias for report figures (GIST/RHMC comparisons)
################################################################################

# Plot comparing GIST vs fixed L and RHMC vs fixed step size (Gaussian KE, div, neck)
png(filename='./out/graphics/rolling_bias/k1_div_neck.png')
plot(y_rolling_bias$hmc_fixedE_k1_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k1_div_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k1_div_neck, col='green')
lines(y_rolling_bias$gist_randomE_k1_div_neck, col='red')
legend('topright', legend=c('Fixed path length and step size',
                            'Fixed path length, random step size',
                            'GIST, fixed step size',
                            'GIST, randomised step size'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing GIST vs fixed L and RHMC vs fixed step size (Gaussian KE, nodiv, neck)
png(filename='./out/graphics/rolling_bias/k1_nodiv_neck.png')
plot(y_rolling_bias$hmc_fixedE_k1_nodiv_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k1_nodiv_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k1_nodiv_neck, col='green')
lines(y_rolling_bias$gist_randomE_k1_nodiv_neck, col='red')
legend('topright', legend=c('Fixed path length and step size',
                            'Fixed path length, random step size',
                            'GIST, fixed step size',
                            'GIST, randomised step size'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing GIST vs fixed L and RHMC vs fixed step size (Laplace KE)
png(filename='./out/graphics/rolling_bias/k2_div_neck.png')
plot(y_rolling_bias$hmc_fixedE_k2_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k2_div_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k2_div_neck, col='green')
lines(y_rolling_bias$gist_randomE_k2_div_neck, col='red')
legend('topright', legend=c('Fixed path length and step size',
                            'Fixed path length, random step size',
                            'GIST, fixed step size',
                            'GIST, randomised step size'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing GIST vs fixed L and RHMC vs fixed step size (Hyperbolic KE)
png(filename='./out/graphics/rolling_bias/k3_div_neck.png')
plot(y_rolling_bias$hmc_fixedE_k3_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k3_div_neck, col='blue')
lines(y_rolling_bias$gist_fixedE_k3_div_neck, col='green')
lines(y_rolling_bias$gist_randomE_k3_div_neck, col='red')
legend('topright', legend=c('Fixed path length and step size',
                            'Fixed path length, random step size',
                            'GIST, fixed step size',
                            'GIST, randomised step size'),
       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
dev.off()

################################################################################
# Graph rolling bias for report figures (KE comparisons)
################################################################################

# Plot comparing KEs (hmc_fixedE_div_neck)
png(filename='./out/graphics/rolling_bias/hmc_fixedE_div_neck.png')
plot(y_rolling_bias$hmc_fixedE_k1_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_fixedE_k2_div_neck, col='green')
lines(y_rolling_bias$hmc_fixedE_k3_div_neck, col='red')
legend('topright', legend=c('Gaussian',
                            'Laplace',
                            'Hyperbolic'),
       col=c('black', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing KEs (hmc_randomE_div_neck)
png(filename='./out/graphics/rolling_bias/hmc_randomE_div_neck.png')
plot(y_rolling_bias$hmc_randomE_k1_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$hmc_randomE_k2_div_neck, col='green')
lines(y_rolling_bias$hmc_randomE_k3_div_neck, col='red')
legend('topright', legend=c('Gaussian',
                            'Laplace',
                            'Hyperbolic'),
       col=c('black', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing KEs (gist_fixedE_div_neck)
png(filename='./out/graphics/rolling_bias/gist_fixedE_div_neck.png')
plot(y_rolling_bias$gist_fixedE_k1_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_fixedE_k2_div_neck, col='green')
lines(y_rolling_bias$gist_fixedE_k3_div_neck, col='red')
legend('topright', legend=c('Gaussian',
                            'Laplace',
                            'Hyperbolic'),
       col=c('black', 'green', 'red'), pch=16, cex=1.5)
dev.off()

# Plot comparing KEs (gist_randomE_div_neck)
png(filename='./out/graphics/rolling_bias/gist_randomE_div_neck.png')
plot(y_rolling_bias$gist_randomE_k1_div_neck, type='l', ylim=c(0,3),
     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
     cex.axis=1.3, cex.lab=1.3)
lines(y_rolling_bias$gist_randomE_k2_div_neck, col='green')
lines(y_rolling_bias$gist_randomE_k3_div_neck, col='red')
legend('topright', legend=c('Gaussian',
                            'Laplace',
                            'Hyperbolic'),
       col=c('black', 'green', 'red'), pch=16, cex=1.5)
dev.off()

################################################################################
# Graph rolling bias for all samplers
################################################################################

# NUTS/Stan
png(filename='./out/graphics/rolling_bias/stan.png')
plot(stan_rolling_bias, type='l', xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators', ylim=c(0, 3),
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
#png(filename='./out/graphics/rolling_bias_sq/hmc_fixedE_k1.png')
#plot(y_rolling_bias_sq$hmc_fixedE_k1_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$hmc_fixedE_k1_div_neck, col='blue')
#lines(y_rolling_bias_sq$hmc_fixedE_k1_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$hmc_fixedE_k1_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for HMC, randomE, k1
#png(filename='./out/graphics/rolling_bias_sq/hmc_randomE_k1.png')
#plot(y_rolling_bias_sq$hmc_randomE_k1_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$hmc_randomE_k1_div_neck, col='blue')
#lines(y_rolling_bias_sq$hmc_randomE_k1_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$hmc_randomE_k1_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for HMC, fixedE, k2
#png(filename='./out/graphics/rolling_bias_sq/hmc_fixedE_k2.png')
#plot(y_rolling_bias_sq$hmc_fixedE_k2_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$hmc_fixedE_k2_div_neck, col='blue')
#lines(y_rolling_bias_sq$hmc_fixedE_k2_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$hmc_fixedE_k2_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for HMC, randomE, k2
#png(filename='./out/graphics/rolling_bias_sq/hmc_randomE_k2.png')
#plot(y_rolling_bias_sq$hmc_randomE_k2_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$hmc_randomE_k2_div_neck, col='blue')
#lines(y_rolling_bias_sq$hmc_randomE_k2_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$hmc_randomE_k2_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for GIST, fixedE, k1
#png(filename='./out/graphics/rolling_bias_sq/gist_fixedE_k1.png')
#plot(y_rolling_bias_sq$gist_fixedE_k1_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$gist_fixedE_k1_div_neck, col='blue')
#lines(y_rolling_bias_sq$gist_fixedE_k1_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$gist_fixedE_k1_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for GIST, randomE, k1
#png(filename='./out/graphics/rolling_bias_sq/gist_randomE_k1.png')
#plot(y_rolling_bias_sq$gist_randomE_k1_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$gist_randomE_k1_div_neck, col='blue')
#lines(y_rolling_bias_sq$gist_randomE_k1_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$gist_randomE_k1_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for GIST, fixedE, k2
#png(filename='./out/graphics/rolling_bias_sq/gist_fixedE_k2.png')
#plot(y_rolling_bias_sq$gist_fixedE_k2_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$gist_fixedE_k2_div_neck, col='blue')
#lines(y_rolling_bias_sq$gist_fixedE_k2_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$gist_fixedE_k2_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

# Plots for GIST, randomE, k2
#png(filename='./out/graphics/rolling_bias_sq/gist_randomE_k2.png')
#plot(y_rolling_bias_sq$gist_randomE_k2_div_mouth, type='l', ylim=c(0,5),
#     xlab='Number of Iterations', ylab='Absolute Error of MCMC Estimators',
#     cex.axis=1.3, cex.lab=1.3)
#lines(y_rolling_bias_sq$gist_randomE_k2_div_neck, col='blue')
#lines(y_rolling_bias_sq$gist_randomE_k2_nodiv_mouth, col='green')
#lines(y_rolling_bias_sq$gist_randomE_k2_nodiv_neck, col='red')
#legend('topright', legend=c('Div limit, start in mouth',
#                            'Div limit, start in neck',
#                            'No div limit, start in mouth',
#                            'No div limit, start in neck'),
#       col=c('black', 'blue', 'green', 'red'), pch=16, cex=1.5)
#dev.off()

################################################################################
# Print acceptance divergence statistics
################################################################################
for(name in names(sampling_data)){
  print(name)
  accept_diverge_stats(sampling_data[[name]])
}

# Consider amount of times L_max was hit during a trajectory
sampling_data_gist <- sampling_data[grepl("gist", names(sampling_data))]
for(name in names(sampling_data_gist)){
  print(name)
  print(sum(sampling_data[[name]][,13]==1000))
}
