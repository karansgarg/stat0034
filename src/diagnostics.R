# Title: diagnostics.R
# Description: R file containing functions to evaluate output of MCMC algorithms

################################################################################
# Load relevant packages
################################################################################
library(coda)

################################################################################
# Convert chain to MCMC, with option to discard based on warmup
################################################################################
make_mcmc <- function(chain, warmup){
  # Discard portion of chain based on warmup specified
  # Make chain an MCMC object
}

################################################################################
# Produce summary diagnostics of MCMC chain
################################################################################
#results_mcmc