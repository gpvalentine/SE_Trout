# SE_Trout Nmix simulation
# G. Valentine and Y. Kanno
# Feb. 2023

# The purpose of this script is to verify the identifiability of the Bayesian hierarchical models used for the SE_Trout project
# This script generates trout observed trout counts, abundance, and an environmental covariate

library(MCMCvis)
library(data.table)
library(tidyverse)

set.seed(1234)

## settings
nSites <- 100
nYears <- 20
nPasses <- 3
x <- matrix(rnorm(nSites*nYears), nrow = nSites, ncol = nYears)  # climate cov (mean standardized)
alpha0 <- runif(nSites, min = 3, max = 4.5)  # site specific intercept
alpha1 <- rnorm(nSites, mean = -0.3, sd = 0.3)  # spatially heterogeneous climate cov effect
eps <- rnorm(nYears, mean = 0, sd = 0.4)  # time random effect
gamma <- matrix(rnorm(nSites*nYears, mean = 0, sd = 0.2), nrow = nSites, ncol = nYears) # time-site random effect
p <- 0.6  # capture prob per pass
N <- lam <- array(dim = c(nSites, nYears))
y <- array(dim = c(nSites, nYears, nPasses))

## simulate abundance
for (i in 1:nSites){
  for (t in 1:nYears){
    lam[i,t] <- exp(alpha0[i] + alpha1[i]*x[i,t] + eps[t] + gamma[i,t])
    N[i,t] <- rpois(1, lam[i,t])
  }
}

## simulate data
for (i in 1:nSites){
  for (t in 1:nYears){
    y[i,t,1] <- rbinom(1, N[i,t], p)
    y[i,t,2] <- rbinom(1, N[i,t] - y[i,t,1], p)
    y[i,t,3] <- rbinom(1, N[i,t] - y[i,t,1] - y[i,t,2], p)
  }
}

# Write model
sink("Analysis/nMix_JAGS_files/nMix_sim.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nSites){
    omega[i] ~ dunif(0, 5)
  }
  
  ## Betas/Slopes
  # Beta for environmental covariate (site-specific)
    # mu.beta.cov - mean parameter for betas
    mu.beta ~ dnorm(0, 0.01)

    # tau.beta.cov - precision parameter for betas
    # sd parameter for tau.beta
    sd.beta ~ dunif(0, 10)
    tau.beta <- 1/(sd.beta^2)
    s2.beta <- sd.beta^2
      
    for (i in 1:nSites){
      # beta.i - Site-specific coefficients for the covariate
      beta[i] ~ dnorm(mu.beta, tau.beta)
    }
    
  ## Random Effects
  # tau.epsilon: precision parameter for epsilon.t
  sd.eps ~ dunif(0, 10)
  tau.eps <- 1/(sd.eps^2)
    
  # epsilon.t - first order random effect
  for (t in 1:nYears) { 
    eps[t] ~ dnorm(0, tau.eps)
  }
  
    
  # gamma.it - second order random effect
  for (i in 1:nSites) {
    
    # tau.gamma: precision parameter for gamma.it
    sd.gam[i] ~ dunif(0, 10)
    tau.gam[i] <- 1/(sd.gam[i]^2)
    
    for (t in 1:nYears) {
      gam[i,t] ~ dnorm(0, tau.gam[i])
    }
  }
  
  ## Detection probability
  # p - detection probability for each source
  p ~ dbeta(((0.6^2 - 0.6^3 - (0.6 * 0.1^2))/0.1^2), ((0.6 - (2*0.6^2) + 0.6^3 - 0.1^2 + (0.6 * 0.1^2))/0.1^2)) # moment matching for mean 0.6 and variance 0.1
  
  ## Process
  for (i in 1:nSites) {
    for (t in 1:nYears) {
      # Data
      N[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- omega[i] + beta[i] * x[i,t] + eps[t] + gam[i,t]
    }
  }
  
  ## Observation
  for (i in 1:nSites) {
    for (t in 1:nYears) {
    # Pass 1
    y[i,t,1] ~ dbin(p, N[i,t])
    # Pass 2
    y[i,t,2] ~ dbin(p, (N[i,t] - y[i,t,1]))
    # Pass 3
    y[i,t,3] ~ dbin(p, (N[i,t] - y[i,t,1] - y[i,t,2]))
    }
  }
  
  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nSites){
    s2.gam[i] <- 1/tau.gam[i]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nSites = nSites, 
                  nYears = nYears,
                  x = x,
                  y = y)


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "sd.beta", "p", "sd.eps", "s2.eps", "sd.gam", "s2.gam")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.inits <- array(numeric(), dim = c(nSites, nYears))
for (i in 1:nSites) {
  for (t in 1:nYears) {
    N.inits[i,t] <- round(as.numeric(ifelse(is.na((y[i,t,1] + y[i,t,2] + y[i,t,3])),
                                                rpois(1, lambda = 200),
                                                (y[i,t,1] + y[i,t,2] + y[i,t,3] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = runif(nSites, 0, 5),
                             sd.beta = runif(1, 0, 10),
                             mu.beta = rnorm(1, -0.5, 0.01),
                             p = 0.6,
                             N = N.inits)


# MCMC settings
ni <- 20000
nc <- 3
nb <- 5000
nt <- 1

set.seed(1234)

# Fit Model
nMix_sim <- jagsUI::jags(data = jags_data,
                                   parameters.to.save = jags_params,
                                   model.file = "Analysis/nMix_JAGS_files/nMix_sim.jags",
                                   n.chains = nc,
                                   n.iter = ni,
                                   n.burnin = nb,
                                   n.thin = nt,
                                   parallel = T,
                                   inits = init_vals)

nMix_sim_params <- MCMCsummary(nMix_sim, HPD = T) %>% 
  rownames_to_column("param")

MCMCtrace(nMix_sim, params = "omega", pdf = F)

# Export
fwrite(nMix_sim_params, "Analysis/nMix_sim_output.csv")

# Plot posteriors of estimated values vs original values
sim_plot_data <- data.frame(sample_val = rbind(as.matrix(nMix_sim$sims.list$mu.beta),
                                               as.matrix(nMix_sim$sims.list$sd.beta),
                                               as.matrix(nMix_sim$sims.list$sd.eps),
                                               as.matrix(nMix_sim$sims.list$sd.gam[,1]),
                                               as.matrix(nMix_sim$sims.list$p)),
                            parameter = rep(c("mu.beta", "sd.beta", "sd.eps", "sd.gam", "p"),
                                            each = length(nMix_sim$sims.list$mu.beta)))

# specify order of the parameters to plot
sim_plot_data$parameter <- factor(sim_plot_data$parameter, c("mu.beta", "sd.beta", "sd.eps", "sd.gam", "p"))

hlines <- data.frame(y = c(-0.3, 0.3, 0.4, 0.2, 0.6), # these are the generating values for the dataset
                     x = c(1,2,3,4,5))

CombinedModel_Sim.plot <- ggplot() + 
  geom_violin(data = sim_plot_data,
              aes(x = parameter,
                  y = sample_val)) +
  geom_segment(data = hlines,
               aes(x = x - 0.5, xend = x + 0.5,
                   y = y, yend = y),
               color = "red",
               size = 0.8,
               linetype = "dashed") + 
  labs(x = "Parameter",
       y = "Value") +
  theme_classic()
