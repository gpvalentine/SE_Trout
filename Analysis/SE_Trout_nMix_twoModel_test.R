# George Valentine
# Synchrony model testing: does it work to split the full model into two? One with covariates and one with the two random effects
# Jan. 2023


####
# Model with just env. covariates
sink("Analysis/nMix_JAGS_files/YOY_climateEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }
  
  ## Betas/Slopes
    
  # Betas for environmental covariates (site-specific)
  for (m in 1:3) {
    
    # mu.beta.cov - mean parameter for beta.covs
    mu.beta.cov[m] ~ dnorm(0, 0.01)

    # tau.beta.cov - precision parameter for beta.covs
    # sd parameter for tau.beta.covs
    sd.cov[m] ~ dunif(0, 10)
    tau.beta.cov[m] <- 1/(sd.cov[m]^2)
    s2.beta.cov[m] <- sd.cov[m]^2
      
    for (i in 1:nReps){
      
      # beta.cov.i - Site-specific coefficients for mean max summer temp, max 0.9Q winter flow, max 0.9Q spring flow
      beta.cov[m,i] ~ dnorm(mu.beta.cov[m], tau.beta.cov[m])
      
    }
  }
  
  ## Detection probability
  # p.j - detection probability for each source
  for (j in 1:nSources) {
      p[j] ~ dbeta(((0.5^2 - 0.5^3 - (0.5 * 0.1^2))/0.1^2), ((0.5 - (2*0.5^2) + 0.5^3 - 0.1^2 + (0.5 * 0.1^2))/0.1^2)) # moment matching for mean 0.5 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta.cov[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta.cov[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_YOY[i,t] ~ dbin(p[Sources[i]], N.YOY[i,t])
    # Pass 2
    p2_YOY[i,t] ~ dbin(p[Sources[i]], (N.YOY[i,t] - p1_YOY[i,t]))
    # Pass 3
    p3_YOY[i,t] ~ dbin(p[Sources[i]], (N.YOY[i,t] - p1_YOY[i,t] - p2_YOY[i,t]))
    }
  }
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_YOY.new[i,t] ~ dbin(p[Sources[i]], N.YOY[i,t])
    }
  }
  
  # Get means, CVs of data and predicted data
  mean_p1_YOY <- mean(p1_YOY[,1:nYears])
  CV_p1_YOY <- sd(p1_YOY[,1:nYears])/mean(p1_YOY[,1:nYears])
  mean_p1_YOY.new <- mean(p1_YOY.new)
  CV_p1_YOY.new <- sd(p1_YOY.new)/mean(p1_YOY.new)
  
  # Calculate p values
  pval.mean_p1 <- step(mean_p1_YOY.new - mean_p1_YOY)
  pval.CV_p1 <- step(CV_p1_YOY.new - CV_p1_YOY)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = sample_areas_wide,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_YOY = p1_YOY,
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)


# Parameters to save
jags_params <- c("beta.cov", "mu.beta.cov", "tau.beta.cov", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

set.seed(1234)

# Fit Model
YOY_climateEffects <- jagsUI::jags(data = jags_data,
                                  parameters.to.save = jags_params,
                                  model.file = "Analysis/nMix_JAGS_files/YOY_climateEffects.jags",
                                  n.chains = nc,
                                  n.iter = ni,
                                  n.burnin = nb,
                                  n.thin = nt,
                                  parallel = T,
                                  inits = init_vals)

YOY_climateEffects_params <- MCMCsummary(YOY_climateEffects, HPD = T)

fwrite(YOY_climateEffects_params, "C:/Users/georgepv/OneDrive - Colostate/Lab PC Backup/Desktop/YOY_climateEffects_params.csv")

# What was the range of covariate effects on YOY abundance?
YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta.cov\\[3")) %>% 
  .[,2] %>% 
  range()

####
# Model with just random effects
sink("Analysis/nMix_JAGS_files/YOY_randomEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
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
  for (i in 1:nReps) {
    
    # tau.gamma: precision parameter for gamma.it
    sd.gam[i] ~ dunif(0, 10)
    tau.gam[i] <- 1/(sd.gam[i]^2)
    
    for (t in 1:nYears) {
      gam[i,t] ~ dnorm(0, tau.gam[i])
    }
  }
  
  ## Detection probability
  # p.j - detection probability for each source
  for (j in 1:nSources) {
      p[j] ~ dbeta(((0.5^2 - 0.5^3 - (0.5 * 0.1^2))/0.1^2), ((0.5 - (2*0.5^2) + 0.5^3 - 0.1^2 + (0.5 * 0.1^2))/0.1^2)) # moment matching for mean 0.5 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_YOY[i,t] ~ dbin(p[Sources[i]], N.YOY[i,t])
    # Pass 2
    p2_YOY[i,t] ~ dbin(p[Sources[i]], (N.YOY[i,t] - p1_YOY[i,t]))
    # Pass 3
    p3_YOY[i,t] ~ dbin(p[Sources[i]], (N.YOY[i,t] - p1_YOY[i,t] - p2_YOY[i,t]))
    }
  }

  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nReps){
    s2.gam[i] <- 1/tau.gam[i]
  }
    
  ## ICC
  for (i in 1:nReps) {
    ICC.YOY[i] <- s2.eps/(s2.eps + s2.gam[i])
  }
  mean_ICC.YOY <- mean(ICC.YOY)
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_YOY.new[i,t] ~ dbin(p[Sources[i]], N.YOY[i,t])
    }
  }
  
  # Get means, CVs of data and predicted data
  mean_p1_YOY <- mean(p1_YOY[,1:nYears])
  CV_p1_YOY <- sd(p1_YOY[,1:nYears])/mean(p1_YOY[,1:nYears])
  mean_p1_YOY.new <- mean(p1_YOY.new)
  CV_p1_YOY.new <- sd(p1_YOY.new)/mean(p1_YOY.new)
  
  # Calculate p values
  pval.mean_p1 <- step(mean_p1_YOY.new - mean_p1_YOY)
  pval.CV_p1 <- step(CV_p1_YOY.new - CV_p1_YOY)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = sample_areas_wide,
                  p1_YOY = p1_YOY,
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)


# Parameters to save
jags_params <- c("p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

set.seed(1234)

# Fit Model
YOY_randomEffects <- jagsUI::jags(data = jags_data,
                                  parameters.to.save = jags_params,
                                  model.file = "Analysis/nMix_JAGS_files/YOY_randomEffects.jags",
                                  n.chains = nc,
                                  n.iter = ni,
                                  n.burnin = nb,
                                  n.thin = nt,
                                  parallel = T,
                                  inits = init_vals)

YOY_randomEffects_params <- MCMCsummary(YOY_BKT_nMix_full, HPD = T)

# calculate average ICC
YOY_randomEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^ICC")) %>% 
  .[,2] %>% 
  mean()

fwrite(YOY_randomEffects_params, "C:/Users/georgepv/OneDrive - Colostate/Lab PC Backup/Desktop/YOY_randomEffects_params.csv")
