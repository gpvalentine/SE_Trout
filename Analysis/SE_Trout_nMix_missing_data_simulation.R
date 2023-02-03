## SE_Trout nMix missing data simulation
# George V

set.seed(1234)

# How much missing data do I have?
mean(is.na(p1_YOY[,1:34]))
mean(is.na(p2_adult[,1:34]))

mean(as.matrix(p1_YOY[,1:34]), na.rm = T)

# Simulate three passes of complete data with the same missing values as p1, p2, p3
p1_0 <- matrix(rpois(lambda = mean(as.matrix(p1_YOY[,1:34]), na.rm = T), nrow(p1)*34),
               nrow = nrow(p1))
p2_0 <- matrix(rpois(lambda = mean(as.matrix(p2_YOY[,1:34]), na.rm = T), nrow(p2)*34),
               nrow = nrow(p2))
p3_0 <- matrix(rpois(lambda = mean(as.matrix(p3_YOY[,1:34]), na.rm = T), nrow(p3)*34),
               nrow = nrow(p3))

p1_sim <- p1_YOY
p1_sim[,1:34][p1_sim[,1:34] >= 0] <- 1
p1_sim[,1:34] <- p1_sim[,1:34] * p1_0

p2_sim <- p2_YOY
p2_sim[,1:34][p2_sim[,1:34] >= 0] <- 1
p2_sim[,1:34] <- p2_sim[,1:34] * p2_0

p3_sim <- p3_YOY
p3_sim[,1:34][p3_sim[,1:34] >= 0] <- 1
p3_sim[,1:34] <- p3_sim[,1:34] * p3_0


######################################
## Run model
# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = sample_areas_wide,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_YOY = p1_sim,
                  p2_YOY = p2_sim, 
                  p3_YOY = p3_sim,
                  Sources = p1_sim$Source)


# Parameters to save
jags_params <- c("beta.cov", "s2.beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", 
                 "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_sim[i,t] + p2_sim[i,t] + p3_sim[i,t])),
                                                rpois(1, lambda = 300),
                                                (p1_sim[i,t] + p2_sim[i,t] + p3_sim[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, -0.5, 0.01),
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
YOY_nMix_full_sim <- jagsUI::jags(data = jags_data,
                                  parameters.to.save = jags_params,
                                  model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_full.jags",
                                  n.chains = nc,
                                  n.iter = ni,
                                  n.burnin = nb,
                                  n.thin = nt,
                                  parallel = T,
                                  inits = init_vals)

YOY_nMix_full_sim_params <- MCMCsummary(YOY_nMix_full_sim,
                                        HPD = T)

# What was the range of covariate effects on YOY abundance?
YOY_nMix_full_sim_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta.cov\\[3")) %>% 
  .[,2] %>% 
  range()
