# George Valentine
# Bayesian Spatial Synchrony Analysis v13

# v13 updates:
  # let tau.beta.covs differ by covariate m
  # remove hyperparameters mu.beta[b]
  # add random spatial effect
  # replace gamma priors for taus with uniform priors for SDs, then calculate variance and then precision
  # Add posterior predictive check (PPC) to full models

# Load required packages
library(tidyverse)
library(data.table)
library(jagsUI)
library(lubridate)
library(truncnorm)  # For creating truncated normal distributions
library(RColorBrewer) # colors for plotting
library(MCMCvis)      # For summarizing MCMC outputs
library(HDInterval)   # For calculating 95% credible intervals
library(gridExtra)    # For combining plots


# Load data

# Spatiotemporal Covariate Data
SE_COMID_temp_covars <- fread("SE_COMID_temp_covars.csv")
SE_COMID_flow_covars <- fread("SE_COMID_flow_covars.csv")

# Trout data
# Set working directory for trout data files
filepath <- "/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/Compiled Trout Data (All sources)/"

SE_Ind_Final <- fread(file = paste0(filepath, "/SE_Ind_Final.csv"))
SE_Sample_Final <- fread(file = paste0(filepath, "/SE_Sample_Final.csv"))
SE_Site_Final <- fread(file = paste0(filepath, "/SE_Site_Final.csv"))

# Some sites don't have any BKT. Make a dataframe of just the sites that have >10% BKT. This is Matt Kulp's suggestion rather than using sites with just any BKT
BKT_Sites <- SE_Ind_Final %>% 
  group_by(SiteID,
           SPP) %>% 
  dplyr::summarise(count = n()) %>% 
  ungroup(SPP) %>% 
  mutate(pct_BKT = prop.table(count) * 100) %>% 
  filter(SPP == "BKT",
         pct_BKT >= 10)
  

########################################################

# COMID/Site data for model
COMID_data <- SE_Site_Final %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         Lat <= 43, # filter for just sites south of the north border of PA
         COMID %in% SE_COMID_flow_covars$COMID, # filter for COMIDs for which we have flow data - there are several sites in Maine with no flow
         COMID %in% SE_COMID_temp_covars$COMID) %>% # same goes for temperature - there's one site in Maine with no temp
  group_by(COMID) %>% 
  summarise(COMID_Area = sum(Length_m * Width_m), # Calculates the sum of site areas within the stream segment
            Lat = mean(Lat),
            Long = mean(Long)) %>% 
  na.omit() %>%  # remove rows for sites that don't have area, lat, or long measurements
  mutate(Lat_Scaled = c(scale(Lat)), # use default mean 0, SD 1
         Long_Scaled = c(scale(Long)))
  
# Make a dataframe of the COMIDs and years when sampling happened
SampleYears <- SE_Sample_Final %>% 
  filter(SiteID %in% BKT_Sites$SiteID,) %>%  # Filter to sites with records of BKT
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  dplyr::select(COMID, Year, Source) %>% 
  unique() %>% 
  filter(COMID %in% COMID_data$COMID, # Filter because we can only use data from COMIDs with area and coordinates
         !is.na(COMID),
         Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015)  # Filter to the latest year that we have flow data

## YOY
# Tally YOY counts by pass at each COMID, year combination
YOY_BKT_passCounts <- SE_Ind_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         COMID %in% COMID_data$COMID, # Filter because we can only use data from COMIDs with all the covariates
         Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015) %>% # Filter to the latest year that we have flow data
  group_by(COMID,
           Year,
           Source) %>% 
  summarise(P1_Count_YOY_BKT = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 1), # Filter here for YOY BKT
            P2_Count_YOY_BKT = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 2),
            P3_Count_YOY_BKT = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 3)) %>% 
  .[!duplicated(.[,c(1:2,4:6)]),] # remove any duplicate rows not already cleaned from the data

# Join sample years where YOY BKT were collected to all sample years
  # This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had not fish
YOY_BKT_passCounts_SampleYears <- SampleYears %>% 
  left_join(YOY_BKT_passCounts)

# Matt Kulp sent a file that includes a list of sites where exotic salmonids were removed and BKT were stocked to restore the stream.
# We don't want
GSMNP_Restored_Sites <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/NPS-GSMNP Trout Data (Matt Kulp)/GSMNP Restored Sites from M Kulp 4.22.csv")

GSMNP_Restored_Sites <- GSMNP_Restored_Sites %>% 
  filter(!Kulp_Notes == "")

# see if those restored sites overlap with what we're about to run the synchrony analysis on
YOY_BKT_passCounts_SampleYears %>% 
  filter(COMID %in% GSMNP_Restored_Sites$COMID) %>% 
  distinct(COMID)
# 3 of the restored sites are in the dataset. Filter them out
YOY_BKT_passCounts_SampleYears <- YOY_BKT_passCounts_SampleYears %>% 
  filter(!COMID %in% GSMNP_Restored_Sites$COMID)

# YOY abundance data for model
# Now make a separate, wide data frame for each pass
p1_YOY <- YOY_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P2_Count_YOY_BKT,
                -P3_Count_YOY_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_YOY_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

# save a vector of the names of the sources
sources <- unique(p1_YOY$Source)

# and change the Source column to numeric 
p1_YOY <- p1_YOY %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% 
  relocate(COMID, Source, .after = last_col()) # move the info columns to the end of the df to allow subsetting by year in the model
  #column_to_rownames(var = "COMID")

p2_YOY <- YOY_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P1_Count_YOY_BKT,
                -P3_Count_YOY_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_YOY_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())
  #column_to_rownames(var = "COMID")

p3_YOY <- YOY_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P1_Count_YOY_BKT,
                -P2_Count_YOY_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_YOY_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())
 #column_to_rownames(var = "COMID")

## Adults
# Tally adult counts by pass at each COMID, year combination
Adult_BKT_passCounts <- SE_Ind_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         COMID %in% COMID_data$COMID, # Filter because we can only use data from COMIDs with all the covariates
         Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015) %>% # Filter to the latest year that we have flow data
  group_by(COMID,
           Year,
           Source) %>% 
  summarise(P1_Count_adult_BKT = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 1), # Filter here for adult BKT
            P2_Count_adult_BKT = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 2),
            P3_Count_adult_BKT = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 3)) %>% 
  .[!duplicated(.[,c(1:2,4:6)]),] # remove any duplicate rows not already cleaned from the data

# Join sample years where adult BKT were collected to all sample years
# This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had not fish
Adult_BKT_passCounts_SampleYears <- SampleYears %>% 
  left_join(Adult_BKT_passCounts)

# Matt Kulp sent a file that includes a list of sites where exotic salmonids were removed and BKT were stocked to restore the stream.
# Obviously we don't want
# see if those restored sites overlap with what we're about to run the synchrony analysis on
Adult_BKT_passCounts_SampleYears %>% 
  filter(COMID %in% GSMNP_Restored_Sites$COMID) %>% 
  distinct(COMID)
# 3 of the restored sites are in the dataset. Filter them out
Adult_BKT_passCounts_SampleYears <- Adult_BKT_passCounts_SampleYears %>% 
  filter(!COMID %in% GSMNP_Restored_Sites$COMID)

# adult abundance data for model
# Now make a separate, wide data frame for each pass
p1_adult <- Adult_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P2_Count_adult_BKT,
                -P3_Count_adult_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_adult_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())
  #column_to_rownames(var = "COMID")

p2_adult <- Adult_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P1_Count_adult_BKT,
                -P3_Count_adult_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_adult_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())
  #column_to_rownames(var = "COMID")

p3_adult <- Adult_BKT_passCounts_SampleYears %>% 
  dplyr::select(-P1_Count_adult_BKT,
                -P2_Count_adult_BKT) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_adult_BKT) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,-1]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())
  #column_to_rownames(var = "COMID")

## COVARIATES
# Make wide dataframes of spatiotemporal covariates
# Temperature
Mean_Max_Summer_Temp_Scaled <- SE_COMID_temp_covars %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         # Year >= 1980,  # Constrain to the latest year for which we have flow data
         # Year != 1982) %>% # Filter out data from 1982, for which there is no trout data
         Year %in% (YOY_BKT_passCounts_SampleYears$Year - 1)) %>% # filter for years which we have trout data, minus one year b/c we are using temp from the prior year
  pivot_wider(names_from = Year,
              values_from = Mean_Max_Summer_Temp_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Mean_Max_Summer_Temp_Scaled <- Mean_Max_Summer_Temp_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Winter flow
Max_0.9Q_WinterFlow_Scaled <- SE_COMID_flow_covars %>%
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         # Year >= 1980,  # Constrain to the latest year for which we have flow data
         # Year != 1982) %>% # Filter out data from 1982, for which there is no trout data
         Year %in% YOY_BKT_passCounts_SampleYears$Year) %>% # filter for years which we have trout data
  dplyr::select(-Max_0.9Q_SpringFlow_Scaled) %>% 
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_WinterFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_WinterFlow_Scaled <- Max_0.9Q_WinterFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Spring flow
Max_0.9Q_SpringFlow_Scaled <- SE_COMID_flow_covars %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         # Year >= 1980,  # Constrain to the latest year for which we have flow data
         # Year != 1982) %>% # Filter out data from 1982, for which there is no trout data
         Year %in% YOY_BKT_passCounts_SampleYears$Year) %>% # filter for years which we have trout data
  dplyr::select(-Max_0.9Q_WinterFlow_Scaled) %>% 
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_SpringFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_SpringFlow_Scaled <- Max_0.9Q_SpringFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Filter the COMID data for just those with fish observations in the time frame of interest
# use a right join to get duplicates for the COMIDs with multiple sources of data
COMID_data <- COMID_data %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Set sample sizes
nReps <- nrow(p1_YOY)
nYears <- ncol(p1_YOY) - 3 # subtract 3 b/c the last three columns have info, not counts
nSources <- length(unique(p1_YOY$Source))

########################################
# Specify "full" models (includes all environmental covariates) for YOY and adults

## YOY
sink("YOY_BKT_nMix_v13_full.jags")
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
      
    for (i in 1:nReps){
      
      # beta.cov.i - Site-specific coefficients for mean max summer temp, max 0.9Q winter flow, max 0.9Q spring flow
      beta.cov[m,i] ~ dnorm(mu.beta.cov[m], tau.beta.cov[m])
    }
  }
  
  
  ## Random Effects
  # tau.epsilon: precision parameter for epsilon.t
  sd.eps ~ dunif(0, 10)
  tau.eps <- 1/(sd.eps^2)
  #tau.eps ~ dgamma(2, 2)
    
  # epsilon.t - first order random effect
  for (t in 1:nYears) { 
    eps[t] ~ dnorm(0, tau.eps)
  }
  
    
  # gamma.it - second order random effect
  for (i in 1:nReps) {
  
    # tau.gamma: precision parameter for gamma.it
    sd.gam[i] ~ dunif(0, 10)
    tau.gam[i] <- 1/(sd.gam[i]^2)
    #tau.gam[i] ~ dgamma(2, 2)
    
    for (t in 1:nYears) {
      gam[i,t] ~ dnorm(0, tau.gam[i])
    }
  }
  
  ## Detection probability
  # p.j - detection probability for each source
  for (j in 1:nSources) {
      p[j] ~ dbeta(12, 12)
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta.cov[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta.cov[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t] + eps[t] + gam[i,t]
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
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_YOY_PPC[i,t] ~ dbin(p[Sources[i]], N.YOY[i,t])
    }
  }
  
  # Get means, CVs of data and predicted data
  mean_p1_YOY <- mean(p1_YOY[,1:nYears])
  CV_p1_YOY <- sd(p1_YOY[,1:nYears])/mean(p1_YOY[,1:nYears])
  mean_p1_YOY_PPC <- mean(p1_YOY_PPC[,1:nYears])
  CV_p1_YOY_PPC <- sd(p1_YOY_PPC[,1:nYears])/mean(p1_YOY_PPC[,1:nYears])
  
  # Calculate p values
  pval.mean_p1 <- step(mean_p1_YOY_PPC - mean_p1_YOY)
  pval.CV_p1 <- step(CV_p1_YOY_PPC - CV_p1_YOY)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_YOY = p1_YOY,
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.YOY", 
                 "pval.mean_p1", "pval.CV_p1")

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
                             mu.beta.cov = rnorm(3, 0, 0.01),
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

# Fit Model
YOY_BKT_nMix_v13_full <- jagsUI::jags(data = jags_data,
                                parameters.to.save = jags_params,
                                model.file = "YOY_BKT_nMix_v13_full.jags",
                                n.chains = nc,
                                n.iter = ni,
                                n.burnin = nb,
                                n.thin = nt,
                                parallel = T,
                                inits = init_vals)

v13_YOY_full_params <- as.data.frame(YOY_BKT_nMix_v13_full$summary)
view(v13_YOY_full_params)
#traceplot(YOY_BKT_nMix_v13_full, parameters = "mu.beta.cov[2]")
nrow(YOY_BKT_nMix_v13_full$sims.list$p)


# Highest density intervals for parameters
p1_95 <- as.data.frame(MCMCpstr(YOY_BKT_nMix_v13_full, params = "p", func = function(x)hdi(x, 0.95)))


## Adults
sink("Adult_BKT_nMix_v13_full.jags")
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
    
    for (i in 1:nReps){
      
      # beta.cov.i - Site-specific coefficients for mean max summer temp, max 0.9Q winter flow, max 0.9Q spring flow
      beta.cov[m,i] ~ dnorm(mu.beta.cov[m], tau.beta.cov[m])
    }
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
    p[j] ~ dbeta(18, 10) # gives mean detection of ~0.65
  }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      
      # Data
      N.adult[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta.cov[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta.cov[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_adult[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
      # Pass 2
      p2_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t]))
      # Pass 3
      p3_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
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
    ICC.adult[i] <- s2.eps/(s2.eps + s2.gam[i])
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_adult = p1_adult, 
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.adult")

# create and populate an array of initial values for N.Adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                          rpois(1, lambda = 300),
                                          (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 1.54)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.65, times = nSources),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# Fit Model
Adult_BKT_nMix_v13_full <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Adult_BKT_nMix_v13_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_Adult_full_params <- as.data.frame(Adult_BKT_nMix_v13_full$summary)
view(v13_Adult_full_params)
#traceplot(Adult_BKT_nMix_v13_full, parameters = "")
nrow(Adult_BKT_nMix_v13_full$sims.list$s2.gam)

################################
# Full Models for north vs south regions of SE US

# Define these regions: we draw the line at Radford, VA in the New River Valley. Per the suggestion of Jake Rash
N_Sites <- SE_Site_Final %>% 
  filter(Lat > 37.13) %>% 
  select(COMID)
S_Sites <- SE_Site_Final %>% 
  filter(Lat <= 37.13) %>% 
  select(COMID)

# filter data to these regions
p1_YOY_N <- p1_YOY %>% 
  filter(COMID %in% N_Sites$COMID)
p2_YOY_N <- p2_YOY %>% 
  filter(COMID %in% N_Sites$COMID)
p3_YOY_N <- p3_YOY %>% 
  filter(COMID %in% N_Sites$COMID)

p1_adult_N <- p1_adult %>% 
  filter(COMID %in% N_Sites$COMID)
p2_adult_N <- p2_adult %>% 
  filter(COMID %in% N_Sites$COMID)
p3_adult_N <- p3_adult %>% 
  filter(COMID %in% N_Sites$COMID)

p1_YOY_S <- p1_YOY %>% 
  filter(COMID %in% S_Sites$COMID)
p2_YOY_S <- p2_YOY %>% 
  filter(COMID %in% S_Sites$COMID)
p3_YOY_S <- p3_YOY %>% 
  filter(COMID %in% S_Sites$COMID)

p1_adult_S <- p1_adult %>% 
  filter(COMID %in% S_Sites$COMID)
p2_adult_S <- p2_adult %>% 
  filter(COMID %in% S_Sites$COMID)
p3_adult_S <- p3_adult %>% 
  filter(COMID %in% S_Sites$COMID)

COMID_data_N <- COMID_data %>% 
  filter(COMID %in% N_Sites$COMID)
COMID_data_S <- COMID_data %>% 
  filter(COMID %in% S_Sites$COMID)

Mean_Max_Summer_Temp_Scaled_N <- Mean_Max_Summer_Temp_Scaled %>% 
  filter(COMID %in% N_Sites$COMID)
Max_0.9Q_WinterFlow_Scaled_N <- Max_0.9Q_WinterFlow_Scaled %>% 
  filter(COMID %in% N_Sites$COMID)
Max_0.9Q_SpringFlow_Scaled_N <- Max_0.9Q_SpringFlow_Scaled %>% 
  filter(COMID %in% N_Sites$COMID)

Mean_Max_Summer_Temp_Scaled_S <- Mean_Max_Summer_Temp_Scaled %>% 
  filter(COMID %in% S_Sites$COMID)
Max_0.9Q_WinterFlow_Scaled_S <- Max_0.9Q_WinterFlow_Scaled %>% 
  filter(COMID %in% S_Sites$COMID)
Max_0.9Q_SpringFlow_Scaled_S <- Max_0.9Q_SpringFlow_Scaled %>% 
  filter(COMID %in% S_Sites$COMID)

# Set sample sizes
nReps_N <- nrow(p1_YOY_N)
nReps_S <- nrow(p1_YOY_S)
nSources_N <- length(unique(p1_YOY_N$Source))
nSources_S <- length(unique(p1_YOY_S$Source))

## YOY analysis
# Models can stay the same, we just subset the data

# NORTH
# Bundle data
jags_data <- list(nReps = nReps_N, 
                  nYears = nYears,
                  nSources = nSources_N,
                  Area = COMID_data_N$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_YOY = p1_YOY_N,
                  p2_YOY = p2_YOY_N, 
                  p3_YOY = p3_YOY_N,
                  Sources = as.numeric(as.factor(p1_YOY_N$Source)))


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.YOY")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_N, nYears))
for (i in 1:nReps_N) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_YOY_N[i,t] + p2_YOY_N[i,t] + p3_YOY_N[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_YOY_N[i,t] + p2_YOY_N[i,t] + p3_YOY_N[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps_N, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps_N, 0, 10),
                             gam = array(rnorm(nReps_N * nYears, 0, 10), dim = c(nReps_N, nYears)),
                             p = rep(0.5, times = nSources_N),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 50000
nt <- 1

# Fit Model
YOY_BKT_nMix_v13_full_N <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "YOY_BKT_nMix_v13_full.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

v13_YOY_full_params_N <- as.data.frame(YOY_BKT_nMix_v13_full_N$summary)
traceplot(YOY_BKT_nMix_v13_full_N, parameters = "s2.gam[86]")
MCMCtrace(YOY_BKT_nMix_v13_full_N, params = c("s2.gam[86]", "alpha[91]"), ISB = F, pdf = F, Rhat = T)
view(v13_YOY_full_params_N)

# SOUTH
# Bundle data
jags_data <- list(nReps = nReps_S, 
                  nYears = nYears,
                  nSources = nSources_S,
                  Area = COMID_data_S$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_YOY = p1_YOY_S,
                  p2_YOY = p2_YOY_S, 
                  p3_YOY = p3_YOY_S,
                  Sources = as.numeric(as.factor(p1_YOY_S$Source)))


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.YOY")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_S, nYears))
for (i in 1:nReps_S) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_YOY_S[i,t] + p2_YOY_S[i,t] + p3_YOY_S[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_YOY_S[i,t] + p2_YOY_S[i,t] + p3_YOY_S[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps_S, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps_S, 0, 10),
                             gam = array(rnorm(nReps_S * nYears, 0, 10), dim = c(nReps_S, nYears)),
                             p = rep(0.5, times = nSources_S),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 50000
nt <- 1

# Fit Model
YOY_BKT_NMix_v13_full_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "YOY_BKT_NMix_v13_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_YOY_full_params_S <- as.data.frame(YOY_BKT_NMix_v13_full_S$summary)
view(v13_YOY_full_params_S)
traceplot(YOY_BKT_NMix_v13_full_S, "alpha[16]")

## Adult analysis
# Models can stay the same, we just subset the data

# NORTH
# Bundle data
jags_data <- list(nReps = nReps_N, 
                  nYears = nYears,
                  nSources = nSources_N,
                  Area = COMID_data_N$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_adult = p1_adult_N,
                  p2_adult = p2_adult_N, 
                  p3_adult = p3_adult_N,
                  Sources = as.numeric(as.factor(p1_adult_N$Source)))


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.adult")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_N, nYears))
for (i in 1:nReps_N) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_adult_N[i,t] + p2_adult_N[i,t] + p3_adult_N[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_adult_N[i,t] + p2_adult_N[i,t] + p3_adult_N[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps_N, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps_N, 0, 10),
                             gam = array(rnorm(nReps_N * nYears, 0, 10), dim = c(nReps_N, nYears)),
                             p = rep(0.5, times = nSources_N),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 50000
nt <- 1

# Fit Model
Adult_BKT_nMix_v13_full_N <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Adult_BKT_nMix_v13_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_adult_full_params_N <- MCMCsummary(Adult_BKT_nMix_v13_full_N,
                                      HPD = T)

# SOUTH
# Bundle data
jags_data <- list(nReps = nReps_S, 
                  nYears = nYears,
                  nSources = nSources_S,
                  Area = COMID_data_S$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_adult = p1_adult_S,
                  p2_adult = p2_adult_S, 
                  p3_adult = p3_adult_S,
                  Sources = as.numeric(as.factor(p1_adult_S$Source)))


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam",  "ICC.adult")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_S, nYears))
for (i in 1:nReps_S) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_adult_S[i,t] + p2_adult_S[i,t] + p3_adult_S[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_adult_S[i,t] + p2_adult_S[i,t] + p3_adult_S[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps_S, 0, 0.001),
                             sd.cov = runif(3, 0, 10),
                             mu.beta.cov = rnorm(3, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps_S, 0, 10),
                             gam = array(rnorm(nReps_S * nYears, 0, 10), dim = c(nReps_S, nYears)),
                             p = rep(0.5, times = nSources_S),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 100000
nc <- 3
nb <- 50000
nt <- 1

# Fit Model
Adult_BKT_nMix_v13_full_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Adult_BKT_nMix_v13_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_adult_full_params_S <- MCMCsummary(Adult_BKT_nMix_v13_full_S,
                                       HPD = T)
################################
# Partial model with just summer temperature as environmental covariate
## YOY
sink("YOY_BKT_nMix_v13_partialSummTemp.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(12, 12)
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Mean_Max_Summer_Temp_Scaled[i,t] + eps[t] + gam[i,t]
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
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  p1_YOY = p1_YOY,
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                          rpois(1, lambda = 200),
                                          (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(1, 0, 10),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)

# Fit Model
YOY_BKT_nMix_v13_partialSummTemp <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "YOY_BKT_nMix_v13_partialSummTemp.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_YOY_partialSummTemp_params <- as.data.frame(YOY_BKT_nMix_v13_partialSummTemp$summary)
view(v13_YOY_partialSummTemp_params)
#traceplot(YOY_BKT_nMix_v13_partialSummTemp, parameters = "beta[1]")
nrow(YOY_BKT_nMix_v13_partialSummTemp$sims.list$s2.gam)

## Adult
sink("Adult_BKT_nMix_v13_partialSummTemp.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(18, 10) # gives mean detection of ~0.65
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Mean_Max_Summer_Temp_Scaled[i,t] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
    }
  }

  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nReps){
    s2.gam[i] <- 1/tau.gam[i]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  p1_adult = p1_adult,
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_adult$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                          rpois(1, lambda = 200),
                                          (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(1, 0, 10),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)

# Fit Model
Adult_BKT_nMix_v13_partialSummTemp <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "Adult_BKT_nMix_v13_partialSummTemp.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

v13_Adult_partialSummTemp_params <- as.data.frame(Adult_BKT_nMix_v13_partialSummTemp$summary)
view(v13_Adult_partialSummTemp_params)
#traceplot(Adult_BKT_nMix_v13_partialSummTemp, parameters = "beta[1]")
nrow(Adult_BKT_nMix_v13_partialSummTemp$sims.list$s2.gam)


################################
# Partial model with just winter flow as environmental covariate

## YOY
sink("YOY_BKT_nMix_v13_partialWintFlow.jags")
cat("
model{
  
  ### Priors ###
  
    ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(12, 12)
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Max_0.9Q_WinterFlow_Scaled[i,t] + eps[t] + gam[i,t]
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
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  p1_YOY = p1_YOY,
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)


# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                          rpois(1, lambda = 200),
                                          (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(1, 0, 10),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)

# Fit Model
YOY_BKT_nMix_v13_partialWintFlow <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "YOY_BKT_nMix_v13_partialWintFlow.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

v13_YOY_partialWintFlow_params <- as.data.frame(YOY_BKT_nMix_v13_partialWintFlow$summary)
view(v13_YOY_partialWintFlow_params)
#traceplot(YOY_BKT_nMix_v13_partialWintFlow, parameters = "beta[1]")
nrow(YOY_BKT_nMix_v13_partialWintFlow$sims.list$s2.gam)

## Adult
sink("Adult_BKT_nMix_v13_partialWintFlow.jags")
cat("
model{
  
  ### Priors ###
  
    ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(18, 10) # gives mean detection of ~0.65
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Max_0.9Q_WinterFlow_Scaled[i,t] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
    }
  }

  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nReps){
    s2.gam[i] <- 1/tau.gam[i]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  p1_adult = p1_adult, 
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                            rpois(1, lambda = 200),
                                            (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov= runif(1, 0, 10),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)

# Fit Model
Adult_BKT_nMix_v13_partialWintFlow <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "Adult_BKT_nMix_v13_partialWintFlow.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

v13_Adult_partialWintFlow_params <- as.data.frame(Adult_BKT_nMix_v13_partialWintFlow$summary)
view(v13_Adult_partialWintFlow_params)
#traceplot(Adult_BKT_nMix_v13_partialWintFlow, parameters = "beta[1]")

################################
# Partial model with just spring flow as environmental covariate
sink("YOY_BKT_nMix_v13_partialSprFlow.jags")
cat("
model{
  
  ### Priors ###
  
    ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(12, 12)
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Max_0.9Q_SpringFlow_Scaled[i,t] + eps[t] + gam[i,t]
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
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_YOY = p1_YOY, 
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                          rpois(1, lambda = 200),
                                          (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(1, 0, 10),
                             mu.beta = rnorm(1, -3.5, 0.001),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)

# Fit Model
YOY_BKT_nMix_v13_partialSprFlow <- jagsUI::jags(data = jags_data,
                                                parameters.to.save = jags_params,
                                                model.file = "YOY_BKT_nMix_v13_partialSprFlow.jags",
                                                n.chains = nc,
                                                n.iter = ni,
                                                n.burnin = nb,
                                                n.thin = nt,
                                                parallel = T,
                                                inits = init_vals)

v13_YOY_partialSprFlow_params <- as.data.frame(YOY_BKT_nMix_v13_partialSprFlow$summary)
view(v13_YOY_partialSprFlow_params)
#traceplot(YOY_BKT_nMix_v13_partialSprFlow, parameters = "beta[1]")

## Adult
sink("Adult_BKT_nMix_v13_partialSprFlow.jags")
cat("
model{
  
  ### Priors ###
  
    ## Site fixed effect
  for (i in 1:nReps){
    alpha[i] ~ dnorm(0, 0.001)
  }

  # Beta for summer temperature (site-specific)
  
  # mu.beta.cov - mean parameter for beta.covs
  mu.beta.cov ~ dnorm(0, 0.01)
  
  # tau.beta.cov - precision parameter for beta.covs
  # sd parameter for tau.beta.covs
  sd.cov ~ dunif(0, 10)
  tau.beta.cov <- 1/(sd.cov^2)
    
  for (i in 1:nReps){
    # beta.cov.i - Site-specific coefficient for mean max summer temp
    beta.cov[i] ~ dnorm(mu.beta.cov, tau.beta.cov)
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
      p[j] ~ dbeta(18, 10) # gives mean detection of ~0.65
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta.cov[i] * Max_0.9Q_SpringFlow_Scaled[i,t] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
    }
  }

  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nReps){
    s2.gam[i] <- 1/tau.gam[i]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_adult = p1_adult, 
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "beta.cov", "mu.beta.cov", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.Adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                            rpois(1, lambda = 200),
                                            (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.cov = runif(1, 0, 10),
                             mu.beta.cov = rnorm(1, 0, 0.01),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)

# Fit Model
Adult_BKT_nMix_v13_partialSprFlow <- jagsUI::jags(data = jags_data,
                                                  parameters.to.save = jags_params,
                                                  model.file = "Adult_BKT_nMix_v13_partialSprFlow.jags",
                                                  n.chains = nc,
                                                  n.iter = ni,
                                                  n.burnin = nb,
                                                  n.thin = nt,
                                                  parallel = T,
                                                  inits = init_vals)

v13_Adult_partialSprFlow_params <- as.data.frame(Adult_BKT_nMix_v13_partialSprFlow$summary)
view(v13_Adult_partialSprFlow_params)
#traceplot(Adult_BKT_nMix_v13_partialSprFlow, parameters = "s2.eps")

################################
# Specify "null" model (no environmental covariates)
sink("YOY_BKT_nMix_v13_null.jags")
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
      p[j] ~ dbeta(12, 12)
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
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
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  p1_YOY = p1_YOY, 
                  p2_YOY = p2_YOY, 
                  p3_YOY = p3_YOY,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.YOY.inits[i,t] <- as.numeric(ifelse(is.na((p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t])),
                                          rpois(1, lambda = 200),
                                          (p1_YOY[i,t] + p2_YOY[i,t] + p3_YOY[i,t] + 1) * 2))
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

# Fit Model
YOY_BKT_nMix_v13_null <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "YOY_BKT_nMix_v13_null.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

v13_YOY_null_params <- as.data.frame(YOY_BKT_nMix_v13_null$summary)
view(v13_YOY_null_params)
#traceplot(YOY_BKT_nMix_v13_null, parameters = "beta[1]")


## Adults
sink("Adult_BKT_nMix_v13_null.jags")
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
      p[j] ~ dbeta(18, 10) # gives mean detection of ~0.65
    }
  
  ## Process
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Sources[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
    }
  }

  ### Derived quantities ###
  
  ## sigma^2.epsilon: random effect 1 variance
  s2.eps <- 1/tau.eps
  
  ## sigma^2gamma: random effect 2 variance
  for (i in 1:nReps){
    s2.gam[i] <- 1/tau.gam[i]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = COMID_data$COMID_Area,
                  p1_adult = p1_adult, 
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_YOY$Source)

# Parameters to save
jags_params <- c("alpha", "p", "s2.eps",  "s2.gam")

# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

# create and populate an array of initial values for N.Adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                            rpois(1, lambda = 200),
                                            (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears, 0, 10),
                             sd.gam = runif(nReps, 0, 10),
                             gam = array(rnorm(nReps * nYears, 0, 10), dim = c(nReps, nYears)),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)

# Fit Model
Adult_BKT_nMix_v13_null <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Adult_BKT_nMix_v13_null.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

v13_Adult_null_params <- as.data.frame(Adult_BKT_nMix_v13_null$summary)
view(v13_Adult_null_params)
#traceplot(Adult_BKT_nMix_v13_null, parameters = "beta[1]")


###################
# Save model posterior outputs
fwrite(v13_YOY_full_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_full_params.csv")
fwrite(v13_Adult_full_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_Adult_full_params.csv")
fwrite(v13_YOY_partialSummTemp_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_partialSummTemp_params.csv")
fwrite(v13_Adult_partialSummTemp_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_Adult_partialSummTemp_params.csv")
fwrite(v13_YOY_partialWintFlow_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_partialWintFlow_params.csv")
fwrite(v13_Adult_partialWintFlow_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_Adult_partialWintFlow_params.csv")
fwrite(v13_YOY_partialSprFlow_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_partialSprFlow_params.csv")
fwrite(v13_Adult_partialSprFlow_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_Adult_partialSprFlow_params.csv")
fwrite(v13_YOY_null_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_null_params.csv")
fwrite(v13_Adult_null_params, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/(v13_Adult_null_params.csv")

###################
## ICC Values
# Join site data to ICC values
YOY_ICCs <- v13_YOY_full_params %>% 
  rownames_to_column(., "param") %>% 
  .[858:1026,] %>% 
  mutate(i = row_number()) %>% 
  cbind(COMID = COMID_data$COMID) %>%
  left_join(SE_Site_Final)

Adult_ICCs <- v13_Adult_full_params %>% 
  rownames_to_column(., "param") %>% 
  .[858:1026,] %>% 
  mutate(i = row_number()) %>% 
  cbind(COMID = COMID_data$COMID) %>%
  left_join(SE_Site_Final)

fwrite(YOY_ICCs, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/YOY_ICCs_v13.csv")

# Plot ICC values on a map
US_states <- map_data("state")

# YOY
YOY_ICC_map <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_ICCs, 
             aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  scale_color_viridis_c() +
  theme_classic()

ggsave("BKT_Nmix_YOY_ICCs_Map.jpg",
       plot = YOY_ICC_map,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")

# Adult
Adult_ICC_map <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = Adult_ICCs, 
             aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       title = "Posterior ICC Means for Adult BKT",
       color = "ICC") +
  scale_color_viridis_c() +
  theme_classic()

ggsave("BKT_Nmix_Adult_ICCs_Map.jpg",
       plot = Adult_ICC_map,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")

# Is ICC correlated with any site-level variables?
library(corrplot)
YOY_corrPlot <- corrplot(cor(YOY_ICCs[,c(2,18:22)], method="spearman", use="pairwise.complete.obs"))
# get values
YOY_corrPlot$corrPos
# mean ICC is not really correlated with any of these site-level covars

Adult_corrPlot <- corrplot(cor(Adult_ICCs[,c(2,18:22)], method="spearman", use="pairwise.complete.obs"))
# get values
Adult_corrPlot$corrPos

v13_YOY_partialWintFlow_params %>% 
  rownames_to_column(., "param") %>% 
  view()
###########################
# Calculate C values for all environmental covariates using posterior means from full and null models

# Global C: posterior mean values
Global_Cs <- data.frame(Covariate = c(rep("All", 2),
                                      rep("Summtemp", 2),
                                      rep("WintFlow", 2),
                                      rep("SprFlow", 2)),
                        Age_Class = rep(c("YOY", "Adult"), 4),
                        Global_C = c((1 - (v13_YOY_full_params['s2.eps',1]/v13_YOY_null_params['s2.eps',1])),
                                     (1 - (v13_Adult_full_params['s2.eps',1]/v13_Adult_null_params['s2.eps',1])),
                                     (1 - (v13_YOY_partialSummTemp_params['s2.eps',1]/v13_YOY_null_params['s2.eps',1])),
                                     (1 - (v13_Adult_partialSummTemp_params['s2.eps',1]/v13_Adult_null_params['s2.eps',1])),
                                     (1 - (v13_YOY_partialWintFlow_params['s2.eps',1]/v13_YOY_null_params['s2.eps',1])),
                                     (1 - (v13_Adult_partialWintFlow_params['s2.eps',1]/v13_Adult_null_params['s2.eps',1])),
                                     (1 - (v13_YOY_partialSprFlow_params['s2.eps',1]/v13_YOY_null_params['s2.eps',1])),
                                     (1 - (v13_Adult_partialSprFlow_params['s2.eps',1]/v13_Adult_null_params['s2.eps',1]))))

## Calculate 95% HDPIs on global C values
# # Get number of iterations in posterior distributions
# n.it <- length(YOY_BKT_nMix_v13_full$sims.list$s2.eps)
# 
# # Create an empty dataframe that stores posterior distributions for overall C values
# C_Vals <- data.frame(C_YOY_allCovs = numeric(),
#                      C_Adult_allCovs = numeric(),
#                      C_YOY_partialSummTemp = numeric(),
#                      C_Adult_partialSummTemp = numeric(),
#                      C_YOY_partialWintFlow = numeric(),
#                      C_Adult_partialWintFlow = numeric(),
#                      C_YOY_partialSprFlow = numeric(),
#                      C_Adult_partialSprFlow = numeric())
# 
# # for() loop that calculates posterior distributions for overall C values 
# for(i in 1:n.it) {
#   C_Vals[i,"C_YOY_allCovs"] <- 1 - (YOY_BKT_nMix_v13_full$sims.list$s2.eps[i]/YOY_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_Adult_allCovs"] <- 1 - (Adult_BKT_nMix_v13_full$sims.list$s2.eps[i]/Adult_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_YOY_partialSummTemp"] <- 1 - (YOY_BKT_nMix_v13_partialSummTemp$sims.list$s2.eps[i]/YOY_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_Adult_partialSummTemp"] <- 1 - (Adult_BKT_nMix_v13_partialSummTemp$sims.list$s2.eps[i]/Adult_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_YOY_partialWintFlow"] <- 1 - (YOY_BKT_nMix_v13_partialWintFlow$sims.list$s2.eps[i]/YOY_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_Adult_partialWintFlow"] <- 1 - (Adult_BKT_nMix_v13_partialWintFlow$sims.list$s2.eps[i]/Adult_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_YOY_partialSprFlow"] <- 1 - (YOY_BKT_nMix_v13_partialSprFlow$sims.list$s2.eps[i]/YOY_BKT_nMix_v13_null$sims.list$s2.eps[i])
#   C_Vals[i,"C_Adult_partialSprFlow"] <- 1 - (Adult_BKT_nMix_v13_partialSprFlow$sims.list$s2.eps[i]/Adult_BKT_nMix_v13_null$sims.list$s2.eps[i])
# }
# 
# # C values for covariates for adults are all <0 b/c those covariates have nonsignificant effects on log density
# # Visualize C values
# C_Vals_Summary <- data.frame(C_Value = rep(c("All Covariates",
#                                              "Mean 0.9Q Summer Air Temperature (Year t-1)",
#                                              "Max 0.9Q Winter Flow (Year t)",
#                                              "Max 0.9Q Spring Flow (Year t)"), 2),
#                                  Age_Class = c(rep("YOY", 4), rep("Adult", 4)),
#                                  Pstr_Mean = c(mean(C_Vals$C_YOY_allCovs),
#                                                mean(C_Vals$C_YOY_partialSummTemp),
#                                                mean(C_Vals$C_YOY_partialWintFlow),
#                                                mean(C_Vals$C_YOY_partialSprFlow),
#                                                mean(C_Vals$C_Adult_allCovs),
#                                                mean(C_Vals$C_Adult_partialSummTemp),
#                                                mean(C_Vals$C_Adult_partialWintFlow),
#                                                mean(C_Vals$C_Adult_partialSprFlow)),
#                                  Pstr_0.025_HDI = c(as.numeric(hdi(C_Vals$C_YOY_allCovs, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialSummTemp, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialWintFlow, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialSprFlow, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_allCovs, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialSummTemp, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialWintFlow, 0.95)[1]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialSprFlow, 0.95)[1])),
#                                  Pstr_0.975_HDI = c(as.numeric(hdi(C_Vals$C_YOY_allCovs, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialSummTemp, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialWintFlow, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_YOY_partialSprFlow, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_allCovs, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialSummTemp, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialWintFlow, 0.95)[2]),
#                                                     as.numeric(hdi(C_Vals$C_Adult_partialSprFlow, 0.95)[2])))
# C_Vals_plot <- ggplot(data = C_Vals_Summary) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
#   geom_linerange(aes(x = C_Value,
#                      ymin = Pstr_0.025_HDI, #includes 95% highest density intervals
#                      ymax = Pstr_0.975_HDI,
#                      color = Age_Class),
#                  position = position_dodge(.25),
#                  size = 0.75) +
#   geom_point(aes(x = C_Value,
#                  y = Pstr_Mean,
#                  color = Age_Class),
#              position = position_dodge(.25)) +
#   labs(title = "Covariate Effects on Synchrony in\n Log Density of BKT",
#     x = element_blank(),
#     y = "Posterior C Value",
#     color = "Age Class") +
#   scale_color_brewer(palette = "Dark2") +
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 11)) +
#   theme_classic()

# # Save plot
# ggsave("BKT_Nmix_C_Vals.jpg",
#        plot = C_Vals_plot,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")

# Site-specific Cs
# create an empty array to store the values
# we use an array here because this essentially adds a third dimension (unique site) to the for() loop above
C_Vals_Site <- data.frame(C_YOY_allCovs = numeric(),
                          C_Adult_allCovs = numeric(),
                          C_YOY_partialSummTemp = numeric(),
                          C_Adult_partialSummTemp = numeric(),
                          C_YOY_partialWintFlow = numeric(),
                          C_Adult_partialWintFlow = numeric(),
                          C_YOY_partialSprFlow = numeric(),
                          C_Adult_partialSprFlow = numeric())


for (i in 1:nReps){
  print(i)
  # Calculate C.gams for the given site
  C_Vals_Site[i, "C_YOY_allCovs"] <- 1 - (v13_YOY_full_params[688 + i, 1]/v13_YOY_null_params[178 + i, 1])
  C_Vals_Site[i, "C_Adult_allCovs"] <- 1 - (v13_Adult_full_params[688 + i, 1]/v13_Adult_null_params[178 + i, 1])
  C_Vals_Site[i, "C_YOY_partialSummTemp"] <- 1 - (v13_YOY_partialSummTemp_params[348 + i, 1]/v13_YOY_null_params[178 + i, 1])
  C_Vals_Site[i, "C_Adult_partialSummTemp"] <- 1 - (v13_Adult_partialSummTemp_params[348 + i, 1]/v13_Adult_null_params[178 + i, 1])
  C_Vals_Site[i, "C_YOY_partialWintFlow"] <- 1 - (v13_YOY_partialWintFlow_params[348 + i, 1]/v13_YOY_null_params[178 + i, 1])
  C_Vals_Site[i, "C_Adult_partialWintFlow"] <- 1 - (v13_Adult_partialWintFlow_params[348 + i, 1]/v13_Adult_null_params[178 + i, 1])
  C_Vals_Site[i, "C_YOY_partialSprFlow"] <- 1 - (v13_YOY_partialSprFlow_params[348 + i, 1]/v13_YOY_null_params[178 + i, 1])
  C_Vals_Site[i, "C_Adult_partialSprFlow"] <- 1 - (v13_Adult_partialSprFlow_params[348 + i, 1]/v13_YOY_null_params[178 + i, 1])
}

# bind COMIDs to C.gam values
C_Vals_Site <- data.frame(COMID = COMID_data$COMID,
                     Lat = COMID_data$Lat,
                     Long = COMID_data$Long) %>% 
  cbind(C_Vals_Site)

# Create maps to show posterior means of site-specific C values 
# summer temp
C_Val_YOY_SummTemp_map <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = C_Vals_Site, 
             aes(x = Long, y = Lat, color = C_YOY_partialSummTemp), 
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       title = "Mean 0.9Q Summer Air\n Temperature (Year t-1)",
       color = "C Value") +
  scale_color_viridis_c(limits = c(0,0.8)) +
  theme_classic()

# winter flow
C_Val_YOY_WintFlow_map <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = C_Vals_Site, 
             aes(x = Long, y = Lat, color = C_YOY_partialWintFlow), 
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       title = "Max 0.9Q Winter Flow (Year t)",
       color = "C Value") +
  scale_color_viridis_c(limits = c(0,0.8)) +
  theme_classic()

# spring flow
C_Val_YOY_SprFlow_map <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = C_Vals_Site, 
             aes(x = Long, y = Lat, color = C_YOY_partialSprFlow), 
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       title = "Max 0.9Q Spring Flow (Year t)",
       color = "C Value") +
  scale_color_viridis_c(limits = c(0,0.8)) +
  theme_classic()

# Combine the three maps into one
# Compound_C_Vals_map <- grid.arrange(C_Val_YOY_SummTemp_map, C_Val_YOY_WintFlow_map, C_Val_YOY_SprFlow_map,
#                                     nrow = 1)

# save plots
ggsave("BKT_Nmix_Compound_C_Vals_SummTemp_Map.jpg",
       plot = C_Val_YOY_SummTemp_map,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")
ggsave("BKT_Nmix_Compound_C_Vals_WintFlow_Map.jpg",
       plot = C_Val_YOY_WintFlow_map,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")
ggsave("BKT_Nmix_Compound_C_Vals_SprFlow_Map.jpg",
       plot = C_Val_YOY_SprFlow_map,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")


# C_Vals_Site <- array(NA, dim = c(n.it, nReps, ncol(C_Vals)))
# 
# for (i in 1:n.it){
#   for (j in 1:nReps){
#     C_Vals_Site[i,j,1] <- 1 - (YOY_BKT_nMix_v13_full$sims.list$s2.gam[i]/YOY_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,2] <- 1 - (Adult_BKT_nMix_v13_full$sims.list$s2.gam[i]/Adult_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,3] <- 1 - (YOY_BKT_nMix_v13_partialSummTemp$sims.list$s2.gam[i]/YOY_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,4] <- 1 - (Adult_BKT_nMix_v13_partialSummTemp$sims.list$s2.gam[i]/Adult_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,5] <- 1 - (YOY_BKT_nMix_v13_partialWintFlow$sims.list$s2.gam[i]/YOY_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,6] <- 1 - (Adult_BKT_nMix_v13_partialWintFlow$sims.list$s2.gam[i]/Adult_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,7] <- 1 - (YOY_BKT_nMix_v13_partialSprFlow$sims.list$s2.gam[i]/YOY_BKT_nMix_v13_null$sims.list$s2.gam[i])
#     C_Vals_Site[i,j,8] <- 1 - (Adult_BKT_nMix_v13_partialSprFlow$sims.list$s2.gam[i]/Adult_BKT_nMix_v13_null$sims.list$s2.gam[i])
#   }
# }

################
# Summarize covariate effects in table
Cov_Effects <- data.frame(
  Covariate = rep(c("Mean 0.9Q Summer Air Temperature (Year t-1)",
                "Max 0.9Q Winter Flow (Year t)",
                "Max 0.9Q Spring Flow (Year t)"), times = 2),
  Age_Class = c(rep("YOY", times = 3),
                    rep("Adult", times = 3))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(v13_YOY_full_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),],
              v13_Adult_full_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),])) %>% 
  # and 95% highest density intervals
  cbind(rbind(as.data.frame(MCMCpstr(YOY_BKT_nMix_v13_full, params = "mu.beta.cov", func = function(x)hdi(x, 0.95))),
              as.data.frame(MCMCpstr(Adult_BKT_nMix_v13_full, params = "mu.beta.cov", func = function(x)hdi(x, 0.95)))))

# and make a plot to visualize
cov_effects_plot <- ggplot(data = Cov_Effects) +
  geom_linerange(aes(x = Covariate,
                     ymin = mu.beta.cov.lower, #includes 95% highest density intervals
                     ymax = mu.beta.cov.upper,
                     color = Age_Class),
                 position = position_dodge(.25),
                 size = 0.75) +
  geom_point(aes(x = Covariate,
                 y = mean,
                 color = Age_Class),
             position = position_dodge(.25)) +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  labs(#title = "Covariate Effects on Log Density of BKT",
       color = "Age Class",
       x = element_blank(),
       y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 11)) +
  theme_classic()

# save plot
ggsave("BKT_Nmix_EnvCov_PstrEffects.jpg",
       plot = cov_effects_plot,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")


################
# Summarize detection probability in table
# Get string of sources and their locations
sources2 <- data.frame(Source = sources,
                      Location = c("GA DNR", "MD DNR", "SNP", "VA DWR", "NCWRC",  "USGS", "TWRA", "GSMNP"))

Detect_probs <- data.frame(
  Location = rep(sources2$Location, times = 2),
  Age_Class = c(rep("YOY", times = 8),
                rep("Adult", times = 8))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(v13_YOY_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),],
              v13_Adult_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),])) %>% 
  # and 95% highest density intervals
  cbind(rbind(as.data.frame(MCMCpstr(YOY_BKT_nMix_v13_full, params = "p", func = function(x)hdi(x, 0.95))),
              as.data.frame(MCMCpstr(Adult_BKT_nMix_v13_full, params = "p", func = function(x)hdi(x, 0.95)))))

# and make a plot to visualize
Detect_probs_plot <- ggplot(data = Detect_probs) +
  geom_linerange(aes(x = Location,
                     ymin = p.lower, #includes 95% highest density intervals
                     ymax = p.upper,
                     color = Age_Class),
                 position = position_dodge(.25),
                 size = 1) +
  geom_point(aes(x = Location,
                 y = mean,
                 color = Age_Class),
             position = position_dodge(.25)) +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "Covariate Effects on Log Density of BKT",
    color = "Age Class",
    x = element_blank(),
    y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme_classic()

# save plot
ggsave("BKT_Nmix_SourcePs.jpg",
       plot = Detect_probs_plot,
       path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
       width = 500,
       height = 350,
       units = "mm",
       scale = 0.25,
       dpi = "retina")
