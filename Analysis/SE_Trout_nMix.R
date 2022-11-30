# George Valentine
# Bayesian Spatial Synchrony Analysis

# Updates from v13:
  # Remove versioning of model files - now git can take care of it (this would have been v14)
  # Move into organized project file for writing and version control
  # save all model objects and results to Analysis and Results folders

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
library(here)         # for smart file paths - eliminates exact paths


# Load data

# Spatiotemporal Covariate Data
SE_COMID_temp_covars <- fread(here("Data", "SE_Trout_COMID_temp_covars.csv"))
SE_COMID_flow_covars <- fread(here("Data", "SE_Trout_COMID_flow_covars.csv"))


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

# create a crosswalk file to associate the agency with the source of the data
source_agency_crosswalk <- data.frame(Source = unique(SE_Site_Final$Source),
                                      Agency = c("GA DNR",
                                                 "NPS GSMNP",
                                                 "ME DFW",
                                                 "MD DNR",
                                                 "NCWRC",
                                                 "NPS Shen",
                                                 "USGS",
                                                 "Clemson U",
                                                 "TWRA",
                                                 "TWRA",
                                                 "USFWS",
                                                 "VT FWD",
                                                 "VA DWR",
                                                 "WV DNR"))
  

########################################################

# COMID/Site data for model
COMID_data <- SE_Site_Final %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         Lat <= 39.716667, # filter for just sites south of the Mason-Dixon Line
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
  filter(SiteID %in% BKT_Sites$SiteID) %>%  # Filter to sites with records of BKT
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  dplyr::select(COMID, Year, Source) %>% 
  unique() %>% 
  filter(COMID %in% COMID_data$COMID, # Filter because we can only use data from COMIDs with area and coordinates
         !is.na(COMID),
         Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015)  # Filter to the latest year that we have flow data

# What was the percentage of single- vs multipass electrofishing?
passes_pcts.table <- SE_Sample_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  filter(COMID %in% COMID_data$COMID,  # Filter because we can only use data from COMIDs with area and coordinates
         !is.na(NumPasses)) %>%
  mutate(Multipass = ifelse(NumPasses > 1, 1, 0)) %>% 
  group_by(Multipass) %>% 
  summarise(Count = n()) %>% 
  mutate(Percent = Count/sum(Count) * 100) %>% 
  as.data.frame()

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
sources <- sort(unique(p1_YOY$Source))

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
  mutate(Mean_Max_Summer_Temp_Scaled = c(scale(Mean_Max_Summer_Temp))) %>% # center and scale the covariate
  dplyr::select(-Mean_Max_Summer_Temp) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% (YOY_BKT_passCounts_SampleYears$Year - 1)) %>% # filter for years which we have trout data, minus one year b/c we are using temp from the prior year
  pivot_wider(names_from = Year,
              values_from = Mean_Max_Summer_Temp_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Mean_Max_Summer_Temp_Scaled <- Mean_Max_Summer_Temp_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Winter flow
Max_0.9Q_WinterFlow_Scaled <- SE_COMID_flow_covars %>%
  mutate(Max_0.9Q_WinterFlow_Scaled = c(scale(Max_0.9Q_WinterFlow))) %>% # center and scale the covariate
  dplyr::select(-Max_0.9Q_WinterFlow) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_SampleYears$Year) %>% # filter for years which we have trout data
  dplyr::select(-Max_0.9Q_SpringFlow) %>% 
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_WinterFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_WinterFlow_Scaled <- Max_0.9Q_WinterFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Spring flow
Max_0.9Q_SpringFlow_Scaled <- SE_COMID_flow_covars %>% 
  mutate(Max_0.9Q_SpringFlow_Scaled = c(scale(Max_0.9Q_SpringFlow))) %>% # center and scale the covariate
  dplyr::select(-Max_0.9Q_SpringFlow) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_SampleYears$Year) %>% # filter for years which we have trout data
  dplyr::select(-Max_0.9Q_WinterFlow) %>% 
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
sink("Analysis/nMix_JAGS_files/YOY_BKT_nMix_full.jags")
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
YOY_BKT_nMix_full <- jagsUI::jags(data = jags_data,
                                parameters.to.save = jags_params,
                                model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_full.jags",
                                n.chains = nc,
                                n.iter = ni,
                                n.burnin = nb,
                                n.thin = nt,
                                parallel = T,
                                inits = init_vals)

#YOY_BKT_nMix_full_params <- as.data.frame(YOY_BKT_nMix_full$summary)
YOY_BKT_nMix_full_params <- MCMCsummary(YOY_BKT_nMix_full,
                                        HPD = T)


## Adults
sink("Analysis/nMix_JAGS_files/Adult_BKT_nMix_full.jags")
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
Adult_BKT_nMix_full <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

#Adult_BKT_nMix_full_params <- as.data.frame(Adult_BKT_nMix_full$summary)
Adult_BKT_nMix_full_params <- MCMCsummary(Adult_BKT_nMix_full,
                                        HPD = T)

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
YOY_BKT_nMix_full_N <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_full.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

#YOY_BKT_nMix_full_N_params <- as.data.frame(YOY_BKT_nMix_full_N$summary)
#MCMCtrace(YOY_BKT_nMix_full_N_params, params = c("s2.gam[86]", "alpha[91]"), ISB = F, pdf = F, Rhat = T)
YOY_BKT_nMix_full_N_params <- MCMCsummary(YOY_BKT_nMix_full_N,
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
YOY_BKT_nMix_full_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

#YOY_BKT_NMix_full_S_params <- as.data.frame(YOY_BKT_NMix_full_S$summary)
YOY_BKT_nMix_full_S_params <- MCMCsummary(YOY_BKT_nMix_full_S,
                                          HPD = T)

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
Adult_BKT_nMix_full_N <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

Adult_BKT_nMix_full_N_params <- MCMCsummary(Adult_BKT_nMix_full_N,
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
Adult_BKT_nMix_full_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_full.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

Adult_BKT_nMix_full_S_params <- MCMCsummary(Adult_BKT_nMix_full_S,
                                       HPD = T)
################################
# Partial model with just summer temperature as environmental covariate
## YOY
sink("Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialSummTemp.jags")
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
YOY_BKT_nMix_partialSummTemp <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialSummTemp.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

YOY_BKT_nMix_partialSummTemp_params <- MCMCsummary(YOY_BKT_nMix_partialSummTemp,
                                            HPD = T)

## Adult
sink("Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialSummTemp.jags")
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
Adult_BKT_nMix_partialSummTemp <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialSummTemp.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

Adult_BKT_nMix_partialSummTemp_params <- MCMCsummary(Adult_BKT_nMix_partialSummTemp,
                                                     HPD = T)

################################
# Partial model with just winter flow as environmental covariate

## YOY
sink("Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialWintFlow.jags")
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
YOY_BKT_nMix_partialWintFlow <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialWintFlow.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

YOY_BKT_nMix_partialWintFlow_params <- MCMCsummary(YOY_BKT_nMix_partialWintFlow, HPD = T)

## Adult
sink("Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialWintFlow.jags")
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
Adult_BKT_nMix_partialWintFlow <- jagsUI::jags(data = jags_data,
                                                   parameters.to.save = jags_params,
                                                   model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialWintFlow.jags",
                                                   n.chains = nc,
                                                   n.iter = ni,
                                                   n.burnin = nb,
                                                   n.thin = nt,
                                                   parallel = T,
                                                   inits = init_vals)

Adult_BKT_nMix_partialWintFlow_params <- MCMCsummary(Adult_BKT_nMix_partialWintFlow, HPD = T)

################################
# Partial model with just spring flow as environmental covariate
sink("Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialSprFlow.jags")
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
YOY_BKT_nMix_partialSprFlow <- jagsUI::jags(data = jags_data,
                                                parameters.to.save = jags_params,
                                                model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_partialSprFlow.jags",
                                                n.chains = nc,
                                                n.iter = ni,
                                                n.burnin = nb,
                                                n.thin = nt,
                                                parallel = T,
                                                inits = init_vals)

YOY_BKT_nMix_partialSprFlow_params <- MCMCsummary(YOY_BKT_nMix_partialSprFlow,
                                                  HPD = T)

## Adult
sink("Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialSprFlow.jags")
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
Adult_BKT_nMix_partialSprFlow <- jagsUI::jags(data = jags_data,
                                                  parameters.to.save = jags_params,
                                                  model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_partialSprFlow.jags",
                                                  n.chains = nc,
                                                  n.iter = ni,
                                                  n.burnin = nb,
                                                  n.thin = nt,
                                                  parallel = T,
                                                  inits = init_vals)

Adult_BKT_nMix_partialSprFlow_params <- MCMCsummary(Adult_BKT_nMix_partialSprFlow,
                                                    HPD = T)

################################
# Specify "null" model (no environmental covariates)
sink("Analysis/nMix_JAGS_files/YOY_BKT_nMix_null.jags")
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
YOY_BKT_nMix_null <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "Analysis/nMix_JAGS_files/YOY_BKT_nMix_null.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

YOY_BKT_nMix_null_params <- MCMCsummary(YOY_BKT_nMix_null,
                                        HPD = T)

## Adults
sink("Analysis/nMix_JAGS_files/Adult_BKT_nMix_null.jags")
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
Adult_BKT_nMix_null <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_BKT_nMix_null.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

Adult_BKT_nMix_null_params <- MCMCsummary(Adult_BKT_nMix_null,
                                          HPD = T)

######################################
### Post-Hoc
######################################
# Create a map of sites
US_states <- map_data("state")

sites_map_data <- SE_Site_Final %>% 
  filter(COMID %in% p1_YOY$COMID) %>% 
  left_join(source_agency_crosswalk)

sites_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = sites_map_data, 
             aes(x = Long, y = Lat, color = Agency), alpha = 0.5) +
  coord_map("albers",
            parameters = c(29.5, 45.5),
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       color = "Source") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

######################################
# Create a table of source data

# Get string of sources and their locations
sources2 <- data.frame(Source = sources) %>% 
  left_join(source_agency_crosswalk)

sources.table <- p1_YOY %>% 
  .[,c(1:34, 37)] %>% 
  pivot_longer(cols = 1:34,
               names_to = "Year",
               values_to = "Count") %>% 
  group_by(Source, Year) %>% 
  summarise(Count = sum(Count, na.rm = T)) %>% 
  ungroup(Year) %>% 
  summarise(Data_Range = paste(Year[which.min(Count)], "-", Year[which.max(Count)]),
            NYears_Data = sum(Count > 0)) %>% 
  cbind(Agency = sources2$Agency) %>% 
  dplyr::select(-Source) %>% 
  .[,c("Agency", "Data_Range", "NYears_Data")]

###################
## ICC Values
# Join site data to ICC values
YOY_ICCs <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data[,c(1,3,4)])

# Export for Shiny app
fwrite(YOY_ICCs, "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Synchrony_App/YOY_ICCS.csv")

YOY_ICCs_N <- YOY_BKT_nMix_full_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data_N[,c(1,3,4)])

YOY_ICCs_S <- YOY_BKT_nMix_full_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data_S[,c(1,3,4)])

Adult_ICCs <- Adult_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data[,c(1,3,4)])

Adult_ICCs_N <- Adult_BKT_nMix_full_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data_N[,c(1,3,4)])

Adult_ICCs_S <- Adult_BKT_nMix_full_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC")) %>% 
  cbind(COMID_data_S[,c(1,3,4)])

# Do the most synchronous/asynchronous sites from the N/S show up in the the overall ICCs?
YOY_highestICCs_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = head(arrange(YOY_ICCs_N, desc(mean))), 
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(YOY_ICCs_S, desc(mean))), 
              aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(YOY_ICCs, desc(mean))), 
              aes(x = Long, y = Lat), shape = 3) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

YOY_lowestICCs_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = head(arrange(YOY_ICCs_N, mean)),
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(YOY_ICCs_S, mean)),
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(YOY_ICCs, mean)), 
             aes(x = Long, y = Lat), shape = 3) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

Adult_highestICCs_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = head(arrange(Adult_ICCs_N, desc(mean))), 
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(Adult_ICCs_S, desc(mean))), 
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(Adult_ICCs, desc(mean))), 
             aes(x = Long, y = Lat), shape = 3) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for Adult BKT",
       color = "ICC") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

Adult_lowestICCs_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = head(arrange(Adult_ICCs_N, mean)), 
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(Adult_ICCs_S, mean)), 
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(Adult_ICCs, mean)), 
             aes(x = Long, y = Lat), shape = 3) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 43)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for Adult BKT",
       color = "ICC") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))


# Plot ICC values on a map
# YOY
YOY_ICC_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_ICCs, 
             aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  scale_color_viridis_c() +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

# ggsave("BKT_Nmix_YOY_ICCs_Map.jpg",
#        plot = YOY_ICC_map,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")

# Adult
Adult_ICC_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = Adult_ICCs, 
             aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for Adult BKT",
       color = "ICC") +
  scale_color_viridis_c() +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

# ggsave("BKT_Nmix_Adult_ICCs_Map.jpg",
#        plot = Adult_ICC_map,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")

# Is ICC correlated with any site-level variables?
library(corrplot)
ICC_corr_data <- YOY_ICCs %>% 
  left_join(SE_Site_Final, by = "COMID")

YOY_corrPlot <- corrplot(cor(ICC_corr_data[,c(2,14:18)], method="spearman", use="pairwise.complete.obs"))
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
                        Global_C = c((1 - (YOY_BKT_nMix_full_params['s2.eps',1]/YOY_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (Adult_BKT_nMix_full_params['s2.eps',1]/Adult_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (YOY_BKT_nMix_partialSummTemp_params['s2.eps',1]/YOY_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (Adult_BKT_nMix_partialSummTemp_params['s2.eps',1]/Adult_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (YOY_BKT_nMix_partialWintFlow_params['s2.eps',1]/YOY_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (Adult_BKT_nMix_partialWintFlow_params['s2.eps',1]/Adult_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (YOY_BKT_nMix_partialSprFlow_params['s2.eps',1]/YOY_BKT_nMix_null_params['s2.eps',1])),
                                     (1 - (Adult_BKT_nMix_partialSprFlow_params['s2.eps',1]/Adult_BKT_nMix_null_params['s2.eps',1]))))



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
  C_Vals_Site[i, "C_YOY_allCovs"] <- 1 - (YOY_BKT_nMix_full_params[648 + i, 1]/YOY_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_Adult_allCovs"] <- 1 - (Adult_BKT_nMix_full_params[648 + i, 1]/Adult_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_YOY_partialSummTemp"] <- 1 - (YOY_BKT_nMix_partialSummTemp_params[328 + i, 1]/YOY_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_Adult_partialSummTemp"] <- 1 - (Adult_BKT_nMix_partialSummTemp_params[328 + i, 1]/Adult_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_YOY_partialWintFlow"] <- 1 - (YOY_BKT_nMix_partialWintFlow_params[328 + i, 1]/YOY_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_Adult_partialWintFlow"] <- 1 - (Adult_BKT_nMix_partialWintFlow_params[328 + i, 1]/Adult_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_YOY_partialSprFlow"] <- 1 - (YOY_BKT_nMix_partialSprFlow_params[328 + i, 1]/YOY_BKT_nMix_null_params[168 + i, 1])
  C_Vals_Site[i, "C_Adult_partialSprFlow"] <- 1 - (Adult_BKT_nMix_partialSprFlow_params[328 + i, 1]/Adult_BKT_nMix_null_params[168 + i, 1])
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
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
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
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
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
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
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
# ggsave("BKT_Nmix_Compound_C_Vals_SummTemp_Map.jpg",
#        plot = C_Val_YOY_SummTemp_map,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")
# ggsave("BKT_Nmix_Compound_C_Vals_WintFlow_Map.jpg",
#        plot = C_Val_YOY_WintFlow_map,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")
# ggsave("BKT_Nmix_Compound_C_Vals_SprFlow_Map.jpg",
#        plot = C_Val_YOY_SprFlow_map,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")


################
# Summarize covariate effects in table
Cov_Effects <- data.frame(
  Covariate = rep(c("Mean 0.9Q Summer Air Temperature (Year t-1)",
                "Max 0.9Q Winter Flow (Year t)",
                "Max 0.9Q Spring Flow (Year t)"), times = 6),
  Life_Stage = rep(c(rep("YOY", times = 3),
                    rep("Adult", times = 3)), times = 3),
  Subregion = c(rep("N+S", times = 6),
                rep("N", times = 6),
                rep("S", times = 6))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(YOY_BKT_nMix_full_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4],
              Adult_BKT_nMix_full_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4],
              YOY_BKT_nMix_full_N_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4],
              Adult_BKT_nMix_full_N_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4],
              YOY_BKT_nMix_full_S_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4],
              Adult_BKT_nMix_full_S_params[c("mu.beta.cov[1]", "mu.beta.cov[2]", "mu.beta.cov[3]"),1:4]))

# Reorder the subregions so the the N+S region plots first
Cov_Effects$Subregion <- factor(Cov_Effects$Subregion, c("N+S", "N", "S"))

# and make a plot to visualize
cov_effects.plot <- ggplot(data = Cov_Effects) +
  geom_pointrange(aes(x = Covariate,
                      y = mean,
                     ymin = `95%_HPDL`, #includes 95% highest density intervals
                     ymax = `95%_HPDU`,
                     #color = Life_Stage,
                     linetype = Subregion),
                 position = position_dodge(.35),
                 size = 0.5,
                 fatten = 1) +
  facet_wrap(~ Life_Stage) +
  # scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  labs(#title = "Covariate Effects on Log Density of BKT",
       #color = "Life Stage",
       x = element_blank(),
       y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 11)) +
  theme_classic() + 
  theme(text = element_text(family =  "serif"),
        legend.key.size = unit(2,"line"))

# # save plot
# ggsave("BKT_Nmix_EnvCov_PstrEffects.jpg",
#        plot = cov_effects_plot,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")


################
# Summarize detection probability in table

Detect_probs <- data.frame(
  Agency = rep(sources2$Agency, times = 2),
  Life_Stage = c(rep("YOY", times = 8),
                rep("Adult", times = 8))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(YOY_BKT_nMix_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),1:4],
              Adult_BKT_nMix_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),1:4]))

# and make a plot to visualize
Detect_probs.plot <- ggplot(data = Detect_probs) +
  geom_linerange(aes(x = Agency,
                     ymin = `95%_HPDL`, #includes 95% highest density intervals
                     ymax = `95%_HPDU`,
                     color = Life_Stage),
                 position = position_dodge(.25),
                 size = 0.5) +
  geom_point(aes(x = Agency,
                 y = mean,
                 color = Life_Stage),
             position = position_dodge(.25)) +
  scale_color_brewer(palette = "Dark2") +
  labs(#title = "Covariate Effects on Log Density of BKT",
    color = "Life Stage",
    x = element_blank(),
    y = element_blank()) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

# save plot
# ggsave("BKT_Nmix_SourcePs.jpg",
#        plot = Detect_probs_plot,
#        path = "C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/Writing/Figures/",
#        width = 500,
#        height = 350,
#        units = "mm",
#        scale = 0.25,
#        dpi = "retina")


########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here("results", "v1.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
