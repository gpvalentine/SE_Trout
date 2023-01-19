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
library(PerformanceAnalytics) # For conducting correlation tests


# Load data

# Spatiotemporal Covariate Data
SE_COMID_temp_covars <- fread(here("Data", "SE_Trout_COMID_temp_covars.csv"))
SE_COMID_flow_covars <- fread(here("Data", "SE_Trout_COMID_flow_covars.csv"))


# Trout data
SE_Ind_Final <- fread(here("Data", "SE_Ind_Final.csv"))
SE_Sample_Final <- fread(here("Data", "SE_Sample_Final.csv"))
SE_Site_Final <- fread(here("Data", "SE_Site_Final.csv"))

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
                                      Agency = c("GA Dept. Natural Resources",
                                                 "Great Smoky Mountain Nat'l Park",
                                                 "ME Dept. Inland Fisheries & Wildlife",
                                                 "MD Dept. Natural Resources",
                                                 "NC Wildlife Resources Commission",
                                                 "Shenandoah Nat'l Park",
                                                 "US Geological Survey",
                                                 "Clemson Univ.",
                                                 "TN Wildlife Resources Agency",
                                                 "TN Wildlife Resources Agency",
                                                 "US Fish & Wildlife Service",
                                                 "VT Fish & Wildlife Dept.",
                                                 "VA Dept. Wildlife Resources",
                                                 "WV Dept. Natural Resources"))
  

########################################################

# COMID/Site data for model
sample_areas <- SE_Sample_Final %>% 
  left_join(SE_Site_Final[,c("SiteID", "COMID", "Lat", "Long", "Length_m", "Width_m")], by = "SiteID") %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         Lat <= 39.716667, # filter for just sites south of the Mason-Dixon Line
         COMID %in% SE_COMID_flow_covars$COMID, # filter for COMIDs for which we have flow data - there are several sites in Maine with no flow
         COMID %in% SE_COMID_temp_covars$COMID) %>% # same goes for temperature - there's one site in Maine with no temp
  group_by(COMID,
           Year = year(Date),
           Source) %>% 
  summarise(Area_Sampled = sum(Length_m * Width_m)) %>%  # Calculates the sum of site areas within the stream segment that were sampled that year by that source
  filter(Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015,    # Filter to the latest year that we have flow data
         !is.na(Area_Sampled)) %>% 
  ungroup()

# Make a dataframe of the COMIDs and years when sampling happened
# SampleYears <- SE_Sample_Final %>% 
#   filter(SiteID %in% BKT_Sites$SiteID) %>%  # Filter to sites with records of BKT
#   left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
#   mutate(Year = year(Date)) %>% 
#   dplyr::select(COMID, Year, Source) %>% 
#   unique() %>% 
#   filter(COMID %in% sample_areas$COMID, # Filter because we can only use data from COMIDs with area and coordinates
#          !is.na(COMID),
#          Year >= 1981, # Filter to one year after the earliest year that we have temperature data
#          Year <= 2015)  # Filter to the latest year that we have flow data

# What was the percentage of single- vs multipass electrofishing?
passes_pcts.table <- SE_Sample_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% # create a year column
  inner_join(sample_areas) %>%  # use an inner join to filter for only segment-source-year combos that have areas
  filter(!is.na(NumPasses)) %>%
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
  right_join(sample_areas) %>%  # use a right join to filter for only segment-source-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
  group_by(COMID,
           Year,
           Source) %>% 
  summarise(P1_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 1), # Filter here for YOY BKT
            P2_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 2),
            P3_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 3)) %>% 
  ungroup()
  #.[!duplicated(.[,c(1:2,4:6)]),] # remove any duplicate rows not already cleaned from the data


      # what's the proportion of segments that have multiple sites in them?
      sites_in_segments <- SE_Ind_Final %>% 
        left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
        mutate(Year = year(Date)) %>% 
        right_join(sample_areas) %>%  # use a right join to filter for only segment-source-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
        group_by(COMID,
                 Year) %>% 
        summarise(nSites = length(unique(SiteID)))
      
      summary(sites_in_segments$nSites)
      sum(sites_in_segments$nSites > 1)/nrow(sites_in_segments)

# Matt Kulp sent a file that includes a list of sites where exotic salmonids were removed and BKT were stocked to restore the stream.
# We don't want
GSMNP_Restored_Sites <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/NPS-GSMNP Trout Data (Matt Kulp)/GSMNP Restored Sites from M Kulp 4.22.csv")

GSMNP_Restored_Sites <- GSMNP_Restored_Sites %>% 
  filter(!Kulp_Notes == "")

# see if those restored sites overlap with what we're about to run the synchrony analysis on
YOY_BKT_passCounts%>% 
  filter(COMID %in% GSMNP_Restored_Sites$COMID) %>% 
  distinct(COMID)
# 2 of the restored sites are in the dataset. Filter them out
YOY_BKT_passCounts <- YOY_BKT_passCounts %>%
  filter(!COMID %in% GSMNP_Restored_Sites$COMID)

# YOY abundance data for model
# Now make a separate, wide data frame for each pass
p1_YOY <- YOY_BKT_passCounts %>% 
  dplyr::select(-P2_Count_YOY,
                -P3_Count_YOY) %>%
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_YOY) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

# save a vector of the names of the sources
sources <- sort(unique(p1_YOY$Source))

p2_YOY <- YOY_BKT_passCounts %>% 
  dplyr::select(-P1_Count_YOY,
                -P3_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_YOY) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())

p3_YOY <- YOY_BKT_passCounts %>% 
  dplyr::select(-P1_Count_YOY,
                -P2_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_YOY) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())

## Adults
# Tally adult counts by pass at each COMID, year combination
Adult_BKT_passCounts <- SE_Ind_Final %>% 
  left_join(SE_Site_Final[,c("SiteID", "COMID")]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  right_join(sample_areas) %>%  # use a right join to filter for only segment-source-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
  group_by(COMID,
           Year,
           Source) %>% 
  summarise(P1_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 1), # Filter here for adult BKT
            P2_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 2),
            P3_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 3)) %>% 
  ungroup()
  #.[!duplicated(.[,c(1:2,4:6)]),] # remove any duplicate rows not already cleaned from the data

# Matt Kulp sent a file that includes a list of sites where exotic salmonids were removed and BKT were stocked to restore the stream.
# Obviously we don't want
# see if those restored sites overlap with what we're about to run the synchrony analysis on
Adult_BKT_passCounts %>% 
  filter(COMID %in% GSMNP_Restored_Sites$COMID) %>% 
  distinct(COMID)
# 2 of the restored sites are in the dataset. Filter them out
Adult_BKT_passCounts<- Adult_BKT_passCounts %>% 
  filter(!COMID %in% GSMNP_Restored_Sites$COMID)

# adult abundance data for model
# Now make a separate, wide data frame for each pass
p1_adult <- Adult_BKT_passCounts %>% 
  dplyr::select(-P2_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_adult) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())

p2_adult <- Adult_BKT_passCounts %>% 
  dplyr::select(-P1_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_adult) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())

p3_adult <- Adult_BKT_passCounts %>% 
  dplyr::select(-P1_Count_adult,
                -P2_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_adult) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>% 
  filter(nYears_data >= 5) %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% # change the Source column to numeric 
  relocate(COMID, Source, .after = last_col())

## COVARIATES
# Make wide dataframes of spatiotemporal covariates
# Temperature
Mean_Max_Summer_Temp_Scaled <- SE_COMID_temp_covars %>% 
  group_by(COMID) %>% 
  mutate(Mean_Max_Summer_Temp_Scaled = c(scale(Mean_Max_Summer_Temp))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Mean_Max_Summer_Temp_Scaled) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% (YOY_BKT_passCounts$Year - 1)) %>% # filter for years which we have trout data, minus one year b/c we are using temp from the prior year
  pivot_wider(names_from = Year,
              values_from = Mean_Max_Summer_Temp_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Mean_Max_Summer_Temp_Scaled <- Mean_Max_Summer_Temp_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Winter flow
Max_0.9Q_WinterFlow_Scaled <- SE_COMID_flow_covars %>%
  group_by(COMID) %>% 
  mutate(Max_0.9Q_WinterFlow_Scaled = c(scale(Max_0.9Q_WinterFlow))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Max_0.9Q_WinterFlow_Scaled) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts$Year) %>% # filter for years which we have trout data
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_WinterFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_WinterFlow_Scaled <- Max_0.9Q_WinterFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Spring flow
Max_0.9Q_SpringFlow_Scaled <- SE_COMID_flow_covars %>%
  group_by(COMID) %>% 
  mutate(Max_0.9Q_SpringFlow_Scaled = c(scale(Max_0.9Q_SpringFlow))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Max_0.9Q_SpringFlow_Scaled) %>% 
  filter(COMID %in% p1_YOY$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts$Year) %>% # filter for years which we have trout data
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_SpringFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

    # Are winter and spring flows correlated?
    # flow_scaled <- SE_COMID_flow_covars %>%
    #   group_by(COMID) %>% 
    #   mutate(wintFlow_scaled = c(scale(Max_0.9Q_WinterFlow)),
    #          sprFlow_scaled = c(scale(Max_0.9Q_SpringFlow)))
    # chart.Correlation(flow_scaled[,6:7], method = "spearman")
    # no - they're not highly correlated. Let's keep spring flow in there.

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_SpringFlow_Scaled <- Max_0.9Q_SpringFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY$COMID))

# Filter the sample areas data for just those with fish observations in the time frame of interest
# use a right join to get duplicates for the COMIDs with multiple sources of data
# Then pivot wider to include it in the model
sample_areas_wide <- sample_areas %>% 
  pivot_wider(names_from = Year, 
              values_from = Area_Sampled) %>% 
  arrange(COMID) %>% 
  relocate(COMID, .after = last_col()) %>%  # move the info columns to the end of the df to allow subsetting by year in the model
  right_join(p1_YOY[,c("COMID", "Source")]) %>% # use a right join to filter for areas of the years AND sources where fish were caught
  dplyr::select(-Source)

# the sample areas cannot have missing values, 
# so we impute any segment-year that is missing a sample area as the average sample area for that segment.
# calculate the average sample area of each segment
mean_areas <- (rowSums(sample_areas_wide[,1:ncol(sample_areas_wide)-1], na.rm = T)/apply(!is.na(sample_areas_wide[,1:ncol(sample_areas_wide)-1]), MARGIN = 1, sum))

# then replace the NAs with the average sample area at that segment
for (i in 1:nrow(sample_areas_wide)) {
  sample_areas_wide[i,][is.na(sample_areas_wide[i,])] <- mean_areas[i]
  print(i)
}

# and change the Source column to numeric 
p1_YOY <- p1_YOY %>% 
  mutate(Source = as.numeric(as.factor(Source))) %>% 
  relocate(COMID, Source, .after = last_col()) # move the info columns to the end of the df to allow subsetting by year in the model

# Set sample sizes
nReps <- nrow(p1_YOY)
nYears <- ncol(p1_YOY) - 3 # subtract 3 b/c the last three columns have info, not counts
nSources <- length(unique(p1_YOY$Source))

########################################
# two models each for YOY and adults: one with just climate effects and one with just the two random effects

# YOY model with just env. covariates
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
    
    # mu.beta.cov - mean parameter for betas
    mu.beta[m] ~ dnorm(0, 0.01)

    # tau.beta.cov - precision parameter for betas
    # sd parameter for tau.beta
    sd.beta[m] ~ dunif(0, 10)
    tau.beta[m] <- 1/(sd.beta[m]^2)
    s2.beta[m] <- sd.beta[m]^2
      
    for (i in 1:nReps){
      
      # beta.i - Site-specific coefficients for mean max summer temp, max 0.9Q winter flow, max 0.9Q spring flow
      beta[m,i] ~ dnorm(mu.beta[m], tau.beta[m])
      
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
      
      log(lambda[i,t]) <- alpha[i] + beta[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t]
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
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

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
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 50000
nc <- 3
nb <- 10000
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

# What was the range of covariate effects on YOY abundance?
YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta\\[3")) %>% 
  .[,2] %>% 
  range()

## YOY model with just random effects
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
jags_params <- c("alpha", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1")

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
ni <- 50000
nc <- 3
nb <- 10000
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

YOY_randomEffects_params <- MCMCsummary(YOY_randomEffects, HPD = T)

#MCMCtrace(YOY_randomEffects_params, params = , pdf = F)

### Adults
# Adult model with just env. covariates
sink("Analysis/nMix_JAGS_files/Adult_climateEffects.jags")
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
    
    # mu.beta - mean parameter for betas
    mu.beta[m] ~ dnorm(0, 0.01)

    # tau.beta - precision parameter for betas
    # sd parameter for tau.betas
    sd.beta[m] ~ dunif(0, 10)
    tau.beta[m] <- 1/(sd.beta[m]^2)
    s2.beta[m] <- sd.beta[m]^2
      
    for (i in 1:nReps){
      
      # beta.i - Site-specific coefficients for mean max summer temp, max 0.9Q winter flow, max 0.9Q spring flow
      beta[m,i] ~ dnorm(mu.beta[m], tau.beta[m])
      
    }
  }
  
  ## Detection probability
  # p.j - detection probability for each source
  for (j in 1:nSources) {
      p[j] ~ dbeta(((0.65^2 - 0.65^3 - (0.65 * 0.1^2))/0.1^2), ((0.65 - (2*0.65^2) + 0.65^3 - 0.1^2 + (0.65 * 0.1^2))/0.1^2)) # moment matching for mean 0.65 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- alpha[i] + beta[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t]
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
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_adult.new[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    }
  }
  
  # Get means, CVs of data and predicted data
  mean_p1_adult <- mean(p1_adult[,1:nYears])
  CV_p1_adult <- sd(p1_adult[,1:nYears])/mean(p1_adult[,1:nYears])
  mean_p1_adult.new <- mean(p1_adult.new)
  CV_p1_adult.new <- sd(p1_adult.new)/mean(p1_adult.new)
  
  # Calculate p values
  pval.mean_p1 <- step(mean_p1_adult.new - mean_p1_adult)
  pval.CV_p1 <- step(CV_p1_adult.new - CV_p1_adult)
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
                  p1_adult = p1_adult,
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_adult$Source)


# Parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2)))
  }
}

# Set initial values
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 50000
nc <- 3
nb <- 10000
nt <- 1

set.seed(1234)

# Fit Model
Adult_climateEffects <- jagsUI::jags(data = jags_data,
                                   parameters.to.save = jags_params,
                                   model.file = "Analysis/nMix_JAGS_files/Adult_climateEffects.jags",
                                   n.chains = nc,
                                   n.iter = ni,
                                   n.burnin = nb,
                                   n.thin = nt,
                                   parallel = T,
                                   inits = init_vals)

Adult_climateEffects_params <- MCMCsummary(Adult_climateEffects, HPD = T)

# What was the range of covariate effects on adult abundance?
Adult_climateEffects %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta.cov\\[1")) %>% 
  .[,2] %>% 
  range()

## Adult model with just random effects
sink("Analysis/nMix_JAGS_files/Adult_randomEffects.jags")
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
      p[j] ~ dbeta(((0.65^2 - 0.65^3 - (0.65 * 0.1^2))/0.1^2), ((0.65 - (2*0.65^2) + 0.65^3 - 0.1^2 + (0.65 * 0.1^2))/0.1^2)) # moment matching for mean 0.65 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
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
    
  ## ICC
  for (i in 1:nReps) {
    ICC.adult[i] <- s2.eps/(s2.eps + s2.gam[i])
  }
  mean_ICC.adult <- mean(ICC.adult)
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_adult.new[i,t] ~ dbin(p[Sources[i]], N.adult[i,t])
    }
  }
  
  # Get means, CVs of data and predicted data
  mean_p1_adult <- mean(p1_adult[,1:nYears])
  CV_p1_adult <- sd(p1_adult[,1:nYears])/mean(p1_adult[,1:nYears])
  mean_p1_adult.new <- mean(p1_adult.new)
  CV_p1_adult.new <- sd(p1_adult.new)/mean(p1_adult.new)
  
  # Calculate p values
  pval.mean_p1 <- step(mean_p1_adult.new - mean_p1_adult)
  pval.CV_p1 <- step(CV_p1_adult.new - CV_p1_adult)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nReps = nReps, 
                  nYears = nYears,
                  nSources = nSources,
                  Area = sample_areas_wide,
                  p1_adult = p1_adult,
                  p2_adult = p2_adult, 
                  p3_adult = p3_adult,
                  Sources = p1_adult$Source)


# Parameters to save
jags_params <- c("alpha", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps, nYears))
for (i in 1:nReps) {
  for (t in 1:nYears) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na((p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t])),
                                                rpois(1, lambda = 200),
                                                (p1_adult[i,t] + p2_adult[i,t] + p3_adult[i,t] + 1) * 2)))
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


# MCMC settings
ni <- 100000
nc <- 3
nb <- 25000
nt <- 1

set.seed(1234)

# Fit Model
Adult_randomEffects <- jagsUI::jags(data = jags_data,
                                  parameters.to.save = jags_params,
                                  model.file = "Analysis/nMix_JAGS_files/Adult_randomEffects.jags",
                                  n.chains = nc,
                                  n.iter = ni,
                                  n.burnin = nb,
                                  n.thin = nt,
                                  parallel = T,
                                  inits = init_vals)

Adult_randomEffects_params <- MCMCsummary(Adult_randomEffects, HPD = T)

#MCMCtrace(Adult_randomEffects_params, params = , pdf = F)


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

sample_areas_wide_N <- sample_areas_wide %>% 
  filter(COMID %in% N_Sites$COMID)
sample_areas_wide_S <- sample_areas_wide %>% 
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

# NORTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_N, 
                  nYears = nYears,
                  nSources = nSources_N,
                  Area = sample_areas_wide_N,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_YOY = p1_YOY_N,
                  p2_YOY = p2_YOY_N, 
                  p3_YOY = p3_YOY_N,
                  Sources = as.numeric(as.factor(p1_YOY_N$Source)))


# Parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

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
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 25000
nc <- 3
nb <- 10000
nt <- 1

# Fit Model
YOY_climateEffects_N <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "Analysis/nMix_JAGS_files/YOY_climateEffects.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

YOY_climateEffects_N_params <- MCMCsummary(YOY_climateEffects_N, HPD = T)

# SOUTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_S, 
                  nYears = nYears,
                  nSources = nSources_S,
                  Area = sample_areas_wide_S,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_YOY = p1_YOY_S,
                  p2_YOY = p2_YOY_S, 
                  p3_YOY = p3_YOY_S,
                  Sources = as.numeric(as.factor(p1_YOY_S$Source)))


# Parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

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
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 100000
nc <- 3
nb <- 50000
nt <- 1

# Fit Model
YOY_climateEffects_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/YOY_climateEffects.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

YOY_climateEffects_S_params <- MCMCsummary(YOY_climateEffects_S, HPD = T)

## Adult analysis
# Models can stay the same, we just subset the data

# NORTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_N, 
                  nYears = nYears,
                  nSources = nSources_N,
                  Area = sample_areas_wide_N,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_adult = p1_adult_N,
                  p2_adult = p2_adult_N, 
                  p3_adult = p3_adult_N,
                  Sources = as.numeric(as.factor(p1_adult_N$Source)))


# Parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

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
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)

# MCMC settings
ni <- 25000
nc <- 3
nb <- 10000
nt <- 1

# Fit Model
Adult_climateEffects_N <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_climateEffects.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

Adult_climateEffects_N_params <- MCMCsummary(Adult_climateEffects_N,  HPD = T)

# SOUTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_S, 
                  nYears = nYears,
                  nSources = nSources_S,
                  Area = sample_areas_wide_S,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_adult = p1_adult_S,
                  p2_adult = p2_adult_S, 
                  p3_adult = p3_adult_S,
                  Sources = as.numeric(as.factor(p1_adult_S$Source)))


# Parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

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
init_vals <- function() list(alpha = rnorm(nReps, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nSources),
                             N.adult = N.adult.inits)
# MCMC settings
ni <- 25000
nc <- 3
nb <- 10000
nt <- 1

# Fit Model
Adult_climateEffects_S <- jagsUI::jags(data = jags_data,
                                        parameters.to.save = jags_params,
                                        model.file = "Analysis/nMix_JAGS_files/Adult_climateEffects.jags",
                                        n.chains = nc,
                                        n.iter = ni,
                                        n.burnin = nb,
                                        n.thin = nt,
                                        parallel = T,
                                        inits = init_vals)

Adult_climateEffects_S_params <- MCMCsummary(Adult_climateEffects_S,  HPD = T)


######################################
### Post-Hoc
######################################

# make a dataframe of segments with coordinates for mapping
COMIDs_source <- SE_Site_Final %>% 
  filter(COMID %in% p1_YOY$COMID) %>% 
  group_by(COMID,
           Source) %>% 
  summarize(State = first(State),
            Lat = first(Lat),
            Long = first(Long),
            Elev_m = first(Elev_m))

segment_data <- data.frame(Source = sources) %>% 
  mutate(sourceNum = row_number()) %>% 
  right_join(p1_YOY[,c("COMID", "Source")], by = c("sourceNum" = "Source")) %>% 
  left_join(source_agency_crosswalk) %>% 
  left_join(COMIDs_source, by = c("COMID", "Source")) %>% 
  dplyr::select(-sourceNum)

# Export to make the map in QGIS
#fwrite(segment_data, here("Data", "SE_Trout_Segments_filtered.csv"))

# Create a map of segments
US_states <- map_data("state")

segments_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = segment_data, 
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

# make north and south versions too
segment_data_N <- segment_data %>% 
  filter(COMID %in% N_Sites$COMID)

segment_data_S <- segment_data %>% 
  filter(COMID %in% S_Sites$COMID)

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

#####################################
# Create a table of segment summaries

# Get NHDPlus data for the segments considered for the analysis
SE_segments_NHDPlus <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")
SE_segments_NHDPlus <- SE_segments_NHDPlus %>%
  filter(COMID %in% sample_areas$COMID) %>% 
  dplyr::select(COMID,
                SLOPE,
                LENGTHKM,
                AreaSqKM,
                MINELEVSMO,
                MAXELEVSMO,
                AreaSqKM,
                StreamOrde)

# filter for N and S subregions
SE_segments_NHDPlus_N <- SE_segments_NHDPlus %>% 
  filter(COMID %in% N_Sites$COMID)

SE_segments_NHDPlus_S <- SE_segments_NHDPlus %>% 
  filter(COMID %in% S_Sites$COMID)

# Get widths of sites within the segments considered for the analysis
SE_segments <- SE_Site_Final %>% 
  filter(COMID %in% sample_areas$COMID) %>% 
  select(COMID,
         Width_m,
         Elev_m)

# filter for N and S subregions
SE_segments_N <- SE_segments %>% 
  filter(COMID %in% N_Sites$COMID)

SE_segments_S <- SE_segments %>% 
  filter(COMID %in% S_Sites$COMID)

Segment_Summary.table <- data.frame(Variable = c("Channel Slope (%)",
                                                 "Length (km)",
                                                 "Catchment area (km^2)",
                                                 "Elevation (m)",
                                                 "Stream order",
                                                 "Wetted width (m)"),
                                    North_Mean = c(mean(SE_segments_NHDPlus_N$SLOPE, na.rm = T)*100,
                                             mean(SE_segments_NHDPlus_N$LENGTHKM, na.rm = T),
                                             mean(SE_segments_NHDPlus_N$AreaSqKM, na.rm = T),
                                             mean(rowMeans(SE_segments_NHDPlus_N[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                             mean(SE_segments_NHDPlus_N$StreamOrde, na.rm = T),
                                             mean(SE_segments_N$Width_m, na.rm = T)),
                                    North_SD = c(sd(SE_segments_NHDPlus_N$SLOPE, na.rm = T)*100,
                                                sd(SE_segments_NHDPlus_N$LENGTHKM, na.rm = T),
                                                sd(SE_segments_NHDPlus_N$AreaSqKM, na.rm = T),
                                                sd(rowMeans(SE_segments_NHDPlus_N[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                sd(SE_segments_NHDPlus_N$StreamOrde, na.rm = T),
                                                sd(SE_segments_N$Width_m, na.rm = T)),
                                    South_Mean = c(mean(SE_segments_NHDPlus_S$SLOPE, na.rm = T)*100,
                                                   mean(SE_segments_NHDPlus_S$LENGTHKM, na.rm = T),
                                                   mean(SE_segments_NHDPlus_S$AreaSqKM, na.rm = T),
                                                   mean(rowMeans(SE_segments_NHDPlus_S[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                   mean(SE_segments_NHDPlus_S$StreamOrde, na.rm = T),
                                                   mean(SE_segments_S$Width_m, na.rm = T)),
                                    South_SD = c(sd(SE_segments_NHDPlus_S$SLOPE, na.rm = T)*100,
                                                 sd(SE_segments_NHDPlus_S$LENGTHKM, na.rm = T),
                                                 sd(SE_segments_NHDPlus_S$AreaSqKM, na.rm = T),
                                                 sd(rowMeans(SE_segments_NHDPlus_S[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                 sd(SE_segments_NHDPlus_S$StreamOrde, na.rm = T),
                                                 sd(SE_segments_S$Width_m, na.rm = T)))

#####################################
# Create a table of covariate summaries

# Filter for temp data at sites and years where we have trout data
temp_summary_data <- SE_COMID_temp_covars %>% 
  filter(COMID %in% p1_YOY$COMID,
         Year %in% 1981:2015)

flow_summary_data <- SE_COMID_flow_covars %>% 
  filter(COMID %in% p1_YOY$COMID,
         Year %in% 1982:2015)

nMix_Covar_Summary.table <- data.frame(Covar = c("Mean 90th percentile summer temperature (C)",
                                       "Maximum 90th percentile winter streamflow (ft^3/s)", 
                                       "Maximum 90th percentile spring streamflow (ft^3/s)"),
                                       Mean = c(mean(temp_summary_data$Mean_Max_Summer_Temp),
                                                mean(flow_summary_data$Max_0.9Q_WinterFlow),
                                                mean(flow_summary_data$Max_0.9Q_SpringFlow)),
                                       sd = c(sd(temp_summary_data$Mean_Max_Summer_Temp),
                                              sd(flow_summary_data$Max_0.9Q_WinterFlow),
                                              sd(flow_summary_data$Max_0.9Q_SpringFlow)))

###################
# Stock-recruitment curves
stock_recruit_data <- SE_Ind_Final %>%
  group_by(SampleID) %>% 
  summarise(Year = first(year(Date)),
            SiteID = first(SiteID),
            n_YOY = sum(TL_mm <= 90),
            n_Adults = sum(TL_mm > 90)) %>% 
  ungroup() %>% 
  left_join(SE_Site_Final[,c("SiteID", "COMID")]) %>%
  filter(!is.na(COMID),
         COMID %in% p1_YOY$COMID) %>% 
  group_by(COMID,
          Year) %>% 
  summarise(log_YOY_t = log(mean(n_YOY)+1),
            log_Adults_t = log(mean(n_Adults))+1) %>% 
  mutate(log_YOY_tp1 = lead(log_YOY_t),
         log_Adults_tp1 = lead(log_Adults_t))

SSR_1.plot <- stock_recruit_data %>% 
  ggplot(aes(x = log_YOY_t,
             y = log_Adults_tp1)) +
  #geom_point(size = 0.5) +
  geom_smooth(aes(group = as.factor(COMID)), 
              method = "lm", 
              color = "black",
              size = 0.33,
              se = FALSE) +
  geom_smooth(method = "lm", 
              color = "blue",
              se = FALSE) +
  theme_classic() +
  labs(x = "Log (mean YOY abundance) in year t",
       y = "Log (mean Adult abundance) in year t+1")

SSR_2.plot <- stock_recruit_data %>% 
  ggplot(aes(x = log_Adults_t,
             y = log_YOY_tp1)) +
  #geom_point(size = 0.5) +
  geom_smooth(aes(group = as.factor(COMID)), 
              method = "lm", 
              color = "black",
              size = 0.33,
              se = FALSE) +
  geom_smooth(method = "lm", 
              color = "blue",
              se = FALSE) +
  theme_classic() +
  labs(x = "Log (mean Adult abundance) in year t",
       y = "Log (mean YOY abundance) in year t+1")

###################
# PPC p-values

YOY_full_PPCs <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(mean)

Adult_full_PPCs <- Adult_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(mean)

PPC_pvals.table <- data.frame(Life_Stage = c(rep("YOY", 2),
                                             rep("Adult", 2)),
                              stat = rep(c("mean", "CV"), 2),
                              pVal = rbind(YOY_full_PPCs,
                                           Adult_full_PPCs))

###################
## ICC Values
# Join site data to ICC values
YOY_ICCs.table <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data)

Adult_ICCs.table <- Adult_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data)

YOY_ICCs_N.table <- YOY_BKT_nMix_full_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data_N)

YOY_ICCs_S.table <- YOY_BKT_nMix_full_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data_S)

Adult_ICCs_N.table <- Adult_BKT_nMix_full_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data_N)

Adult_ICCs_S.table <- Adult_BKT_nMix_full_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data_S)

# Do the most synchronous/asynchronous sites from the N/S show up in the the overall ICCs?
YOY_highestICCs_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = head(arrange(YOY_ICCs_N.table, desc(mean))), 
             aes(x = Long, y = Lat), shape = 5) + # diamonds
  geom_point(data = head(arrange(YOY_ICCs_S.table, desc(mean))), 
              aes(x = Long, y = Lat), shape = 1) + # circles
  geom_point(data = head(arrange(YOY_ICCs.table, desc(mean))), 
              aes(x = Long, y = Lat), shape = 3) + # crosses
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
  geom_point(data = head(arrange(YOY_ICCs_N.table, mean)),
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(YOY_ICCs_S.table, mean)),
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(YOY_ICCs.table, mean)), 
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
  geom_point(data = head(arrange(Adult_ICCs_N.table, desc(mean))), 
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(Adult_ICCs_S.table, desc(mean))), 
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(Adult_ICCs.table, desc(mean))), 
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
  geom_point(data = head(arrange(Adult_ICCs_N.table, mean)), 
             aes(x = Long, y = Lat), shape = 5) +
  geom_point(data = head(arrange(Adult_ICCs_S.table, mean)), 
             aes(x = Long, y = Lat), shape = 1) +
  geom_point(data = head(arrange(Adult_ICCs.table, mean)), 
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
  geom_point(data = YOY_ICCs.table, 
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

# Adult
Adult_ICC_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = Adult_ICCs.table, 
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

# Is ICC correlated with any site-level variables?
library(corrplot)

# load NHDplus data
NHDplus_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")
# and filter to our sites
NHDplus_data <- NHDplus_data %>% 
  filter(COMID %in% segment_data$COMID)

# load streamcat data
StreamCat_Data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/StreamCat Data/StreamCat_Covars_Combined.csv")
# and filter to our sites
StreamCat_Data <- StreamCat_Data %>% 
  filter(COMID %in% segment_data$COMID)

ICC_corr_data <- YOY_ICCs.table[,c(2,9,12:14)] %>% 
  left_join(SE_Site_Final, by = "COMID") %>% 
  left_join(NHDplus_data) %>% 
  left_join(StreamCat_Data)

YOY_ICC_corrPlot <- corrplot(cor(select_if(ICC_corr_data, is.numeric) , method="spearman", use="pairwise.complete.obs"))
# get values
YOY_ICC_corrs.table <- as.data.frame(YOY_corrPlot$corrPos) %>% 
  filter(xName == "mean")
# mean ICC is not really correlated with any of these site-level covars


################
# Covariate effects/betas
# Summarize covariate effects in table
# plot data
mu.beta_samples.table <- rbind(data.frame(sample_val = rbind(as.matrix(YOY_BKT_nMix_full$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(YOY_BKT_nMix_full$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(YOY_BKT_nMix_full$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_BKT_nMix_full$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("N+S")),
                               data.frame(sample_val = rbind(as.matrix(Adult_BKT_nMix_full$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(Adult_BKT_nMix_full$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(Adult_BKT_nMix_full$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_BKT_nMix_full$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("N+S")),
                               data.frame(sample_val = rbind(as.matrix(YOY_BKT_nMix_full_N$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(YOY_BKT_nMix_full_N$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(YOY_BKT_nMix_full_N$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_BKT_nMix_full_N$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("N")),
                               data.frame(sample_val = rbind(as.matrix(Adult_BKT_nMix_full_N$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(Adult_BKT_nMix_full_N$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(Adult_BKT_nMix_full_N$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_BKT_nMix_full_N$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("N")),
                               data.frame(sample_val = rbind(as.matrix(YOY_BKT_nMix_full_S$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(YOY_BKT_nMix_full_S$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(YOY_BKT_nMix_full_S$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_BKT_nMix_full_S$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("S")),
                               data.frame(sample_val = rbind(as.matrix(Adult_BKT_nMix_full_S$sims.list$mu.beta.cov[,1]),
                                                             as.matrix(Adult_BKT_nMix_full_S$sims.list$mu.beta.cov[,2]),
                                                             as.matrix(Adult_BKT_nMix_full_S$sims.list$mu.beta.cov[,3])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_BKT_nMix_full_S$sims.list$mu.beta.cov[,1])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("S")))

# Reorder the subregions so the N+S region plots first
mu.beta_samples.table$subregion <- factor(mu.beta_samples.table$subregion, c("N+S", "N", "S"))
# Reorder the covariates so that summer temperature plots first
mu.beta_samples.table$covar <- factor(mu.beta_samples.table$covar, c("Summer Air Temperature", "Winter Flow", "Spring Flow"))

# plot
cov_effects.plot <- ggplot(mu.beta_samples.table) +
  geom_violin(aes(x = sample_val,
                  y = fct_rev(covar),
                  fill = fct_rev(subregion)),
              trim = F,
              alpha = 0.75,
              color = NA) +
  facet_grid(life_stage ~ .) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  theme_classic() +
  labs(x = "Value",
       y = element_blank(),
       fill = "Subregion") +
  scale_fill_brewer(palette = "Dark2",
                    limits = c("N+S", "N", "S")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 11))

# map betas (segment-specific covariate effects) in space
YOY_climate_effects.table <- rbind(YOY_BKT_nMix_full_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta.cov\\[1,")) %>% 
                                     select(mean),
                                   YOY_BKT_nMix_full_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta.cov\\[2,")) %>% 
                                     select(mean),
                                   YOY_BKT_nMix_full_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta.cov\\[3,")) %>% 
                                     select(mean)) %>% 
  cbind(segment_data,
        rbind(data.frame(covar = rep("Summer Temperature", 159)),
              data.frame(covar = rep("Winter Flow", 159)),
              data.frame(covar = rep("Spring Flow", 159))))

# reorder the covariates so that summer temperature plots first
YOY_climate_effects.table$covar <- factor(YOY_climate_effects.table$covar, c("Summer Temperature", "Winter Flow", "Spring Flow"))

library(latex2exp)

# YOY_climate_effects_map.plot <- ggplot() +
#   geom_polygon(data = US_states, 
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = YOY_climate_effects.table, 
#              aes(x = Long, y = Lat, color = mean), 
#              alpha = 0.5) +
#   coord_map("bonne",
#             lat0 = 40,
#             xlim = c(-85, -76),
#             ylim = c(34.5, 40)) +
#   labs(x = "Long",
#        y = "Lat",
#        color = TeX(r'($\beta$ Value)')) +
#   scale_color_viridis_c() +
#   theme_classic() +
#   facet_grid(. ~ covar)



YOY_SummTemp_betas.plot <- ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Summer Temperature",],
             aes(x = Long, y = Lat, color = mean),
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)'),
       title = "a)") +
  scale_color_viridis_c() +
  theme_classic()

YOY_WintFlow_betas.plot <- ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Winter Flow",],
             aes(x = Long, y = Lat, color = mean),
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)'),
       title = "b)") +
  scale_color_viridis_c() +
  theme_classic()

YOY_SprFlow_betas.plot <- ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Spring Flow",],
             aes(x = Long, y = Lat, color = mean),
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)'),
       title = "c)") +
  scale_color_viridis_c() +
  theme_classic()

YOY_climate_effects_map.plot <- grid.arrange(YOY_SummTemp_betas.plot,
                                             YOY_WintFlow_betas.plot,
                                             YOY_SprFlow_betas.plot)

####
# plots of first-pass trout count versus each climate covariate where each segment = each line
count_vs_climate_data <- pivot_longer(Mean_Max_Summer_Temp_Scaled, cols = 1:34, names_to = "Year", values_to = "SummTemp") %>% 
  left_join(pivot_longer(Max_0.9Q_WinterFlow_Scaled, cols = 1:34, names_to = "Year", values_to = "WintFlow")) %>%
  left_join(pivot_longer(Max_0.9Q_SpringFlow_Scaled, cols = 1:34, names_to = "Year", values_to = "SprFlow")) %>% 
  left_join(pivot_longer(p1_YOY, cols = 1:34, names_to = "Year", values_to = "count_YOY"))

count_vs_SummTemp.plot <- ggplot(count_vs_climate_data) +
  geom_line(aes(x = SummTemp,
                y = log(count_YOY),
                group = COMID)) +
  theme_classic()

count_vs_WintFlow.plot <- ggplot(count_vs_climate_data) +
  geom_line(aes(x = WintFlow,
                y = log(count_YOY),
                group = COMID)) +
  #lims(x = c(-1,0)) +
  theme_classic()

count_vs_SprFlow.plot <- ggplot(count_vs_climate_data) +
  geom_line(aes(x = SprFlow,
                y = count_YOY,
                group = COMID)) +
  theme_classic()


####
# are covariate effects (betas) correlated geographically?
chart.Correlation(pivot_wider(YOY_climate_effects.table, names_from = covar, values_from = mean)[,5:10], method = "spearman")

beta_corr_data <- YOY_climate_effects.table %>% 
  pivot_wider(names_from = covar, values_from = mean) %>% 
  left_join(NHDplus_data) %>% 
  left_join(StreamCat_Data)

YOY_beta_corrPlot <- corrplot(cor(select_if(beta_corr_data, is.numeric) , method="spearman", use="pairwise.complete.obs"))
# get values
YOY_beta_corrs.table <- as.data.frame(YOY_beta_corrPlot$corrPos) %>% 
  filter(xName %in% c("Summer Temperature", "Winter Flow", "Spring Flow"))
# betas are not really correlated with any of these site-level covars

####
# What proportion of beta.cov values are significant at 95% HDPIs for YOY?
# Create a table to store these values
pct_effects_negative.table <- data.frame(Covar = c("Summ Temp", "Wint Flow", "Spr Flow"),
                                   pct_sgn_neg = rep(NA, times = 3))

# Summer temp
YOY_summTemp_betas <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov\\[1"))

pct_effects_negative.table[1,2] <- sum(YOY_summTemp_betas[,"95%_HPDU"] < 0)/nrow(YOY_summTemp_betas) * 100

# Winter Flow
YOY_wintFlow_betas <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov\\[2"))

pct_effects_negative.table[2,2] <- sum(YOY_wintFlow_betas[,"95%_HPDU"] < 0)/nrow(YOY_wintFlow_betas) * 100

# Spring Flow
YOY_sprFlow_betas <- YOY_BKT_nMix_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov\\[3"))

pct_effects_negative.table[3,2] <- sum(YOY_sprFlow_betas[,"95%_HPDU"] < 0)/nrow(YOY_sprFlow_betas) * 100

####
# What was the range of covariate effects on YOY abundance?
cov_effect_ranges.table <- data.frame(Covar = c("Summ Temp", "Wint Flow", "Spr Flow"),
                                      min_mean = rbind(min(YOY_summTemp_betas$mean),
                                                  min(YOY_wintFlow_betas$mean),
                                                  min(YOY_sprFlow_betas$mean)),
                                      min_95HDPL = rbind(min(YOY_summTemp_betas$`95%_HPDL`),
                                                         min(YOY_wintFlow_betas$`95%_HPDL`),
                                                         min(YOY_sprFlow_betas$`95%_HPDL`)),
                                      min_95HDPU = rbind(min(YOY_summTemp_betas$`95%_HPDU`),
                                                         min(YOY_wintFlow_betas$`95%_HPDU`),
                                                         min(YOY_sprFlow_betas$`95%_HPDU`)),
                                      max_mean = rbind(max(YOY_summTemp_betas$mean),
                                                  max(YOY_wintFlow_betas$mean),
                                                  max(YOY_sprFlow_betas$mean)),
                                      max_95HDPL = rbind(max(YOY_summTemp_betas$`95%_HPDL`),
                                                         max(YOY_wintFlow_betas$`95%_HPDL`),
                                                         max(YOY_sprFlow_betas$`95%_HPDL`)),
                                      max_95HDPU = rbind(max(YOY_summTemp_betas$`95%_HPDU`),
                                                         max(YOY_wintFlow_betas$`95%_HPDU`),
                                                         max(YOY_sprFlow_betas$`95%_HPDU`)))

################
# Summarize detection probability in table

Detect_probs.table <- data.frame(
  Agency = rep(sources2$Agency, times = 2),
  Life_Stage = c(rep("YOY", times = 8),
                rep("Adult", times = 8))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(YOY_BKT_nMix_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),1:4],
              Adult_BKT_nMix_full_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]"),1:4]))

# and make a plot to visualize
Detect_probs.plot <- ggplot(data = Detect_probs.table) +
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
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  theme_classic() + 
  theme(text = element_text(family =  "serif"),
        axis.text.x = element_text(angle = -45, hjust=0))


########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here("results", "v1.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
