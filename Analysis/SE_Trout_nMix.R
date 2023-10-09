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
library(viridis)      # Colorblind-friendly plotting palettes


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
source_agency_crosswalk <- data.frame(Agency = unique(SE_Site_Final$Agency),
                                      Agency = c("GA Dept. Natural Resources",
                                                 "Great Smoky Mountain Nat'l Park",
                                                 "ME Dept. Inland Fisheries & Wildlife",
                                                 "MD Dept. Natural Resources",
                                                 "NC Wildlife Resources Commission",
                                                 "Shenandoah Nat'l Park",
                                                 "GA Dept. Natural Resources",
                                                 "MD Dept. Natural Resources",
                                                 "SC Dept. Natural Resources",
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
  left_join(source_agency_crosswalk) %>% 
  left_join(SE_Site_Final[,c("SiteID", "COMID", "Lat", "Long", "Length_m", "Width_m")], by = "SiteID") %>% 
  filter(SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         Lat <= 39.716667, # filter for just sites south of the Mason-Dixon Line
         COMID %in% SE_COMID_flow_covars$COMID, # filter for COMIDs for which we have flow data - there are several sites in Maine with no flow
         COMID %in% SE_COMID_temp_covars$COMID) %>% # same goes for temperature - there's one site in Maine with no temp
  group_by(COMID,
           Year = year(Date),
           Agency) %>% 
  summarise(Area_Sampled = sum(Length_m * Width_m)) %>%  # Calculates the sum of site areas within the stream segment that were sampled that year by that agency
  filter(Year >= 1981, # Filter to one year after the earliest year that we have temperature data
         Year <= 2015,    # Filter to the latest year that we have flow data
         !is.na(Area_Sampled)) %>% 
  ungroup()

# Make a dataframe of the COMIDs and years when sampling happened
# SampleYears <- SE_Sample_Final %>% 
#   filter(SiteID %in% BKT_Sites$SiteID) %>%  # Filter to sites with records of BKT
#   left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
#   mutate(Year = year(Date)) %>% 
#   dplyr::select(COMID, Year, Agency) %>% 
#   unique() %>% 
#   filter(COMID %in% sample_areas$COMID, # Filter because we can only use data from COMIDs with area and coordinates
#          !is.na(COMID),
#          Year >= 1981, # Filter to one year after the earliest year that we have temperature data
#          Year <= 2015)  # Filter to the latest year that we have flow data

# What was the percentage of single- vs multipass electrofishing?
passes_pcts.table <- SE_Sample_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  left_join(source_agency_crosswalk) %>% 
  mutate(Year = year(Date)) %>% # create a year column
  inner_join(sample_areas) %>%  # use an inner join to filter for only segment-agency-year combos that have areas
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
  left_join(source_agency_crosswalk) %>% # join in agency name
  mutate(Year = year(Date)) %>% 
  right_join(sample_areas) %>%  # use a right join to filter for only segment-agency-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
  left_join(SE_Sample_Final[,c("SampleID", "NumPasses")]) %>% # join in the number of passes taken at each sample
  filter(!is.na(NumPasses)) %>% 
  group_by(SampleID) %>% 
  summarize(COMID = first(COMID),
            Year = first(Year),
            Agency = first(Agency),
            NumPasses = first(NumPasses),
            P1_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 1), # Filter here for YOY BKT
            P2_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 2),
            P3_Count_YOY = sum(SPP == "BKT" & TL_mm <= 90 & PassNo == 3)) %>% 
  ungroup()

# where only one pass was conducted, set counts for p2 and p3 to NA
YOY_BKT_passCounts[which(YOY_BKT_passCounts$NumPasses == 1) ,c("P2_Count_YOY", "P3_Count_YOY")] <- NA

# Now summarize by year, agency, and COMID
YOY_BKT_passCounts <- YOY_BKT_passCounts %>% 
  group_by(COMID,
           Year,
           Agency) %>% 
  summarise(P1_Count_YOY = sum(P1_Count_YOY), 
            P2_Count_YOY = sum(P2_Count_YOY),
            P3_Count_YOY = sum(P3_Count_YOY)) %>% 
  ungroup()
  #.[!duplicated(.[,c(1:2,4:6)]),] # remove any duplicate rows not already cleaned from the data



      # what's the proportion of segments that have multiple sites in them?
      sites_in_segments <- SE_Sample_Final %>% 
        left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
        left_join(source_agency_crosswalk) %>% 
        mutate(Year = year(Date)) %>% 
        right_join(sample_areas) %>%  # use a right join to filter for only segment-agency-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
        group_by(COMID,
                 Year,
                 Agency) %>% 
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

# Make separate dataframes of count data for climate models and for random effects models
YOY_BKT_passCounts_CE <- YOY_BKT_passCounts
YOY_BKT_passCounts_RE <- YOY_BKT_passCounts %>% filter(!is.na(P2_Count_YOY))

# YOY abundance data for model
# Now make a separate, wide data frame for each pass
p1_YOY_CE <- YOY_BKT_passCounts_CE %>% 
  dplyr::select(-P2_Count_YOY,
                -P3_Count_YOY) %>%
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_YOY) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

# save a vector of the names of the agencies
agencies_CE <- sort(unique(p1_YOY_CE$Agency))

p1_YOY_RE <- YOY_BKT_passCounts_RE %>% 
  dplyr::select(-P2_Count_YOY,
                -P3_Count_YOY) %>%
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_YOY) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

# save a vector of the names of the agencies
agencies_RE <- sort(unique(p1_YOY_RE$Agency))

p2_YOY_CE <- YOY_BKT_passCounts_CE %>% 
  dplyr::select(-P1_Count_YOY,
                -P3_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_YOY) %>% 
  arrange(COMID) %>% 
  right_join(p1_YOY_CE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p2_YOY_RE <- YOY_BKT_passCounts_RE %>% 
  dplyr::select(-P1_Count_YOY,
                -P3_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_YOY) %>% 
  arrange(COMID) %>% 
  right_join(p1_YOY_RE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p3_YOY_CE <- YOY_BKT_passCounts_CE %>% 
  dplyr::select(-P1_Count_YOY,
                -P2_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_YOY) %>% 
  arrange(COMID) %>% 
  right_join(p1_YOY_CE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p3_YOY_RE <- YOY_BKT_passCounts_RE %>% 
  dplyr::select(-P1_Count_YOY,
                -P2_Count_YOY) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_YOY) %>% 
  arrange(COMID) %>% 
  right_join(p1_YOY_RE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

## Adults
# Tally adult counts by pass at each COMID, year combination
Adult_BKT_passCounts <- SE_Ind_Final %>% 
  left_join(SE_Site_Final[,c("SiteID", "COMID")]) %>% # Join in COMIDs
  left_join(source_agency_crosswalk) %>% # join in agency name
  mutate(Year = year(Date)) %>% 
  right_join(sample_areas) %>%  # use a right join to filter for only segment-agency-year combos that have areas. This allows for the possibility that there were samples not accounted for in SE_Ind_Final because they were taken but had no fish
  left_join(SE_Sample_Final[,c("SampleID", "NumPasses")]) %>% # join in the number of passes taken at each sample
  filter(!is.na(NumPasses)) %>% 
  group_by(SampleID) %>% 
  summarize(COMID = first(COMID),
            Year = first(Year),
            Agency = first(Agency),
            NumPasses = first(NumPasses),
            P1_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 1), # Filter here for YOY BKT
            P2_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 2),
            P3_Count_adult = sum(SPP == "BKT" & TL_mm > 90 & PassNo == 3)) %>% 
  ungroup()

# where only one pass was conducted, set counts for p2 and p3 to NA
Adult_BKT_passCounts[which(Adult_BKT_passCounts$NumPasses == 1) ,c("P2_Count_adult", "P3_Count_adult")] <- NA

# Now summarize by year, agency, and COMID
Adult_BKT_passCounts <- Adult_BKT_passCounts %>% 
  group_by(COMID,
           Year,
           Agency) %>% 
  summarise(P1_Count_adult = sum(P1_Count_adult), 
            P2_Count_adult = sum(P2_Count_adult),
            P3_Count_adult = sum(P3_Count_adult)) %>% 
  ungroup()

# Matt Kulp sent a file that includes a list of sites where exotic salmonids were removed and BKT were stocked to restore the stream.
# Obviously we don't want
# see if those restored sites overlap with what we're about to run the synchrony analysis on
Adult_BKT_passCounts %>% 
  filter(COMID %in% GSMNP_Restored_Sites$COMID) %>% 
  distinct(COMID)
# 2 of the restored sites are in the dataset. Filter them out
Adult_BKT_passCounts<- Adult_BKT_passCounts %>% 
  filter(!COMID %in% GSMNP_Restored_Sites$COMID)

# Make separate dataframes of count data for climate models and for random effects models
Adult_BKT_passCounts_CE <- Adult_BKT_passCounts
Adult_BKT_passCounts_RE <- Adult_BKT_passCounts %>% filter(!is.na(P2_Count_adult))

# adult abundance data for model
# Now make a separate, wide data frame for each pass
p1_adult_CE <- Adult_BKT_passCounts_CE %>% 
  dplyr::select(-P2_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_adult) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

p1_adult_RE <- Adult_BKT_passCounts_RE %>% 
  dplyr::select(-P2_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P1_Count_adult) %>% 
  arrange(COMID) %>% 
  mutate(nYears_data = rowSums(!is.na(.[,3:ncol(.)]))) %>%  # count the number of years data at that site and pass
  filter(nYears_data >= 5)

p2_adult_CE <- Adult_BKT_passCounts_CE %>% 
  dplyr::select(-P1_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_adult) %>% 
  arrange(COMID) %>% 
  right_join(p1_adult_CE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p2_adult_RE <- Adult_BKT_passCounts_RE %>% 
  dplyr::select(-P1_Count_adult,
                -P3_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P2_Count_adult) %>% 
  arrange(COMID) %>% 
  right_join(p1_adult_RE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p3_adult_CE <- Adult_BKT_passCounts_CE %>% 
  dplyr::select(-P1_Count_adult,
                -P2_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_adult) %>% 
  arrange(COMID) %>% 
  right_join(p1_adult_CE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p3_adult_RE <- Adult_BKT_passCounts_RE %>% 
  dplyr::select(-P1_Count_adult,
                -P2_Count_adult) %>% 
  arrange(Year) %>% 
  pivot_wider(names_from = Year, 
              values_from = P3_Count_adult) %>% 
  arrange(COMID) %>% 
  right_join(p1_adult_RE[,c("COMID", "Agency")]) %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p1_adult_CE <- p1_adult_CE %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

p1_adult_RE <- p1_adult_RE %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% # change the Agency column to numeric 
  relocate(COMID, Agency, .after = last_col())

## COVARIATES
# Make wide dataframes of spatiotemporal covariates
# Temperature
Mean_Max_Summer_Temp_Scaled <- SE_COMID_temp_covars %>% 
  group_by(COMID) %>% 
  mutate(Mean_Max_Summer_Temp_Scaled = c(scale(Mean_Max_Summer_Temp))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Mean_Max_Summer_Temp_Scaled) %>% 
  filter(COMID %in% p1_YOY_CE$COMID, # Filter to just data for which we have samples
         Year %in% (YOY_BKT_passCounts$Year - 1)) %>% # filter for years which we have trout data, minus one year b/c we are using temp from the prior year
  pivot_wider(names_from = Year,
              values_from = Mean_Max_Summer_Temp_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Mean_Max_Summer_Temp_Scaled <- Mean_Max_Summer_Temp_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY_CE$COMID))

# Winter flow
Max_0.9Q_WinterFlow_Scaled <- SE_COMID_flow_covars %>%
  group_by(COMID) %>% 
  mutate(Max_0.9Q_WinterFlow_Scaled = c(scale(Max_0.9Q_WinterFlow))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Max_0.9Q_WinterFlow_Scaled) %>% 
  filter(COMID %in% p1_YOY_CE$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts$Year) %>% # filter for years which we have trout data
  pivot_wider(names_from = Year,
              values_from = Max_0.9Q_WinterFlow_Scaled) %>% 
  relocate(COMID, .after = last_col())

# use a right join to get duplicates for the COMIDs with multiple sources of data
Max_0.9Q_WinterFlow_Scaled <- Max_0.9Q_WinterFlow_Scaled %>% 
  right_join(data.frame(COMID = p1_YOY_CE$COMID))

# Spring flow
Max_0.9Q_SpringFlow_Scaled <- SE_COMID_flow_covars %>%
  group_by(COMID) %>% 
  mutate(Max_0.9Q_SpringFlow_Scaled = c(scale(Max_0.9Q_SpringFlow))) %>% # center and scale the covariate
  dplyr::select(COMID,
                Year,
                Max_0.9Q_SpringFlow_Scaled) %>% 
  filter(COMID %in% p1_YOY_CE$COMID, # Filter to just data for which we have samples
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
  right_join(data.frame(COMID = p1_YOY_CE$COMID))

# Filter the sample areas data for just those with fish observations in the time frame of interest
# use a right join to get duplicates for the COMIDs with multiple sources of data
# Then pivot wider to include it in the model
sample_areas_wide_CE <- sample_areas %>% 
  pivot_wider(names_from = Year, 
              values_from = Area_Sampled) %>% 
  arrange(COMID) %>% 
  relocate(COMID, .after = last_col()) %>%  # move the info columns to the end of the df to allow subsetting by year in the model
  right_join(p1_YOY_CE[,c("COMID", "Agency")]) %>% # use a right join to filter for areas of the years AND agencies where fish were caught
  dplyr::select(-Agency) %>% 
  .[,sort(colnames(.))] # arrange the columns (years)

sample_areas_wide_RE <- sample_areas %>% 
  pivot_wider(names_from = Year, 
              values_from = Area_Sampled) %>% 
  arrange(COMID) %>% 
  relocate(COMID, .after = last_col()) %>%  # move the info columns to the end of the df to allow subsetting by year in the model
  right_join(p1_YOY_RE[,c("COMID", "Agency")]) %>% # use a right join to filter for areas of the years AND sources where fish were caught
  dplyr::select(-Agency) %>% 
  .[,sort(colnames(.))] # arrange the columns (years)

# the sample areas cannot have missing values, 
# so we impute any segment-year that is missing a sample area as the average sample area for that segment.
# calculate the average sample area of each segment
mean_areas_CE <- (rowSums(sample_areas_wide_CE[,1:ncol(sample_areas_wide_CE)-1], na.rm = T)/apply(!is.na(sample_areas_wide_CE[,1:ncol(sample_areas_wide_CE)-1]), MARGIN = 1, sum))
mean_areas_RE <- (rowSums(sample_areas_wide_RE[,1:ncol(sample_areas_wide_RE)-1], na.rm = T)/apply(!is.na(sample_areas_wide_RE[,1:ncol(sample_areas_wide_RE)-1]), MARGIN = 1, sum))

# then replace the NAs with the average sample area at that segment
for (i in 1:nrow(sample_areas_wide_CE)) {
  sample_areas_wide_CE[i,][is.na(sample_areas_wide_CE[i,])] <- mean_areas_CE[i]
  print(i)
}
for (i in 1:nrow(sample_areas_wide_RE)) {
  sample_areas_wide_RE[i,][is.na(sample_areas_wide_RE[i,])] <- mean_areas_RE[i]
  print(i)
}

# and change the Agency column to numeric 
p1_YOY_CE <- p1_YOY_CE %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% 
  relocate(COMID, Agency, .after = last_col()) # move the info columns to the end of the df to allow subsetting by year in the model
p1_YOY_RE <- p1_YOY_RE %>% 
  mutate(Agency = as.numeric(as.factor(Agency))) %>% 
  relocate(COMID, Agency, .after = last_col()) # move the info columns to the end of the df to allow subsetting by year in the model


# Set sample sizes
nReps_CE <- nrow(p1_YOY_CE)
nReps_RE <- nrow(p1_YOY_RE)
nYears_CE <- ncol(p1_YOY_CE) - 3 # subtract 3 b/c the last three columns have info, not counts
nYears_RE <- ncol(p1_YOY_RE) - 3 # subtract 3 b/c the last three columns have info, not counts
nAgencies_CE <- length(unique(p1_YOY_CE$Agency))
nAgencies_RE <- length(unique(p1_YOY_RE$Agency))

########################################
# two models each for YOY and adults: one with just climate effects and one with just the two random effects

# YOY model with just env. covariates
sink("Analysis/nMix_JAGS_files/YOY_climateEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    omega[i] ~ dnorm(0, 0.001)
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
  # p.j - detection probability for each agency
  for (j in 1:nAgencies) {
      p[j] ~ dbeta(((0.5^2 - 0.5^3 - (0.5 * 0.1^2))/0.1^2), ((0.5 - (2*0.5^2) + 0.5^3 - 0.1^2 + (0.5 * 0.1^2))/0.1^2)) # moment matching for mean 0.5 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- omega[i] + beta[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_YOY[i,t] ~ dbin(p[Agencies[i]], N.YOY[i,t])
    # Pass 2
    p2_YOY[i,t] ~ dbin(p[Agencies[i]], (N.YOY[i,t] - p1_YOY[i,t]))
    # Pass 3
    p3_YOY[i,t] ~ dbin(p[Agencies[i]], (N.YOY[i,t] - p1_YOY[i,t] - p2_YOY[i,t]))
    }
  }
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_YOY.new[i,t] ~ dbin(p[Agencies[i]], N.YOY[i,t])
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
jags_data <- list(nReps = nReps_CE, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE,
                  Area = sample_areas_wide_CE,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_YOY = p1_YOY_CE,
                  p2_YOY = p2_YOY_CE, 
                  p3_YOY = p3_YOY_CE,
                  Agencies = p1_YOY_CE$Agency)


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1", "N.YOY")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_CE, nYears_CE))
for (i in 1:nReps_CE) {
  for (t in 1:nYears_CE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_CE[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_CE[i,t], p2_YOY_CE[i,t], p3_YOY_CE[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

YOY_climateEffects_params <- MCMCsummary(YOY_climateEffects, HPD = T, excl = "N.YOY")

# What was the range of covariate effects on YOY abundance?
YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta\\[1")) %>% 
  .[,2] %>% 
  range()

#MCMCtrace(YOY_climateEffects, params = "N.YOY", pdf = F)

# get abundance over time
N.YOY_climateEffects <- MCMCsummary(YOY_climateEffects,
                                    func = "median",
                                    HPD = T, params = "N.YOY")

## YOY model with just random effects
sink("Analysis/nMix_JAGS_files/YOY_randomEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    omega[i] ~ dnorm(0, 0.1)
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
  # p.j - detection probability for each agemcy
  for (j in 1:nAgencies) {
      p[j] ~ dbeta(((0.5^2 - 0.5^3 - (0.5 * 0.1^2))/0.1^2), ((0.5 - (2*0.5^2) + 0.5^3 - 0.1^2 + (0.5 * 0.1^2))/0.1^2)) # moment matching for mean 0.5 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.YOY[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- omega[i] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_YOY[i,t] ~ dbin(p[Agencies[i]], N.YOY[i,t])
    # Pass 2
    p2_YOY[i,t] ~ dbin(p[Agencies[i]], (N.YOY[i,t] - p1_YOY[i,t]))
    # Pass 3
    p3_YOY[i,t] ~ dbin(p[Agencies[i]], (N.YOY[i,t] - p1_YOY[i,t] - p2_YOY[i,t]))
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
      p1_YOY.new[i,t] ~ dbin(p[Agencies[i]], N.YOY[i,t])
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
jags_data <- list(nReps = nReps_RE, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE,
                  Area = sample_areas_wide_RE,
                  p1_YOY = p1_YOY_RE,
                  p2_YOY = p2_YOY_RE, 
                  p3_YOY = p3_YOY_RE,
                  Agencies = p1_YOY_RE$Agency)


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1", "N.YOY")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_RE, nYears_RE))
for (i in 1:nReps_RE) {
  for (t in 1:nYears_RE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_RE[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_RE[i,t], p2_YOY_RE[i,t], p3_YOY_RE[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE, 0, 0.1),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE, 0, 10),
                             gam = array(rnorm(nReps_RE * nYears_RE, 0, 10), dim = c(nReps_RE, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE),
                             N.YOY = N.YOY.inits)


# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

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

YOY_randomEffects_params <- MCMCsummary(YOY_randomEffects, HPD = T, excl = "N.YOY")

#MCMCtrace(YOY_randomEffects, params = "omega", pdf = F)

# plot abundance over time
#N.YOY_randomEffects <- MCMCsummary(YOY_randomEffects, HPD = T, params = "N.YOY")

### Adults
# Adult model with just env. covariates
sink("Analysis/nMix_JAGS_files/Adult_climateEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    omega[i] ~ dnorm(0, 0.001)
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
  # p.j - detection probability for each agency
  for (j in 1:nAgencies) {
      p[j] ~ dbeta(((0.65^2 - 0.65^3 - (0.65 * 0.1^2))/0.1^2), ((0.65 - (2*0.65^2) + 0.65^3 - 0.1^2 + (0.65 * 0.1^2))/0.1^2)) # moment matching for mean 0.65 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- omega[i] + beta[1,i] * Mean_Max_Summer_Temp_Scaled[i,t] + beta[2,i] * Max_0.9Q_WinterFlow_Scaled[i,t] + beta[3,i] * Max_0.9Q_SpringFlow_Scaled[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Agencies[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Agencies[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Agencies[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
    }
  }
  
  ### Posterior Predictive Check ###
  # Predict new data
  for (i in 1:nReps) {
    for (t in 1:nYears) {
      # Pass 1
      p1_adult.new[i,t] ~ dbin(p[Agencies[i]], N.adult[i,t])
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
jags_data <- list(nReps = nReps_CE, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE,
                  Area = sample_areas_wide_CE,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled,
                  p1_adult = p1_adult_CE,
                  p2_adult = p2_adult_CE, 
                  p3_adult = p3_adult_CE,
                  Agencies = p1_adult_CE$Agency)


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1", "N.adult")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_CE, nYears_CE))
for (i in 1:nReps_CE) {
  for (t in 1:nYears_CE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_CE[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_adult_CE[i,t], p2_adult_CE[i,t], p3_adult_CE[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

Adult_climateEffects_params <- MCMCsummary(Adult_climateEffects, HPD = T, excl = "N.adult")

# What was the range of covariate effects on adult abundance?
Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta\\[3")) %>% 
  .[,2] %>% 
  range()

N.adult_climateEffects <- MCMCsummary(Adult_climateEffects,
                                    func = "median",
                                    HPD = T, params = "N.adult")

## Adult model with just random effects
sink("Analysis/nMix_JAGS_files/Adult_randomEffects.jags")
cat("
model{
  
  ### Priors ###
  
  ## Site fixed effect
  for (i in 1:nReps){
    omega[i] ~ dnorm(0, 0.001)
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
  # p.j - detection probability for each agency
  for (j in 1:nAgencies) {
      p[j] ~ dbeta(((0.65^2 - 0.65^3 - (0.65 * 0.1^2))/0.1^2), ((0.65 - (2*0.65^2) + 0.65^3 - 0.1^2 + (0.65 * 0.1^2))/0.1^2)) # moment matching for mean 0.65 and variance 0.1
    }
  
  ## Process
  # Full model (all env. covars)
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    
      # Data
      N.adult[i,t] ~ dpois((Area[i,t] / 1000) * lambda[i,t])
      
      log(lambda[i,t]) <- omega[i] + eps[t] + gam[i,t]
    }
  }
  
  
  ## Observation
  for (i in 1:nReps) {
    for (t in 1:nYears) {
    # Pass 1
    p1_adult[i,t] ~ dbin(p[Agencies[i]], N.adult[i,t])
    # Pass 2
    p2_adult[i,t] ~ dbin(p[Agencies[i]], (N.adult[i,t] - p1_adult[i,t]))
    # Pass 3
    p3_adult[i,t] ~ dbin(p[Agencies[i]], (N.adult[i,t] - p1_adult[i,t] - p2_adult[i,t]))
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
      p1_adult.new[i,t] ~ dbin(p[Agencies[i]], N.adult[i,t])
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
jags_data <- list(nReps = nReps_RE, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE,
                  Area = sample_areas_wide_RE,
                  p1_adult = p1_adult_RE,
                  p2_adult = p2_adult_RE, 
                  p3_adult = p3_adult_RE,
                  Agencies = p1_adult_RE$Agency)


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.adult", "mean_ICC.adult", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_RE, nYears_RE))
for (i in 1:nReps_RE) {
  for (t in 1:nYears_RE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_RE[i,t]),
                                                  rpois(1, lambda = 200),
                                                  sum(p1_adult_RE[i,t], p2_adult_RE[i,t], p3_adult_RE[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE, 0, 10),
                             gam = array(rnorm(nReps_RE * nYears_RE, 0, 10), dim = c(nReps_RE, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE),
                             N.adult = N.adult.inits)


# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

set.seed(1234)
save.image()

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

MCMCtrace(Adult_randomEffects, params = "omega", pdf = F)


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
p1_YOY_CE_N <- p1_YOY_CE %>% 
  filter(COMID %in% N_Sites$COMID)
p2_YOY_CE_N <- p2_YOY_CE %>% 
  filter(COMID %in% N_Sites$COMID)
p3_YOY_CE_N <- p3_YOY_CE %>% 
  filter(COMID %in% N_Sites$COMID)

p1_YOY_RE_N <- p1_YOY_RE %>% 
  filter(COMID %in% N_Sites$COMID)
p2_YOY_RE_N <- p2_YOY_RE %>% 
  filter(COMID %in% N_Sites$COMID)
p3_YOY_RE_N <- p3_YOY_RE %>% 
  filter(COMID %in% N_Sites$COMID)

p1_adult_CE_N <- p1_adult_CE %>% 
  filter(COMID %in% N_Sites$COMID)
p2_adult_CE_N <- p2_adult_CE %>% 
  filter(COMID %in% N_Sites$COMID)
p3_adult_CE_N <- p3_adult_CE %>% 
  filter(COMID %in% N_Sites$COMID)

p1_adult_RE_N <- p1_adult_RE %>% 
  filter(COMID %in% N_Sites$COMID)
p2_adult_RE_N <- p2_adult_RE %>% 
  filter(COMID %in% N_Sites$COMID)
p3_adult_RE_N <- p3_adult_RE %>% 
  filter(COMID %in% N_Sites$COMID)

p1_YOY_CE_S <- p1_YOY_CE %>% 
  filter(COMID %in% S_Sites$COMID)
p2_YOY_CE_S <- p2_YOY_CE %>% 
  filter(COMID %in% S_Sites$COMID)
p3_YOY_CE_S <- p3_YOY_CE %>% 
  filter(COMID %in% S_Sites$COMID)

p1_YOY_RE_S <- p1_YOY_RE %>% 
  filter(COMID %in% S_Sites$COMID)
p2_YOY_RE_S <- p2_YOY_RE %>% 
  filter(COMID %in% S_Sites$COMID)
p3_YOY_RE_S <- p3_YOY_RE %>% 
  filter(COMID %in% S_Sites$COMID)

p1_adult_CE_S <- p1_adult_CE %>% 
  filter(COMID %in% S_Sites$COMID)
p2_adult_CE_S <- p2_adult_CE %>% 
  filter(COMID %in% S_Sites$COMID)
p3_adult_CE_S <- p3_adult_CE %>% 
  filter(COMID %in% S_Sites$COMID)

p1_adult_RE_S <- p1_adult_RE %>% 
  filter(COMID %in% S_Sites$COMID)
p2_adult_RE_S <- p2_adult_RE %>% 
  filter(COMID %in% S_Sites$COMID)
p3_adult_RE_S <- p3_adult_RE %>% 
  filter(COMID %in% S_Sites$COMID)

sample_areas_wide_CE_N <- sample_areas_wide_CE %>% 
  filter(COMID %in% N_Sites$COMID)
sample_areas_wide_RE_N <- sample_areas_wide_RE %>% 
  filter(COMID %in% N_Sites$COMID)
sample_areas_wide_CE_S <- sample_areas_wide_CE %>% 
  filter(COMID %in% S_Sites$COMID)
sample_areas_wide_RE_S <- sample_areas_wide_RE %>% 
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
nReps_CE_N <- nrow(p1_YOY_CE_N)
nReps_RE_N <- nrow(p1_YOY_RE_N)
nReps_CE_S <- nrow(p1_YOY_CE_S)
nReps_RE_S <- nrow(p1_YOY_RE_S)
nAgencies_CE_N <- length(unique(p1_YOY_CE_N$Agency))
nAgencies_RE_N <- length(unique(p1_YOY_RE_N$Agency))
nAgencies_CE_S <- length(unique(p1_YOY_CE_S$Agency))
nAgencies_RE_S <- length(unique(p1_YOY_RE_S$Agency))

## YOY analysis
# Models can stay the same, we just subset the data

# NORTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_CE_N, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE_N,
                  Area = sample_areas_wide_CE_N,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_YOY = p1_YOY_CE_N,
                  p2_YOY = p2_YOY_CE_N, 
                  p3_YOY = p3_YOY_CE_N,
                  Agencies = as.numeric(as.factor(p1_YOY_CE_N$Agency)))


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_CE_N, nYears_CE))
for (i in 1:nReps_CE_N) {
  for (t in 1:nYears_CE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_CE_N[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_CE_N[i,t], p2_YOY_CE_N[i,t], p3_YOY_CE_N[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE_N, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE_N),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

# NORTH - random effects model
# Bundle data
jags_data <- list(nReps = nReps_RE_N, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE_N,
                  Area = sample_areas_wide_RE_N,
                  p1_YOY = p1_YOY_RE_N,
                  p2_YOY = p2_YOY_RE_N, 
                  p3_YOY = p3_YOY_RE_N,
                  Agencies = as.numeric(as.factor(p1_YOY_RE_N$Agency)))


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_RE_N, nYears_RE))
for (i in 1:nReps_RE_N) {
  for (t in 1:nYears_RE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_RE_N[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_RE_N[i,t], p2_YOY_RE_N[i,t], p3_YOY_RE_N[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE_N, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE_N, 0, 10),
                             gam = array(rnorm(nReps_RE_N * nYears_RE, 0, 10), dim = c(nReps_RE_N, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE_N),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

# Fit Model
YOY_randomEffects_N <- jagsUI::jags(data = jags_data,
                                     parameters.to.save = jags_params,
                                     model.file = "Analysis/nMix_JAGS_files/YOY_randomEffects.jags",
                                     n.chains = nc,
                                     n.iter = ni,
                                     n.burnin = nb,
                                     n.thin = nt,
                                     parallel = T,
                                     inits = init_vals)

YOY_randomEffects_N_params <- MCMCsummary(YOY_randomEffects_N, HPD = T)

# SOUTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_CE_S, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE_S,
                  Area = sample_areas_wide_CE_S,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_YOY = p1_YOY_CE_S,
                  p2_YOY = p2_YOY_CE_S, 
                  p3_YOY = p3_YOY_CE_S,
                  Agencies = as.numeric(as.factor(p1_YOY_CE_S$Agency)))


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_CE_S, nYears_CE))
for (i in 1:nReps_CE_S) {
  for (t in 1:nYears_CE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_CE_S[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_CE_S[i,t], p2_YOY_CE_S[i,t], p3_YOY_CE_S[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE_S, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE_S),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

# SOUTH - random effects model
# Bundle data
jags_data <- list(nReps = nReps_RE_S, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE_S,
                  Area = sample_areas_wide_RE_S,
                  p1_YOY = p1_YOY_RE_S,
                  p2_YOY = p2_YOY_RE_S, 
                  p3_YOY = p3_YOY_RE_S,
                  Agencies = as.numeric(as.factor(p1_YOY_RE_S$Agency)))


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.YOY", "mean_ICC.YOY", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.YOY. Initial values must all be great than or equal to the sum of observed counts
N.YOY.inits <- array(numeric(), dim = c(nReps_RE_S, nYears_RE))
for (i in 1:nReps_RE_S) {
  for (t in 1:nYears_RE) {
    N.YOY.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_YOY_RE_S[i,t]),
                                                rpois(1, lambda = 200),
                                                sum(p1_YOY_RE_S[i,t], p2_YOY_RE_S[i,t], p3_YOY_RE_S[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE_S, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE_S, 0, 10),
                             gam = array(rnorm(nReps_RE_S * nYears_RE, 0, 10), dim = c(nReps_RE_S, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE_S),
                             N.YOY = N.YOY.inits)

# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

# Fit Model
YOY_randomEffects_S <- jagsUI::jags(data = jags_data,
                                    parameters.to.save = jags_params,
                                    model.file = "Analysis/nMix_JAGS_files/YOY_randomEffects.jags",
                                    n.chains = nc,
                                    n.iter = ni,
                                    n.burnin = nb,
                                    n.thin = nt,
                                    parallel = T,
                                    inits = init_vals)

YOY_randomEffects_S_params <- MCMCsummary(YOY_randomEffects_S, HPD = T)

## Adult analysis
# Models can stay the same, we just subset the data

# NORTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_CE_N, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE_N,
                  Area = sample_areas_wide_CE_N,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_N,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_N,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_N,
                  p1_adult = p1_adult_CE_N,
                  p2_adult = p2_adult_CE_N, 
                  p3_adult = p3_adult_CE_N,
                  Agencies = as.numeric(as.factor(p1_adult_CE_N$Agency)))


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_CE_N, nYears_CE))
for (i in 1:nReps_CE_N) {
  for (t in 1:nYears_CE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_CE_N[i,t]),
                                                  rpois(1, lambda = 200),
                                                  sum(p1_adult_CE_N[i,t], p2_adult_CE_N[i,t], p3_adult_CE_N[i,t], 1, na.rm = T) * 2)))
  }
}


# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE_N, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE_N),
                             N.adult = N.adult.inits)

# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

# NORTH - random effects model
# Bundle data
jags_data <- list(nReps = nReps_RE_N, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE_N,
                  Area = sample_areas_wide_RE_N,
                  p1_adult = p1_adult_RE_N,
                  p2_adult = p2_adult_RE_N, 
                  p3_adult = p3_adult_RE_N,
                  Agencies = as.numeric(as.factor(p1_adult_RE_N$Agency)))


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.adult", "mean_ICC.adult", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_RE_N, nYears_RE))
for (i in 1:nReps_RE_N) {
  for (t in 1:nYears_RE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_RE_N[i,t]),
                                                  rpois(1, lambda = 200),
                                                  sum(p1_adult_RE_N[i,t], p2_adult_RE_N[i,t], p3_adult_RE_N[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE_N, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE_N, 0, 10),
                             gam = array(rnorm(nReps_RE_N * nYears_RE, 0, 10), dim = c(nReps_RE_N, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE_N),
                             N.adult = N.adult.inits)

# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

# Fit Model
Adult_randomEffects_N <- jagsUI::jags(data = jags_data,
                                    parameters.to.save = jags_params,
                                    model.file = "Analysis/nMix_JAGS_files/Adult_randomEffects.jags",
                                    n.chains = nc,
                                    n.iter = ni,
                                    n.burnin = nb,
                                    n.thin = nt,
                                    parallel = T,
                                    inits = init_vals)

Adult_randomEffects_N_params <- MCMCsummary(Adult_randomEffects_N, HPD = T)

# SOUTH - climate model
# Bundle data
jags_data <- list(nReps = nReps_CE_S, 
                  nYears = nYears_CE,
                  nAgencies = nAgencies_CE_S,
                  Area = sample_areas_wide_CE_S,
                  Mean_Max_Summer_Temp_Scaled = Mean_Max_Summer_Temp_Scaled_S,
                  Max_0.9Q_WinterFlow_Scaled = Max_0.9Q_WinterFlow_Scaled_S,
                  Max_0.9Q_SpringFlow_Scaled = Max_0.9Q_SpringFlow_Scaled_S,
                  p1_adult = p1_adult_CE_S,
                  p2_adult = p2_adult_CE_S, 
                  p3_adult = p3_adult_CE_S,
                  Agencies = as.numeric(as.factor(p1_adult_CE_S$Agency)))


# Parameters to save
jags_params <- c("omega", "beta", "mu.beta", "s2.beta", "p", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_CE_S, nYears_CE))
for (i in 1:nReps_CE_S) {
  for (t in 1:nYears_CE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_CE_S[i,t]),
                                                  rpois(1, lambda = 200),
                                                  sum(p1_adult_CE_S[i,t], p2_adult_CE_S[i,t], p3_adult_CE_S[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_CE_S, 0, 0.001),
                             sd.beta = runif(3, 0, 10),
                             mu.beta = rnorm(3, -0.5, 0.01),
                             p = rep(0.5, times = nAgencies_CE_S),
                             N.adult = N.adult.inits)
# MCMC settings
ni <- 25000
nc <- 3
nb <- 5000
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

# SOUTH - random effects model
# Bundle data
jags_data <- list(nReps = nReps_RE_S, 
                  nYears = nYears_RE,
                  nAgencies = nAgencies_RE_S,
                  Area = sample_areas_wide_RE_S,
                  p1_adult = p1_adult_RE_S,
                  p2_adult = p2_adult_RE_S, 
                  p3_adult = p3_adult_RE_S,
                  Agencies = as.numeric(as.factor(p1_adult_RE_S$Agency)))


# Parameters to save
jags_params <- c("omega", "p", "s2.eps",  "s2.gam",  "ICC.adult", "mean_ICC.adult", "pval.mean_p1", "pval.CV_p1")

# create and populate an array of initial values for N.adult. Initial values must all be great than or equal to the sum of observed counts
N.adult.inits <- array(numeric(), dim = c(nReps_RE_S, nYears_RE))
for (i in 1:nReps_RE_S) {
  for (t in 1:nYears_RE) {
    N.adult.inits[i,t] <- round(as.numeric(ifelse(is.na(p1_adult_RE_S[i,t]),
                                                  rpois(1, lambda = 200),
                                                  sum(p1_adult_RE_S[i,t], p2_adult_RE_S[i,t], p3_adult_RE_S[i,t], 1, na.rm = T) * 2)))
  }
}

# Set initial values
init_vals <- function() list(omega = rnorm(nReps_RE_S, 0, 0.001),
                             sd.eps = runif(1, 0, 10),
                             eps = rnorm(nYears_RE, 0, 10),
                             sd.gam = runif(nReps_RE_S, 0, 10),
                             gam = array(rnorm(nReps_RE_S * nYears_RE, 0, 10), dim = c(nReps_RE_S, nYears_RE)),
                             p = rep(0.5, times = nAgencies_RE_S),
                             N.adult = N.adult.inits)

# MCMC settings
ni <- 50000
nc <- 3
nb <- 20000
nt <- 1

# Fit Model
Adult_randomEffects_S <- jagsUI::jags(data = jags_data,
                                      parameters.to.save = jags_params,
                                      model.file = "Analysis/nMix_JAGS_files/Adult_randomEffects.jags",
                                      n.chains = nc,
                                      n.iter = ni,
                                      n.burnin = nb,
                                      n.thin = nt,
                                      parallel = T,
                                      inits = init_vals)

Adult_randomEffects_S_params <- MCMCsummary(Adult_randomEffects_S, HPD = T)

#MCMCtrace(Adult_randomEffects_S, params = "omega", pdf = F)

save.image()
######################################
### Post-Hoc
######################################

# make a dataframe of segments with coordinates for mapping
COMIDs_agency <- SE_Site_Final %>% 
  filter(COMID %in% p1_YOY_CE$COMID) %>% 
  left_join(source_agency_crosswalk) %>% 
  group_by(COMID,
           Agency) %>% 
  summarize(State = first(State),
            Lat = first(Lat),
            Long = first(Long),
            Elev_m = first(Elev_m))

segment_data_CE <- data.frame(Agency = agencies) %>% 
  mutate(agencyNum = row_number()) %>% 
  right_join(p1_YOY_CE[,c("COMID", "Agency")], by = c("agencyNum" = "Agency")) %>% 
  left_join(COMIDs_agency, by = c("COMID", "Agency")) %>% 
  dplyr::select(-agencyNum) %>% 
  arrange(COMID)

segment_data_RE <- data.frame(Agency = agencies_RE) %>% 
  mutate(agencyNum = row_number()) %>% 
  right_join(p1_YOY_RE[,c("COMID", "Agency")], by = c("agencyNum" = "Agency")) %>% 
  left_join(COMIDs_agency, by = c("COMID", "Agency")) %>% 
  dplyr::select(-agencyNum) %>% 
  arrange(COMID)

# Export to make the map in QGIS
fwrite(segment_data_CE, here("Data", "SE_Trout_Segments_filtered.csv"))

# Create a map of segments
US_states <- map_data("state")

segments_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = segment_data_CE, 
             aes(x = Long, y = Lat, color = Agency), alpha = 0.5) +
  coord_map("albers",
            parameters = c(29.5, 45.5),
            xlim = c(-85, -76),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       color = "Agency") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() + 
  theme(text = element_text(family =  "serif"))

# for presentations
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = segment_data,
#              aes(x = Long, y = Lat, color = Agency), alpha = 0.5) +
#   geom_hline(aes(yintercept = 37.13),
#              linetype = "dashed") +
#   coord_map("albers",
#             parameters = c(29.5, 45.5),
#             xlim = c(-83.75, -71),
#             ylim = c(33, 41.7)) +
#   labs(color = "") +
#   scale_color_brewer(palette = "Set1") +
#   theme_void() +
#   theme(legend.position = "none")

# make north and south versions too
segment_data_CE_N <- segment_data_CE %>% 
  filter(COMID %in% N_Sites$COMID)
segment_data_RE_N <- segment_data_RE %>% 
  filter(COMID %in% N_Sites$COMID)

segment_data_CE_S <- segment_data_CE %>% 
  filter(COMID %in% S_Sites$COMID)
segment_data_RE_S <- segment_data_RE %>% 
  filter(COMID %in% S_Sites$COMID)

######################################
# Create a table of agency data
agencies.table <- p1_YOY_CE %>% 
  .[,c(1:34, 37)] %>% 
  pivot_longer(cols = 1:34,
               names_to = "Year",
               values_to = "Count") %>% 
  group_by(agency = Agency, Year) %>% 
  dplyr::summarise(Count = sum(Count, na.rm = T)) %>% 
  ungroup(Year) %>% 
  filter(Count > 0) %>% 
  dplyr::summarise(Data_Range = paste(min(Year), "-", max(Year)),
                   NYears_Data = n()) %>% 
  cbind(Agency = agencies) %>%
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

# Climate covariates
summTemp_N <- SE_COMID_temp_covars %>% 
  filter(COMID %in% N_Sites$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_CE$Year) # filter for years which we have trout data
summTemp_S <- SE_COMID_temp_covars %>% 
  filter(COMID %in% S_Sites$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_CE$Year) # filter for years which we have trout data

streamFlow_N <- SE_COMID_flow_covars %>% 
  filter(COMID %in% N_Sites$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_CE$Year) # filter for years which we have trout data
streamFlow_S <- SE_COMID_flow_covars %>% 
  filter(COMID %in% S_Sites$COMID, # Filter to just data for which we have samples
         Year %in% YOY_BKT_passCounts_CE$Year) # filter for years which we have trout data

Segment_Summary.table <- data.frame(Variable = c("Channel Slope (%)",
                                                 "Length (km)",
                                                 "Catchment area (km^2)",
                                                 "Elevation (m)",
                                                 "Stream order",
                                                 "Wetted width (m)",
                                                 "Max Summer\nTemperature (C)",
                                                 "Max 0.9Q Winter\nStream Flow (CFS)",
                                                 "Max 0.9Q Spring\nStream Flow (CFS)"),
                                    North_Mean = c(mean(SE_segments_NHDPlus_N$SLOPE, na.rm = T)*100,
                                             mean(SE_segments_NHDPlus_N$LENGTHKM, na.rm = T),
                                             mean(SE_segments_NHDPlus_N$AreaSqKM, na.rm = T),
                                             mean(rowMeans(SE_segments_NHDPlus_N[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                             round(mean(SE_segments_NHDPlus_N$StreamOrde, na.rm = T), 0),
                                             mean(SE_segments_N$Width_m, na.rm = T),
                                             mean(summTemp_N$Mean_Max_Summer_Temp, na.rm = T),
                                             mean(streamFlow_N$Max_0.9Q_WinterFlow, na.rm = T),
                                             mean(streamFlow_N$Max_0.9Q_SpringFlow, na.rm = T)),
                                    North_SD = c(sd(SE_segments_NHDPlus_N$SLOPE, na.rm = T)*100,
                                                sd(SE_segments_NHDPlus_N$LENGTHKM, na.rm = T),
                                                sd(SE_segments_NHDPlus_N$AreaSqKM, na.rm = T),
                                                sd(rowMeans(SE_segments_NHDPlus_N[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                round(sd(SE_segments_NHDPlus_N$StreamOrde, na.rm = T), 0),
                                                sd(SE_segments_N$Width_m, na.rm = T),
                                                sd(summTemp_N$Mean_Max_Summer_Temp, na.rm = T),
                                                sd(streamFlow_N$Max_0.9Q_WinterFlow, na.rm = T),
                                                sd(streamFlow_N$Max_0.9Q_SpringFlow, na.rm = T)),
                                    South_Mean = c(mean(SE_segments_NHDPlus_S$SLOPE, na.rm = T)*100,
                                                   mean(SE_segments_NHDPlus_S$LENGTHKM, na.rm = T),
                                                   mean(SE_segments_NHDPlus_S$AreaSqKM, na.rm = T),
                                                   mean(rowMeans(SE_segments_NHDPlus_S[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                   round(mean(SE_segments_NHDPlus_S$StreamOrde, na.rm = T), 0),
                                                   mean(SE_segments_S$Width_m, na.rm = T),
                                                   mean(summTemp_S$Mean_Max_Summer_Temp, na.rm = T),
                                                   mean(streamFlow_S$Max_0.9Q_WinterFlow, na.rm = T),
                                                   mean(streamFlow_S$Max_0.9Q_SpringFlow, na.rm = T)),
                                    South_SD = c(sd(SE_segments_NHDPlus_S$SLOPE, na.rm = T)*100,
                                                 sd(SE_segments_NHDPlus_S$LENGTHKM, na.rm = T),
                                                 sd(SE_segments_NHDPlus_S$AreaSqKM, na.rm = T),
                                                 sd(rowMeans(SE_segments_NHDPlus_S[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                 round(sd(SE_segments_NHDPlus_S$StreamOrde, na.rm = T), 0),
                                                 sd(SE_segments_S$Width_m, na.rm = T),
                                                 sd(summTemp_S$Mean_Max_Summer_Temp, na.rm = T),
                                                 sd(streamFlow_S$Max_0.9Q_WinterFlow, na.rm = T),
                                                 sd(streamFlow_S$Max_0.9Q_SpringFlow, na.rm = T)))

# What is the average length of *sites* that were electrofished?
SE_Site_Final %>% 
  filter(COMID %in% sample_areas$COMID) %>% 
  summarise( mean_len = mean(Length_m, na.rm = T))

#####################################
# Create a table of covariate summaries

# Filter for temp data at sites and years where we have trout data
temp_summary_data <- SE_COMID_temp_covars %>% 
  filter(COMID %in% p1_YOY_CE$COMID,
         Year %in% 1981:2015)

flow_summary_data <- SE_COMID_flow_covars %>% 
  filter(COMID %in% p1_YOY_CE$COMID,
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
         COMID %in% p1_YOY_CE$COMID) %>% 
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
              linewidth = 0.33,
              se = FALSE) +
  geom_smooth(method = "lm", 
              color = "blue",
              se = FALSE) +
  theme_classic() +
  labs(x = "Log (mean YOY abundance) in year t",
       y = "Log (mean adult abundance) in year t+1") +
  theme(aspect.ratio = 1)

SSR_2.plot <- stock_recruit_data %>% 
  ggplot(aes(x = log_Adults_t,
             y = log_YOY_tp1)) +
  #geom_point(size = 0.5) +
  geom_smooth(aes(group = as.factor(COMID)), 
              method = "lm", 
              color = "black",
              linewidth = 0.33,
              se = FALSE) +
  geom_smooth(method = "lm", 
              color = "blue",
              se = FALSE) +
  theme_classic() +
  labs(x = "Log (mean adult abundance) in year t",
       y = "Log (mean YOY abundance) in year t+1") +
  theme(aspect.ratio = 1)

###################
# PPC p-values

YOY_climateEffects_PPCs <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(2)

YOY_randomEffects_PPCs <- YOY_randomEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(2)

Adult_climateEffects_PPCs <- Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(2)

Adult_randomEffects_PPCs <- Adult_randomEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(param %in% c("pval.mean_p1", "pval.CV_p1")) %>% 
  dplyr::select(2)

PPC_pvals.table <- data.frame(Life_Stage = c(rep("YOY", 4),
                                             rep("Adult", 4)),
                              stat = rep(c("mean", "CV"), 4),
                              model = rep(c(rep("climateEffects", 2), rep("randomEffects", 2)), 2),
                              pVal = rbind(YOY_climateEffects_PPCs,
                                           YOY_randomEffects_PPCs,
                                           Adult_climateEffects_PPCs,
                                           Adult_randomEffects_PPCs))

###################
## ICC Values
# Join site data to ICC values
YOY_ICCs.table <- YOY_randomEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data_RE)

Adult_ICCs.table <- Adult_randomEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data_RE)

YOY_ICCs_N.table <- YOY_randomEffects_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data_RE_N)

YOY_ICCs_S.table <- YOY_randomEffects_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.YOY\\[")) %>% 
  cbind(segment_data_RE_S)

Adult_ICCs_N.table <- Adult_randomEffects_N_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data_RE_N)

Adult_ICCs_S.table <- Adult_randomEffects_S_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "ICC.adult\\[")) %>% 
  cbind(segment_data_RE_S)

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
             aes(x = Long, y = Lat), shape = 5) + # diamonds
  geom_point(data = head(arrange(YOY_ICCs_S.table, mean)),
             aes(x = Long, y = Lat), shape = 1) + # circles
  geom_point(data = head(arrange(YOY_ICCs.table, mean)),
             aes(x = Long, y = Lat), shape = 3) + # crosses
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  theme_classic()
  #theme(text = element_text(family =  "serif"))

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
  geom_point(data = head(arrange(YOY_ICCs.table, mean), 5),
             aes(x = Long, y = Lat), shape = 4, size = 2, stroke = 1) + # crosses
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for YOY BKT",
       color = "ICC") +
  scale_color_viridis() +
  theme_classic()  
  #theme(text = element_text(family =  "serif"))

# # plot for presentations
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = YOY_ICCs.table,
#           aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
#   geom_point(data = head(arrange(YOY_ICCs.table, mean), 5),
#              aes(x = Long, y = Lat), shape = 3) + # crosses
#   geom_segment(aes(x = -85, xend = -75,
#                    y = 37.13, yend = 37.13),
#              linetype = "dashed") +
#   coord_map("albers",
#             parameters = c(29.5, 45.5),
#             xlim = c(-83.75, -71),
#             ylim = c(33, 41.7)) +
#   labs(color = "ICC") +
#   scale_color_viridis(option = "H") +
#   theme_void()

# Adult
Adult_ICC_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = Adult_ICCs.table, 
             aes(x = Long, y = Lat, color = mean), alpha = 0.5) +
  geom_point(data = head(arrange(Adult_ICCs.table, mean), 5),
             aes(x = Long, y = Lat), shape = 4, size = 2, stroke = 1) + # crosses
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -74),
            ylim = c(34.5, 40)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior ICC Means for Adult BKT",
       color = "ICC") +
  scale_color_viridis() +
  theme_classic() 
  #theme(text = element_text(family =  "serif"))

# Is ICC correlated with any site-level variables?
library(corrplot)

# load NHDplus data
NHDplus_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")
# and filter to our sites
NHDplus_data <- NHDplus_data %>% 
  filter(COMID %in% segment_data_RE$COMID)

# load streamcat data
StreamCat_Data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/StreamCat Data/StreamCat_Covars_Combined.csv")
# and filter to our sites
StreamCat_Data <- StreamCat_Data %>% 
  filter(COMID %in% segment_data_RE$COMID)

ICC_corr_data <- YOY_ICCs.table[,c(2,9,11:13)] %>%
  left_join(NHDplus_data) %>% 
  left_join(StreamCat_Data) %>% 
  select(where(is.numeric),
         -COMID,
         -fid,
         -OBJECTID,
         -GNIS_ID,
         -REACHCODE,
         -WBAREACOMI,
         -FCODE,
         -FromNode,
         -ToNode,
         -Hydroseq,
         -LevelPathI,
         -TerminalPa,
         -Divergence,
         -StartFlag,
         -TerminalFl,
         -DnLevel,
         -UpLevelPat,
         -DnLevelPat,
         -UpHydroseq,
         -DnMinorHyd,
         -DnDrainCou,
         -DnHydroseq,
         -RtnDiv,
         -VPUIn,
         -VPUOut,
         -Tidal,
         -ToNode,
         -FromNode,
         -LevelPathI,
         -DnLevelPat,
         -TerminalPa,
         -FromMeas,
         -ToMeas)

YOY_ICC_corrPlot <- corrplot(cor(ICC_corr_data , method="spearman", use="pairwise.complete.obs"))
# get values
YOY_ICC_corrs.table <- as.data.frame(YOY_ICC_corrPlot$corrPos) %>% 
  filter(xName == "mean")
  # mean ICC is not really correlated with any of these site-level covars
YOY_ICC_corrs.table <- YOY_ICC_corrs.table %>% 
  arrange(-abs(corr)) %>% 
  .[2:11,c(2,5)]

####
# Gap analysis for asynchronous sites

# Load stream segments that are within protected areas
# These data come from the USGS protected areas database
PA_StreamSegments <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/SE_Trout/Data/SE_Trout_PA_flowlines.csv")

# % BKT segments with 25% best ICCs in protected areas:
nrow(YOY_ICCs.table %>% filter(mean < quantile(mean, 0.25), COMID %in% PA_StreamSegments$COMID))/nrow(YOY_ICCs.table %>% filter(mean < quantile(mean, 0.25))) * 100

# % of top 5 ICCs in protected areas
nrow(YOY_ICCs.table %>% slice_min(mean, n = 5) %>% filter(COMID %in% PA_StreamSegments$COMID))/nrow(YOY_ICCs.table %>% slice_min(mean, n = 5)) * 100

# % of ALL segments in protected areas
nrow(YOY_ICCs.table %>% filter(COMID %in% PA_StreamSegments$COMID))/nrow(YOY_ICCs.table) * 100

# 25% most synchronous segments
nrow(YOY_ICCs.table %>% filter(mean > quantile(mean, 0.75), COMID %in% PA_StreamSegments$COMID))/nrow(YOY_ICCs.table %>% filter(mean > quantile(mean, 0.75))) * 100

# and 5 most synchronous segments
nrow(YOY_ICCs.table %>% slice_max(mean, n = 5) %>% filter(COMID %in% PA_StreamSegments$COMID))/nrow(YOY_ICCs.table %>% slice_max(mean, n = 5)) * 100

################
# Covariate effects/betas
YOY_mu.betas.table <- data.frame(rbind(data.frame(subregion = rep("All", 3)),
                                   data.frame(subregion = rep("North", 3)),
                                   data.frame(subregion = rep("South", 3))),
                             rbind(data.frame(covar = c("SummTemp", "WintFlow", "SprFlow")),
                                   data.frame(covar = c("SummTemp", "WintFlow", "SprFlow")),
                                   data.frame(covar = c("SummTemp", "WintFlow", "SprFlow"))),
                             rbind(YOY_climateEffects_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "mu.beta")) %>% 
                                     select(2,4,5),
                                   YOY_climateEffects_N_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "mu.beta")) %>% 
                                     select(2,4,5),
                                   YOY_climateEffects_S_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "mu.beta")) %>% 
                                     select(2,4,5)))

Adult_mu.betas.table <- data.frame(rbind(data.frame(subregion = rep("All", 3)),
                                       data.frame(subregion = rep("North", 3)),
                                       data.frame(subregion = rep("South", 3))),
                                 rbind(data.frame(covar = c("SummTemp", "WintFlow", "SprFlow")),
                                       data.frame(covar = c("SummTemp", "WintFlow", "SprFlow")),
                                       data.frame(covar = c("SummTemp", "WintFlow", "SprFlow"))),
                                 rbind(Adult_climateEffects_params %>% 
                                         rownames_to_column(., "param") %>% 
                                         filter(str_detect(param, "mu.beta")) %>% 
                                         select(2,4,5),
                                       Adult_climateEffects_N_params %>% 
                                         rownames_to_column(., "param") %>% 
                                         filter(str_detect(param, "mu.beta")) %>% 
                                         select(2,4,5),
                                       Adult_climateEffects_S_params %>% 
                                         rownames_to_column(., "param") %>% 
                                         filter(str_detect(param, "mu.beta")) %>% 
                                         select(2,4,5)))
  
# Summarize samples of covariate effects in table
# plot data
# mu.beta_samples.table <- rbind(data.frame(sample_val = rbind(as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,1]),
#                                                              as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,2]),
#                                                              as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(YOY_climateEffects_N$sims.list$mu.beta[,1])),
#                                           life_stage = rep("YOY"),
#                                           subregion = rep("North")),
#                                data.frame(sample_val = rbind(as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,1]),
#                                                              as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,2]),
#                                                              as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(Adult_climateEffects_N$sims.list$mu.beta[,1])),
#                                           life_stage = rep("Adult"),
#                                           subregion = rep("North")),
#                                data.frame(sample_val = rbind(as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,1]),
#                                                              as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,2]),
#                                                              as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(YOY_climateEffects_S$sims.list$mu.beta[,1])),
#                                           life_stage = rep("YOY"),
#                                           subregion = rep("South")),
#                                data.frame(sample_val = rbind(as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,1]),
#                                                              as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,2]),
#                                                              as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(Adult_climateEffects_S$sims.list$mu.beta[,1])),
#                                           life_stage = rep("Adult"),
#                                           subregion = rep("South")),
#                                data.frame(sample_val = rbind(as.matrix(YOY_climateEffects$sims.list$mu.beta[,1]),
#                                                              as.matrix(YOY_climateEffects$sims.list$mu.beta[,2]),
#                                                              as.matrix(YOY_climateEffects$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(YOY_climateEffects$sims.list$mu.beta[,1])),
#                                           life_stage = rep("YOY"),
#                                           subregion = rep("Both")),
#                                data.frame(sample_val = rbind(as.matrix(Adult_climateEffects$sims.list$mu.beta[,1]),
#                                                              as.matrix(Adult_climateEffects$sims.list$mu.beta[,2]),
#                                                              as.matrix(Adult_climateEffects$sims.list$mu.beta[,3])),
#                                           covar = rep(c("Summer Air Temperature",
#                                                         "Winter Flow",
#                                                         "Spring Flow"),
#                                                       each = length(Adult_climateEffects$sims.list$mu.beta[,1])),
#                                           life_stage = rep("Adult"),
#                                           subregion = rep("Both")))

# Summarize samples of covariate effects in table - ADAPTED TO SUBSET FOR 95% CI of samples
mu.beta_samples.table <- rbind(data.frame(sample_val = rbind(as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,1][YOY_climateEffects_N$sims.list$mu.beta[,1] > quantile(YOY_climateEffects_N$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects_N$sims.list$mu.beta[,1] < quantile(YOY_climateEffects_N$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,2][YOY_climateEffects_N$sims.list$mu.beta[,2] > quantile(YOY_climateEffects_N$sims.list$mu.beta[,2], 0.025) & YOY_climateEffects_N$sims.list$mu.beta[,2] < quantile(YOY_climateEffects_N$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(YOY_climateEffects_N$sims.list$mu.beta[,3][YOY_climateEffects_N$sims.list$mu.beta[,3] > quantile(YOY_climateEffects_N$sims.list$mu.beta[,3], 0.025) & YOY_climateEffects_N$sims.list$mu.beta[,3] < quantile(YOY_climateEffects_N$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_climateEffects_N$sims.list$mu.beta[,1][YOY_climateEffects_N$sims.list$mu.beta[,1] > quantile(YOY_climateEffects_N$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects_N$sims.list$mu.beta[,1] < quantile(YOY_climateEffects_N$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("North")),
                               data.frame(sample_val = rbind(as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,1][Adult_climateEffects_N$sims.list$mu.beta[,1] > quantile(Adult_climateEffects_N$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects_N$sims.list$mu.beta[,1] < quantile(Adult_climateEffects_N$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,2][Adult_climateEffects_N$sims.list$mu.beta[,2] > quantile(Adult_climateEffects_N$sims.list$mu.beta[,2], 0.025) & Adult_climateEffects_N$sims.list$mu.beta[,2] < quantile(Adult_climateEffects_N$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(Adult_climateEffects_N$sims.list$mu.beta[,3][Adult_climateEffects_N$sims.list$mu.beta[,3] > quantile(Adult_climateEffects_N$sims.list$mu.beta[,3], 0.025) & Adult_climateEffects_N$sims.list$mu.beta[,3] < quantile(Adult_climateEffects_N$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_climateEffects_N$sims.list$mu.beta[,1][Adult_climateEffects_N$sims.list$mu.beta[,1] > quantile(Adult_climateEffects_N$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects_N$sims.list$mu.beta[,1] < quantile(Adult_climateEffects_N$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("North")),
                               data.frame(sample_val = rbind(as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,1][YOY_climateEffects_S$sims.list$mu.beta[,1] > quantile(YOY_climateEffects_S$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects_S$sims.list$mu.beta[,1] < quantile(YOY_climateEffects_S$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,2][YOY_climateEffects_S$sims.list$mu.beta[,2] > quantile(YOY_climateEffects_S$sims.list$mu.beta[,2], 0.025) & YOY_climateEffects_S$sims.list$mu.beta[,2] < quantile(YOY_climateEffects_S$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(YOY_climateEffects_S$sims.list$mu.beta[,3][YOY_climateEffects_S$sims.list$mu.beta[,3] > quantile(YOY_climateEffects_S$sims.list$mu.beta[,3], 0.025) & YOY_climateEffects_S$sims.list$mu.beta[,3] < quantile(YOY_climateEffects_S$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_climateEffects_S$sims.list$mu.beta[,1][YOY_climateEffects_S$sims.list$mu.beta[,1] > quantile(YOY_climateEffects_S$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects_S$sims.list$mu.beta[,1] < quantile(YOY_climateEffects_S$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("South")),
                               data.frame(sample_val = rbind(as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,1][Adult_climateEffects_S$sims.list$mu.beta[,1] > quantile(Adult_climateEffects_S$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects_S$sims.list$mu.beta[,1] < quantile(Adult_climateEffects_S$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,2][Adult_climateEffects_S$sims.list$mu.beta[,2] > quantile(Adult_climateEffects_S$sims.list$mu.beta[,2], 0.025) & Adult_climateEffects_S$sims.list$mu.beta[,2] < quantile(Adult_climateEffects_S$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(Adult_climateEffects_S$sims.list$mu.beta[,3][Adult_climateEffects_S$sims.list$mu.beta[,3] > quantile(Adult_climateEffects_S$sims.list$mu.beta[,3], 0.025) & Adult_climateEffects_S$sims.list$mu.beta[,3] < quantile(Adult_climateEffects_S$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_climateEffects_S$sims.list$mu.beta[,1][Adult_climateEffects_S$sims.list$mu.beta[,1] > quantile(Adult_climateEffects_S$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects_S$sims.list$mu.beta[,1] < quantile(Adult_climateEffects_S$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("South")),
                               data.frame(sample_val = rbind(as.matrix(YOY_climateEffects$sims.list$mu.beta[,1][YOY_climateEffects$sims.list$mu.beta[,1] > quantile(YOY_climateEffects$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects$sims.list$mu.beta[,1] < quantile(YOY_climateEffects$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(YOY_climateEffects$sims.list$mu.beta[,2][YOY_climateEffects$sims.list$mu.beta[,2] > quantile(YOY_climateEffects$sims.list$mu.beta[,2], 0.025) & YOY_climateEffects$sims.list$mu.beta[,2] < quantile(YOY_climateEffects$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(YOY_climateEffects$sims.list$mu.beta[,3][YOY_climateEffects$sims.list$mu.beta[,3] > quantile(YOY_climateEffects$sims.list$mu.beta[,3], 0.025) & YOY_climateEffects$sims.list$mu.beta[,3] < quantile(YOY_climateEffects$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(YOY_climateEffects$sims.list$mu.beta[,1][YOY_climateEffects$sims.list$mu.beta[,1] > quantile(YOY_climateEffects$sims.list$mu.beta[,1], 0.025) & YOY_climateEffects$sims.list$mu.beta[,1] < quantile(YOY_climateEffects$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("YOY"),
                                          subregion = rep("Both")),
                               data.frame(sample_val = rbind(as.matrix(Adult_climateEffects$sims.list$mu.beta[,1][Adult_climateEffects$sims.list$mu.beta[,1] > quantile(Adult_climateEffects$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects$sims.list$mu.beta[,1] < quantile(Adult_climateEffects$sims.list$mu.beta[,1], 0.975)]),
                                                             as.matrix(Adult_climateEffects$sims.list$mu.beta[,2][Adult_climateEffects$sims.list$mu.beta[,2] > quantile(Adult_climateEffects$sims.list$mu.beta[,2], 0.025) & Adult_climateEffects$sims.list$mu.beta[,2] < quantile(Adult_climateEffects$sims.list$mu.beta[,2], 0.975)]),
                                                             as.matrix(Adult_climateEffects$sims.list$mu.beta[,3][Adult_climateEffects$sims.list$mu.beta[,3] > quantile(Adult_climateEffects$sims.list$mu.beta[,3], 0.025) & Adult_climateEffects$sims.list$mu.beta[,3] < quantile(Adult_climateEffects$sims.list$mu.beta[,3], 0.975)])),
                                          covar = rep(c("Summer Air Temperature",
                                                        "Winter Flow",
                                                        "Spring Flow"),
                                                      each = length(Adult_climateEffects$sims.list$mu.beta[,1][Adult_climateEffects$sims.list$mu.beta[,1] > quantile(Adult_climateEffects$sims.list$mu.beta[,1], 0.025) & Adult_climateEffects$sims.list$mu.beta[,1] < quantile(Adult_climateEffects$sims.list$mu.beta[,1], 0.975)])),
                                          life_stage = rep("Adult"),
                                          subregion = rep("Both")))


# Reorder the subregions so the N+S region plots first
mu.beta_samples.table$subregion <- factor(mu.beta_samples.table$subregion, c("North", "South", "Both"))
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
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  theme_classic() +
  #theme(panel.grid.major.y = element_line(colour = "grey80")) +
  geom_hline(yintercept=seq(0,length(mu.beta_samples.table[,1]),1)+.5,color="grey80") +
  labs(x = "Value",
       y = element_blank(),
       fill = "Subregion") +
  scale_fill_brewer(palette = "Dark2",
                    limits = c("North", "South", "Both")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 11))

# map betas (segment-specific covariate effects) in space
YOY_climate_effects.table <- rbind(YOY_climateEffects_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta\\[1,")) %>% 
                                     select(mean),
                                   YOY_climateEffects_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta\\[2,")) %>% 
                                     select(mean),
                                   YOY_climateEffects_params %>% 
                                     rownames_to_column(., "param") %>% 
                                     filter(str_detect(param, "beta\\[3,")) %>% 
                                     select(mean)) %>% 
  cbind(segment_data_CE,
        rbind(data.frame(covar = rep("Summer Temperature", 144)),
              data.frame(covar = rep("Winter Flow", 144)),
              data.frame(covar = rep("Spring Flow", 144)))) %>% 
  mutate(bin = cut(mean, breaks = c(Inf, 2, 1, -1, -2, -Inf)))

# reorder the covariates so that summer temperature plots first
YOY_climate_effects.table$covar <- factor(YOY_climate_effects.table$covar, c("Summer Temperature", "Winter Flow", "Spring Flow"))

library(latex2exp)

YOY_betas.plot <- ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_climate_effects.table,
             aes(x = Long, y = Lat, color = mean),
             alpha = 0.75) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.8, 39.7)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)')) +
  #scale_color_viridis(option = "H") +
  scale_color_binned(type = "viridis") +
  theme_classic() +
  facet_grid(covar ~ .)

YOY_SummTemp_betas.plot <- ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Summer Temperature",],
             aes(x = Long, y = Lat, color = mean),
             alpha = 0.75) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.8, 39.7)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)')) +
  #scale_color_viridis(option = "H") +
  scale_color_binned(type = "viridis") +
  theme_classic() 

# plot for presentations
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Summer Temperature",],
#              aes(x = Long, y = Lat, color = mean),
#              alpha = 0.5) +
#   geom_segment(aes(x = -85, xend = -75,
#                    y = 37.13, yend = 37.13),
#              linetype = "dashed") +
#   coord_map("albers",
#             parameters = c(29.5, 45.5),
#             xlim = c(-83.75, -76),
#             ylim = c(33, 41.7)) +
#   labs(color = TeX(r'($\beta$ Value)')) +
#   scale_color_viridis(option = "H") +
#   theme_void()

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
            ylim = c(34.8, 39.7)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)')) +
  #scale_color_viridis(option = "H") +
  scale_color_binned(type = "viridis") +
  theme_classic() 

# # plot for presentations
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Winter Flow",],
#              aes(x = Long, y = Lat, color = mean),
#              alpha = 0.5) +
#   geom_segment(aes(x = -85, xend = -75,
#                    y = 37.13, yend = 37.13),
#              linetype = "dashed") +
#   coord_map("albers",
#             parameters = c(29.5, 45.5),
#             xlim = c(-83.75, -76),
#             ylim = c(33, 41.7)) +
#   labs(color = TeX(r'($\beta$ Value)')) +
#   scale_color_viridis(option = "H") +
#   theme_void()

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
            ylim = c(34.8, 39.7)) +
  labs(x = "Long",
       y = "Lat",
       color = TeX(r'($\beta$ Value)')) +
  #scale_color_viridis(option = "H") +
  scale_color_binned(type = "viridis") +
  theme_classic()

# # plot for presentations
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = YOY_climate_effects.table[YOY_climate_effects.table$covar == "Spring Flow",],
#              aes(x = Long, y = Lat, color = mean),
#              alpha = 0.5) +
#   geom_segment(aes(x = -85, xend = -75,
#                    y = 37.13, yend = 37.13),
#                linetype = "dashed") +
#   coord_map("albers",
#             parameters = c(29.5, 45.5),
#             xlim = c(-83.75, -76),
#             ylim = c(33, 41.7)) +
#   labs(color = TeX(r'($\beta$ Value)')) +
#   scale_color_viridis(option = "H") +
#   theme_void()

####
# plots of first-pass trout count versus each climate covariate where each segment = each line
count_vs_climate_data <- pivot_longer(Mean_Max_Summer_Temp_Scaled, cols = 1:34, names_to = "Year", values_to = "SummTemp") %>% 
  left_join(pivot_longer(Max_0.9Q_WinterFlow_Scaled, cols = 1:34, names_to = "Year", values_to = "WintFlow")) %>%
  left_join(pivot_longer(Max_0.9Q_SpringFlow_Scaled, cols = 1:34, names_to = "Year", values_to = "SprFlow")) %>% 
  left_join(pivot_longer(p1_YOY_CE, cols = 1:34, names_to = "Year", values_to = "count_YOY"))

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
#chart.Correlation(pivot_wider(YOY_climate_effects.table, names_from = covar, values_from = mean)[,4:10], method = "spearman")

beta_corr_data <- YOY_climate_effects.table %>% 
  pivot_wider(names_from = covar, values_from = mean) %>% 
  left_join(NHDplus_data) %>% 
  left_join(StreamCat_Data) %>% 
  select(where(is.numeric),
         -COMID,
         -fid,
         -OBJECTID,
         -GNIS_ID,
         -REACHCODE,
         -WBAREACOMI,
         -FCODE,
         -FromNode,
         -ToNode,
         -Hydroseq,
         -LevelPathI,
         -TerminalPa,
         -Divergence,
         -StartFlag,
         -TerminalFl,
         -DnLevel,
         -UpLevelPat,
         -DnLevelPat,
         -UpHydroseq,
         -DnMinorHyd,
         -DnDrainCou,
         -DnHydroseq,
         -RtnDiv,
         -VPUIn,
         -VPUOut,
         -Tidal)

YOY_beta_corrPlot <- corrplot(cor(beta_corr_data, method="spearman", use="pairwise.complete.obs"))
# get values
YOY_beta_corrs.table <- as.data.frame(YOY_beta_corrPlot$corrPos) %>% 
  filter(xName %in% c("Summer Temperature", "Winter Flow", "Spring Flow"),
         !yName %in% c("Summer Temperature", "Winter Flow", "Spring Flow"))
# betas are not really correlated with any of these site-level covars
YOY_beta_corrs.table <- YOY_beta_corrs.table %>% 
  arrange(-abs(corr)) %>% 
  .[1:10, c(1,2,5)]

####
# What proportion of beta values are significant at 95% HDPIs for YOY?
# Create a table to store these values
pct_effects_negative.table <- data.frame(Covar = c("Summ Temp", "Wint Flow", "Spr Flow"),
                                   pct_sgn_neg = rep(NA, times = 3))

# Summer temp
YOY_summTemp_betas <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[1"))

pct_effects_negative.table[1,2] <- sum(YOY_summTemp_betas[,"95%_HPDU"] < 0)/nrow(YOY_summTemp_betas) * 100

# Winter Flow
YOY_wintFlow_betas <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[2"))

pct_effects_negative.table[2,2] <- sum(YOY_wintFlow_betas[,"95%_HPDU"] < 0)/nrow(YOY_wintFlow_betas) * 100

# Spring Flow
YOY_sprFlow_betas <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[3"))

pct_effects_negative.table[3,2] <- sum(YOY_sprFlow_betas[,"95%_HPDU"] < 0)/nrow(YOY_sprFlow_betas) * 100

# Summer temp
Adult_summTemp_betas <- Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[1"))
# Winter Flow
Adult_wintFlow_betas <- Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[2"))
# Spring Flow
Adult_sprFlow_betas <- Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[3"))

####
# What was the range of covariate effects on YOY abundance?
cov_effect_ranges.table <- data.frame(Life_Stage = c(rep("YOY", 3), rep("Adult", 3)),
                                      Covar = rep(c("Summ Temp", "Wint Flow", "Spr Flow"),2),
                                      min_mean = rbind(min(YOY_summTemp_betas$mean),
                                                  min(YOY_wintFlow_betas$mean),
                                                  min(YOY_sprFlow_betas$mean),
                                                  min(Adult_summTemp_betas$mean),
                                                  min(Adult_wintFlow_betas$mean),
                                                  min(Adult_sprFlow_betas$mean)),
                                      min_95HDPL = rbind(min(YOY_summTemp_betas$`95%_HPDL`),
                                                         min(YOY_wintFlow_betas$`95%_HPDL`),
                                                         min(YOY_sprFlow_betas$`95%_HPDL`),
                                                         min(Adult_summTemp_betas$`95%_HPDL`),
                                                         min(Adult_wintFlow_betas$`95%_HPDL`),
                                                         min(Adult_sprFlow_betas$`95%_HPDL`)),
                                      min_95HDPU = rbind(min(YOY_summTemp_betas$`95%_HPDU`),
                                                         min(YOY_wintFlow_betas$`95%_HPDU`),
                                                         min(YOY_sprFlow_betas$`95%_HPDU`),
                                                         min(Adult_summTemp_betas$`95%_HPDU`),
                                                         min(Adult_wintFlow_betas$`95%_HPDU`),
                                                         min(Adult_sprFlow_betas$`95%_HPDU`)),
                                      max_mean = rbind(max(YOY_summTemp_betas$mean),
                                                  max(YOY_wintFlow_betas$mean),
                                                  max(YOY_sprFlow_betas$mean),
                                                  max(Adult_summTemp_betas$mean),
                                                  max(Adult_wintFlow_betas$mean),
                                                  max(Adult_sprFlow_betas$mean)),
                                      max_95HDPL = rbind(max(YOY_summTemp_betas$`95%_HPDL`),
                                                         max(YOY_wintFlow_betas$`95%_HPDL`),
                                                         max(YOY_sprFlow_betas$`95%_HPDL`),
                                                         max(Adult_summTemp_betas$`95%_HPDL`),
                                                         max(Adult_wintFlow_betas$`95%_HPDL`),
                                                         max(Adult_sprFlow_betas$`95%_HPDL`)),
                                      max_95HDPU = rbind(max(YOY_summTemp_betas$`95%_HPDU`),
                                                         max(YOY_wintFlow_betas$`95%_HPDU`),
                                                         max(YOY_sprFlow_betas$`95%_HPDU`),
                                                         max(Adult_summTemp_betas$`95%_HPDU`),
                                                         max(Adult_wintFlow_betas$`95%_HPDU`),
                                                         max(Adult_sprFlow_betas$`95%_HPDU`))) %>% 
  mutate(HDP_range = max_95HDPU - min_95HDPL,
         mean_range = max_mean - min_mean)

# visualize in a plot
climate_effect_ranges.plot <- ggplot(data = cov_effect_ranges.table) +
  geom_linerange(aes(y = Covar,
                     xmin = min_mean,
                     xmax = max_mean, 
                     color = Life_Stage),
                 position = position_dodge(.25),
                 size = 1) +
  scale_color_brewer(palette = "Dark2", breaks = c("YOY", "Adult")) +
  labs(color = "Life Stage",
    x = "Effect",
    y = element_blank()) +
  theme_classic()

# What were the variances of covariate effects on YOY abundance?
YOY_climate_vars <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^s2.beta")) %>% 
  .[,c(2,4,5)]

Adult_climate_vars <- Adult_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^s2.beta")) %>% 
  .[,c(2,4,5)]

cov_effect_variances.table <- data.frame(Life_Stage = c(rep("YOY", 3), rep("Adult", 3)),
                                         Covar = rep(c("Summ Temp", "Wint Flow", "Spr Flow"), 2),
                                         rbind(YOY_climate_vars,
                                               Adult_climate_vars))

################
# Summarize detection probability in table

Detect_probs <- data.frame(
  Agency = rep(agencies, times = 2),
  Life_Stage = c(rep("YOY", times = 9),
                rep("Adult", times = 9))) %>% 
  # add in jagsUI model summary values
  cbind(rbind(YOY_climateEffects_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]", "p[9]"),1:4],
              Adult_climateEffects_params[c("p[1]", "p[2]", "p[3]","p[4]", "p[5]", "p[6]", "p[7]", "p[8]", "p[9]"),1:4]))



# and make a plot to visualize
Detect_probs.plot <- ggplot(data = Detect_probs) +
  geom_linerange(aes(y = Agency,
                     xmin = `95%_HPDL`, #includes 95% highest density intervals
                     xmax = `95%_HPDU`,
                     color = Life_Stage),
                 position = position_dodge(.25),
                 linewidth = 0.5) +
  geom_point(aes(y = Agency,
                 x = mean,
                 color = Life_Stage),
             position = position_dodge(.25)) +
  scale_color_brewer(palette = "Dark2", breaks = c("YOY", "Adult")) +
  labs(#title = "Covariate Effects on Log Density of BKT",
    color = "Life Stage",
    x = element_blank(),
    y = element_blank()) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 12)) +
  theme_classic() 
  # theme(text = element_text(family =  "serif"),
  #       axis.text.x = element_text(angle = -45, hjust=0))

# create tables of the two life stages for manuscript
Detect_probs_YOY.table <- Detect_probs[Detect_probs$Life_Stage == "YOY",]
Detect_probs_Adult.table <- Detect_probs[Detect_probs$Life_Stage == "Adult",]

################
# Portfolio effect analysis
# How does the range of coefficients of variation (CVs) of observed count among populations compare to the CV of all populations?
# Evidence of portfolio effect: if the overall CV is lower than that of the individual populations

# Calculate CV for each segment
CV.p1 <- rep(NA, nrow(p1_YOY_CE))

for (i in 1:nrow(p1_YOY_CE)) {
  CV.p1[i] <- sd(p1_YOY_CE[i,1:34], na.rm = T)/rowMeans(p1_YOY_CE[i,1:34], na.rm = T)
}

p1_colMean <- rep(NA, 34)

# Calculate CV for all segments
for (t in 1:34) {
  p1_colMean[t] <- colMeans(p1_YOY_CE[,t], na.rm = T)
}

CV.p1_all <- sd(p1_colMean, na.rm = T)/mean(p1_colMean, na.rm = T)

CVs.table <- data.frame(rbind(as.data.frame(CV.p1), CV.p1_all),
                  rbind(as.data.frame(rep("Segment CVs", length(CV.p1))),
                                "Overall CV"))
colnames(CVs.table) <- c("CV", "param")

CVs.plot <- ggplot(CVs.table) +
  geom_boxplot(aes(x = param,
                   y = CV)) +
  labs(x = "Parameter") +
  theme_classic()


########################################################
# Export plots to the results folder
save.image()
# Save the directory to which to save results files
run_dir <- here::here("results/v4.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
