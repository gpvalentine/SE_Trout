# George Valentine
# Script to format covariates for Bayesian Synchrony Analysis

# Desired covariates:
  # Temperature (from DAYMET):
    # - Summer Mean Daily Max
    # - Summer Mean Daily Mean
  # Flow (from D. Carlisle USGS):
    # - Winter (Dec. t-1 to Feb. t) Max Monthly 0.9Q
    # - Spring (Mar. - May) Max Monthly 0.9Q

  # Also download precip (DAYMET) and check correlation with flow

# Load packages
library(tidyverse)
library(data.table)   # Fast file reading and writing
library(lubridate)    # For working with dates

# Load trout data
# Set working directory for data files
filepath <- "/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/Compiled Trout Data (All sources)/"

# Load data
SE_Ind_Final <- fread(file = paste0(filepath, "/SE_Ind_Final.csv"))
SE_Sample_Final <- fread(file = paste0(filepath, "/SE_Sample_Final.csv"))
SE_Site_Final <- fread(file = paste0(filepath, "/SE_Site_Final.csv"))

##############################################

# Make a dataframe of all COMID-Year combinations
SampleYears <- SE_Sample_Final %>% 
  left_join(SE_Site_Final[,c(1,6)]) %>% # Join in COMIDs
  mutate(Year = year(Date)) %>% 
  dplyr::select(COMID, Year) %>% 
  unique() %>% 
  na.omit()

##############################################
# Use the DAYMET tool from Oak Ridge Nat'l Laboratory to download temperature and precip data
# Info at https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html

# load package
library(daymetr)

# daymetr needs a csv file with the following format: site, latitude, longitude
COMIDs_for_daymetr <- SE_Site_Final %>% 
  group_by(COMID) %>% 
  summarise(Lat = mean(Lat),
            Long = mean(Long))

# Export to WD
fwrite(COMIDs_for_daymetr, "COMIDs_for_daymetr.csv", col.names = F, row.names = F)

# get the beginning and end years for which we need data
end <- max(SampleYears$Year)

# now use the daymetr package to download data for those sites
SE_COMIDs_DAYMET_data <- download_daymet_batch(file_location = "COMIDs_for_daymetr.csv",
                                              start = 1980, # the earliest year that DAYMET data are available
                                              end = end,
                                              internal = T,
                                              simplify = T)

# Filter for just the temperature and precipitation measurements that we want
SE_COMIDs_DAYMET_data_wide <- SE_COMIDs_DAYMET_data %>% 
filter(measurement %in% c("tmax..deg.c.",  "tmin..deg.c.", "prcp..mm.day.")) %>% 
  mutate(originDate = paste0(year, "-01-01"), # set a temporary variable as the first day of the year
         Date = as.Date((as.numeric(yday)-1), origin =  originDate),  # and convert the day of the year to a date
         COMID = as.numeric(site),
         value = as.numeric(value)) %>% 
  pivot_wider(id_cols = c(COMID, Date), # Pivot wider to make tmax and tmin their own columns
              names_from = measurement,
              values_from = value,
              names_sort = T)

# Export data
fwrite(SE_COMIDs_DAYMET_data_wide, "/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/SE_COMIDs_DAYMET_data.csv")
# Import data
SE_COMIDs_DAYMET_data <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/SE_COMIDs_DAYMET_data.csv")

# Calculate the mean maximum summer temperature at each COMID and year
SE_COMID_temp_covars <- SE_COMIDs_DAYMET_data %>% 
  mutate(Year = year(Date),
         is.summer = month(Date) %in% c(6,7,8,9)) %>% 
  group_by(COMID, 
           Year,
           is.summer) %>% 
  summarise(Mean_Max_Summer_Temp = mean(tmax..deg.c.)) %>% 
  filter(is.summer == TRUE,
         Year >= 1980,
         Year <= 2015) %>% 
  ungroup() %>% 
  dplyr::select(-is.summer)

# Export temp covars
fwrite(SE_COMID_temp_covars, "Data/SE_Trout_COMID_temp_covars.csv")

#############################################
# Flow data
SE_flow_data <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Site Flow Data/SE_sites_flow_all_6.22.csv")

# check to see if monthly precipitation and flow are correlated - they should be.
# flow <- SE_flow_data %>% 
#   filter(Year >= 1980) %>%
#   dplyr::select(P50_Q)
# 
# precip <- SE_COMIDs_DAYMET_data %>% 
#   filter(year(Date) <= 2015,
#          !is.na(COMID)) %>% 
#   group_by(COMID,
#            year(Date),
#            month(Date)) %>% 
#   summarise(mean_monthly_precip = mean(prcp..mm.day.)) %>% 
#   ungroup() %>% 
#   dplyr::select(-COMID,
#                 -`year(Date)`,
#                 -`month(Date)`)
# 
# cor(x = flow,
#     y = precip)

# Summarize flow by year:
  # Also fall flow?
  # winter max 0.9Q flow: Dec. t-1 to Feb. t
  # spring max 0.9Q flow: Mar. t - May t

# Summarizing across different winters is hard: 
  # create a second year column that assigns each Dec, to the following year so that we can group by year
SE_flow_data <- SE_flow_data %>% 
  mutate(Year2 = ifelse(Month == 12, NA, Year)) %>%
  fill(Year2, .direction = "up") %>% 
  mutate(Year2 = ifelse(Year == 2015 & Month == 12, NA, Year2)) %>% 
  filter(!is.na(Year2))

# Summarize spring flow
spring_flow <- SE_flow_data %>%
  filter(Month %in% c(3,4,5),
         Year >= 1980,
         Year <= 2015) %>% 
  group_by(COMID,
           Year) %>% 
  summarise(Max_0.9Q_SpringFlow = max(P90_Q),
            Year = first(Year))

# Summarize winter flow
winter_flow <- SE_flow_data %>%
  filter(Month %in% c(1,2,12),
         Year >= 1980,
         Year <= 2015) %>% 
  group_by(COMID,
           Year2) %>% 
  summarise(Max_0.9Q_WinterFlow = max(P90_Q))

# Join flow covariates
SE_COMID_flow_covars <- spring_flow %>% 
  left_join(winter_flow, by = c("COMID", "Year" = "Year2"))

# and export
fwrite(SE_COMID_flow_covars, "Data/SE_Trout_COMID_flow_covars.csv")
