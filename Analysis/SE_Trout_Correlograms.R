# George Valentine
# BKT Spatial Synchrony Project
# Synchrony Analysis with Correlograms

# Updates from v3:
  # Remove versioning of model files - now git can take care of it (this would have been v4)
  # Move into organized project file for writing and version control
  # Calculate density as (mean COMID abundance)/ COMID area rather than mean COMID density

# Load Packages
library(data.table) # Allows fast reading and writing of files
library(ggplot2)    # For plotting
library(tidyverse)  # For data processing
library(lubridate)  # For date processing
library(reshape2)   # For reshaping the data into long format
library(ncf)        # For synchrony analysis
library(grid)       # for labeling the gridded plot
library(gridExtra)  # for creating a labeled, gridded plot
library(FSA)        # For making depletion estimates
library(knitr)      # For making tidy tables

# Set working directory for data files
setwd("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/Compiled Trout Data (All sources)")

# Load data
SE_Ind_Final <- fread("SE_Ind_Final.csv")
SE_Sample_Final <- fread("SE_Sample_Final.csv")
SE_Site_Final <- fread("SE_Site_Final.csv")

#######################################################################################################################
#######################################################################################################################
# Visualize the temporal coverage of the different datasets
BKT_Sample_Temporal_Coverage <- SE_Ind_Final %>% 
  filter(LongTerm == 1,
         SPP == "BKT") %>%
  group_by(SampleID, Date, State) %>% 
  summarise(`log(N)` = log(n())) %>% 
  ggplot() +
  geom_point(aes(x = Date, y = State, color = `log(N)`)) +
  theme_bw() +
  theme(text = element_text(size=20))
# seems like the time period roughly 2000 - 2018 has the most coverage

##############################
## DENSITY ##
##############################

######
# SPATIAL SYNCHRONY OF LOG PREDICTED BKT DENSITY BY NHDPlus v2.1 COMID

# get a list of just sites south of the north border of PA
south_sites <- SE_Site_Final %>% 
  filter(Lat <= 39.716667) %>%  # filter for just sites south of the Mason-Dixon Line
  select(SiteID)

# Some sites don't have any BKT. Make a dataframe of just the sites that have >10% BKT. This is Matt Kulp's suggestion rather than using sites with just any BKT
BKT_Sites <- SE_Ind_Final %>% 
  group_by(SiteID,
           SPP) %>% 
  dplyr::summarise(count = n()) %>% 
  ungroup(SPP) %>% 
  mutate(pct_BKT = prop.table(count) * 100) %>% 
  filter(SPP == "BKT",
         pct_BKT >= 10)

# sum BKT abundance at each sample and pass and filter to just multipass data
# multipass data are necessary for predicting abundance
# this will be passed into the removal() function
BKT_sample_sums <- SE_Ind_Final %>% 
  filter(SPP == "BKT",
         SiteID %in% south_sites$SiteID,
         SiteID %in% BKT_Sites$SiteID, # filter for just the sites with records of BKT
         year(Date) %in% c(1995:2015)) %>% # Filter to time frame of interest
  group_by(SampleID,
           PassNo) %>% 
  summarise(SiteID = first(SiteID),
            #State = first(State),
            Year = first(year(Date)),
            Count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = PassNo, values_from = Count)

# Remove samples with only one pass
BKT_sample_sums$count_passes <- rowSums(apply(!is.na(BKT_sample_sums[,4:7]), MARGIN = 2, as.numeric))
BKT_sample_sums <- BKT_sample_sums %>% 
  filter(count_passes > 1) %>% 
  dplyr::select(-count_passes)

# use removal() from the "FSA" package to estimate abundance from multipass data
BKT_AbundEsts <- apply(BKT_sample_sums[,4:7], MARGIN = 1, FUN = removal, 
                       method = "CarleStrub", # The CarleStrub method is most robust
                       just.ests = TRUE)
BKT_AbundEsts <- as.data.frame(t(BKT_AbundEsts))

# Join the abundance estimates to SampleIDs, and summarize by COMID and year
BKT_COMID_Dens <- BKT_sample_sums %>% 
  left_join(SE_Site_Final[,c("SiteID", "Lat", "Long", "Length_m", "Width_m", "COMID")], by = "SiteID") %>%  # left join in site data (we mostly just want lat/long here)
  cbind(Abund_Est = BKT_AbundEsts[,1]) %>% 
  group_by(SiteID, Year) %>% 
  dplyr::summarise(COMID = first(COMID),
                   Lat = first(Lat),
                   Long = first(Long),
                   SiteID_Avg_Abund = mean(Abund_Est), # avg abundance at site in year
                   SiteID_Area = Length_m * Width_m) %>% 
  group_by(COMID, Year) %>% 
  dplyr::summarise(Lat = mean(Lat),
                   Long = mean(Long),
                   Avg_Density = sum(SiteID_Avg_Abund)/(sum(SiteID_Area)/1000)) # avg density/1000m^2/COMID


# Use dcast to "widen" the data by year
BKT_COMID_Dens_wide <- dcast(data = BKT_COMID_Dens, 
                                 formula = COMID + Lat + Long ~ Year,
                                 value.var = "Avg_Density")

# Make a new column with a count of the number of years of nonzero data, and filter out sites that have fewer than five years' of data during this time period
BKT_COMID_Dens_wide <- BKT_COMID_Dens_wide %>% 
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 5,
         !is.na(COMID)) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

# what percentage of observations are missing data?
mean(is.na(BKT_COMID_Dens_wide[,4:ncol(BKT_COMID_Dens_wide)]))

# get polygons for us states to use as a basemap. These are built into ggplot
states_map <- map_data("state")

# make a quick map to visualize what sites we will be using
BKT_COMID_Sites_Map <- ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_point(data = BKT_COMID_Dens_wide, aes(x = Long, y = Lat)) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude") + 
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.4, -72),
            ylim = c(34.5,45)) +
  theme(text = element_text(size=20))

BKT_COMID_Sites_Map

# Save the map
# ggsave("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/BKT_COMID_Sites_Map.jpeg",
#        )

# calculate the average density of fish/year at each site
BKT_COMID_Dens_wide$mean_density <- (rowSums(BKT_COMID_Dens_wide[,4:ncol(BKT_COMID_Dens_wide)], na.rm = T)/apply(!is.na(BKT_COMID_Dens_wide[,4:ncol(BKT_COMID_Dens_wide)]), MARGIN = 1, sum))

# then replace the NAs with the average density of fish per year at that COMID
for (i in 1:nrow(BKT_COMID_Dens_wide)) {
  BKT_COMID_Dens_wide[i,4:ncol(BKT_COMID_Dens_wide)][is.na(BKT_COMID_Dens_wide[i,4:ncol(BKT_COMID_Dens_wide)])] <- BKT_COMID_Dens_wide[i, "mean_density"]
  print(i)
}

# remove the column where we stored mean density
BKT_COMID_Dens_wide <- BKT_COMID_Dens_wide %>% 
  select(-mean_density)

# Now make a new dataframe with log density
# calculating synchrony based on log density helps remove the undue influence of outliers
BKT_COMID_logDens_wide <- BKT_COMID_Dens_wide
BKT_COMID_logDens_wide[,4:ncol(BKT_COMID_logDens_wide)] <- log(BKT_COMID_logDens_wide[,4:ncol(BKT_COMID_logDens_wide)])

# Now, calculate spatial synchrony in the time series using the Sncf() function
BKT_COMID_logDens.corr <- Sncf(x=BKT_COMID_logDens_wide$Long,
                               y=BKT_COMID_logDens_wide$Lat,
                               z=BKT_COMID_logDens_wide[,4:ncol(BKT_COMID_logDens_wide)],
                               resamp = 1000, 
                               latlon = T,
                               xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide$Long, BKT_COMID_logDens_wide$Lat)))))

plot(BKT_COMID_logDens.corr)
summary(BKT_COMID_logDens.corr)

# Plot using ggplot for export
# create a dataframe of the predicted values
BKT.COMID.logDens.corr.pred.df <- data.frame(x = matrix(unlist(BKT_COMID_logDens.corr$real$predicted$x)),
                                             y = matrix(unlist(BKT_COMID_logDens.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
BKT.COMID.logDens.corr.bootValues.df <- as.data.frame(t(BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
BKT.COMID.logDens.corr.boot.df <- data.frame(x = matrix(unlist(BKT_COMID_logDens.corr$boot$boot.summary$predicted$x)),
                                             ymin = BKT.COMID.logDens.corr.bootValues.df$`0.025`,
                                             ymax = BKT.COMID.logDens.corr.bootValues.df$`0.975`)

# make the plot
BKT.COMID.logDens.synch.plot <- ggplot() +
  geom_ribbon(data = BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey") +
  geom_line(data = BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = BKT.COMID.logDens.corr$real$cbar,
             linetype = "dashed") +
  theme_classic() +
  lims(x = c(0, max(BKT.COMID.logDens.corr.pred.df$x)),
       y = c(-0.25, 0.5)) +
  labs(title = "c.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=12))

BKT.COMID.logDens.synch.plot

# Does synchrony differ between the north and south regions of the SE US?
# Define these regions: we draw the line at Radford, VA in the New River Valley. Per the suggestion of Jake Rash

## North
BKT_COMID_logDens_wide_N <- BKT_COMID_logDens_wide %>% 
  filter(Lat > 37.13)

BKT_COMID_logDens.corr_N <- Sncf(x=BKT_COMID_logDens_wide_N$Long,
                                     y=BKT_COMID_logDens_wide_N$Lat,
                                     z=BKT_COMID_logDens_wide_N[,4:ncol(BKT_COMID_logDens_wide_N)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide_N$Long, BKT_COMID_logDens_wide_N$Lat)))))

plot(BKT_COMID_logDens.corr_N)
summary(BKT_COMID_logDens.corr_N)

## South
BKT_COMID_logDens_wide_S <- BKT_COMID_logDens_wide %>% 
  filter(Lat <= 37.13)

BKT_COMID_logDens.corr_S <- Sncf(x=BKT_COMID_logDens_wide_S$Long,
                                     y=BKT_COMID_logDens_wide_S$Lat,
                                     z=BKT_COMID_logDens_wide_S[,4:ncol(BKT_COMID_logDens_wide_S)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide_S$Long, BKT_COMID_logDens_wide_S$Lat)))))

plot(BKT_COMID_logDens.corr_S)
summary(BKT_COMID_logDens.corr_S)

######
# SPATIAL SYNCHRONY OF LOG PREDICTED YOY BKT DENSITY BY NHDPlus v2.1 COMID

# sum YOY BKT abundance at each sample and pass and filter to just multipass data
# multipass data are necessary for predicting abundance
# this will be passed into the removal() function
YOY_BKT_sample_sums <- SE_Ind_Final %>% 
  filter(SPP == "BKT",
         SiteID %in% south_sites$SiteID,
         TL_mm <= 90, # filter here for YOY
         year(Date) %in% c(1995:2015)) %>% # Filter to time frame of interest
  group_by(SampleID,
           PassNo) %>% 
  summarise(SiteID = first(SiteID),
            Year = first(year(Date)),
            Count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = PassNo, values_from = Count)

# Remove samples with only one pass
YOY_BKT_sample_sums$count_passes <- rowSums(apply(!is.na(YOY_BKT_sample_sums[,4:7]), MARGIN = 2, as.numeric))
YOY_BKT_sample_sums <- YOY_BKT_sample_sums %>% 
  filter(count_passes > 1) %>% 
  dplyr::select(-count_passes)

# use removal() from the "FSA" package to estimate YOY abundance from multipass data
YOY_BKT_AbundEsts <- apply(YOY_BKT_sample_sums[,4:7], MARGIN = 1, FUN = removal, 
                           method = "CarleStrub", # The CarleStrub method is most robust
                           just.ests = TRUE)
YOY_BKT_AbundEsts <- as.data.frame(t(YOY_BKT_AbundEsts))

# Join the YOY abundance estimates to SampleIDs, and summarize by COMID and year
YOY_BKT_COMID_Dens <- YOY_BKT_sample_sums %>% 
  left_join(SE_Site_Final[,c("SiteID", "Lat", "Long", "Length_m", "Width_m", "COMID")], by = "SiteID") %>%  # left join in site data (we mostly just want lat/long here)
  cbind(Abund_Est = YOY_BKT_AbundEsts[,1]) %>% 
  group_by(SiteID, Year) %>% 
  dplyr::summarise(COMID = first(COMID),
                   Lat = first(Lat),
                   Long = first(Long),
                   Avg_Abund = mean(Abund_Est),
                   Avg_Density = Avg_Abund/((Length_m * Width_m)/1000)) %>%  # avg density/Site
  group_by(COMID, Year) %>% 
  dplyr::summarise(Lat = mean(Lat),
                   Long = mean(Long),
                   Avg_Density = mean(Avg_Density)) # Avg density/COMID


# Use dcast to "widen" the data by year
YOY_BKT_COMID_Dens_wide <- dcast(data = YOY_BKT_COMID_Dens, 
                                 formula = COMID + Lat + Long ~ Year,
                                 value.var = "Avg_Density")

# Make a new column with a count of the number of years of nonzero data, and filter out sites that have fewer than five years' of data during this time period
YOY_BKT_COMID_Dens_wide <- YOY_BKT_COMID_Dens_wide %>%  
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 5,
         !is.na(COMID)) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

# what percentage of observations are missing data?
mean(is.na(YOY_BKT_COMID_Dens_wide[,4:ncol(YOY_BKT_COMID_Dens_wide)]))

# make a quick map to visualize what sites we will be using
ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_point(data = YOY_BKT_COMID_Dens_wide, aes(x = Long, y = Lat)) +
  theme_classic() +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.4, -72),
            ylim = c(34.5,45))

# calculate the average density of YOY fish/year at each site
YOY_BKT_COMID_Dens_wide$mean_density <- (rowSums(YOY_BKT_COMID_Dens_wide[,4:ncol(YOY_BKT_COMID_Dens_wide)], na.rm = T)/apply(!is.na(YOY_BKT_COMID_Dens_wide[,4:ncol(YOY_BKT_COMID_Dens_wide)]), MARGIN = 1, sum))

# then replace the NAs with the average density of YOY fish per year at that COMID
for (i in 1:nrow(YOY_BKT_COMID_Dens_wide)) {
  YOY_BKT_COMID_Dens_wide[i,4:ncol(YOY_BKT_COMID_Dens_wide)][is.na(YOY_BKT_COMID_Dens_wide[i,4:ncol(YOY_BKT_COMID_Dens_wide)])] <- YOY_BKT_COMID_Dens_wide[i, "mean_density"]
  print(i)
}

# remove the column where we stored mean density
YOY_BKT_COMID_Dens_wide <- YOY_BKT_COMID_Dens_wide %>% 
  select(-mean_density)

# Now make a new dataframe with log density
# calculating synchrony based on log YOY density helps remove the undue influence of outliers
YOY_BKT_COMID_logDens_wide <- YOY_BKT_COMID_Dens_wide
YOY_BKT_COMID_logDens_wide[,4:ncol(YOY_BKT_COMID_logDens_wide)] <- log(YOY_BKT_COMID_logDens_wide[,4:ncol(YOY_BKT_COMID_logDens_wide)])

# Now, calculate spatial synchrony in the time series using the Sncf() function
YOY_BKT_COMID_logDens.corr <- Sncf(x=YOY_BKT_COMID_logDens_wide$Long,
                                   y=YOY_BKT_COMID_logDens_wide$Lat,
                                   z=YOY_BKT_COMID_logDens_wide[,4:ncol(YOY_BKT_COMID_logDens_wide)],
                                   resamp = 1000, 
                                   latlon = T,
                                   xmax = ((2/3)*max(na.omit(gcdist(YOY_BKT_COMID_logDens_wide$Long, YOY_BKT_COMID_logDens_wide$Lat)))))

plot(YOY_BKT_COMID_logDens.corr)
summary(YOY_BKT_COMID_logDens.corr)

# Plot using ggplot for export
# create a dataframe of the predicted values
YOY_BKT.COMID.logDens.corr.pred.df <- data.frame(x = matrix(unlist(YOY_BKT_COMID_logDens.corr$real$predicted$x)),
                                                 y = matrix(unlist(YOY_BKT_COMID_logDens.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
YOY_BKT.COMID.logDens.corr.bootValues.df <- as.data.frame(t(YOY_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
YOY_BKT.COMID.logDens.corr.boot.df <- data.frame(x = matrix(unlist(YOY_BKT_COMID_logDens.corr$boot$boot.summary$predicted$x)),
                                                 ymin = YOY_BKT.COMID.logDens.corr.bootValues.df$`0.025`,
                                                 ymax = YOY_BKT.COMID.logDens.corr.bootValues.df$`0.975`)

# make the plot
YOY_BKT.COMID.logDens.synch.plot <- ggplot() +
  geom_ribbon(data = YOY_BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey") +
  geom_line(data = YOY_BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = YOY_BKT.COMID.logDens.corr$real$cbar,
             linetype = "dashed") +
  theme_classic() +
  lims(x = c(0, max(YOY_BKT.COMID.logDens.corr.pred.df$x)),
       y = c(-0.25, 0.5)) +
  labs(title = "a.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=12))

YOY_BKT.COMID.logDens.synch.plot

# Does synchrony differ between the north and south regions of the SE US?
# Define these regions: we draw the line at Radford, VA in the New River Valley. Per the suggestion of Jake Rash

## YOY
## North
YOY_BKT_COMID_logDens_wide_N <- YOY_BKT_COMID_logDens_wide %>% 
  filter(Lat > 37.13)

YOY_BKT_COMID_logDens.corr_N <- Sncf(x=YOY_BKT_COMID_logDens_wide_N$Long,
                                   y=YOY_BKT_COMID_logDens_wide_N$Lat,
                                   z=YOY_BKT_COMID_logDens_wide_N[,4:ncol(YOY_BKT_COMID_logDens_wide_N)],
                                   resamp = 1000, 
                                   latlon = T,
                                   xmax = ((2/3)*max(na.omit(gcdist(YOY_BKT_COMID_logDens_wide_N$Long, YOY_BKT_COMID_logDens_wide_N$Lat)))))

plot(YOY_BKT_COMID_logDens.corr_N)
summary(YOY_BKT_COMID_logDens.corr_N)

## South
YOY_BKT_COMID_logDens_wide_S <- YOY_BKT_COMID_logDens_wide %>% 
  filter(Lat <= 37.13)

YOY_BKT_COMID_logDens.corr_S <- Sncf(x=YOY_BKT_COMID_logDens_wide_S$Long,
                                     y=YOY_BKT_COMID_logDens_wide_S$Lat,
                                     z=YOY_BKT_COMID_logDens_wide_S[,4:ncol(YOY_BKT_COMID_logDens_wide_S)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(YOY_BKT_COMID_logDens_wide_S$Long, YOY_BKT_COMID_logDens_wide_S$Lat)))))

plot(YOY_BKT_COMID_logDens.corr_S)
summary(YOY_BKT_COMID_logDens.corr_S)

######
# SPATIAL SYNCHRONY OF LOG PREDICTED ADULT BKT DENSITY BY NHDPlus v2.1 COMID

# sum Adult BKT abundance at each sample and pass and filter to just multipass data
# multipass data are necessary for predicting abundance
# this will be passed into the removal() function
Adult_BKT_sample_sums <- SE_Ind_Final %>% 
  filter(SPP == "BKT",
         TL_mm > 90, # filter here for Adult
         year(Date) %in% c(1980:2015)) %>% # Filter to time frame of interest
  group_by(SampleID,
           PassNo) %>% 
  summarise(SiteID = first(SiteID),
            Year = first(year(Date)),
            Count = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = PassNo, values_from = Count)

# Remove samples with only one pass
Adult_BKT_sample_sums$count_passes <- rowSums(apply(!is.na(Adult_BKT_sample_sums[,4:7]), MARGIN = 2, as.numeric))
Adult_BKT_sample_sums <- Adult_BKT_sample_sums %>% 
  filter(count_passes > 1) %>% 
  dplyr::select(-count_passes)

# use removal() from the "FSA" package to estimate Adult abundance from multipass data
Adult_BKT_AbundEsts <- apply(Adult_BKT_sample_sums[,4:7], MARGIN = 1, FUN = removal, 
                             method = "CarleStrub", # The CarleStrub method is most robust
                             just.ests = TRUE)
Adult_BKT_AbundEsts <- as.data.frame(t(Adult_BKT_AbundEsts))

# Join the adult abundance estimates to SampleIDs, and summarize by COMID and year
Adult_BKT_COMID_Dens <- Adult_BKT_sample_sums %>% 
  left_join(SE_Site_Final[,c("SiteID", "Lat", "Long", "Length_m", "Width_m", "COMID")], by = "SiteID") %>%  # left join in site data (we mostly just want lat/long here)
  cbind(Abund_Est = Adult_BKT_AbundEsts[,1]) %>% 
  group_by(SiteID, Year) %>% 
  dplyr::summarise(COMID = first(COMID),
                   Lat = first(Lat),
                   Long = first(Long),
                   Avg_Abund = mean(Abund_Est),
                   Avg_Density = Avg_Abund/((Length_m * Width_m)/1000)) %>%  # avg density/Site
  group_by(COMID, Year) %>% 
  dplyr::summarise(Lat = mean(Lat),
                   Long = mean(Long),
                   Avg_Density = mean(Avg_Density)) # Avg density/COMID


# Use dcast to "widen" the data by year
Adult_BKT_COMID_Dens_wide <- dcast(data = Adult_BKT_COMID_Dens, 
                                   formula = COMID + Lat + Long ~ Year,
                                   value.var = "Avg_Density")

# Make a new column with a count of the number of years of nonzero data, and filter out sites that have fewer than five years' of data during this time period
Adult_BKT_COMID_Dens_wide <- Adult_BKT_COMID_Dens_wide %>%  
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 5,
         !is.na(COMID)) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

# what percentage of observations are missing data?
mean(is.na(Adult_BKT_COMID_Dens_wide[,4:ncol(Adult_BKT_COMID_Dens_wide)]))

# make a quick map to visualize what sites we will be using
ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_point(data = Adult_BKT_COMID_Dens_wide, aes(x = Long, y = Lat)) +
  theme_classic() +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.4, -72),
            ylim = c(34.5,45))

# calculate the average density of Adult fish/year at each site
Adult_BKT_COMID_Dens_wide$mean_density <- (rowSums(Adult_BKT_COMID_Dens_wide[,4:ncol(Adult_BKT_COMID_Dens_wide)], na.rm = T)/apply(!is.na(Adult_BKT_COMID_Dens_wide[,4:ncol(Adult_BKT_COMID_Dens_wide)]), MARGIN = 1, sum))

# then replace the NAs with the average density of adult fish per year at that COMID
for (i in 1:nrow(Adult_BKT_COMID_Dens_wide)) {
  Adult_BKT_COMID_Dens_wide[i,4:ncol(Adult_BKT_COMID_Dens_wide)][is.na(Adult_BKT_COMID_Dens_wide[i,4:ncol(Adult_BKT_COMID_Dens_wide)])] <- Adult_BKT_COMID_Dens_wide[i, "mean_density"]
  print(i)
}

# remove the column where we stored mean density
Adult_BKT_COMID_Dens_wide <- Adult_BKT_COMID_Dens_wide %>% 
  select(-mean_density)

# Now make a new dataframe with log density
# calculating synchrony based on log adult density helps remove the undue influence of outliers
Adult_BKT_COMID_logDens_wide <- Adult_BKT_COMID_Dens_wide
Adult_BKT_COMID_logDens_wide[,4:ncol(Adult_BKT_COMID_logDens_wide)] <- log(Adult_BKT_COMID_logDens_wide[,4:ncol(Adult_BKT_COMID_logDens_wide)])

# Now, calculate spatial synchrony in the time series using the Sncf() function
Adult_BKT_COMID_logDens.corr <- Sncf(x=Adult_BKT_COMID_logDens_wide$Long,
                                     y=Adult_BKT_COMID_logDens_wide$Lat,
                                     z=Adult_BKT_COMID_logDens_wide[,4:ncol(Adult_BKT_COMID_logDens_wide)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(Adult_BKT_COMID_logDens_wide$Long, Adult_BKT_COMID_logDens_wide$Lat)))))

plot(Adult_BKT_COMID_logDens.corr)
summary(Adult_BKT_COMID_logDens.corr)


# Plot using ggplot for export
# create a dataframe of the predicted values
Adult_BKT.COMID.logDens.corr.pred.df <- data.frame(x = matrix(unlist(Adult_BKT_COMID_logDens.corr$real$predicted$x)),
                                                   y = matrix(unlist(Adult_BKT_COMID_logDens.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
Adult_BKT.COMID.logDens.corr.bootValues.df <- as.data.frame(t(Adult_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
Adult_BKT.COMID.logDens.corr.boot.df <- data.frame(x = matrix(unlist(Adult_BKT_COMID_logDens.corr$boot$boot.summary$predicted$x)),
                                                   ymin = Adult_BKT.COMID.logDens.corr.bootValues.df$`0.025`,
                                                   ymax = Adult_BKT.COMID.logDens.corr.bootValues.df$`0.975`)

# make the plot
Adult_BKT.COMID.logDens.synch.plot <- ggplot() +
  geom_ribbon(data = Adult_BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey") +
  geom_line(data = Adult_BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = Adult_BKT.COMID.logDens.corr$real$cbar,
             linetype = "dashed") +
  theme_classic() +
  lims(x = c(0, 1000),
       y = c(-0.25, 0.5)) +
  labs(title = "b.)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=12))

Adult_BKT.COMID.logDens.synch.plot

# Does synchrony differ between the north and south regions of the SE US?
# Define these regions: we draw the line at Radford, VA in the New River Valley. Per the suggestion of Jake Rash

## Adult
## North
Adult_BKT_COMID_logDens_wide_N <- Adult_BKT_COMID_logDens_wide %>% 
  filter(Lat > 37.13)

Adult_BKT_COMID_logDens.corr_N <- Sncf(x=Adult_BKT_COMID_logDens_wide_N$Long,
                                     y=Adult_BKT_COMID_logDens_wide_N$Lat,
                                     z=Adult_BKT_COMID_logDens_wide_N[,4:ncol(Adult_BKT_COMID_logDens_wide_N)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(Adult_BKT_COMID_logDens_wide_N$Long, Adult_BKT_COMID_logDens_wide_N$Lat)))))

plot(Adult_BKT_COMID_logDens.corr_N)
summary(Adult_BKT_COMID_logDens.corr_N)

## South
Adult_BKT_COMID_logDens_wide_S <- Adult_BKT_COMID_logDens_wide %>% 
  filter(Lat <= 37.13)

Adult_BKT_COMID_logDens.corr_S <- Sncf(x=Adult_BKT_COMID_logDens_wide_S$Long,
                                     y=Adult_BKT_COMID_logDens_wide_S$Lat,
                                     z=Adult_BKT_COMID_logDens_wide_S[,4:ncol(Adult_BKT_COMID_logDens_wide_S)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(Adult_BKT_COMID_logDens_wide_S$Long, Adult_BKT_COMID_logDens_wide_S$Lat)))))

plot(Adult_BKT_COMID_logDens.corr_S)
summary(Adult_BKT_COMID_logDens.corr_S)

#######################
## Now combine all density plots into one figure
compound_density_corrgram.data <- data.frame(AgeClass = c(rep("All", 900),
                                                              rep("YOY", 900),
                                                              rep("Adult", 900)),
                                             Subregion = rep(c(rep("All", 300),
                                                                   rep("North", 300),
                                                                   rep("South", 300)), 3),
                                             x = rbind(
                                               # All BKT, all subregions
                                               matrix(unlist(BKT_COMID_logDens.corr$real$predicted$x)),
                                               # All BKT, northern subregion
                                               matrix(unlist(BKT_COMID_logDens.corr_N$real$predicted$x)),
                                               # All BKT, southern subregion
                                               matrix(unlist(BKT_COMID_logDens.corr_S$real$predicted$x)),
                                               # YOY BKT, all subregions
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr$real$predicted$x)),
                                               # YOY BKT, northern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_N$real$predicted$x)),
                                               # YOY BKT, southern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_S$real$predicted$x)),
                                               # Adult BKT, all subregions
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr$real$predicted$x)),
                                               # Adult BKT, northern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_N$real$predicted$x)),
                                               # Adult BKT, southern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_S$real$predicted$x))),
                                             y = rbind(
                                               # All BKT, all subregions
                                               matrix(unlist(BKT_COMID_logDens.corr$real$predicted$y)),
                                               # All BKT, northern subregion
                                               matrix(unlist(BKT_COMID_logDens.corr_N$real$predicted$y)),
                                               # All BKT, southern subregion
                                               matrix(unlist(BKT_COMID_logDens.corr_S$real$predicted$y)),
                                               # YOY BKT, all subregions
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr$real$predicted$y)),
                                               # YOY BKT, northern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_N$real$predicted$y)),
                                               # YOY BKT, southern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_S$real$predicted$y)),
                                               # Adult BKT, all subregions
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr$real$predicted$y)),
                                               # Adult BKT, northern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_N$real$predicted$y)),
                                               # Adult BKT, southern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_S$real$predicted$y))),
                                             ymin = c(
                                               # All BKT, all subregions
                                               as.data.frame(t(BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,2],
                                               # All BKT, northern subregion
                                               as.data.frame(t(BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,2],
                                               # All BKT, southern subregion
                                               as.data.frame(t(BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,2],
                                               # YOY BKT, all subregions
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,2],
                                               # YOY BKT, northern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,2],
                                               # YOY BKT, southern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,2],
                                               # Adult BKT, all subregions
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,2],
                                               # Adult BKT, northern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,2],
                                               # Adult BKT, southern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,2]),
                                             ymax = c(
                                               # All BKT, all subregions
                                               as.data.frame(t(BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,10],
                                               # All BKT, northern subregion
                                               as.data.frame(t(BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,10],
                                               # All BKT, southern subregion
                                               as.data.frame(t(BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,10],
                                               # YOY BKT, all subregions
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,10],
                                               # YOY BKT, northern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,10],
                                               # YOY BKT, southern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,10],
                                               # Adult BKT, all subregions
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr$boot$boot.summary$predicted$y))[,10],
                                               # Adult BKT, northern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,10],
                                               # Adult BKT, southern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,10]),
                                             Avg_Corr = c(
                                               # All BKT, all subregions
                                               rep(BKT_COMID_logDens.corr$real$cbar, 300),
                                               # All BKT, northern subregion
                                               rep(BKT_COMID_logDens.corr_N$real$cbar, 300),
                                               # All BKT, southern subregion
                                               rep(BKT_COMID_logDens.corr_S$real$cbar, 300),
                                               # YOY BKT, all subregions
                                               rep(YOY_BKT_COMID_logDens.corr$real$cbar, 300),
                                               # YOY BKT, northern subregion
                                               rep(YOY_BKT_COMID_logDens.corr_N$real$cbar, 300),
                                               # YOY BKT, southern subregion
                                               rep(YOY_BKT_COMID_logDens.corr_S$real$cbar, 300),
                                               # Adult BKT, all subregions
                                               rep(Adult_BKT_COMID_logDens.corr$real$cbar, 300),
                                               # Adult BKT, northern subregion
                                               rep(Adult_BKT_COMID_logDens.corr_N$real$cbar, 300),
                                               # Adult BKT, southern subregion
                                               rep(Adult_BKT_COMID_logDens.corr_S$real$cbar, 300)))

# Rearrange the AgeClass category
compound_density_corrgram.data$AgeClass <- factor(compound_density_corrgram.data$AgeClass,
                                                  levels = c("All", "Adult", "YOY"))

compound_density_corrgram.plot <- ggplot(compound_density_corrgram.data) + 
  # Estimate
  geom_line(aes(x = x,
                y = y)) +
  # Envelope
  geom_ribbon(aes(x = x,
                  ymin = ymin,
                  ymax = ymax),
              alpha = 0.5) +
  # Average synchrony
  geom_hline(aes(yintercept = Avg_Corr),
             linetype = "dashed") +
  geom_hline(yintercept = 0) +
  lims(x = c(0, 450),
       y = c(-0.25, 1)) +
  theme_classic() +
  labs(x = "Distance Class (km)",
       y = "Correlation") +
  facet_grid(AgeClass ~ Subregion)
#####################
# and combine summary stats into a table
Corr_Stats <- compound_density_corrgram.data %>% 
  group_by(AgeClass,
           Subregion) %>% 
  summarise(Mean_Corr = first(Avg_Corr),
            Max_Corr = max(y))


#write.table(density_summStats, "BKT log Pred Density Summary Stats.csv")

#############################################################
# some explorations
# Compare synchrony between 2000-2018 in three subsets

corr1 <- spline.correlog(x = BKT_COMID_logDens_wide$Long,
                        y = BKT_COMID_logDens_wide$Lat,
                        z = BKT_COMID_logDens_wide[,4:9],
                        resamp = 1000,
                        latlon = T,
                        xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide$Long, BKT_COMID_logDens_wide$Lat)))))
plot(corr1)
summary(corr1)
corr1$real$cbar

corr2 <- spline.correlog(x = BKT_COMID_logDens_wide$Long,
                         y = BKT_COMID_logDens_wide$Lat,
                         z = BKT_COMID_logDens_wide[,11:16],
                         resamp = 1000,
                         latlon = T,
                         xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide$Long, BKT_COMID_logDens_wide$Lat)))))
plot(corr2)

corr3 <- spline.correlog(x = BKT_COMID_logDens_wide$Long,
                         y = BKT_COMID_logDens_wide$Lat,
                         z = BKT_COMID_logDens_wide[,17:22],
                         resamp = 1000,
                         latlon = T,
                         xmax = ((2/3)*max(na.omit(gcdist(BKT_COMID_logDens_wide$Long, BKT_COMID_logDens_wide$Lat)))))
plot(corr3)

#######################################################################
# Is there any synchrony in summer temperatures or winter flows?


## Max Estimated 0.9Q Flow
# Load flow estimates from Daren Carlisle at USGS
#flow_data <- as.data.frame(fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/SE BKT Recruitment Paper/R Files/SE BKT Recruitment Project/flow_preds_all.csv"))
flow_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/SE_COMID_flow_covars.csv")

# Filter these predictions for just the ones at our BKT COMID sites and just for the years of interest and just for winter months (december-february)
flow_data2 <- flow_data %>%
  inner_join(as.data.frame(BKT_COMID_logDens_wide[,1:3]), by = c("COMID" = "COMID"))

# widen the data to stick it in Sncf
flow_data2_wide <- dcast(data = flow_data2, 
                         formula = COMID + Lat + Long ~ Year,
                         value.var = "Max_0.9Q_WinterFlow_Scaled")

# Now, calculate spatial synchrony in the time series using the Sncf() function
Flow.COMID.corr <- Sncf(x = flow_data2_wide$Long,
                        y = flow_data2_wide$Lat,
                        z = flow_data2_wide[,4:39],
                        resamp = 1000, 
                        latlon = T,
                        xmax = ((2/3)*max(na.omit(gcdist(flow_data2_wide$Long, flow_data2_wide$Lat)))))

plot(Flow.COMID.corr, cex.lab=1.5, cex.axis=1.5)


# Plot using ggplot for export
# create a dataframe of the predicted values
Flow.COMID.corr.pred.df <- data.frame(x = matrix(unlist(Flow.COMID.corr$real$predicted$x)),
                                      y = matrix(unlist(Flow.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
Flow.COMID.corr.bootValues.df <- as.data.frame(t(Flow.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
Flow.COMID.corr.boot.df <- data.frame(x = matrix(unlist(Flow.COMID.corr$boot$boot.summary$predicted$x)),
                                      ymin = Flow.COMID.corr.bootValues.df$`0.025`,
                                      ymax = Flow.COMID.corr.bootValues.df$`0.975`)

# make the plot
Flow.COMID.synch.plot <- ggplot() +
  geom_ribbon(data = Flow.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey") +
  geom_line(data = Flow.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = Flow.COMID.corr$real$cbar,
             linetype = "dashed") +
  theme_classic() +
  lims(y = c(-0.5, 1)) +
  labs(x = "Pairwise Distance (km)",
       y = "Correlation") +
  theme(text = element_text(size=20))

# ggsave("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/WintFlow.COMID.synch.plot.jpeg",
#        width = 5, height = 3)


## Mean Estimated 0.9Q Summer Temp
# Summer temps
#load summer temp data (keep in mind this is 1980-2015)
SummTemp_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/SE_COMID_temp_covars.csv")

# Join in COMIDs
# data are already selected for summer and scaled
SummTemp_data2 <- SummTemp_data %>% 
  inner_join(BKT_COMID_logDens_wide[,1:3]) %>% 
  filter(Year %in% seq(from = 1980, to = 2015, by = 1))

# widen the data to stick it in Sncf
SummTemp_data2_wide <- dcast(data = SummTemp_data2, 
                             formula = COMID + Lat + Long ~ Year,
                             value.var = "Mean_Max_Summer_Temp_Scaled")

# Now, calculate spatial synchrony in the time series using the Sncf() function
SummTemp.COMID.corr <- Sncf(x = SummTemp_data2_wide$Long,
                            y = SummTemp_data2_wide$Lat,
                            z = SummTemp_data2_wide[,4:39],
                            resamp = 1000, 
                            latlon = T,
                            xmax = ((2/3)*max(na.omit(gcdist(SummTemp_data2_wide$Long, SummTemp_data2_wide$Lat)))))

plot(SummTemp.COMID.corr, cex.lab=1.5, cex.axis=1.5)

# Plot using ggplot for export
# create a dataframe of the predicted values
SummTemp.COMID.corr.pred.df <- data.frame(x = matrix(unlist(SummTemp.COMID.corr$real$predicted$x)),
                                          y = matrix(unlist(SummTemp.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
SummTemp.COMID.corr.bootValues.df <- as.data.frame(t(SummTemp.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
SummTemp.COMID.corr.boot.df <- data.frame(x = matrix(unlist(SummTemp.COMID.corr$boot$boot.summary$predicted$x)),
                                          ymin = SummTemp.COMID.corr.bootValues.df$`0.025`,
                                          ymax = SummTemp.COMID.corr.bootValues.df$`0.975`)

# make the plot
SummTemp_corrgram.plot <- ggplot() +
  geom_ribbon(data = SummTemp.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey") +
  geom_line(data = SummTemp.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = SummTemp.COMID.corr$real$cbar,
             linetype = "dashed") +
  theme_classic() +
  lims(y = c(-0.5, 1)) +
  labs(x = "Pairwise Distance (km)",
       y = "Correlation") +
  theme(text = element_text(size=20))

# Measured flow and temperature
NS204_daily_temps <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/Temperature Data Working/Dolloff Temperature Data/Final Files/Final Temperature Data/NS204_temps_daily.csv")
NS204_sites <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/Temperature Data Working/Dolloff Temperature Data/Final Files/Final Temperature Data/NS204_sites.csv")
  
obs_temp <- SE_temps %>% 
  left_join(NS204_sites[,c("SiteID", "COMID", "Lat", "Long")]) %>% 
  dplyr::select(COMID,
                Lat,
                Long,
                Date,
                WaterTemp_c_MAX,
                AirTemp_c_MAX) %>% 
  filter(month(Date) %in% c(6,7,8,9)) %>% # filter for just summer
  mutate(Year = year(Date),
         WaterTemp_c_MAX_scaled = scale(WaterTemp_c_MAX), # scale temperatures
         AirTemp_c_MAX_scaled = scale(AirTemp_c_MAX)) %>% 
  group_by(COMID,
           Year) %>% 
  summarise(Lat = first(Lat),
            Long = first(Long),
            Mean_Max_Summer_WaterTemp_Scaled = mean(WaterTemp_c_MAX_scaled, na.rm = T),
            Mean_Max_Summer_AirTemp_Scaled = mean(AirTemp_c_MAX_scaled, na.rm = T))

# replace NaNs with NAs
obs_temp$Mean_Max_Summer_WaterTemp_Scaled[is.nan(obs_temp$Mean_Max_Summer_WaterTemp_Scaled)] <- NA
obs_temp$Mean_Max_Summer_AirTemp_Scaled[is.nan(obs_temp$Mean_Max_Summer_AirTemp_Scaled)] <- NA

# Use dcast to "widen" the data by year
obs_waterTemp_wide <- dcast(data = obs_temp, 
                                   formula = COMID + Lat + Long ~ Year,
                                   value.var = "Mean_Max_Summer_WaterTemp_Scaled")

obs_airTemp_wide <- dcast(data = obs_temp, 
                            formula = COMID + Lat + Long ~ Year,
                            value.var = "Mean_Max_Summer_AirTemp_Scaled")

# calculate the average mean max summer temp at each site
obs_waterTemp_wide$mean_temp <- (rowSums(obs_waterTemp_wide[,4:ncol(obs_waterTemp_wide)], na.rm = T)/apply(!is.na(obs_waterTemp_wide[,4:ncol(obs_waterTemp_wide)]), MARGIN = 1, sum))
obs_airTemp_wide$mean_temp <- (rowSums(obs_airTemp_wide[,4:ncol(obs_airTemp_wide)], na.rm = T)/apply(!is.na(obs_airTemp_wide[,4:ncol(obs_airTemp_wide)]), MARGIN = 1, sum))

# then replace the NAs with the average mean max summer temp at that COMID
for (i in 1:nrow(obs_waterTemp_wide)) {
  obs_waterTemp_wide[i,4:ncol(obs_waterTemp_wide)][is.na(obs_waterTemp_wide[i,4:ncol(obs_waterTemp_wide)])] <- obs_waterTemp_wide[i, "mean_temp"]
  print(i)
}
for (i in 1:nrow(obs_airTemp_wide)) {
  obs_airTemp_wide[i,4:ncol(obs_airTemp_wide)][is.na(obs_airTemp_wide[i,4:ncol(obs_airTemp_wide)])] <- obs_airTemp_wide[i, "mean_temp"]
  print(i)
}

# remove columns with mean max summer temp
obs_waterTemp_wide <- obs_waterTemp_wide %>% 
  dplyr::select(-mean_temp)
obs_airTemp_wide <- obs_airTemp_wide %>% 
  dplyr::select(-mean_temp)

# Now, calculate spatial synchrony in the time series using the Sncf() function
SummWaterTempObs.COMID.corr <- Sncf(x = obs_waterTemp_wide$Long,
                            y = obs_waterTemp_wide$Lat,
                            z = obs_waterTemp_wide[,4:9],
                            resamp = 1000, 
                            latlon = T,
                            xmax = ((2/3)*max(na.omit(gcdist(obs_waterTemp_wide$Long, obs_waterTemp_wide$Lat)))))

plot(SummWaterTempObs.COMID.corr)

# Now, calculate spatial synchrony in the time series using the Sncf() function
SummAirTempObs.COMID.corr <- Sncf(x = obs_airTemp_wide$Long,
                                    y = obs_airTemp_wide$Lat,
                                    z = obs_airTemp_wide[,4:9],
                                    resamp = 1000, 
                                    latlon = T,
                                    xmax = ((2/3)*max(na.omit(gcdist(obs_airTemp_wide$Long, obs_airTemp_wide$Lat)))))

plot(SummAirTempObs.COMID.corr)

# Prep for ggplot for export
# water
# create a dataframe of the predicted values
SummWaterTempObs.COMID.corr.pred.df <- data.frame(x = matrix(unlist(SummWaterTempObs.COMID.corr$real$predicted$x)),
                                                y = matrix(unlist(SummWaterTempObs.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
SummWaterTempObs.COMID.corr.bootValues.df <- as.data.frame(t(SummWaterTempObs.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
SummWaterTempObs.COMID.corr.boot.df <- data.frame(x = matrix(unlist(SummWaterTempObs.COMID.corr$boot$boot.summary$predicted$x)),
                                                ymin = SummWaterTempObs.COMID.corr.bootValues.df$`0.025`,
                                                ymax = SummWaterTempObs.COMID.corr.bootValues.df$`0.975`)
# air
# create a dataframe of the predicted values
SummAirTempObs.COMID.corr.pred.df <- data.frame(x = matrix(unlist(SummAirTempObs.COMID.corr$real$predicted$x)),
                                          y = matrix(unlist(SummAirTempObs.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
SummAirTempObs.COMID.corr.bootValues.df <- as.data.frame(t(SummAirTempObs.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
SummAirTempObs.COMID.corr.boot.df <- data.frame(x = matrix(unlist(SummAirTempObs.COMID.corr$boot$boot.summary$predicted$x)),
                                          ymin = SummAirTempObs.COMID.corr.bootValues.df$`0.025`,
                                          ymax = SummAirTempObs.COMID.corr.bootValues.df$`0.975`)

########################################
# Create correlogams that combine temperature, flow, AND log BKT density

# uses 4 colors from rcolorbrewer's "dark2" pallete
colors1 <- c("Mean Max Observed\n Summer Water Temp" = "#7570b3",
            "Mean Max Observed\n Summer Air Temp" = "#d95f02",
            "Max 0.9Q Winter Flow" = "#1b9e77",
            "BKT Log Density" = "#e7298a")

compound_corrgram.plot <- ggplot() +
  # log BKT density (both age classes)
  geom_ribbon(data = BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "BKT Log Density"),
              alpha = 0.5) +
  geom_line(data = BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = BKT.COMID.logDens.corr$real$cbar,
             linetype = "dashed") +
  # Winter Flow
  geom_ribbon(data = Flow.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Max 0.9Q Winter Flow"),
              alpha = 0.5) +
  geom_line(data = Flow.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = Flow.COMID.corr$real$cbar,
             linetype = "dashed") +
  # Summer Water Temp
  geom_ribbon(data = SummWaterTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Mean Max Observed\n Summer Water Temp"),
              alpha = 0.5) +
  geom_line(data = SummWaterTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummWaterTempObs.COMID.corr.pred.df$real$cbar,
             linetype = "dashed") +
  # Summer Air Temp
  geom_ribbon(data = SummAirTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Mean Max Observed\n Summer Air Temp"),
              alpha = 0.5) +
  geom_line(data = SummAirTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummAirTempObs.COMID.corr.pred.df$real$cbar,
             linetype = "dashed") +
  # Formatting
  geom_hline(yintercept = 0) +
  theme_classic() +
  lims(y = c(-0.5, 1)) +
  labs(x = "Pairwise Distance (km)",
       y = "Correlation",
       fill = "Legend") +
  scale_fill_manual(values = colors1)


colors2 <- c("Mean Max Daily Observed\nSummer Water Temp" = "#7570b3",
             "Mean Max Daily Observed\nSummer Air Temp" = "#d95f02",
             "Max 0.9Q Winter Flow" = "#1b9e77",
             "BKT Log YOY Density" = "#e7298a")

compound_corrgram_YOY.plot <- ggplot() +
  # log YOY density
  geom_ribbon(data = YOY_BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "BKT Log YOY Density"),
              alpha = 0.5) +
  geom_line(data = YOY_BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = YOY_BKT.COMID.logDens.corr$real$cbar,
             linetype = "dashed") +
  # Winter Flow
  geom_ribbon(data = Flow.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Max 0.9Q Winter Flow"),
              alpha = 0.5) +
  geom_line(data = Flow.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = Flow.COMID.corr$real$cbar,
             linetype = "dashed") +
  # Summer Water Temp
  geom_ribbon(data = SummWaterTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Mean Max Daily Observed\nSummer Water Temp"),
              alpha = 0.5) +
  geom_line(data = SummWaterTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummWaterTempObs.COMID.corr.pred.df$real$cbar,
             linetype = "dashed") +
  # Summer Air Temp
  geom_ribbon(data = SummAirTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Mean Max Daily Observed\nSummer Air Temp"),
              alpha = 0.5) +
  geom_line(data = SummAirTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummAirTempObs.COMID.corr.pred.df$real$cbar,
             linetype = "dashed") +
  # Formatting
  geom_hline(yintercept = 0) +
  theme_classic() +
  lims(y = c(-0.5, 1)) +
  labs(x = "Pairwise Distance (km)",
       y = "Correlation",
       fill = "Legend") +
  scale_fill_manual(values = colors2)
