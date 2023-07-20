# George Valentine
# BKT Spatial Synchrony Project
# Synchrony Analysis with Correlograms v2

# This version uses abundance estimates from the N-mixture models as outputs


# Install packages
library(ncf)

sample_areas_long <- pivot_longer(sample_areas_wide_CE,
             cols = 1:34,
             names_to = "Year",
             values_to = "Area_m2") %>% 
  mutate(Year = as.numeric(Year)) %>% 
  arrange(COMID, Year)

## YOY
N.YOY_climateEffects.table <- N.YOY_climateEffects %>%
  rownames_to_column(., "param") %>%
  mutate(SegmentNo = as.numeric(str_extract(param, "[:digit:]+")),
         YearNo = as.numeric(str_extract(param, "(?<=,)[:digit:]+")),
         Year = YearNo + 1981) %>% # calculate the year of sampling
  left_join(cbind(COMID = p1_YOY_CE$COMID, Agency = p1_YOY_CE$Agency, data.frame(SegmentNo = seq(1,144)))) %>% 
  arrange(COMID, Year) %>% 
  cbind(sample_areas_long[,3]) %>% 
  mutate(logDens = ifelse(func > 0, 
                          log(func/(Area_m2/1000)),
                          0))
          

N.YOY_climateEffects_wide <- N.YOY_climateEffects.table %>% 
  pivot_wider(names_from = Year,
              values_from = logDens,
              id_cols = c(COMID, Agency)) %>% 
  cbind(segment_data_CE[,c("Lat", "Long")]) %>% 
  relocate(Lat, Long, .after = COMID) %>% 
  dplyr::select(-Agency)


YOY_BKT_COMID_logDens.corr2 <- Sncf(x=N.YOY_climateEffects_wide$Long,
                                   y=N.YOY_climateEffects_wide$Lat,
                                   z=N.YOY_climateEffects_wide[,4:ncol(N.YOY_climateEffects_wide)],
                                   resamp = 1000, 
                                   latlon = T,
                                   xmax = ((2/3)*max(na.omit(gcdist(N.YOY_climateEffects_wide$Long, N.YOY_climateEffects_wide$Lat)))))

plot(YOY_BKT_COMID_logDens.corr2)
YOY_corrgram_summary.table <- summary(YOY_BKT_COMID_logDens.corr2)

# Plot using ggplot for export
# create a dataframe of the predicted values
YOY_BKT.COMID.logDens.corr.pred.df <- data.frame(x = matrix(unlist(YOY_BKT_COMID_logDens.corr2$real$predicted$x)),
                                                 y = matrix(unlist(YOY_BKT_COMID_logDens.corr2$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
YOY_BKT.COMID.logDens.corr.bootValues.df <- as.data.frame(t(YOY_BKT_COMID_logDens.corr2$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
YOY_BKT.COMID.logDens.corr.boot.df <- data.frame(x = matrix(unlist(YOY_BKT_COMID_logDens.corr2$boot$boot.summary$predicted$x)),
                                                 ymin = YOY_BKT.COMID.logDens.corr.bootValues.df$`0.025`,
                                                 ymax = YOY_BKT.COMID.logDens.corr.bootValues.df$`0.975`)

# Does synchrony differ between the north and south regions of the SE US?
# Define these regions: we draw the line at Radford, VA in the New River Valley. Per the suggestion of Jake Rash

## YOY
## North
N.YOY_climateEffects_wide_N <- N.YOY_climateEffects_wide %>% 
  filter(Lat > 37.13)

YOY_BKT_COMID_logDens.corr_N <- Sncf(x=N.YOY_climateEffects_wide_N$Long,
                                     y=N.YOY_climateEffects_wide_N$Lat,
                                     z=N.YOY_climateEffects_wide_N[,4:ncol(N.YOY_climateEffects_wide_N)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(N.YOY_climateEffects_wide_N$Long, N.YOY_climateEffects_wide_N$Lat)))))

plot(YOY_BKT_COMID_logDens.corr_N)
summary(YOY_BKT_COMID_logDens.corr_N)

## South
N.YOY_climateEffects_wide_S <- N.YOY_climateEffects_wide %>% 
  filter(Lat <= 37.13)

YOY_BKT_COMID_logDens.corr_S <- Sncf(x=N.YOY_climateEffects_wide_S$Long,
                                     y=N.YOY_climateEffects_wide_S$Lat,
                                     z=N.YOY_climateEffects_wide_S[,4:ncol(N.YOY_climateEffects_wide_S)],
                                     resamp = 1000, 
                                     latlon = T,
                                     xmax = ((2/3)*max(na.omit(gcdist(N.YOY_climateEffects_wide_S$Long, N.YOY_climateEffects_wide_S$Lat)))))

plot(YOY_BKT_COMID_logDens.corr_S)
summary(YOY_BKT_COMID_logDens.corr_S)

## Adults
N.Adult_climateEffects.table <- N.adult_climateEffects %>%
  rownames_to_column(., "param") %>%
  mutate(SegmentNo = as.numeric(str_extract(param, "[:digit:]+")),
         YearNo = as.numeric(str_extract(param, "(?<=,)[:digit:]+")),
         Year = YearNo + 1981) %>% # calculate the year of sampling
  left_join(cbind(COMID = p1_adult_CE$COMID, Agency = p1_adult_CE$Agency, data.frame(SegmentNo = seq(1,144)))) %>% 
  arrange(COMID, Year) %>% 
  cbind(sample_areas_long[,3]) %>% 
  mutate(logDens = ifelse(func > 0, 
                          log(func/(Area_m2/1000)),
                          0))


N.Adult_climateEffects_wide <- N.Adult_climateEffects.table %>% 
  pivot_wider(names_from = Year,
              values_from = logDens,
              id_cols = c(COMID, Agency)) %>% 
  cbind(segment_data_CE[,c("Lat", "Long")]) %>% 
  relocate(Lat, Long, .after = COMID) %>% 
  dplyr::select(-Agency)


Adult_BKT_COMID_logDens.corr2 <- Sncf(x=N.Adult_climateEffects_wide$Long,
                                    y=N.Adult_climateEffects_wide$Lat,
                                    z=N.Adult_climateEffects_wide[,4:ncol(N.Adult_climateEffects_wide)],
                                    resamp = 1000, 
                                    latlon = T,
                                    xmax = ((2/3)*max(na.omit(gcdist(N.Adult_climateEffects_wide$Long, N.Adult_climateEffects_wide$Lat)))))

plot(Adult_BKT_COMID_logDens.corr2)
Adult_corrgram_summary.table <- summary(Adult_BKT_COMID_logDens.corr2)


# Plot using ggplot for export
# create a dataframe of the predicted values
Adult_BKT.COMID.logDens.corr.pred.df <- data.frame(x = matrix(unlist(Adult_BKT_COMID_logDens.corr2$real$predicted$x)),
                                                   y = matrix(unlist(Adult_BKT_COMID_logDens.corr2$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
Adult_BKT.COMID.logDens.corr.bootValues.df <- as.data.frame(t(Adult_BKT_COMID_logDens.corr2$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
Adult_BKT.COMID.logDens.corr.boot.df <- data.frame(x = matrix(unlist(Adult_BKT_COMID_logDens.corr2$boot$boot.summary$predicted$x)),
                                                   ymin = Adult_BKT.COMID.logDens.corr.bootValues.df$`0.025`,
                                                   ymax = Adult_BKT.COMID.logDens.corr.bootValues.df$`0.975`)

## Adult
## North
N.Adult_climateEffects_wide_N <- N.Adult_climateEffects_wide %>% 
  filter(Lat > 37.13)

Adult_BKT_COMID_logDens.corr_N <- Sncf(x=N.Adult_climateEffects_wide_N$Long,
                                       y=N.Adult_climateEffects_wide_N$Lat,
                                       z=N.Adult_climateEffects_wide_N[,4:ncol(N.Adult_climateEffects_wide_N)],
                                       resamp = 1000, 
                                       latlon = T,
                                       xmax = ((2/3)*max(na.omit(gcdist(N.Adult_climateEffects_wide_N$Long, N.Adult_climateEffects_wide_N$Lat)))))

plot(Adult_BKT_COMID_logDens.corr_N)
summary(Adult_BKT_COMID_logDens.corr_N)

## South
N.Adult_climateEffects_wide_S <- N.Adult_climateEffects_wide %>% 
  filter(Lat <= 37.13)

Adult_BKT_COMID_logDens.corr_S <- Sncf(x=N.Adult_climateEffects_wide_S$Long,
                                       y=N.Adult_climateEffects_wide_S$Lat,
                                       z=N.Adult_climateEffects_wide_S[,4:ncol(N.Adult_climateEffects_wide_S)],
                                       resamp = 1000, 
                                       latlon = T,
                                       xmax = ((2/3)*max(na.omit(gcdist(N.Adult_climateEffects_wide_S$Long, N.Adult_climateEffects_wide_S$Lat)))))

plot(Adult_BKT_COMID_logDens.corr_S)
summary(Adult_BKT_COMID_logDens.corr_S)

#######################
## Now combine all density plots into one figure
compound_density_corrgram.data <- data.frame(AgeClass = c(rep("YOY", 600),
                                                          rep("Adult", 600)),
                                             Subregion = rep(c(rep("North", 300),
                                                               rep("South", 300)), 2),
                                             x = rbind(
                                               # YOY BKT, northern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_N$real$predicted$x)),
                                               # YOY BKT, southern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_S$real$predicted$x)),
                                               # Adult BKT, northern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_N$real$predicted$x)),
                                               # Adult BKT, southern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_S$real$predicted$x))),
                                             y = rbind(
                                               # YOY BKT, northern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_N$real$predicted$y)),
                                               # YOY BKT, southern subregion
                                               matrix(unlist(YOY_BKT_COMID_logDens.corr_S$real$predicted$y)),
                                               # Adult BKT, northern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_N$real$predicted$y)),
                                               # Adult BKT, southern subregion
                                               matrix(unlist(Adult_BKT_COMID_logDens.corr_S$real$predicted$y))),
                                             ymin = c(
                                               # YOY BKT, northern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,2],
                                               # YOY BKT, southern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,2],
                                               # Adult BKT, northern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,2],
                                               # Adult BKT, southern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,2]),
                                             ymax = c(
                                               # YOY BKT, northern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,10],
                                               # YOY BKT, southern subregion
                                               as.data.frame(t(YOY_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,10],
                                               # Adult BKT, northern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_N$boot$boot.summary$predicted$y))[,10],
                                               # Adult BKT, southern subregion
                                               as.data.frame(t(Adult_BKT_COMID_logDens.corr_S$boot$boot.summary$predicted$y))[,10]),
                                             Avg_Corr = c(
                                               # YOY BKT, northern subregion
                                               rep(YOY_BKT_COMID_logDens.corr_N$real$cbar, 300),
                                               # YOY BKT, southern subregion
                                               rep(YOY_BKT_COMID_logDens.corr_S$real$cbar, 300),
                                               # Adult BKT, northern subregion
                                               rep(Adult_BKT_COMID_logDens.corr_N$real$cbar, 300),
                                               # Adult BKT, southern subregion
                                               rep(Adult_BKT_COMID_logDens.corr_S$real$cbar, 300)))

# Rearrange the AgeClass category
compound_density_corrgram.data$AgeClass <- factor(compound_density_corrgram.data$AgeClass,
                                                  levels = c( "Adult", "YOY"))

compound_corrgram_BKT.plot <- ggplot(compound_density_corrgram.data) + 
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
  lims(x = c(0, 250),
       y = c(-0.25, 0.5)) +
  theme_classic() +
  labs(x = "Distance Class (km)",
       y = "Correlation") +
  facet_grid(AgeClass ~ Subregion)
#################################################
## Climate Covariates

N.YOY_climateEffects_wide_unique <- N.YOY_climateEffects_wide %>% 
  distinct(COMID, .keep_all = T)

## Max Estimated 0.9Q Flow
# Load flow estimates from Daren Carlisle at USGS
#flow_data <- as.data.frame(fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/SE BKT Recruitment Paper/R Files/SE BKT Recruitment Project/flow_preds_all.csv"))
flow_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/SE_Trout/Data/SE_Trout_COMID_flow_covars.csv")

# Filter these predictions for just the ones at our BKT COMID sites and just for the years of interest and just for winter months (december-february)
# also center and scale the flow
flow_data2 <- flow_data %>%
  filter(Year %in% seq(1982,2015)) %>% 
  right_join(N.YOY_climateEffects_wide_unique[,1:3], by = c("COMID" = "COMID")) 


# widen the data to stick it in Sncf
flow_data2_wide <- dcast(data = flow_data2, 
                         formula = COMID + Lat + Long ~ Year,
                         #value.var = "Max_0.9Q_WinterFlow_Scaled")
                         value.var = "Mean_EstQ_WinterFlow")

# Now, calculate spatial synchrony in the time series using the Sncf() function
Flow.COMID.corr <- Sncf(x = flow_data2_wide$Long,
                        y = flow_data2_wide$Lat,
                        z = flow_data2_wide[,4:37],
                        resamp = 1000, 
                        latlon = T,
                        xmax = ((2/3)*max(na.omit(gcdist(flow_data2_wide$Long, flow_data2_wide$Lat)))))

plot(Flow.COMID.corr, cex.lab=1.5, cex.axis=1.5)
FlowEst_corrgram_summary.table <- summary(Flow.COMID.corr)

# create a dataframe of the predicted values
Flow.COMID.corr.pred.df <- data.frame(x = matrix(unlist(Flow.COMID.corr$real$predicted$x)),
                                      y = matrix(unlist(Flow.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
Flow.COMID.corr.bootValues.df <- as.data.frame(t(Flow.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
Flow.COMID.corr.boot.df <- data.frame(x = matrix(unlist(Flow.COMID.corr$boot$boot.summary$predicted$x)),
                                      ymin = Flow.COMID.corr.bootValues.df$`0.025`,
                                      ymax = Flow.COMID.corr.bootValues.df$`0.975`)


## Mean Estimated Daily Max Summer Temp
# Summer temps
#load summer temp data (keep in mind this is 1980-2015)
SummTemp_data <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/SE_COMID_temp_covars.csv")

# Join in COMIDs
# data are already selected for summer and scaled
SummTemp_data2 <- SummTemp_data %>% 
  right_join(N.YOY_climateEffects_wide_unique[,1:3]) %>% 
  filter(Year %in% seq(from = 1982, to = 2015, by = 1))

# widen the data to stick it in Sncf
SummTemp_data2_wide <- dcast(data = SummTemp_data2, 
                             formula = COMID + Lat + Long ~ Year,
                             value.var = "Mean_Max_Summer_Temp_Scaled")

# Now, calculate spatial synchrony in the time series using the Sncf() function
SummTemp.COMID.corr <- Sncf(x = SummTemp_data2_wide$Long,
                            y = SummTemp_data2_wide$Lat,
                            z = SummTemp_data2_wide[,4:37],
                            resamp = 1000, 
                            latlon = T,
                            xmax = ((2/3)*max(na.omit(gcdist(SummTemp_data2_wide$Long, SummTemp_data2_wide$Lat)))))

plot(SummTemp.COMID.corr, cex.lab=1.5, cex.axis=1.5)
SummTempEst_corrgram_summary.table <- data.frame(avg_corr = SummTemp.COMID.corr$real$cbar,
                                                 corr_length = SummTemp.COMID.corr$real$x.intercept)

# create a dataframe of the predicted values
SummTemp.COMID.corr.pred.df <- data.frame(x = matrix(unlist(SummTemp.COMID.corr$real$predicted$x)),
                                          y = matrix(unlist(SummTemp.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
SummTemp.COMID.corr.bootValues.df <- as.data.frame(t(SummTemp.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
SummTemp.COMID.corr.boot.df <- data.frame(x = matrix(unlist(SummTemp.COMID.corr$boot$boot.summary$predicted$x)),
                                          ymin = SummTemp.COMID.corr.bootValues.df$`0.025`,
                                          ymax = SummTemp.COMID.corr.bootValues.df$`0.975`)


## Measured precipitation ##
# hourly measured precip from 1/1/2008 to 12/31/2013 for roughly the same geographic extent as the trout sites
# downloaded 11/30/22 from https://www.ncei.noaa.gov/maps/hourly/?layers=001
# data info available in "Data" folder in "PRECIP_HLY_documentation.pdf"

# import file of observed precip (100ths of a mm)
SE_precip_hourly <- fread(here("Data", "SE_precip.csv"))

# set missing data to NA and filter out low elevation sites (avg trout site elevation is ~630m)
# warning about "NAs introduced by coercion" is okay here because lat,long columns contain the word "unknown", which cannot be converted to numeric
SE_precip_hourly[SE_precip_hourly == 25399.75] <- NA # remove the NA value

SE_precip_hourly <- SE_precip_hourly %>% 
  #group_by(STATION) %>% 
  mutate(Date = ymd_hm(DATE),
         Elev_m = as.numeric(ELEVATION),
         Lat = as.numeric(LATITUDE),
         Long = as.numeric(LONGITUDE)) %>% 
  filter(#Elev_m > 400,
    !is.na(ELEVATION))

# summarise to get total precip for each winter
# reuse code from "SE_Trout_covariate_formatting.R" to make a second year column to use in summarizing winter flow
SE_precip_winter <- SE_precip_hourly %>% 
  mutate(Month = month(Date),
         Year = year(Date),
         Year2 = ifelse(Month == 12, NA, Year)) %>%
  arrange(Year, Month) %>% 
  group_by(STATION) %>% 
  tidyr::fill(Year2, .direction = "up") %>% 
  mutate(Year2 = ifelse(Year == 2013 & Month == 12, NA, Year2)) %>% 
  filter(!is.na(Year2),
         !is.na(HPCP),
         Month %in% c(12,1,2)) %>% 
  group_by(STATION,
           Year2) %>% 
  summarise(Total_Winter_Precip = sum(HPCP),
            Lat = first(Lat),
            Long = first(Long))
# mutate(Total_Winter_Precip_Scaled = c(scale(Total_Winter_Precip)))

# Use dcast to "widen" the data by year
SE_precip_winter_wide <- dcast(data = SE_precip_winter, 
                               formula = STATION + Lat + Long ~ Year2,
                               value.var = "Total_Winter_Precip")

# plot sites
precip_sites_map.plot <- ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  geom_point(data = SE_precip_winter_wide, aes(x = Long, y = Lat)) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude") + 
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-85, -76),
            ylim = c(34.5, 40))

# Make a new column with a count of the number of years of nonzero data, and filter out sites that have fewer than five years' of data during this time period
SE_precip_winter_wide <- SE_precip_winter_wide %>%  
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 3) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

# calculate the average max 90th percentile precip at each site
SE_precip_winter_wide$mean_precip <- (rowSums(SE_precip_winter_wide[,4:ncol(SE_precip_winter_wide)], na.rm = T)/apply(!is.na(SE_precip_winter_wide[,4:ncol(SE_precip_winter_wide)]), MARGIN = 1, sum))

# then replace the NAs with the average max 90th percentile precip at that site
for (i in 1:nrow(SE_precip_winter_wide)) {
  SE_precip_winter_wide[i,4:ncol(SE_precip_winter_wide)][is.na(SE_precip_winter_wide[i,4:ncol(SE_precip_winter_wide)])] <- SE_precip_winter_wide[i, "mean_precip"]
  print(i)
}

# remove column with  average max 90th percentile precip
SE_precip_winter_wide <- SE_precip_winter_wide %>% 
  dplyr::select(-mean_precip)

# Now, calculate spatial synchrony in the time series using the Sncf() function
WinterPrecipObs.COMID.corr <- Sncf(x = SE_precip_winter_wide$Long,
                                   y = SE_precip_winter_wide$Lat,
                                   z = SE_precip_winter_wide[,4:9],
                                   resamp = 1000, 
                                   latlon = T,
                                   xmax = ((2/3)*max(na.omit(gcdist(SE_precip_winter_wide$Long, SE_precip_winter_wide$Lat)))))

plot(WinterPrecipObs.COMID.corr)
WinterPrecipObs_corrgram_summary.table <- data.frame(avg_corr = WinterPrecipObs.COMID.corr$real$cbar,
                                                     corr_length = WinterPrecipObs.COMID.corr$real$x.intercept)

# Prep for ggplot for export
# create a dataframe of the predicted values
WinterPrecipObs.COMID.corr.pred.df <- data.frame(x = matrix(unlist(WinterPrecipObs.COMID.corr$real$predicted$x)),
                                                 y = matrix(unlist(WinterPrecipObs.COMID.corr$real$predicted$y)))

# ...and a dataframe of the y values for the bootstrap prediction
WinterPrecipObs.COMID.corr.bootValues.df <- as.data.frame(t(WinterPrecipObs.COMID.corr$boot$boot.summary$predicted$y))

# merge x and y values (min and max) for the bootstrap confidence intervals into a dataframe
WinterPrecipObs.COMID.corr.boot.df <- data.frame(x = matrix(unlist(WinterPrecipObs.COMID.corr$boot$boot.summary$predicted$x)),
                                                 ymin = WinterPrecipObs.COMID.corr.bootValues.df$`0.025`,
                                                 ymax = WinterPrecipObs.COMID.corr.bootValues.df$`0.975`)


## Measured temperature ##
NS204_daily_temps <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/Temperature Data Working/Dolloff Temperature Data/Final Files/Final Temperature Data/NS204_temps_daily.csv")
NS204_sites <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Temperature/Temperature Data Working/Dolloff Temperature Data/Final Files/Final Temperature Data/NS204_sites.csv")


# are daily air and water temperatures correlated?
library(corrplot)
temp_corr <- NS204_daily_temps %>% 
  select(AirTemp_c_MEAN,
         WaterTemp_c_MEAN)
daily_airwater_corr.val <- cor(temp_corr, method="pearson", use="pairwise.complete.obs")[2,1]

obs_temp <- NS204_daily_temps %>% 
  left_join(NS204_sites[,c("SiteID", "COMID", "Lat", "Long")]) %>% 
  dplyr::select(COMID,
                Lat,
                Long,
                Date,
                WaterTemp_c_MEAN,
                AirTemp_c_MEAN) %>% 
  filter(month(Date) %in% c(6,7,8,9)) %>% # filter for just summer
  mutate(Year = year(Date)) %>% 
  group_by(COMID,
           Year) %>% 
  summarise(Lat = first(Lat),
            Long = first(Long),
            Mean_Summer_WaterTemp = mean(WaterTemp_c_MEAN, na.rm = T),
            Mean_Summer_AirTemp = mean(AirTemp_c_MEAN, na.rm = T))
# mutate(Mean_Summer_WaterTemp_Scaled = c(scale(Mean_Summer_WaterTemp)),
#        Mean_Summer_AirTemp_Scaled = c(scale(Mean_Summer_AirTemp)))

# replace NaNs with NAs
obs_temp$Mean_Summer_WaterTemp[is.nan(obs_temp$Mean_Summer_WaterTemp)] <- NA
obs_temp$Mean_Summer_AirTemp[is.nan(obs_temp$Mean_Summer_AirTemp)] <- NA

# Use dcast to "widen" the data by year
obs_waterTemp_wide <- dcast(data = obs_temp, 
                            formula = COMID + Lat + Long ~ Year,
                            value.var = "Mean_Summer_WaterTemp")

obs_airTemp_wide <- dcast(data = obs_temp, 
                          formula = COMID + Lat + Long ~ Year,
                          value.var = "Mean_Summer_AirTemp")

# Make a new column with a count of the number of years of nonzero data, and filter out sites that have fewer than five years' of data during this time period
obs_waterTemp_wide <- obs_waterTemp_wide %>%  
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 3) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

obs_airTemp_wide <- obs_airTemp_wide %>%  
  mutate(nYears_data = rowSums(.[,4:ncol(.)] != 0, na.rm = T)) %>% 
  filter(nYears_data >= 3) %>% 
  dplyr::select(-nYears_data) # remove the column because we don't need it anymore

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
SummWaterTempObs_corrgram_summary.table <- summary(SummWaterTempObs.COMID.corr)

# Now, calculate spatial synchrony in the time series using the Sncf() function
SummAirTempObs.COMID.corr <- Sncf(x = obs_airTemp_wide$Long,
                                  y = obs_airTemp_wide$Lat,
                                  z = obs_airTemp_wide[,4:9],
                                  resamp = 1000, 
                                  latlon = T,
                                  xmax = ((2/3)*max(na.omit(gcdist(obs_airTemp_wide$Long, obs_airTemp_wide$Lat)))))

plot(SummAirTempObs.COMID.corr)
SummAirTempObs_corrgram_summary.table <- summary(SummAirTempObs.COMID.corr)

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

#################################################
## Plotting
colors2 <- c("Winter Flow" = "#1b9e77",
             "Summer Air Temp" = "#d95f02",
             "Summer Water Temp" = "#e6ab02",
             "Winter Precip" = "#7570b3",
             "YOY Density" = "#e7298a",
             "Adult Density" = "#66a61e")

compound_corrgram2.plot <- ggplot() +
  # Winter Flow
  geom_ribbon(data = Flow.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Winter Flow"),
              alpha = 0.5) +
  geom_line(data = Flow.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = Flow.COMID.corr$real$cbar,
             linetype = "dashed",
             color = "#1b9e77",
             alpha = 0.5) +
  # Winter precip
  geom_ribbon(data = WinterPrecipObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Winter Precip"),
              alpha = 0.5) +
  geom_line(data = WinterPrecipObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = WinterPrecipObs.COMID.corr$real$cbar,
             linetype = "dashed",
             color = "#7570b3",
             alpha = 0.5) +
  # Summer Water Temp
  geom_ribbon(data = SummWaterTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Summer Water Temp"),
              alpha = 0.5) +
  geom_line(data = SummWaterTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummWaterTempObs.COMID.corr$real$cbar,
             linetype = "dashed",
             color = "#e6ab02",
             alpha = 0.5) +
  # Summer Air Temp
  geom_ribbon(data = SummAirTempObs.COMID.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Summer Air Temp"),
              alpha = 0.5) +
  geom_line(data = SummAirTempObs.COMID.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = SummAirTempObs.COMID.corr$real$cbar,
             linetype = "dashed",
             color = "#d95f02",
             alpha = 0.5) +
  # log YOY density
  geom_ribbon(data = YOY_BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "YOY Density"),
              alpha = 0.5) +
  geom_line(data = YOY_BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = YOY_BKT_COMID_logDens.corr2$real$cbar,
             linetype = "dashed",
             color = "#e7298a",
             alpha = 0.5) +
  # log adult density
  geom_ribbon(data = Adult_BKT.COMID.logDens.corr.boot.df,
              aes(x = x, ymin = ymin, ymax = ymax, fill = "Adult Density"),
              alpha = 0.5) +
  geom_line(data = Adult_BKT.COMID.logDens.corr.pred.df,
            aes(x = x, y = y)) +
  geom_hline(yintercept = Adult_BKT_COMID_logDens.corr2$real$cbar,
             linetype = "dashed",
             color = "#66a61e",
             alpha = 0.5) +
  # Formatting
  geom_hline(yintercept = 0) +
  theme_classic() +
  lims(y = c(-0.25, 1)) +
  labs(x = "Pairwise Distance (km)",
       y = "Correlation",
       fill = "Legend") +
  scale_fill_manual(values = colors2) +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.5)))

########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here("results", "v3.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)