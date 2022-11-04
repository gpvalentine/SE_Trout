## Variogram analysis of Bayesian N-mixture model parameters
## George Valentine
## Oct. 2022

## Load packages
library(tidyverse)
library(data.table)
library(gstat)      # for making variograms
library(geoR)       # for making variograms
library(sp)

## Load data

# Model parameters summary
v13_YOY_full_params <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Spatial Synchrony in Trout Project/R Files/Spatial Synchrony in Trout/Bayesian synchrony model v13 outputs/v13_YOY_full_params.csv")

# COMIDs from model fitting
COMID_data

# Site data
SE_Site_Final

## ICC Values
# Join site data to ICC values
YOY_ICCs <- v13_YOY_full_params %>% 
  rownames_to_column(., "param") %>% 
  .[858:1026,] %>% 
  cbind(COMID_data[,c(1,3:4)])

YOY_ICCs.sp <- YOY_ICCs
coordinates(YOY_ICCs.sp) <- ~Long + Lat
proj4string(YOY_ICCs.sp) = "+proj=longlat"

# create variogram
YOY_ICC_variogram <- gstat::variogram(data = YOY_ICCs.sp,
                                      object = mean~1)

#plot
ggplot(YOY_ICC_variogram) +
  geom_point(aes(x = dist,
                 y = gamma)) +
  labs(x = "Distance Class (km)",
       y = "Variance",
       title = "Variogram of YOY ICCs") +
  theme_classic()

# Semivariogram
# There will be a warning about co-located data. That's okay because we have multiple ICCs at the same sites (different sources). Same goes for some of the other measurements.
YOY_ICCs.geo <- as.geodata(YOY_ICCs[,c(14:15, 2)])
dup.coords(YOY_ICCs.geo)
YOY_ICC_semivariogram <- geoR::variog(YOY_ICCs.geo)
plot(YOY_ICC_semivariogram,
     main = "Semivariogram of YOY ICCs")


## Beta values: Mean Max Summer Temp
# Join site data to ICC values
YOY_Betas_SummTemp <- v13_YOY_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov") & str_detect(param, "1,")) %>% # this filters out the beta[1,] parameters
  cbind(COMID_data[,c(1,3:4)])

YOY_Betas_SummTemp.sp <- YOY_Betas_SummTemp
coordinates(YOY_Betas_SummTemp.sp) <- ~Long + Lat
proj4string(YOY_Betas_SummTemp.sp) = "+proj=longlat"


# create variogram
YOY_Betas_SummTemp_variogram <- gstat::variogram(data = YOY_Betas_SummTemp.sp,
                                      object = mean~1)

#plot
ggplot(YOY_Betas_SummTemp_variogram) +
  geom_point(aes(x = dist,
                 y = gamma)) +
  labs(x = "Distance Class (km)",
       y = "Variance",
       title = "Variogram of YOY SummTemp Betas") +
  theme_classic()

# Semivariogram
YOY_Betas_SummTemp.geo <- as.geodata(YOY_Betas_SummTemp[,c(14:15, 2)])
YOY_Betas_SummTemp_semivariogram <- geoR::variog(YOY_Betas_SummTemp.geo)
plot(YOY_Betas_SummTemp_semivariogram,
     main = "Semivariogram of YOY SummTemp Betas")

## Beta values: Mean Max Winter Flow
# Join site data to ICC values
YOY_Betas_WintFlow <- v13_YOY_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov") & str_detect(param, "2,")) %>% # this filters out the beta[1,] parameters
  cbind(COMID_data[,c(1,3:4)])

YOY_Betas_WintFlow.sp <- YOY_Betas_WintFlow
coordinates(YOY_Betas_WintFlow.sp) <- ~Long + Lat
proj4string(YOY_Betas_WintFlow.sp) = "+proj=longlat"

# create variogram
YOY_Betas_WintFlow_variogram <- gstat::variogram(data = YOY_Betas_WintFlow.sp,
                                                 object = mean~1)

#plot
ggplot(YOY_Betas_WintFlow_variogram) +
  geom_point(aes(x = dist,
                 y = gamma)) +
  labs(x = "Distance Class (km)",
       y = "Variance",
       title = "Variogram of YOY WintFlow Betas") +
  theme_classic()

# Semivariogram
YOY_Betas_WintFlow.geo <- as.geodata(YOY_Betas_WintFlow[,c(14:15, 2)])
dup.coords(YOY_Betas_WintFlow.geo)
YOY_Betas_WintFlow_semivariogram <- geoR::variog(YOY_Betas_WintFlow.geo)
plot(YOY_Betas_WintFlow_semivariogram,
     main = "Semivariogram of YOY WintFlow Betas")

## Beta values: Mean Max Spring Flow
# Join site data to ICC values
YOY_Betas_SprFlow <- v13_YOY_full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta.cov") & str_detect(param, "3,")) %>% # this filters out the beta[1,] parameters
  cbind(COMID_data[,c(1,3:4)])

YOY_Betas_SprFlow.sp <- YOY_Betas_SprFlow
coordinates(YOY_Betas_SprFlow.sp) <- ~Long + Lat
proj4string(YOY_Betas_SprFlow.sp) = "+proj=longlat"

# create variogram
YOY_Betas_SprFlow_variogram <- gstat::variogram(data = YOY_Betas_SprFlow.sp,
                                                 object = mean~1)

#plot
ggplot(YOY_Betas_SprFlow_variogram) +
  geom_point(aes(x = dist))+
  labs(x = "Distance Class (km)",
       y = "Variance",
       title = "Variogram of YOY SprFlow Betas") +
  theme_classic()

# Semivariogram
YOY_Betas_SprFlow.geo <- as.geodata(YOY_Betas_SprFlow[,c(14:15, 2)])
# use sp.transform
dup.coords(YOY_Betas_SprFlow.geo)
YOY_Betas_SprFlow_semivariogram <- geoR::variog(YOY_Betas_SprFlow.geo)
plot(YOY_Betas_SprFlow_semivariogram,
     main = "Semivariogram of YOY SprFlow Betas")
