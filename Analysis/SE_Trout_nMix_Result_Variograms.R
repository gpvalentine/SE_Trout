## Variogram analysis of Bayesian N-mixture model parameters
## George Valentine
## Oct. 2022

## Load packages
library(tidyverse)
library(data.table)
library(gstat)      # for making variograms
library(geoR)       # for making semivariograms
library(sp)

## This script only works with objects created by "SE_Trout_Nmix.R"

## ICC Values
# Join site data to ICC values
YOY_ICCs <- YOY_randomEffects_params %>%
  rownames_to_column(., "param") %>%
  filter(str_detect(param, "ICC.YOY\\[")) %>%
  cbind(segment_data_RE)

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
YOY_ICCs.geo <- as.geodata(YOY_ICCs[,c(12:13, 2)])
dup.coords(YOY_ICCs.geo)
YOY_ICC_semivariogram <- geoR::variog(YOY_ICCs.geo)
YOY_ICCs_semivariogram.plot <- ggplot(data = as.data.frame(YOY_ICC_semivariogram[1:2])) +
  geom_point(aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic()


## Beta values: Mean Max Summer Temp
# Join site data to ICC values
YOY_Betas_SummTemp <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[1,")) %>% # this filters out the beta[1,] parameters
  cbind(segment_data_CE)

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
YOY_Betas_SummTemp.geo <- as.geodata(YOY_Betas_SummTemp[,c(12:13, 2)])
YOY_Betas_SummTemp_semivariogram <- geoR::variog(YOY_Betas_SummTemp.geo)
YOY_SummTempBetas_semivariogram.plot <- ggplot(data = as.data.frame(YOY_Betas_SummTemp_semivariogram[1:2])) +
  geom_point(aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic()

## Beta values: Mean Max Winter Flow
# Join site data to ICC values
YOY_Betas_WintFlow <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[2,")) %>% # this filters out the beta[1,] parameters
  cbind(segment_data_CE)

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
YOY_Betas_WintFlow.geo <- as.geodata(YOY_Betas_WintFlow[,c(12:13, 2)], coords.col = c(1,2))
dup.coords(YOY_Betas_WintFlow.geo)
YOY_Betas_WintFlow_semivariogram <- geoR::variog(YOY_Betas_WintFlow.geo)
YOY_WintFlowBetas_semivariogram.plot <- ggplot(data = as.data.frame(YOY_Betas_WintFlow_semivariogram[1:2])) +
  geom_point(aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic()

## Beta values: Mean Max Spring Flow
# Join site data to ICC values
YOY_Betas_SprFlow <- YOY_climateEffects_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta\\[3,")) %>% # this filters out the beta[1,] parameters
  cbind(segment_data_CE)

YOY_Betas_SprFlow.sp <- YOY_Betas_SprFlow
coordinates(YOY_Betas_SprFlow.sp) <- ~Long + Lat
proj4string(YOY_Betas_SprFlow.sp) = "+proj=longlat"

# create variogram
YOY_Betas_SprFlow_variogram <- gstat::variogram(data = YOY_Betas_SprFlow.sp,
                                                 object = mean~1)

#plot
# ggplot(YOY_Betas_SprFlow_variogram) +
#   geom_point(aes(x = dist))+
#   labs(x = "Distance Class (km)",
#        y = "Variance",
#        title = "Variogram of YOY SprFlow Betas") +
#   theme_classic()

# Semivariogram
YOY_Betas_SprFlow.geo <- as.geodata(YOY_Betas_SprFlow[,c(12:13, 2)])
# use sp.transform
dup.coords(YOY_Betas_SprFlow.geo)
YOY_Betas_SprFlow_semivariogram <- geoR::variog(YOY_Betas_SprFlow.geo)
YOY_SprFlowBetas_semivariogram.plot <- ggplot(data = as.data.frame(YOY_Betas_SprFlow_semivariogram[1:2])) +
  geom_point(aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic()

########################
# make a compound semivariogram
YOY_Betas_compound_semivariogram.table <- rbind(as.data.frame(YOY_Betas_SummTemp_semivariogram[1:2]),
                                                as.data.frame(YOY_Betas_WintFlow_semivariogram[1:2]),
                                                as.data.frame(YOY_Betas_SprFlow_semivariogram[1:2])) %>% 
  cbind(rbind(data.frame(covar = rep("Summer Temperature", 13)),
              data.frame(covar = rep("Winter Flow", 13)),
              data.frame(covar = rep("Spring Flow", 13))))

# reorder the covariates so that summer temperature plots first
YOY_Betas_compound_semivariogram.table$covar <- factor(YOY_Betas_compound_semivariogram.table$covar, c("Summer Temperature", "Winter Flow", "Spring Flow"))


YOY_climate_effects_semivariogram.plot <- ggplot() +
  geom_point(data = YOY_Betas_compound_semivariogram.table,
             aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic() +
  facet_grid(. ~ covar)

########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here("results", "v3.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
