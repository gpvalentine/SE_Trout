# George Valentine
# Script to get segment summary data

# Load Spatiotemporal Covariate Data
SE_COMID_temp_covars <- fread(here("Data", "SE_Trout_COMID_temp_covars.csv"))
SE_COMID_flow_covars <- fread(here("Data", "SE_Trout_COMID_flow_covars.csv"))

# Also used in this script:
  # COMID_data - Segment data from "SE_Trout_nMix.R"

# Get NHDPlus data for the segments considered for the analysis
SE_segments_NHDPlus <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")
SE_segments_NHDPlus <- SE_segments_NHDPlus %>%
  filter(COMID %in% COMID_data$COMID) %>% 
  dplyr::select(SLOPE,
                LENGTHKM,
                AreaSqKM,
                MINELEVSMO,
                MAXELEVSMO,
                AreaSqKM,
                StreamOrde)

# Get widths of sites within the segments considered for the analysis
SE_segments <- SE_Site_Final %>% 
  filter(COMID %in% COMID_data$COMID) %>% 
  select(Width_m,
         Elev_m)

Segment_Summary.table <- data.frame(Variable = c("Mean 90th percentile summer temperature (C)",
                                           "Maximum 90th percentile winter streamflow (ft^3/s)",
                                           "Maximum 90th percentile spring streamflow (ft^3/s)",
                                           "Channel Slope (%)",
                                           "Length (km)",
                                           "Catchment area (km^2)",
                                           "Elevation (m)",
                                           "Stream order",
                                           "Wetted width (m)"),
                                    Mean = c(mean(SE_COMID_temp_covars$Mean_Max_Summer_Temp),
                                             mean(SE_COMID_flow_covars$Max_0.9Q_WinterFlow),
                                             mean(SE_COMID_flow_covars$Max_0.9Q_SpringFlow),
                                             mean(SE_segments_NHDPlus$SLOPE, na.rm = T)*100,
                                             mean(SE_segments_NHDPlus$LENGTHKM, na.rm = T),
                                             mean(SE_segments_NHDPlus$AreaSqKM, na.rm = T),
                                             mean(rowMeans(SE_segments_NHDPlus[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                             mean(SE_segments_NHDPlus$StreamOrde, na.rm = T),
                                             mean(SE_segments$Width_m, na.rm = T)),
                                    Median = c(median(SE_COMID_temp_covars$Mean_Max_Summer_Temp),
                                               median(SE_COMID_flow_covars$Max_0.9Q_WinterFlow),
                                               median(SE_COMID_flow_covars$Max_0.9Q_SpringFlow),
                                               median(SE_segments_NHDPlus$SLOPE, na.rm = T)*100,
                                               median(SE_segments_NHDPlus$LENGTHKM, na.rm = T),
                                               median(SE_segments_NHDPlus$AreaSqKM, na.rm = T),
                                               median(rowMeans(SE_segments_NHDPlus[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                               median(SE_segments_NHDPlus$StreamOrde, na.rm = T),
                                               median(SE_segments$Width_m, na.rm = T)),
                                    Maximum = c(max(SE_COMID_temp_covars$Mean_Max_Summer_Temp),
                                                max(SE_COMID_flow_covars$Max_0.9Q_WinterFlow),
                                                max(SE_COMID_flow_covars$Max_0.9Q_SpringFlow),
                                                max(SE_segments_NHDPlus$SLOPE, na.rm = T)*100,
                                                max(SE_segments_NHDPlus$LENGTHKM, na.rm = T),
                                                max(SE_segments_NHDPlus$AreaSqKM, na.rm = T),
                                                max(rowMeans(SE_segments_NHDPlus[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                max(SE_segments_NHDPlus$StreamOrde, na.rm = T),
                                                max(SE_segments$Width_m, na.rm = T)),
                                    Minimum = c(min(SE_COMID_temp_covars$Mean_Max_Summer_Temp),
                                                min(SE_COMID_flow_covars$Max_0.9Q_WinterFlow),
                                                min(SE_COMID_flow_covars$Max_0.9Q_SpringFlow),
                                                min(SE_segments_NHDPlus$SLOPE, na.rm = T)*100,
                                                min(SE_segments_NHDPlus$LENGTHKM, na.rm = T),
                                                min(SE_segments_NHDPlus$AreaSqKM, na.rm = T),
                                                min(rowMeans(SE_segments_NHDPlus[,c("MAXELEVSMO", "MINELEVSMO")], na.rm = T), na.rm = T)/100,
                                                min(SE_segments_NHDPlus$StreamOrde, na.rm = T),
                                                min(SE_segments$Width_m, na.rm = T)))

