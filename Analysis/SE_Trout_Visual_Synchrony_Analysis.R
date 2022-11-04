###############################################
# Do populations visually show synchrony?
# Make the "spaghetti plot" that shows time series of adundance at individual COMIDs over time

# the "BKT_COMID_Dens" comes from the "SE_Trout_Correlograms.R" file
timeSeries <- BKT_COMID_Dens %>% 
  mutate(Subregion = ifelse(Lat < 37.13, "South", "North")) %>% 
  left_join(SE_Site_Final[,c(6,10)])


YOY_pass1_spaghetti.plot <- ggplot(data = timeSeries) +
  geom_line(aes(x = Year,
                y = Avg_Density,
                group = COMID,
                color = Source),
            size = 0.35,
            alpha = 0.5) +
  labs(y = "Average BKT Density (fish/1000m^2)") +
  theme_classic() +
  facet_wrap(~Subregion)
