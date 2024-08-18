##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Plot station locations by vessel for the Gulf of Alaska
##                Bottom Trawl Survey from 2015 - 2023
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

## Import libraries, connect to Oracle
library(akgfmaps)
library(gapindex)

chl <- gapindex::get_connected(check_access = F)

## Import GOA spatial layers, station locations
goa_base_layers <- akgfmaps::get_base_layers(select.region = "goa")
goa_data <- gapindex::get_data(year_set = c(2023, 2021, 2019, 2017, 2015),
                               survey_set = "GOA", channel = chl)

## Loop over years and plot stations by vessel

pdf(file = "analysis_scripts/GOA/goa_restratification_2025/previous_haul_locations/previous_haul_locations.pdf",
    width = 6, height = 5, onefile = TRUE)
for (cruiseid in c(202301, 202101, 201901, 201701, 201501)) {
  plot(sf::st_geometry(obj = goa_base_layers$survey.area), main = cruiseid )
  plot(sf::st_geometry(obj = goa_base_layers$akland),
       add = TRUE, col = "tan", border = F )
  points(START_LATITUDE ~ START_LONGITUDE,
         data = goa_data$haul,
         subset = CRUISE == cruiseid,
         pch = 16, cex = 0.5,
         col = c("143" = "blue", "148" = "black",
                 "176" = "red", "178" = "green")[paste(VESSEL)])
  legend("bottom", legend = c("143", "148", "176", "178"),
         col = c("143" = "blue", "148" = "black",
                 "176" = "red", "178" = "green"),
         pch = 16, title = "GAP Vessel Code", bty = "n")
}
dev.off()
