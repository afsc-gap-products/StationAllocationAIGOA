##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Station Allocation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For year 2025
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages, laning areas around Kodiak
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::build()
# install.packages("C:/Users/zack.oyafuso/Work/GitHub/StationAllocationAIGOA_0.1.0.tar.gz")
library(StationAllocationAIGOA)
library(terra)
library(akgfmaps)

laning_area <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_laning_area.shp")

goa_stations <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_stations_2025.gpkg")
goa_stations[, c("x", "y")] <-
  terra::crds(x = terra::centroids(x = goa_stations,
                                   inside = TRUE))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate a 550 station allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shallow_boat <- 176 # Alaska Provider
deep_boat <- 148    # Ocean Explorer

goa_stn_allocation <- goa_allocate_stations(
  n = 550,
  min_n_per_stratum = 4,
  species = c(
    "arrowtooth flounder", ## Atherestes stomias
    "Pacific cod", ## Gadus macrocephalus
    "walleye pollock", ## Gadus chalcogrammus
    "rex sole", ## Glyptocephalus zachirus
    "flathead sole", ## Hippoglossoides elassodon
    "Pacific halibut", ## Hippoglossus stenolepis
    "southern rock sole", ## Lepidopsetta bilineata
    "northern rock sole", ## Lepidopsetta polyxystra
    "Pacific ocean perch", ## Sebastes alutus
    "silvergray rockfish", ## Sebastes brevispinis
    "northern rockfish", ## Sebastes polyspinis
    "dusky rockfish", ## Sebastes variabilis
    "REBS rockfish", ## Sebastes aleutianus and S. melanostictus
    "Dover sole", ## Microstomus pacificus
    "shortspine thornyhead" ## Sebastolobus alascanus
  ),
  max_iter = 5000,
  trawl = c("Y", "N", "UNK")[c(1, 3)],
  survey_year = 2025,
  planning_years = c(1996, 1999, seq(from = 2003, to = 2023, by = 2))
)

stn_allocation <- goa_stn_allocation$drawn_stations

stn_allocation <- merge(x = stn_allocation,
                        y = goa_stations[, c("STATION", "x", "y")],
                        by = c("STATION"))

stns_in_lane <- terra::relate(x = terra::vect(x = stn_allocation,
                                              geom = c("x", "y")),
                              y = laning_area,
                              relation = "intersects",
                              pairs = T)

stn_allocation[stns_in_lane[, "id.x"], "VESSEL"] <-
  c(deep_boat, shallow_boat)[stns_in_lane[, "id.y"]]

diff_in_stns <- abs(x = diff(x = table(stn_allocation$VESSEL)))

goa_strata <- StationAllocationAIGOA::goa_stratum_boundaries

for (istratum in goa_strata$STRATUM[goa_strata$USED]) {
  n_in_stratum <- sum(stn_allocation$STRATUM == istratum &
                        is.na(x = stn_allocation$VESSEL))

  vessel_assignments <- vector(length = n_in_stratum)
  vessel_assignments[1:n_in_stratum] <- sample(x = c(deep_boat, shallow_boat))

  stn_allocation[stn_allocation$STRATUM == istratum &
                   is.na(x = stn_allocation$VESSEL), "VESSEL"] <-
    vessel_assignments

}


## Quick Plot
goa_base_layers <-
  akgfmaps::get_base_layers(select.region = "goa", set.crs = "EPSG:3338")

pdf(file = "analysis_scripts/GOA/goa_restratification_2025/previous_haul_locations/GOA_2025_stations.pdf",
    width = 6, height = 5, onefile = TRUE)
plot(sf::st_geometry(obj = goa_base_layers$survey.area))
plot(sf::st_geometry(obj = goa_base_layers$akland),
     add = TRUE, col = "tan", border = F )
points(y ~ x, data = stn_allocation, cex = 0.5, pch = 16,
       col = c("176" = "black",
               "148" = "red")[paste(stn_allocation$VESSEL)])
legend("bottom", legend = c("148", "176"),
       col = c("148" = "black", "176" = "red"),
       pch = 16, title = "GAP Vessel Code", bty = "n")

plot(sf::st_geometry(obj = goa_base_layers$survey.area))
plot(sf::st_geometry(obj = goa_base_layers$akland),
     add = TRUE, col = "tan", border = F )
points(y ~ x, data = stn_allocation, cex = 0.5, pch = 16,
       col = c("176" = "black",
               "148" = "red")[paste(stn_allocation$VESSEL)])
plot(laning_area, add = TRUE, border = c("red", "black"), lwd = 2)
legend("bottom", legend = c("148", "176"),
       col = c("148" = "black", "176" = "red"),
       pch = 16, title = "GAP Vessel Code", bty = "n")
dev.off()
