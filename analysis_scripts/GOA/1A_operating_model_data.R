###############################################################################
## Project:         GOA Groundfish CPUE Data Synthesis
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors:    Lewis Barnett (lewis.barnett@noaa.gov)
## Description:     Create CPUE dataset used for VAST for species of interest
###############################################################################
rm(list = ls())

##################################################
#### Import gapindex, Connect to Oracle
##################################################
library(gapindex)
chl <- gapindex::get_connected(check_access = FALSE)

##################################################
#### Set up haul-level CPUE survey, update with most recent survey year
##################################################
spp_set <- read.csv(file = "data/GOA/species_list/species_list.csv")
planning_years <- c(1996, 1999, seq(from = 2003, to = 2023, by = 2))

## Pull Data from RACEBASE
gapindex_data <- gapindex::get_data(survey_set = "GOA",
                                    year_set = planning_years,
                                    spp_codes = spp_set[,-3],
                                    pull_lengths = FALSE,
                                    channel = chl)

## Calculate and zero-fill CPUE
gapindex_cpue <- gapindex::calc_cpue(gapdata = gapindex_data)

## Format catch and effort data for VAST. Note: The `AreaSwept_km2` field set
## to 1 when using CPUE in the `Catch_KG` field.
data_geostat <-
  with(gapindex_cpue,
       data.frame(Hauljoin = HAULJOIN,
                  Species = SPECIES_CODE,
                  Year = YEAR,
                  Vessel = "missing",
                  AreaSwept_km2 = 1,
                  Catch_KG = CPUE_KGKM2,
                  Lat = LATITUDE_DD_START,
                  Lon = LONGITUDE_DD_START,
                  Depth_m = DEPTH_M,
                  LOG10_DEPTH_M = log10(x = DEPTH_M),
                  LOG10_DEPTH_M_CEN = scale(x = log10(x = DEPTH_M)),
                  Pass = 0) )

##################################################
####   Save
##################################################
if (!dir.exists(paths = "data/GOA/vast_data/"))
  dir.create(path = "data/GOA/vast_data/")
write.csv(x = data_geostat,
          file = "data/GOA/vast_data/goa_data_geostat.csv",
          row.names = F)
saveRDS(object = data_geostat,
        file = "data/GOA/vast_data/goa_data_geostat.RDS")
