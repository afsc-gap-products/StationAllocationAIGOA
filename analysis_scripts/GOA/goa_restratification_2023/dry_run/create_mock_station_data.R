##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create "mock" 2023 stations
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Take the 2021 GOA station data and reclassify to the 2025
##                 restratification
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import AFSC.GAP.DBE R package and terra
##   Connect to Oracle. Make sure you are connected to the VPN
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github("zoyafuso-NOAA/design-based-indices")
library(terra)
library(lubridate)
library(AFSC.GAP.DBE)
sql_channel <- AFSC.GAP.DBE::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get RACEBASE data for species_of_interest in 2021
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
racebase_data <- AFSC.GAP.DBE::get_data(
  year_set = 2019,
  survey_set = "GOA",
  # spp_codes = species_of_interest,
  spp_codes = data.frame("SPECIES_CODE" = c(21720, 21740),
                         "GROUP" = c(21720, 21740)),
  haul_type = 3,
  abundance_haul = c("Y"),
  pull_lengths = TRUE,
  sql_channel = sql_channel)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_2025 <- terra::vect("data/GOA/processed_shapefiles/goa_strata_2025.shp")
names(strata_2025) <- c("SURVEY", "STRATUM", "AREA", "PERIMETER",
                        "INPFC_AREA", "MIN_DEPTH", "MAX_DEPTH",
                        "DESCRIPTION", "SUMMARY_AREA", "SUMMARY_DEPTH",
                        "SUMMARY_AREA_DEPTH", "REGULATORY_AREA_NAME",
                        "STRATUM_TYPE")
# strata_2025 <- strata_2025[, c("SURVEY", "STRATUM", "AREA", "DESCRIPTION")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Reclassify haul station locations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
haul_locs_lat <-
  terra::vect(x = racebase_data$haul[, c("START_LONGITUDE", "START_LATITUDE")],
              geom = c("START_LONGITUDE", "START_LATITUDE"),
              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

haul_locs_aea <-
  terra::project(x = haul_locs_lat,
                 terra::crs(x = strata_2025))

haul_locs_2025 <- terra::extract(x = strata_2025, y = haul_locs_aea)

# library(StationAllocationAIGOA)
# opt_allocation <- StationAllocationAIGOA::goa_allocate_stations(n = 541)
# opt_allocation$ms_allocation
# table(haul_locs_2025$STRATUM)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Reclassify racebase_data based on new stratification
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
racebase_data_mock <- racebase_data
racebase_data_mock$survey$DESIGN_YEAR <- 2023

racebase_data_mock$cruise$YEAR <- racebase_data_mock$cruise$DESIGN_YEAR <- 2023
racebase_data_mock$cruise$CRUISE <- 202301
racebase_data_mock$cruise$CRUISEJOIN <- racebase_data_mock$cruise$CRUISEJOIN * 1000


racebase_data_mock$haul[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$haul[, c("CRUISEJOIN", "HAULJOIN" )] * 1000
racebase_data_mock$haul$CRUISE <- 202301
racebase_data_mock$haul$STRATUM <- haul_locs_2025$STRATUM
racebase_data_mock$haul$START_TIME <-
  racebase_data_mock$haul$START_TIME + years(4)

racebase_data_mock$strata <- cbind(SURVEY = "GOA",
                                   subset(x = AFSC.GAP.DBE::new_stratum_table,
                                          subset = TYPE == "STRATUM" &
                                            SURVEY_DEFINITION_ID == 47 &
                                            DESIGN_YEAR == 2023))
names(racebase_data_mock$strata)[names(racebase_data_mock$strata) ==
                                   "AREA_ID"] <- "STRATUM"

racebase_data_mock$catch$HAULJOIN <- racebase_data$catch$HAULJOIN * 1000

racebase_data_mock$size[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$size[, c("CRUISEJOIN", "HAULJOIN")] * 1000

racebase_data_mock$size$CRUISE <- 202301

racebase_data_mock$specimen[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$specimen[, c("CRUISEJOIN", "HAULJOIN")] * 1000
racebase_data_mock$specimen$CRUISE <- 202301

## Calculate CPUE
racebase_cpue <- AFSC.GAP.DBE::calc_cpue(racebase_tables = racebase_data_mock)
# racebase_cpue$YEAR <- 2023

## Calculate Stratum Biomass
racebase_biomass_stratum <-
  AFSC.GAP.DBE::calc_biomass_stratum(racebase_tables = racebase_data_mock,
                                     cpue = racebase_cpue,
                                     vulnerability = 1)

racebase_biomass_subareas <-
  AFSC.GAP.DBE::calc_agg_biomass(racebase_tables = racebase_data_mock,
                                 biomass_strata = racebase_biomass_stratum)


combined_biomass <-
  rbind(data.frame(
    AREA_ID = racebase_biomass_stratum$STRATUM,
    racebase_biomass_stratum[, c("SURVEY_DEFINITION_ID", "GROUP",
                                 "YEAR",
                                 "BIOMASS_MT", "BIOMASS_VAR",
                                 "POPULATION_COUNT", "POPULATION_VAR")]),
    racebase_biomass_subareas)

combined_biomass <- combined_biomass[order(x = combined_biomass$AREA_ID,
                                           decreasing = T), ]

write.csv(x = combined_biomass,
          file = "C:/Users/zack.oyafuso/Desktop/combined_biomass_GOA_2023.csv")

sizecomp_AIGOA <-
  AFSC.GAP.DBE::calc_size_stratum_AIGOA(
    racebase_tables = racebase_data_mock,
    racebase_cpue = racebase_cpue,
    racebase_stratum_popn = racebase_biomass_stratum)

test_agg_sizecomp <-
  rbind(AFSC.GAP.DBE::calc_agg_size_comp(racebase_tables = racebase_data_mock,
                                         size_comps = sizecomp_AIGOA))

combined_sizecomps <-
  rbind(test_agg_sizecomp,
        data.frame(AREA_ID = sizecomp_AIGOA$STRATUM,
                   subset(x = sizecomp_AIGOA,
                          select = c(SURVEY_DEFINITION_ID, SPECIES_CODE, YEAR,
                                     LENGTH_MM, MALES, FEMALES, UNSEXED))))

write.csv(x = combined_sizecomps,
          file = "C:/Users/zack.oyafuso/Desktop/combined_sizecomps_GOA_2023.csv",
          row.names = F)

test_agecomps <-
  AFSC.GAP.DBE::calc_age_comp(racebase_tables = racebase_data_mock,
                              size_comp = sizecomp_AIGOA)

write.csv(x = test_agecomps,
          file = "C:/Users/zack.oyafuso/Desktop/combined_agecomps_GOA_2023.csv",
          row.names = F)
