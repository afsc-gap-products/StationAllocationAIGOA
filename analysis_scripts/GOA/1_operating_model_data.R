##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Groundfish CPUE Data Pull
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create CPUE dataset used for sdmTMB for species of interest
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to Oracle. Make sure you're on the VPN/network.
##   Set up constants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gapindex)
library(sf)
chl <- gapindex::get_connected(check_access = FALSE)

spp_set <- read.csv(file = "data/GOA/species_list/species_list.csv")
planning_years <- c(1996, 1999, seq(from = 2003, to = 2023, by = 2))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull catch and effort data and calculate CPUE (zero-filled)
##   Standard data: ABUNDANCE_HAUL = "Y" hauls
##
##   Non-standard data: hauls in NMFS areas 519 (Unimak Pass) and 659 (SE
##   Inside) that were previously ABUNDANCE_HAUL = "Y" but are now outside the
##   survey footprint when the survey design was updated in 2025. These hauls
##   are now coded as HAUL_TYPE 24
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Pull Data from RACEBASE
goa_standard_data <-
  gapindex::get_data(survey_set = "GOA",
                     year_set = planning_years,
                     spp_codes = spp_set[,-3],
                     pull_lengths = FALSE,
                     channel = chl)

## Calculate and zero-fill CPUE
goa_standard_cpue <-
  gapindex::calc_cpue(gapdata = goa_standard_data) |> as.data.frame()

goa_nonstandard_data <-
  gapindex::get_data(survey_set = "GOA",
                     year_set = planning_years,
                     spp_codes = spp_set[,-3],
                     abundance_haul = "N", haul_type = 24,
                     pull_lengths = FALSE,
                     taxonomic_source = "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION",
                     channel = chl)
goa_nonstandard_cpue <-
  as.data.frame(x = gapindex::calc_cpue(gapdata = goa_nonstandard_data))

## Combine standard and non-standard CPUE
goa_all_cpue <- rbind(goa_standard_cpue, goa_nonstandard_cpue)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a lat/lon spatial object of the station locations
##   Transform station spatial object to UTM (zone 5)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_sf <- sf::st_as_sf(x = goa_all_cpue,
                        coords = c("LONGITUDE_DD_START", "LATITUDE_DD_START"),
                        crs = "+proj=longlat +datum=WGS84")
goa_all_cpue[, c("E_km_z5", "N_km_z5")] <-
  sf::st_coordinates(sf::st_transform(x = cpue_sf,
                                      crs = "+proj=utm +zone=5 +units=km"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format CPUE data for sdmTMB() and save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_allspp <- with(goa_all_cpue,
                   data.frame(hauljoin = HAULJOIN,
                              year = as.integer(YEAR),
                              lon = LONGITUDE_DD_START,
                              lat = LATITUDE_DD_START,
                              X = E_km_z5,
                              Y = N_km_z5,
                              species = SPECIES_CODE,
                              catch_kg = WEIGHT_KG,
                              catch_n = COUNT,
                              effort_km2 = AREA_SWEPT_KM2,
                              cpue_kg_km2 = CPUE_KGKM2,
                              cpue_n_km2 = CPUE_NOKM2))

## Save catch and effort data
if (!dir.exists(paths = "data/GOA/sdmtmb_data/"))
  dir.create(path = "data/GOA/sdmtmb_data/")
write.csv(x = dat_allspp,
          file = "data/GOA/sdmtmb_data/goa_data_geostat.csv",
          row.names = F)
saveRDS(object = dat_allspp,
        file = "data/GOA/sdmtmb_data/goa_data_geostat.RDS")
