##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create "mock" 2023 stations
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Take the 2021 GOA station data and reclassify to the 2025
##                 restratification
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import gapindex R package and terra
##   Connect to Oracle. Make sure you are connected to the VPN
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github("zoyafuso-NOAA/design-based-indices")
library(terra)
library(lubridate)
library(gapindex)
sql_channel <- gapindex::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Get RACEBASE data for biomass/abundance taxa for year 2019
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
racebase_data <- gapindex::get_data(
  year_set = 2019,
  survey_set = "GOA",
  spp_codes = NULL,
  haul_type = 3,
  abundance_haul = c("Y"),
  pull_lengths = TRUE,
  sql_channel = sql_channel)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_2025 <-
  terra::vect(x = "data/GOA/processed_shapefiles/goa_strata_2025.shp")
names(x = strata_2025) <- c("SURVEY", "STRATUM", "AREA", "PERIMETER",
                            "INPFC_AREA", "MIN_DEPTH", "MAX_DEPTH",
                            "DESCRIPTION", "SUMMARY_AREA", "SUMMARY_DEPTH",
                            "SUMMARY_AREA_DEPTH", "REGULATORY_AREA_NAME",
                            "STRATUM_TYPE")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Reclassify haul station locations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
haul_locs_lat <-
  terra::vect(x = racebase_data$haul[, c("START_LONGITUDE", "START_LATITUDE")],
              geom = c("START_LONGITUDE", "START_LATITUDE"),
              crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

haul_locs_aea <- terra::project(x = haul_locs_lat,
                                terra::crs(x = strata_2025))

haul_locs_2025 <- terra::extract(x = strata_2025, y = haul_locs_aea)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Reclassify racebase_data based on new stratification
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Change year to 2023 in the survey cruise data
racebase_data_mock <- racebase_data
racebase_data_mock$survey$DESIGN_YEAR <- 2025
racebase_data_mock$survey_design$YEAR <-
  racebase_data_mock$survey_design$DESIGN_YEAR <- 2025

## create new (unique) CRUISE, CRUISEJOIN, and HAULJOIN values in the cruise
## and haul data. Modify the START_TIME to 2025 by advancing 6 years
racebase_data_mock$cruise$CRUISE <- 202501
racebase_data_mock$cruise$CRUISEJOIN <-
  racebase_data_mock$cruise$CRUISEJOIN * 1000
racebase_data_mock$cruise$YEAR <-
  racebase_data_mock$cruise$DESIGN_YEAR <- 2025

racebase_data_mock$haul[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$haul[, c("CRUISEJOIN", "HAULJOIN" )] * 1000
racebase_data_mock$haul$CRUISE <- 202501
racebase_data_mock$haul$STRATUM <- haul_locs_2025$STRATUM
racebase_data_mock$haul$START_TIME <-
  racebase_data_mock$haul$START_TIME + lubridate::years(6)

## Set strata table as the 2023 strata
racebase_data_mock$strata <-
  RODBC::sqlQuery(channel = sql_channel,
                  query = paste0("SELECT * FROM GAP_PRODUCTS.AREA ",
                                 "WHERE SURVEY_DEFINITION_ID = 47 ",
                                 "AND DESIGN_YEAR = 2025 AND TYPE = 'STRATUM'"))
names(racebase_data_mock$strata)[names(racebase_data_mock$strata) ==
                                   "AREA_ID"] <- "STRATUM"
racebase_data_mock$strata$SURVEY <- "GOA"

racebase_data_mock$subarea <-
  RODBC::sqlQuery(channel = sql_channel,
                  query = paste0("SELECT * FROM GAP_PRODUCTS.AREA ",
                                 "WHERE SURVEY_DEFINITION_ID = 47 ",
                                 "AND DESIGN_YEAR = 2025 AND TYPE != 'STRATUM'"))

racebase_data_mock$stratum_groups <-
  RODBC::sqlQuery(channel = sql_channel,
                  query = paste0("SELECT * FROM GAP_PRODUCTS.STRATUM_GROUPS ",
                                 "WHERE SURVEY_DEFINITION_ID = 47 ",
                                 "AND DESIGN_YEAR = 2025"))

racebase_data_mock$catch$HAULJOIN <- racebase_data$catch$HAULJOIN * 1000

racebase_data_mock$size[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$size[, c("CRUISEJOIN", "HAULJOIN")] * 1000

racebase_data_mock$size$CRUISE <- 202501

racebase_data_mock$specimen[, c("CRUISEJOIN", "HAULJOIN")] <-
  racebase_data_mock$specimen[, c("CRUISEJOIN", "HAULJOIN")] * 1000
racebase_data_mock$specimen$CRUISE <- 202501

## Calculate CPUE
racebase_cpue <- gapindex::calc_cpue(racebase_tables = racebase_data_mock)

## Calculate biomass/abundance/mean and var CPUE across strata
racebase_biomass_stratum <-
  gapindex::calc_biomass_stratum(racebase_tables = racebase_data_mock,
                                 cpue = racebase_cpue)

## Aggregate biomass/abundance/mean and var CPUE across subareas
racebase_biomass_subareas <-
  gapindex::calc_biomass_subarea(racebase_tables = racebase_data_mock,
                                 biomass_strata = racebase_biomass_stratum)

## Calculate size composition across strata
racebase_sizecomp_stratum <-
  gapindex::calc_sizecomp_stratum(
    racebase_tables = racebase_data_mock,
    racebase_cpue = racebase_cpue,
    racebase_stratum_popn = racebase_biomass_stratum,
    spatial_level = "stratum",
    fill_NA_method = "AIGOA")

## Aggregate size composition across subareas
racebase_sizecomp_subareas <- gapindex::calc_sizecomp_subareas(
  racebase_tables = racebase_data_mock,
  size_comps = racebase_sizecomp_stratum)

## Calculate age-length key
racebase_alk <- gapindex::calc_ALK(racebase_tables = racebase_data_mock,
                                   unsex = "all",
                                   global = F)

## Calculate age composition across strata
racebase_agecomp_stratum <- gapindex::calc_agecomp_stratum(
  racebase_tables = racebase_data_mock,
  alk = racebase_alk,
  size_comp = racebase_sizecomp_stratum)

## Aggregate age composition for the entire region
racebase_agecomp_subarea <- gapindex::calc_agecomp_region(
  racebase_tables = racebase_data_mock,
  age_comps_stratum = racebase_agecomp_stratum)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Combine stratum and subarea tables
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove extra columns from `production_cpue`
racebase_cpue <- subset(x = racebase_cpue,
                        select = c(HAULJOIN, SPECIES_CODE,
                                   WEIGHT_KG, COUNT, AREA_SWEPT_KM2,
                                   CPUE_KGKM2, CPUE_NOKM2) )

names(x = racebase_biomass_stratum)[
  names(x = racebase_biomass_stratum) == "STRATUM"
] <- "AREA_ID"

racebase_biomass <-
  rbind(racebase_biomass_stratum[, names(racebase_biomass_subareas)],
        racebase_biomass_subareas)
racebase_biomass <- subset(x = racebase_biomass, select = -SURVEY)

names(x = racebase_sizecomp_stratum)[
  names(x = racebase_sizecomp_stratum) == "STRATUM"
] <- "AREA_ID"

racebase_sizecomp <-
  rbind(racebase_sizecomp_subareas,
        racebase_sizecomp_stratum[, names(racebase_sizecomp_subareas)])


names(x = racebase_agecomp_stratum$age_comp)[
  names(x = racebase_agecomp_stratum$age_comp) == "STRATUM"
] <- "AREA_ID"

racebase_agecomp <-
  rbind(racebase_agecomp_subarea,
        racebase_agecomp_stratum$age_comp[, names(racebase_agecomp_subarea)])
racebase_agecomp <- subset(x = racebase_agecomp, select = -SURVEY)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull table metadata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
main_metadata_columns <-
  RODBC::sqlQuery(channel = sql_channel,
                  query = "SELECT * FROM GAP_PRODUCTS.METADATA_COLUMN")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Append to existing tables in Oracle (schema GAP_PRODUCTS)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (idata in c("agecomp",  "sizecomp", "biomass",
                 "cpue")[-4]) {

  match_idx <-
    match(x = names(x = get(x = paste0("racebase_", idata))),
          table = toupper(x = main_metadata_columns$METADATA_COLNAME))

  metadata_columns <-
    with(main_metadata_columns,
         data.frame(colname = toupper(x = METADATA_COLNAME[match_idx]),
                    colname_long = METADATA_COLNAME_LONG[match_idx],
                    units = METADATA_UNITS[match_idx],
                    datatype = METADATA_DATATYPE[match_idx],
                    colname_desc = METADATA_COLNAME_DESC[match_idx]))

  RODBC::sqlSave(channel = sql_channel,
                 dat = get(x = paste0("racebase_", idata)),
                 tablename = paste0("GAP_PRODUCTS.", toupper(x = idata)),
                 append = TRUE, rownames = F)

}
