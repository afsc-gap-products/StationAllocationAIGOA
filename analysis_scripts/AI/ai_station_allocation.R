##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Aleutian Island BTS Station Allocation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##
## Description:   Station Allocation Procedure for the 2026 Aleutian Islands
##                bottom trawl survey. Make sure you are connected to the VPN
##                or network to connect to Oracle and to access the price data
##                spreadsheet in the G: drive.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Constants, update for a new survey
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options(scipen = 999999)
current_year <- 2026
n_allocation <- 400
n_bonus <- 20
n_new <- 8
n_cut <- c(-20, -40)
vessel_ids <- c(176, 148)
output_dir <- paste0("Y:/RACE_GF/ALEUTIAN/AI ", current_year,
                     "/Station Allocation/")
n_allocaton_with_bonus <- n_allocation + n_bonus
set.seed(seed = current_year)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages, connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library(devtools)
# devtools::build()
# devtools::install_github("afsc-gap-products/StationAllocationAIGOA")

library(StationAllocationAIGOA)
library(xlsx)
library(gapindex)
library(akgfmaps)
channel <- gapindex::get_connected(check_access = F)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query AI Strata information,
##   Based on the total stratum area and edge (
##   i.e., perimeter), determine whether a stratum is either
##   "Large": Edge:Area ratio < mean(Edge:Area) ratio across strata AND
##                  Total Area > 1000 km^2
##   OR
##   "Thin": otherwise
##   This information will be used to partition the allocation of new stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata <- akgfmaps::get_base_layers(select.region = "ai",
                                    design.year = 2024)$survey.strata
strata <- data.frame(stratum = strata$STRATUM,
                     area_km2 = strata$AREA / 1000 / 1000,
                     perim_km = sf::st_perimeter(x = strata) )
strata$perim_km <- units::set_units(x = strata$perim_km, "km")

strata$edge_to_area <- with(strata, perim_km / 2 / sqrt(area_km2 * pi) )

strata$strata_type <-
  ifelse(test = strata$edge_to_area < mean(x = strata$edge_to_area) &
           strata$area_km2 > 1000,
         yes = "large",
         no = "thin")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull candidate AI stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
candidate_stations <- StationAllocationAIGOA::get.ai.stations(channel = channel)

survey_stations <- akgfmaps::get_base_layers(select.region = "ai",
                                             set.crs = "EPSG:3338",
                                             design.year = 1991)$survey.grid
## output station grid with updated trawlability information
ai_grid_gis <-
  merge(x = survey_stations,
        y = RODBC::sqlQuery(
          channel = channel,
          query = "
                            SELECT STATIONID AS STATION, STRATUM, TRAWLABLE
                FROM AI.AIGRID_GIS"),
        by = c("STATION", "STRATUM"))

sf::write_sf(ai_grid_gis, dsn = paste0(output_dir, "aigridgis.gpkg"))

ai_hauls <- RODBC::sqlQuery(channel = channel,
                            query = paste0("SELECT * FROM RACEBASE.HAUL
                WHERE REGION = 'AI' "))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import price data from COAR. Because there is a species complex, it is
##   easier to run gapindex to calculate mean and sd across strata for each
##   species/species complex.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
price_data <- readxl::read_excel(
  path = paste0(output_dir,
                "Ex vessel prices/GOA_AI_PLANNING_SPECIES_01152026.xlsx"),
  sheet = "ex_vessel_prices"
)

## Create df
planning_species <-
  readxl::read_excel(
    path = paste0(output_dir,
                  "Ex vessel prices/GOA_AI_PLANNING_SPECIES_01152026.xlsx"),
    sheet = "species_groupings"
  )

## Get survey from Oracle
planning_species_data <- gapindex::get_data(year_set = 1991:2024,
                                            survey_set = "AI",
                                            spp_codes = planning_species,
                                            channel = channel)

## Calculate CPUE
planning_species_cpue <-
  gapindex::calc_cpue(gapdata = planning_species_data)

## Calculate CPUE and Biomass
planning_species_biomass <-
  gapindex::calc_biomass_stratum(gapdata = planning_species_data,
                                 cpue = planning_species_cpue)
planning_species_biomass_region <-
  gapindex::calc_biomass_subarea(gapdata = planning_species_data,
                                 biomass_stratum = planning_species_biomass) |>
  subset(subset = AREA_ID == 99904)

spp_start_year <-
  RODBC::sqlQuery(channel = channel,
                  query = "SELECT * FROM GAP_PRODUCTS.SPECIES_YEAR")

for (ispp in 1:nrow(x = spp_start_year)) {
  removed_idx <-
    which(planning_species_biomass_region$SPECIES_CODE ==
            spp_start_year$SPECIES_CODE[ispp] &
            planning_species_biomass_region$YEAR <
            spp_start_year$YEAR_STARTED[ispp])

  if (length(x = removed_idx) > 0) {
    planning_species_biomass_region <-
      planning_species_biomass_region[-removed_idx, ]
  }

  removed_idx <-
    which(planning_species_biomass$SPECIES_CODE ==
            spp_start_year$SPECIES_CODE[ispp] &
            planning_species_biomass$YEAR <
            spp_start_year$YEAR_STARTED[ispp])

  if (length(x = removed_idx) > 0) {
    planning_species_biomass <-
      planning_species_biomass[-removed_idx, ]
  }

}

## Calculate the stratum standard deviation of CPUE from the
## standard error of mean-CPUE.
planning_species_biomass$CPUE_KGKM2_SD <-
  sqrt(planning_species_biomass$CPUE_KGKM2_VAR *
         planning_species_biomass$N_HAUL)

stratum_stats <- subset(x = planning_species_biomass,
                        select = c(STRATUM, SPECIES_CODE, YEAR,
                                   CPUE_KGKM2_MEAN, CPUE_KGKM2_SD))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Economic value is used as the species weightings used in effort allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate the mean biomass across years for each planning taxon.
species_weightings <-
  stats::aggregate(BIOMASS_MT ~ SPECIES_CODE,
                   data = planning_species_biomass_region,
                   FUN = mean)

## Merge price data to species_weighting using SPECIES_CODE as key
species_weightings <-
  merge(x = species_weightings,
        y = subset(x = price_data,
                   select = c(common_name_coar, species_code,
                              ex_vessel_price_per_mt)),
        by.x = "SPECIES_CODE",
        by.y = "species_code")

## Calculate total economic value by multiplying by ex-vessel price.
species_weightings$econ_value <- with(species_weightings,
                                      BIOMASS_MT * ex_vessel_price_per_mt)

species_weightings <- species_weightings[order(species_weightings$econ_value,
                                               decreasing = TRUE), ]
species_weightings$econ_value_wgt <- (species_weightings$econ_value /
                                        sum(species_weightings$econ_value)) |>
  round(digits = 3)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate optimal allocation at 400 and 420 total stations
##   Minimum number of stations per stratum: 4 stations
##   Note: the n_total arguments are adjusted due to rounding errors in the
##   allocation procedure
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allocation_w_bonus <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = n_allocaton_with_bonus,
  min_n_stratum = 4)
sum(allocation_w_bonus$ms_allocation)

allocation <- allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = n_allocation + 1,
  min_n_stratum = 4)
sum(allocation$ms_allocation)

allocation_cut1 <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = n_allocation - 21,
  min_n_stratum = 4)
sum(allocation_cut1$ms_allocation)

allocation_cut2 <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = n_allocation - 37,
  min_n_stratum = 4)
sum(allocation_cut2$ms_allocation)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Bind the station effort across strata under 400 and 420 total stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station_allocation <-
  data.frame(
    stratum = names(x = allocation_w_bonus$ms_allocation),
    region = sapply(X = substr(x = names(x = allocation_w_bonus$ms_allocation),
                               start = 1,
                               stop = 1),
                    FUN = function(x)
                      switch(x,
                             "2" = "Western",
                             "3" = "Central", "4" = "Central",
                             "5" = "Eastern", "6" = "Eastern",
                             "7" = "SBS")),
    avail = tapply(X = candidate_stations$STATIONID ,
                   INDEX = candidate_stations$STRATUM,
                   FUN = length),
    allocation = as.numeric(x = allocation$ms_allocation),
    allocation_w_bonus = as.numeric(x = allocation_w_bonus$ms_allocation)
  )

## Merge strata information to station_allocation using "stratum" as the key
station_allocation <- merge(x = station_allocation,
                            y = strata,
                            by = "stratum")

## Bonus stations: the difference between the allocation at 420 total stations
## versus the allocation at 400 total stations.
station_allocation$bonus_stn <- with(station_allocation,
                                     allocation_w_bonus - allocation)
station_allocation$new_stn <- 0

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   For the other regions, we will randomly choose a large and a thin stratum
##   from which one allocated station from the chosen stratum will be converted
##   to a "new" station to be searched.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (iregion in c("Western", "Central",
                  "Eastern", "SBS")) { ## Loop over regions -- start
  chosen_large_stratum <-
    sample(x = which(station_allocation$region == iregion &
                       station_allocation$strata_type == "large"),
           size = 1)
  station_allocation$new_stn[chosen_large_stratum] <- 1
  station_allocation$allocation[chosen_large_stratum] <-
    station_allocation$allocation[chosen_large_stratum] - 1

  chosen_thin_stratum <-
    sample(x = which(station_allocation$region == iregion &
                       station_allocation$strata_type == "thin"),
           size = 1)
  station_allocation$new_stn[chosen_thin_stratum] <- 1
  station_allocation$allocation[chosen_thin_stratum] <-
    station_allocation$allocation[chosen_thin_stratum] - 1
} ## Loop over regions -- end

station_allocation <-
  with(station_allocation,
       data.frame(stratum,
                  strata_type,
                  region,
                  prescribed_stn = allocation,
                  new_stn,
                  bonus_stn))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Randomly draw stations based on the allocation, assign to Vessels
##   Turn the centroids of the stations into a spatial object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
picked_stns_pts <-
  StationAllocationAIGOA::draw_strs_stations(
    stn_allocation =
      data.frame(stratum = names(x = allocation$ms_allocation),
                 ms_allocation = as.numeric(x = allocation$ms_allocation)),
    survey_year = 2026,
    stn_df = cbind(candidate_stations,
                   TRAWLABLE = "Y"),
    trawl = "Y",
    min_area_km2 = 0
  ) |>
  subset(select = c(STATIONID, STRATUM)) |>
  StationAllocationAIGOA::assign_stations_to_vessels() |>
  cbind(STATION_TYPE = "assigned") |>
  merge(x = survey_stations,
        by.x = c("STATION", "STRATUM"),
        by.y = c("STATIONID", "STRATUM"))  %>%
  dplyr::group_by(STATION, STRATUM, VESSEL) %>%
  dplyr::summarise() %>%
  sf::st_centroid(of_largest_polygon = TRUE)

picked_stns_pts[, c("LONGITUDE", "LATITUDE")] <-
  picked_stns_pts |> sf::st_transform(crs = "EPSG:4326") |> sf::st_coordinates()

# picked_stns_pts$LONGITUDE <-
#   ifelse(test = picked_stns_pts$LONGITUDE > 0,
#          yes = picked_stns_pts$LONGITUDE - 360,
#          no = picked_stns_pts$LONGITUDE)
#
# plot(LATITUDE ~ LONGITUDE, data = picked_stns_pts,
#      pch = 16,
#      col = c("176" = "red", "148" = "black")[paste(picked_stns_pts$VESSEL)])
# plot(st_geometry(ai_land), add = TRUE)
#
# ai_land <- akgfmaps::get_base_layers(select.region = 'ai',
#                                      set.crs = "EPSG:4326",
#                                      split.land.at.180 = TRUE)$akland
# st_geometry(ai_land)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Create shapefile of prescribed stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
picked_stns <- rbind(
  data.frame(st_drop_geometry(x = picked_stns_pts),
             STATION_TYPE = "assigned"),
  ## Append New Stations
  apply(X = station_allocation |> subset(subset = new_stn > 0),
        MARGIN = 1,
        FUN = function(x)
          data.frame(STATION = NA,
                     STRATUM = rep(x["stratum"],
                                   x["new_stn"]),
                     STRATUM_TYPE = x["strata_type"])

  ) |>
    do.call(what = "rbind") |>
    StationAllocationAIGOA::assign_stations_to_vessels(order_by = "STRATUM_TYPE") |>
    cbind(LONGITUDE = NA, LATITUDE = NA, STATION_TYPE = "new") |>
    subset(select = -STRATUM_TYPE),

  ## Append Bonus Stations
  apply(X = station_allocation |> subset(subset = bonus_stn > 0),
        MARGIN = 1,
        FUN = function(x)
          data.frame(STATION = NA,
                     STRATUM = rep(x["stratum"],
                                   x["bonus_stn"]),
                     VESSEL = NA,
                     LONGITUDE = NA,
                     LATITUDE = NA,
                     STATION_TYPE = "bonus")
  ) |>
    do.call(what = "rbind")
)

## Randomly replace allocated stations with new stations
removed_stns <- c()
for (istn in which(picked_stns$STATION_TYPE == "new")) {
  temp_vessel <- picked_stns$VESSEL[istn]
  temp_stratum <- picked_stns$STRATUM[istn]

  removed_stns <- c(removed_stns,
                    sample(x = which(picked_stns$VESSEL == temp_vessel &
                                       picked_stns$STRATUM == temp_stratum &
                                       picked_stns$STATION_TYPE == "assigned"),
                           size = 1))

}

picked_stns <- picked_stns[-removed_stns, ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Write Results to G: Drive
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## First write the Station Allocation
xlsx::write.xlsx(
  x = picked_stns,
  file = paste0(output_dir, "ai_", current_year, "_station_allocation_",
                n_allocation, "stn.xlsx"),
  sheetName = "Station Allocation",
  row.names = FALSE,
  showNA = FALSE
)

## Next, write the Stratum Allocation by vessel and stratum as a new tab
stratum_allocation <-
  table(picked_stns$STRATUM, picked_stns$VESSEL) |>
  as.data.frame.matrix()

xlsx::write.xlsx(
  x = stratum_allocation,
  file = paste0(output_dir, "ai_", current_year, "_station_allocation_",
                n_allocation, "stn.xlsx"),
  sheetName = "Stratum Allocation",
  row.names = TRUE,
  showNA = FALSE,
  append = TRUE
)

## Next, write the general plan for dropping stations by stratum as a new tab
xlsx::write.xlsx(
  x = data.frame(stratum = names(x = allocation$ms_allocation) |> as.numeric(),
                 "Drop 20 stns" = allocation$ms_allocation - allocation_cut1$ms_allocation,
                 "Drop 40 stns" = allocation$ms_allocation - allocation_cut2$ms_allocation,
                 check.names = F),
  file = paste0(output_dir, "ai_", current_year, "_station_allocation_",
                n_allocation, "stn.xlsx"),
  sheetName = "Priority Strata for Stn Drops",
  row.names = TRUE,
  showNA = FALSE,
  append = TRUE
)

## Next, output the location of the assigned stations as a geopackage for
## charts
sf::write_sf(
  obj = picked_stns_pts,
  dsn = paste0(output_dir, "ai_", current_year, "_station_allocation_",
               n_allocation, "stn.gpkg")
)


merge(x = ai_hauls |>
        subset(select = c(STATIONID, CRUISE, STRATUM, PERFORMANCE,
                          HAUL_TYPE, ACCESSORIES, ABUNDANCE_HAUL)),
      by.x = c( "STRATUM", "STATIONID"),
      y = picked_stns |> subset(subset = STATION_TYPE == "assigned"),
      by.y = c( "STRATUM", "STATION")) |>
  merge(y = ai_grid_gis |>
          subset(subset = is.na(x = TRAWLABLE),
                 select = c(STATION, STRATUM, TRAWLABLE)),
        by.x = c( "STRATUM", "STATIONID"),
        by.y = c( "STRATUM", "STATION"))


merge(x = ai_hauls |>
        subset(select = c(STATIONID, CRUISE, STRATUM, PERFORMANCE,
                          HAUL_TYPE, ACCESSORIES, ABUNDANCE_HAUL)),
      by.x = c( "STRATUM", "STATIONID"),
      y = picked_stns |> subset(subset = STATION_TYPE == "assigned"),
      by.y = c( "STRATUM", "STATION")) |>
  merge(y = ai_grid_gis |>
          subset(subset = is.na(x = TRAWLABLE),
                 select = c(STATION, STRATUM, TRAWLABLE)),
        by.x = c( "STRATUM", "STATIONID"),
        by.y = c( "STRATUM", "STATION")) |> subset(select = c( "STRATUM", "STATIONID", "VESSEL")) |> unique()
