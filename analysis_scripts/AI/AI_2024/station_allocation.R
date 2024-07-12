##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Aleutian Island BTS Station Allocation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##
## Description:   Station Allocation Procedure for the 2024 Aleutian Islands
##                bottom trawl survey. Make sure you are connected to the VPN
##                or network to connect to Oracle and to access the price data
##                spreadsheet in the G: drive.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

options(scipen = 999999)
current_year <- 2024
set.seed(seed = current_year)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages, connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
devtools::build()

install.packages("C:/Users/zack.oyafuso/Work/GitHub/StationAllocationAIGOA_0.1.0.tar.gz",
                 repos = NULL,
                 type="source")

library(StationAllocationAIGOA)
library(xlsx)
library(gapindex)
channel <- gapindex::get_connected()

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
strata <- terra::vect(x = "data/AI/strata_polygon/ai_strata.shp")
strata <- data.frame(stratum = sort(x = unique(x = strata$STRATUM)),
                     area_km2 = tapply(X = strata$AREA / 1000 / 1000,
                                       INDEX = strata$STRATUM,
                                       FUN = sum),
                     perim_km = tapply(X = strata$PERIMETER / 1000,
                                       INDEX = strata$STRATUM,
                                       FUN = sum))
strata$edge_to_area <-
  with(strata, perim_km / 2 / sqrt(area_km2 * pi) )

strata$strata_type <-
  ifelse(test = strata$edge_to_area < mean(x = strata$edge_to_area) &
           strata$area_km2 > 1000,
         yes = "large",
         no = "thin")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import price data from COAR. Because there is a species complex, it is
##   easier to run gapindex to calculate mean and sd across strata for each
##   species/species complex.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
price_data <- readxl::read_excel(
  path = "G:/ALEUTIAN/AI 2024/Station Allocation/ai_exvessel_2022.xlsx")

## Create df
planning_species <-
  data.frame(
    GROUP_CODE = c(price_data$species_code[price_data$species_code != 10260],
                   10260, 10260, 10260),
    SPECIES_CODE = c(price_data$species_code[price_data$species_code != 10260],
                     10260, 10261, 10262)
  )

## Get survey from Oracle
planning_species_data <- gapindex::get_data(year_set = 1991:2022,
                                            survey_set = "AI",
                                            spp_codes = planning_species,
                                            sql_channel = channel)

## Calculate CPUE
planning_species_cpue <-
  gapindex::calc_cpue(racebase_tables = planning_species_data)

## Calculate CPUE and Biomass
planning_species_biomass <-
  gapindex::calc_biomass_stratum(racebase_tables = planning_species_data,
                                 cpue = planning_species_cpue)

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
                   data = planning_species_biomass,
                   FUN = mean)

## Merge price data to species_weighting using SPECIES_CODE as key
species_weightings <-
  merge(x = species_weightings,
        y = subset(x = price_data,
                   select = c(species_code, ex_vessel_price_per_mt)),
        by.x = "SPECIES_CODE", by.y = "species_code")

## Calculate total economic value by multiplying by ex-vessel price.
species_weightings$econ_value <- with(species_weightings,
                                      BIOMASS_MT * ex_vessel_price_per_mt)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Pull candidate AI stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
candidate_stations <- StationAllocationAIGOA::get.ai.stations(channel = channel)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate optimal allocation at 400 and 420 total stations
##   Minimum number of stations per stratum: 2 stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_320 <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = 321,
  min_n_stratum = 2)
sum(n_320)

n_340 <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = 340,
  min_n_stratum = 2)
sum(n_340)

n_341 <- StationAllocationAIGOA::allocate.effort(
  channel = channel,
  species_weightings = species_weightings,
  stratum_stats = stratum_stats,
  n_total = 341,
  min_n_stratum = 2)
sum(n_341)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Bind the station effort across strata under 400 and 420 total stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station_allocation <-
  data.frame(stratum = names(x = n_320),
             region = sapply(X = substr(x = names(x = n_320),
                                        start = 1,
                                        stop = 1),
                             FUN = function(x)
                               switch(x,
                                      "2" = "Western",
                                      "3" = "Central", "4" = "Central",
                                      "5" = "Eastern", "6" = "Eastern",
                                      "7" = "SBS")),
             avail = tapply(X = candidate_stations$stationid,
                            INDEX = candidate_stations$stratum,
                            FUN = length),
             n_320,
             n_340)

## Merge strata information to station_allocation using "stratum" as the key
station_allocation <- merge(x = station_allocation,
                            y = strata,
                            by = "stratum")

## Bonus stations: the difference between the allocation at 340 total stations
## versus the allocation at 320 total stations.
station_allocation$bonus_stn <- with(station_allocation, n_340 - n_320)
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
  station_allocation$n_320[chosen_large_stratum] <-
    station_allocation$n_320[chosen_large_stratum] - 1

  chosen_thin_stratum <-
    sample(x = which(station_allocation$region == iregion &
                       station_allocation$strata_type == "thin"),
           size = 1)
  station_allocation$new_stn[chosen_thin_stratum] <- 1
  station_allocation$n_320[chosen_thin_stratum] <-
    station_allocation$n_320[chosen_thin_stratum] - 1
} ## Loop over regions -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Randomly draw stations based on the allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station_allocation <-
  with(station_allocation,
       data.frame(stratum, region,
                  prescribed_320 = n_320,
                  new_stn,
                  bonus_stn))

picked_stns <- data.frame()

for (istratum in 1:length(x = station_allocation$stratum))
{ ## Loop over strata -- start

  ## Append randomly drawn prescribed stations
  picked_stns <-
    rbind(
      picked_stns,
      data.frame( candidate_stations[
        sample(x = which(x = candidate_stations$stratum ==
                           station_allocation$stratum[istratum]),
               size = station_allocation$prescribed_320[istratum]),
      ],
      station_type = "prescribed")
    )

  ## If there are new stations in the stratum ("new_stn"), append blank records
  ## that just have the stratum information
  if (station_allocation$new_stn[istratum] > 0)
    picked_stns <-
      rbind(picked_stns,
            data.frame(stationid = NA,
                       stratum = rep(station_allocation$stratum[istratum],
                                     station_allocation$new_stn[istratum]),
                       longitude = NA,
                       latitude = NA,
                       station_type = "new_stn") )

  ## If there are bonus stations in the stratum ("bonus_stn"), append blank
  ## records that just have the stratum information
  if (station_allocation$bonus_stn[istratum] > 0)
    picked_stns <-
      rbind(picked_stns,
            data.frame(stationid = NA,
                       stratum = rep(station_allocation$stratum[istratum],
                                     station_allocation$bonus_stn[istratum]),
                       longitude = NA,
                       latitude = NA,
                       station_type = "bonus_stn"))

} ## Loop over strata -- end
names(x = picked_stns) <- toupper(x = names(x = picked_stns))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Write Results to G: Drive
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xlsx::write.xlsx(
  x = subset(x = picked_stns,
             select = c(STRATUM, STATION_TYPE, STATIONID, LONGITUDE, LATITUDE)),
  file = paste0("G:/ALEUTIAN/AI ", current_year, "/Station Allocation/",
                "ai_", current_year, "_station_allocation_320stn.xlsx"),
  row.names = F, showNA = F
)
