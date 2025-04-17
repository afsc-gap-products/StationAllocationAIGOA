##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Calculate Total distance travelled and duration of
##                allocations at different effort levels per vessel
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(terra)
library(gapindex)
library(tidyr)
library(navmaps)
library(StationAllocationAIGOA)

#' Solve a traveling Salesman problem (TSP) for a survey station alaocation
#'
#' This function optimizes the sampling order of stations by solving a Traveling Salesman Problem (TSP) using the nearest insertion method to minimize the distance traveled.
#'
#' @param x An `sf` object representing the survey stations.
#'
#' @return A list containing:
#' \item{distance_nodes}{An `sf` object of survey stations ordered based on
#' the optimal TSP solution, with an added column for inter-station distances in kilometers.}
#' \item{tsp}{The TSP solution object.}
#'
#' @importFrom sf st_distance st_transform
#' @importFrom TSP as.ATSP TSP solve_TSP
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' # Load and transform survey stations
#' x <- system.file("extdata", "goa_station_allocation_520.shp", package = "navmaps") |>
#' sf::st_read() |>
#'   sf::st_transform(crs = "EPSG:32606") |> # UTM zone 2
#'   dplyr::filter(VESSEL == 148) # Ocean Explorer
#'
#' # Solve TSP for station order
#' tsp_out <- planning_solve_station_tsp(x = x)
#'
#' # Estimate sampling days
#' survey_days <- planning_calc_survey_days(
#'   station_nodes = tsp_out$distance_nodes,
#'   max_daily_hr = 12,
#'   processing_time_hr = 1.5,
#'   max_daily_stn = 6,
#'   transit_speed_kmh = 1.852*7, # Converted from knots to km/h
#'   set_retrieve_hr = 0.5,
#'   set_on_arrival = FALSE
#' )
#'
#' print(survey_days)
#' }
#'
#' @export
#'

planning_solve_station_tsp <-
  function(x, hamilton = TRUE) {

    tsp_data <- x |>
      sf::st_distance() |>
      as.matrix() |>
      matrix(
        nrow = nrow(x),
        ncol = nrow(x)
      ) |>
      TSP::TSP()

    if (hamilton) {
      # Insert dummy city to hack the solver to calculate the hamiltonian path
      tsp_data <- TSP::insert_dummy(tsp_data, label = "cut")
    }

    tsp_solution <- TSP::solve_TSP(x = tsp_data,
                                   method = "nearest_insertion")

    if (hamilton) {
      # Insert dummy city to hack the solver to calculate the hamiltonian path
      tsp_solution <- TSP::cut_tour(tsp_solution, "cut")
    }

    x$node <- as.numeric(attr(tsp_solution, "names"))
    x <- x[as.numeric(attr(tsp_solution, "names")), ]
    x$order <- 1:nrow(x)
    x$distance <- c(0,
                    sf::st_distance(x = x[1:(nrow(x)-1), ],
                                    y = x[2:nrow(x), ],
                                    by_element = TRUE)
    )

    x$distance <- as.numeric(x$distance / 1000)

    return(list(distance_nodes = x,
                tsp = tsp_solution)
    )

  }

channel <- gapindex::get_connected(check_access = FALSE)

inpfc <- RODBC::sqlQuery(
  channel = channel,
  query = "select STRATUM, AREA_ID as INFPC from gap_products.stratum_groups
 where survey_definition_id = 47
and area_id IN (919, 929, 939, 949, 959)")
nmfs <- RODBC::sqlQuery(
  channel = channel,
  query = "select STRATUM, AREA_ID as NMFS from gap_products.stratum_groups
 where survey_definition_id = 47
and area_id IN (610, 620, 630, 640, 650)")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import polygons used in the allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
laning_area <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_laning_area.shp") |>
  terra::project("EPSG:3338")


goa_stations <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_stations_2025.gpkg")
goa_stations[, c("x", "y")] <-
  terra::crds(x = terra::centroids(x = goa_stations,
                                   inside = TRUE))
goa_stations[, c("LONGITUDE", "LATITUDE")] <-
  terra::crds(x = terra::project(terra::centroids(x = goa_stations,
                                                  inside = TRUE),
                                 "EPSG:4326"))

goa_strata <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_strata_2025.gpkg")

## `goa_base` are basic shape layers from the akgfmaps package
goa_base <- akgfmaps::get_base_layers(select.region = "goa",
                                      set.crs = "EPSG:3338")

shallow_boat <- 176 # Alaska Provider
deep_boat <- 148    # Ocean Explorer

# Setup parameters
max_daily_hr = 12 # Maximum hours worked in a day
processing_time_hr = 1.25 # Minimum processing time in hours
max_daily_stn = 6 # Maximum number of stations sampled in a day
# transit_speed_kmh = 1.852*9 # Transit speed between stations in kilometers/hour
set_retrieve_hr = 0.5 # Time to set and retrieve the gear in hours
set_on_arrival = FALSE # Set on arrival at the station or set based on the minimum processing time?

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_2023 <-
  sf::st_read(dsn = "G:/GOA/GOA 2023/Station Allocation/Stations_2023.shp")
goa_2023$VESSEL <- ifelse(goa_2023$vessel == "AKP", 176, 148)

goa_2021 <-
  sf::st_read(dsn = "G:/GOA/GOA 2021/Station allocation/GOA_2021_Stns.shp")
goa_2021$VESSEL <- ifelse(goa_2021$vessel == "OEX", 148, 176)

hist_2023_res <- hist_2023_distance <- data.frame()

for (iyear in c(2021, 2023)) {
  for (ivessel in c(shallow_boat, deep_boat)) {
    vessel_tsp <- sf::st_as_sf(x = subset(x = get(paste0("goa_", iyear)),
                                          subset = VESSEL == ivessel),
                               coords = c("x", "y"),
                               crs = "EPGS:3338") |>
      planning_solve_station_tsp()

    vessel_tsp$distance_nodes <-
      merge(x = vessel_tsp$distance_nodes, by.x = "stratum",
            y = inpfc, by.y = "STRATUM")

    median_distance_by_inpfc <-
      tapply(X = vessel_tsp$distance_nodes$distance,
             INDEX = vessel_tsp$distance_nodes$INFPC,
             FUN = median)

    mean_distance_by_inpfc <-
      tapply(X = vessel_tsp$distance_nodes$distance,
             INDEX = vessel_tsp$distance_nodes$INFPC,
             FUN = mean)

    hist_2023_distance <-
      rbind(hist_2023_distance,
            data.frame(year = iyear,
                       region = c(919,929,939,949,959,99903),
                       vessel = ivessel,
                       mean_dist = c(mean_distance_by_inpfc,
                                     mean(vessel_tsp$distance_nodes$distance)),
                       median_dist = c(median_distance_by_inpfc,
                                       median(vessel_tsp$distance_nodes$distance))))
    station_cluster = vessel_tsp$distance_nodes

    for (ispeed in c(8, 8.5, 9)) {
      vessel_days <-
        planning_calc_survey_days(
          station_nodes = vessel_tsp$distance_nodes,
          max_daily_hr = max_daily_hr,
          processing_time_hr = processing_time_hr,
          max_daily_stn = max_daily_stn,
          transit_speed_kmh = 1.852*ispeed,
          set_retrieve_hr = set_retrieve_hr,
          set_on_arrival = TRUE
        )

      vessel_days$station_nodes |>
        dplyr::group_by(day) |>
        dplyr::summarise(n_hauls = n()) |>
        dplyr::group_by(n_hauls) |>
        dplyr::summarise(n = n())

      total_duration <- vessel_days$total_days

      hist_2023_res  <-
        rbind(hist_2023_res,
              data.frame(year = iyear,
                         vessel = ivessel,
                         total_n =  nrow(vessel_days$station_nodes),
                         speed = ispeed,
                         total_days = total_duration,
                         stn_rate = nrow(vessel_days$station_nodes) /
                           total_duration))
    }
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

result_df <- result_distance <- data.frame()

stn_per_day <-
  array(data = 0,
        dim = c(length(seq(from = 520, to = 350, by = -10)), 6, 3, 2),
         dimnames = list(seq(from = 520, to = 350, by = -10),
                         NULL,
                         c(8, 8.5, 9),
                         c(148, 176)))
expected_cv <-
  matrix(nrow = length(seq(from = 520, to = 350, by = -10)),
         ncol = 15,
         dimnames = list(seq(from = 520, to = 350, by = -10),
                         c("arrowtooth flounder", "Pacific cod", "walleye pollock", "rex sole",
                           "flathead sole", "Pacific halibut", "southern rock sole", "northern rock sole",
                           "Pacific ocean perch", "silvergray rockfish", "northern rockfish", "dusky rockfish",
                           "REBS rockfish", "Dover sole", "shortspine thornyhead")))

for (ieffort in seq(from = 520, to = 350, by = -10)) {
  current_year <- 2025
  set.seed(seed = current_year)

  goa_stn_allocation <-
    StationAllocationAIGOA::goa_allocate_stations(n = ieffort)
  stn_allocation <- goa_stn_allocation$drawn_stations
  n_strata <- nrow(x = goa_stn_allocation$ms_allocation)

  expected_cv[paste(ieffort), ] <-
    goa_stn_allocation$expected_cv$expected_ms_cv

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Assign stations to vessels w/o any laning first
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Attach easting and northings to the allocated stations
  stn_allocation <-
    merge(x = stn_allocation,
          y = goa_stations[, c("STATION", "x", "y", "LONGITUDE", "LATITUDE")],
          by = c("STATION"))
  stn_allocation$VESSEL <- NA
  stn_allocation$LANE <- F

  for (istratum in 1:n_strata) { ## Loop over strata -- start
    temp_stratum <- goa_stn_allocation$ms_allocation$stratum[istratum]
    nh <- goa_stn_allocation$ms_allocation$ms_allocation[istratum]
    stn_allocation$VESSEL[stn_allocation$STRATUM == temp_stratum] <-
      sample(c(shallow_boat, deep_boat))[(1:nh)%%2 + 1]
  }

  table(stn_allocation$VESSEL)
  table(stn_allocation$STRATUM, stn_allocation$VESSEL)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Identify stations in the lane above Kodiak, assign to the shallow boat
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stns_in_lane <-
    terra::relate(x = terra::vect(x = stn_allocation,
                                  geom = c("x", "y")),
                  y = laning_area[laning_area$Name == "N_of_Kodiak"],
                  relation = "intersects",
                  pairs = T)
  strata_in_lane <- unique(stn_allocation[stns_in_lane[, "id.x"],]$STRATUM)
  strata_out_lane <- unique(stn_allocation[-stns_in_lane[, "id.x"],]$STRATUM)

  stn_allocation[stns_in_lane[, "id.x"], "VESSEL"] <- shallow_boat
  stn_allocation[stns_in_lane[, "id.x"], "LANE"] <- T
  table(stn_allocation$VESSEL)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Tabulate how many stations need to be rebalanced across strata so
  ##   that the vessels have the same number of assigned stations
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_stns_to_rebalance <- ceiling(diff(table(stn_allocation$VESSEL)) / 2)
  stn_allocation_nonlane <- table(stn_allocation$STRATUM[!stn_allocation$LANE])
  stn_allocation_nonlane <- stn_allocation_nonlane[stn_allocation_nonlane > 4]

  ## Randomly choose strata to rebalance stations, with probabilities
  ## proportional to the station allocation.
  strata_to_switch <- table(sample(x = names(stn_allocation_nonlane),
                                   size = n_stns_to_rebalance, replace = TRUE,
                                   prob = stn_allocation_nonlane))

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Rebalance stations between vessels
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (istratum in 1:length(x = strata_to_switch)) {
    stn_to_switch <- sample(
      x = which(stn_allocation$VESSEL == shallow_boat &
                  stn_allocation$STRATUM == names(strata_to_switch)[istratum] &
                  !stn_allocation$LANE),
      size = strata_to_switch[istratum]
    )

    stn_allocation$VESSEL[stn_to_switch] <- deep_boat
  }

  table(stn_allocation$VESSEL) ## Should be equal between vessels
  table(stn_allocation$STRATUM, stn_allocation$VESSEL)

  for (ivessel in c(shallow_boat, deep_boat)) {
    vessel_tsp <- sf::st_as_sf(x = subset(x = stn_allocation,
                                          subset = VESSEL == ivessel),
                               coords = c("x", "y"),
                               crs = "EPGS:3338") |>
      planning_solve_station_tsp()

    result_distance <-
      rbind(
        result_distance,
        data.frame(
          vessel = ivessel,
          total_n =  nrow(vessel_tsp$distance_nodes),
          region = c(610, 620, 630, 640, 650, 99903),
          mean_dist = c(tapply(X = vessel_tsp$distance_nodes$distance,
                               INDEX = vessel_tsp$distance_nodes$REP_AREA,
                               FUN = mean),
                        mean(vessel_tsp$distance_nodes$distance)),
          median_dist = c(tapply(X = vessel_tsp$distance_nodes$distance,
                                 INDEX = vessel_tsp$distance_nodes$REP_AREA,
                                 FUN = median),
                          median(vessel_tsp$distance_nodes$distance))
        )
      )

    station_cluster = vessel_tsp$distance_nodes

    for (ispeed in c(8, 8.5, 9)) {
      vessel_days <-
        planning_calc_survey_days(
          station_nodes = vessel_tsp$distance_nodes,
          max_daily_hr = max_daily_hr,
          processing_time_hr = processing_time_hr,
          max_daily_stn = max_daily_stn,
          transit_speed_kmh = 1.852*ispeed,
          set_retrieve_hr = set_retrieve_hr,
          set_on_arrival = TRUE
        )

      stn_tab <- vessel_days$station_nodes |>
        dplyr::group_by(day) |>
        dplyr::summarise(n_hauls = n()) |>
        dplyr::group_by(n_hauls) |>
        dplyr::summarise(n = n())

      stn_per_day[paste(ieffort),
                  stn_tab$n_hauls,
                  paste(ispeed),
                  paste(ivessel)] <-
        stn_tab$n

      (total_duration <- vessel_days$total_days)

      result_df  <- rbind(result_df,
                          data.frame(vessel = ivessel,
                                     total_n = ieffort,
                                     speed = ispeed,
                                     total_days = total_duration,
                                     ideal_rate = ieffort / 2 / total_duration))
    }
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
matplot(rownames(expected_cv), expected_cv[, order(expected_cv[1, ],
                                                 decreasing = FALSE)],
        las = 1, col = "black",
        xlab = c("Total Effort"), ylab = "Expected CV")
abline(h = c(0.2, 0.25, 0.3), lty = c(1, 2, 1))

par(mar = c(1, 1, 1, 1))
plot(1, type = "n", axes = F, ann = F)
legend("center",
       legend = colnames(expected_cv)[order(expected_cv[1, ],
                                            decreasing = TRUE)],
       bty = "n",
       pch = rev(c(1:9, 0, letters[1:5])), cex = 1.5)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


par(mfcol = c(2, 2), mar = c(5, 5, 1, 1))
for (ivessel in c(shallow_boat, deep_boat)) {
  plot(1,
       type = "n",
       xlab = "Vessel Effort", ylab = "Total Duration", main = ivessel,
       xlim = c(350, 520) / 2,
       ylim = c(43, 65), las = 1)

  for (ispeed in c(8, 8.5, 9)) {
    lines(total_days ~ I(total_n/2),
          data = result_df,
          subset = vessel == ivessel & speed == ispeed,
          col = c("8" = "black", "8.5" = "red", "9" = "blue")[paste(ispeed)])
    points(total_days ~ I(total_n/2),
           data = result_df,
           subset = vessel == ivessel & speed == ispeed,
           pch = 16,
           col = c("8" = "black", "8.5" = "red", "9" = "blue")[paste(ispeed)])
  }

  abline(h = seq(0,65,5), col = "grey", lty = "dashed")
  legend("bottomright",
         legend = paste(c(8, 8.5, 9), "knots"),
         col = c("black", "red", "blue"),
         pch = 16, lty = 1)

  # points(x = hist_2023_res$total_n[hist_2023_res$vessel == ivessel &
  #                                    hist_2023_res$year == 2023],
  #        y = hist_2023_res$total_days[hist_2023_res$vessel == ivessel &
  #                                       hist_2023_res$year == 2023],
  #        pch = 1, col = c("black", "red", "blue"), cex = 2)
  #
  # points(x = hist_2023_res$total_n[hist_2023_res$vessel == ivessel &
  #                                    hist_2023_res$year == 2021],
  #        y = hist_2023_res$total_days[hist_2023_res$vessel == ivessel &
  #                                       hist_2023_res$year == 2021],
  #        pch = 1, col = c("black", "red", "blue"), cex = 2)
}


for (ivessel in c(shallow_boat, deep_boat)) {
  plot(1,
       type = "n",
       xlab = "Vessel Effort", ylab = "Average Station Rate",
       main = ivessel,
       xlim = c(170, 260),
       ylim = c(3, 4.75), las = 1)
  abline(h = seq(2,5,0.25), col = "grey", lty = "dashed")

  for (ispeed in c(8, 8.5, 9)) {
    lines(ideal_rate ~ I(total_n/2),
          data = result_df,
          subset = vessel == ivessel & speed == ispeed,
          col = c("8" = "black", "8.5" = "red", "9" = "blue")[paste(ispeed)])
    points(ideal_rate ~ I(total_n/2),
           data = result_df,
           subset = vessel == ivessel & speed == ispeed,
           pch = 16,
           col = c("8" = "black", "8.5" = "red", "9" = "blue")[paste(ispeed)])
  }
  legend("bottomright",
         legend = paste(c(9, 8.5, 8), "knots"),
         col = c("blue", "red", "black"),
         pch = 16, lty = 1)
}

plot(1,
     type = "n",
     xlab = "Vessel Effort", ylab = "Mean Distance",
     main = ivessel,
     xlim = c(170, 265),
     ylim = c(27, 40), las = 1)
abline(h = seq(2,5,0.25), col = "grey", lty = "dashed")
legend("topright",
       legend = paste(c(148, 176)),
       col = c("black", "red"),
       pch = 16, lty = 1)

for (ivessel in c(shallow_boat, deep_boat)) {
  lines(mean_dist ~ total_n,
        data = result_distance,
        subset = vessel == ivessel & region  == 99903,
        col = c("148" = "black", "176" = "red")[paste(ivessel)])
  points(mean_dist ~ total_n,
         data = result_distance,
         subset = vessel == ivessel & region  == 99903,
         pch = 16,
         col = c("148" = "black", "176" = "red")[paste(ivessel)])
}


plot(1,
     type = "n",
     xlab = "Vessel Effort", ylab = "Median Distance",
     main = ivessel,
     xlim = c(170, 260),
     ylim = c(22, 36), las = 1)
abline(h = seq(2,5,0.25), col = "grey", lty = "dashed")
legend("topright",
       legend = paste(c(148, 176)),
       col = c("black", "red"),
       pch = 16, lty = 1)

for (ivessel in c(shallow_boat, deep_boat)) {
  lines(median_dist ~ total_n,
        data = result_distance,
        subset = vessel == ivessel & region == 99903,
        col = c("148" = "black", "176" = "red")[paste(ivessel)])
  points(median_dist ~ total_n,
         data = result_distance,
         subset = vessel == ivessel & region == 99903,
         pch = 16,
         col = c("148" = "black", "176" = "red")[paste(ivessel)])
}

points(median_dist ~ total_n,
       data = hist_2023_distance,
       subset = region == 99903,
       pch = 16, col = c("red", "black"), cex = 2)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stns_2009 <-
  terra::vect(x = "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/Near/GOA2009_825Stns/") |>
  merge(by.x = "stratum", y = inpfc, by.y = "STRATUM")
stns_2013 <-
  terra::vect(x = "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/Near/GOA2013_550Stns/") |>
  merge(by.x = "stratum", y = inpfc, by.y = "STRATUM")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages, laning areas around Kodiak
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github(repo = "afsc-gap-products/StationAllocationAIGOA")
library(StationAllocationAIGOA)
library(terra)
library(akgfmaps)
library(RColorBrewer)
library(xlsx)



proposed_stn_dist <- data.frame()

for (ieffort in c(550, 520, 400)) {
  current_year <- 2025
  set.seed(seed = current_year)

  goa_stn_allocation <- StationAllocationAIGOA::goa_allocate_stations(n = ieffort)
  stn_allocation <- goa_stn_allocation$drawn_stations
  n_strata <- nrow(x = goa_stn_allocation$ms_allocation)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Assign stations to vessels w/o any laning first
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Attach easting and northings to the allocated stations
  stn_allocation <-
    merge(x = stn_allocation,
          y = goa_stations[, c("STATION", "x", "y", "LONGITUDE", "LATITUDE")],
          by = c("STATION"))
  stn_allocation$VESSEL <- NA
  stn_allocation$LANE <- F

  for (istratum in 1:n_strata) { ## Loop over strata -- start
    temp_stratum <- goa_stn_allocation$ms_allocation$stratum[istratum]
    nh <- goa_stn_allocation$ms_allocation$ms_allocation[istratum]
    stn_allocation$VESSEL[stn_allocation$STRATUM == temp_stratum] <-
      sample(c(shallow_boat, deep_boat))[(1:nh)%%2 + 1]
  }

  table(stn_allocation$VESSEL)
  table(stn_allocation$STRATUM, stn_allocation$VESSEL)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Identify stations in the lane above Kodiak, assign to the shallow boat
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stns_in_lane <-
    terra::relate(x = terra::vect(x = stn_allocation,
                                  geom = c("x", "y")),
                  y = laning_area[laning_area$Name == "N_of_Kodiak"],
                  relation = "intersects",
                  pairs = T)
  strata_in_lane <- unique(stn_allocation[stns_in_lane[, "id.x"],]$STRATUM)
  strata_out_lane <- unique(stn_allocation[-stns_in_lane[, "id.x"],]$STRATUM)

  stn_allocation[stns_in_lane[, "id.x"], "VESSEL"] <- shallow_boat
  stn_allocation[stns_in_lane[, "id.x"], "LANE"] <- T
  table(stn_allocation$VESSEL)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Tabulate how many stations need to be rebalanced across strata so
  ##   that the vessels have the same number of assigned stations
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_stns_to_rebalance <- ceiling(diff(table(stn_allocation$VESSEL)) / 2)
  stn_allocation_nonlane <- table(stn_allocation$STRATUM[!stn_allocation$LANE])
  stn_allocation_nonlane <- stn_allocation_nonlane[stn_allocation_nonlane > 4]

  ## Randomly choose strata to rebalance stations, with probabilities
  ## proportional to the station allocation.
  strata_to_switch <- table(sample(x = names(stn_allocation_nonlane),
                                   size = n_stns_to_rebalance, replace = TRUE,
                                   prob = stn_allocation_nonlane))

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Rebalance stations between vessels
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (istratum in 1:length(x = strata_to_switch)) {
    stn_to_switch <- sample(
      x = which(stn_allocation$VESSEL == shallow_boat &
                  stn_allocation$STRATUM == names(strata_to_switch)[istratum] &
                  !stn_allocation$LANE),
      size = strata_to_switch[istratum]
    )

    stn_allocation$VESSEL[stn_to_switch] <- deep_boat
  }

  table(stn_allocation$VESSEL) ## Should be equal between vessels
  table(stn_allocation$STRATUM, stn_allocation$VESSEL)

  goa_stn_sp <-
    stn_allocation |>
    terra::vect(geom = c("LONGITUDE", "LATITUDE"),
                crs = "EPSG:4269") |>
    terra::project("EPSG:3338")

  plot(goa_stn_sp,
       col = c("148" = "red",
               "176" = "black")[paste(goa_stn_sp$VESSEL)] )

  goa_hauls_dist <-
    do.call(what = rbind,
            args = lapply(X = split(x = goa_stn_sp, f = goa_stn_sp$VESSEL),
                          FUN = function(df) {
                            df$distance <- terra::nearest(x = df)$distance * 1e-3
                            df <- merge(x = as.data.frame(df), y = nmfs)
                            return(df)
                          }))

  proposed_stn_dist <-
    rbind(proposed_stn_dist,
          data.frame(effort = ieffort,
                     aggregate(distance ~ VESSEL + NMFS,
                               FUN = function(df) round(x = summary(df),
                                                        digits = 2),
                               data = goa_hauls_dist)))
}


par(mfcol = c(1,2), mar = c(4, 4, 2, 1), oma = c(1, 1, 0, 0))
for (iarea in 1:5) {
  plot(distance ~ HAULJOIN,
       data = historical_data,
       main = c(919, 929, 939, 949, 959)[iarea],
       subset = INFPC == c(919, 929, 939, 949, 959)[iarea],
       pch = 16, cex = 2,
       las = 1, ylim = c(8, 30), xlim = c(2,5),
       col = c("148" = "red",
               "176" = "black",
               "143" = "blue")[paste(historical_data$VESSEL)],
       xlab = "Mean Station Rate", ylab = "Mean Min. Distance")

  legend("bottom", ncol = 3, legend = c(148, 176, 143),
         col = c("red", 'black', "blue"),
         pch = 16)

  plot(distance ~ effort,
       data = proposed_stn_dist,
       main = c(610,620,630,640,650)[iarea],
       subset = NMFS == c(610,620,630,640,650)[iarea], pch = 16, cex = 2,
       las = 1, ylim = c(8, 30),
       xlab = "Total Stations", ylab = "",
       col = c("148" = "red",
               "176" = "black")[paste(proposed_stn_dist$VESSEL)])

  legend("bottom", ncol = 2, legend = c(148, 176),
         col = c("red", 'black'),
         pch = 16)
}

