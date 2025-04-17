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
# devtools::install_github(repo = "afsc-gap-products/StationAllocationAIGOA")
library(StationAllocationAIGOA)
library(terra)
library(akgfmaps)
library(RColorBrewer)
library(xlsx)

laning_area <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_laning_area.shp")

ca_border <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/CanadaBorder.shp")

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

output_dir <- "G:/GOA/GOA 2025/Station Allocation/"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate a 520 station allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shallow_boat <- 176 # Alaska Provider
deep_boat <- 148    # Ocean Explorer
current_year <- 2025
total_n <- 450

goa_stn_allocation <- StationAllocationAIGOA::goa_allocate_stations(
  n = total_n,
  min_n_per_stratum = 4,
  survey_year = current_year
)

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
                y = laning_area,
                relation = "intersects",
                pairs = T)
strata_in_lane <- unique(stn_allocation[stns_in_lane[, "id.x"],]$STRATUM)
strata_out_lane <- unique(stn_allocation[-stns_in_lane[, "id.x"],]$STRATUM)

stn_allocation[stns_in_lane[stns_in_lane[, "id.y"] == 1, "id.x" ] ,
               "VESSEL"] <- deep_boat
stn_allocation[stns_in_lane[stns_in_lane[, "id.y"] == 2, "id.x" ] ,
               "VESSEL"] <- shallow_boat
stn_allocation[stns_in_lane[, "id.x"],
               "LANE"] <- T
table(stn_allocation$VESSEL)

plot(y ~ x, data = stn_allocation, pch = 16,
     col = c("148" = "black", "176" = "red")[paste(stn_allocation$VESSEL)])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Tabulate how many stations need to be rebalanced across strata so
##   that the vessels have the same number of assigned stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_stns_to_rebalance <- abs(ceiling(diff(table(stn_allocation$VESSEL)) / 2))
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
lesser_boat <- as.numeric(x = names(x = which.min(x = table(stn_allocation$VESSEL))))
greater_boat <- as.numeric(x = names(x = which.max(x = table(stn_allocation$VESSEL))))

for (istratum in 1:length(x = strata_to_switch)) {
  stn_to_switch <- sample(
    x = which(stn_allocation$VESSEL == greater_boat &
                stn_allocation$STRATUM == names(strata_to_switch)[istratum] &
                !stn_allocation$LANE),
    size = strata_to_switch[istratum]
  )

  stn_allocation$VESSEL[stn_to_switch] <- lesser_boat
}

table(stn_allocation$VESSEL) ## Should be equal between vessels
table(stn_allocation$STRATUM, stn_allocation$VESSEL)

# {
#   ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##  Assign stations to vessels w/o any laning first
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.seed(seed = current_year)
# ## Attach easting and northings to the allocated stations
# stn_allocation <-
#   merge(x = stn_allocation,
#         y = goa_stations[, c("STATION", "x", "y", "LONGITUDE", "LATITUDE")],
#         by = c("STATION"))
# stn_allocation$VESSEL <- NA
# stn_allocation$LANE <- F
#
# for (istratum in 1:n_strata) { ## Loop over strata -- start
#   temp_stratum <- goa_stn_allocation$ms_allocation$stratum[istratum]
#   nh <- goa_stn_allocation$ms_allocation$ms_allocation[istratum]
#   stn_allocation$VESSEL[stn_allocation$STRATUM == temp_stratum] <-
#     sample(c(shallow_boat, deep_boat))[(1:nh)%%2 + 1]
# }
#
# table(stn_allocation$VESSEL)
# table(stn_allocation$STRATUM, stn_allocation$VESSEL)
#
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##  Identify stations in the lane above Kodiak, assign to the shallow boat
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# stns_in_lane <-
#   terra::relate(x = terra::vect(x = stn_allocation,
#                                 geom = c("x", "y")),
#                 y = laning_area[laning_area$Name == "N_of_Kodiak"],
#                 relation = "intersects",
#                 pairs = T)
# strata_in_lane <- unique(stn_allocation[stns_in_lane[, "id.x"],]$STRATUM)
# strata_out_lane <- unique(stn_allocation[-stns_in_lane[, "id.x"],]$STRATUM)
#
# stn_allocation[stns_in_lane[, "id.x"], "VESSEL"] <- shallow_boat
# stn_allocation[stns_in_lane[, "id.x"], "LANE"] <- T
# table(stn_allocation$VESSEL)
#
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##   Tabulate how many stations need to be rebalanced across strata so
# ##   that the vessels have the same number of assigned stations
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# n_stns_to_rebalance <- ceiling(diff(table(stn_allocation$VESSEL)) / 2)
# stn_allocation_nonlane <- table(stn_allocation$STRATUM[!stn_allocation$LANE])
# stn_allocation_nonlane <- stn_allocation_nonlane[stn_allocation_nonlane > 4]
#
# ## Randomly choose strata to rebalance stations, with probabilities
# ## proportional to the station allocation.
# strata_to_switch <- table(sample(x = names(stn_allocation_nonlane),
#                                  size = n_stns_to_rebalance, replace = TRUE,
#                                  prob = stn_allocation_nonlane))
#
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##   Rebalance stations between vessels
# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for (istratum in 1:length(x = strata_to_switch)) {
#   stn_to_switch <- sample(
#     x = which(stn_allocation$VESSEL == shallow_boat &
#                 stn_allocation$STRATUM == names(strata_to_switch)[istratum] &
#                 !stn_allocation$LANE),
#     size = strata_to_switch[istratum]
#   )
#
#   stn_allocation$VESSEL[stn_to_switch] <- deep_boat
# }
#
# table(stn_allocation$VESSEL) ## Should be equal between vessels
# table(stn_allocation$STRATUM, stn_allocation$VESSEL)
# }

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Replace 2025 stations currently in the disputed Canadian Border
##   From Stratum 352: 739-76-352 749-83-352
##   From Stratum 252: 743-80-252
##   Currently we don't have shapefile to filter out station in disputed
##   Canadian waters
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ca_stns <- terra::intersect(x = goa_stations[goa_stations$STATION %in%
                                    stn_allocation$STATION, ],
                 y = ca_border)

ca_stns <- data.frame(
  CA_STATION = ca_stns$STATION,
  STRATUM = ca_stns$STRATUM,
  ALT_STATION = sample(x = goa_stations$STATION[
    goa_stations$STRATUM == 352 &
      !goa_stations$STATION %in%
      terra::intersect(goa_stations, ca_border)$STATION &
      goa_stations$AREA_KM2 >= 5 &
      goa_stations$TRAWLABLE != "N"
  ],
  size = 1),
  VESSEL = c(176),
  stringsAsFactors = FALSE)

updated_ca_stns <- cbind(as.data.frame(
  x = goa_stations[match(x = ca_stns$ALT_STATION,
                         table = goa_stations$STATION),
                   c("STATION", "GRIDID", "NMFS_AREA", "REP_AREA",
                     "STRATUM", "TRAWLABLE", "AREA_KM2", "x", "y",
                     "LONGITUDE", "LATITUDE")]),
  VESSEL = ca_stns$VESSEL, LANE = F)

stn_allocation <-
  rbind(stn_allocation[-match(x = ca_stns$CA_STATION,
                              table = stn_allocation$STATION), ],
        updated_ca_stns)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot map of drawn stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = paste0(output_dir, "goa_2025_station_allocation_",
                  total_n ,"_map.pdf"),
    width = 8, height = 6, onefile = TRUE)
for (iarea in c(610, 620, 630, 640, 650)) { ## Loop over area -- start

  ## temporary objects
  temp_strata <- goa_strata[goa_strata$REP_AREA == iarea]
  n_strata <- nrow(x = temp_strata)
  stratum_cols <- c(RColorBrewer::brewer.pal(name = "Set1",
                                             n = n_strata - 1), "cyan")

  ## Plot blank
  plot(temp_strata, border = F, axes = F, col = "white", lwd = 0.1)

  ## add land
  plot(goa_base$akland, add = TRUE, col = "tan", border = T, lwd = 0.1)

  ## add strata
  plot(temp_strata, border = TRUE, axes = F, col = stratum_cols,
       lwd = 0.1, add = TRUE)

  ## Add outlines of stations
  plot(goa_stations[goa_stations$REP_AREA == iarea],
       border = T, lwd = 0.05, add = TRUE )

  ## Add chosen stations
  with(stn_allocation[stn_allocation$REP_AREA == iarea, ],
       points(x, y, pch = c("176" = 15, "148" = 16)[paste(VESSEL)], cex = 0.3)
  )

  ## Add untrawlable areas
  plot(goa_stations[goa_stations$REP_AREA == iarea
                    & goa_stations$TRAWLABLE == "N"],
       border = FALSE, col = "grey", add = TRUE )

  ## Legend
  legend_labels <- with(as.data.frame(temp_strata),
                        paste0("Stratum ", STRATUM, ": ",
                               DEPTH_MIN_M , " - ", DEPTH_MAX_M, " m"))

  legend(c("610" = "topleft", "620" = "topleft", "630" = "bottomright",
           "640" = "bottom", "650" = "bottomleft")[paste(iarea)],
         legend = c(legend_labels, "Untrawlable"),
         title = "Stratum Legend",
         fill = c(stratum_cols, "grey"),
         xpd = NA)

  legend(c("610" = "bottom", "620" = "bottom", "630" = "bottom",
           "640" = "bottomleft", "650" = "bottom")[paste(iarea)],
         legend = c(176, 148), pch = 15:16, ncol = 2, xpd = NA)

  plot(ca_border, add = TRUE)

  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}  ## Loop over area -- end
dev.off()

## Write drawn stations to geopackage
stn_centroids_aea <-
  terra::centroids(x = merge(x = goa_stations,
                             y = stn_allocation[, c("STATION", "VESSEL")],
                             by = "STATION"),
                   inside = TRUE)

writeVector(
  x = stn_centroids_aea,
  file = paste0(output_dir, "goa_2025_station_allocation_",
                total_n, "_aea.gpkg"),
  overwrite = file.exists(paste0(output_dir, "goa_2025_station_allocation_",
                                 total_n, "_aea.gpkg"))
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate priority strata for bonus stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_allocation_bonus <- StationAllocationAIGOA::goa_allocate_stations(
  n = total_n + 30,
  min_n_per_stratum = 4,
  survey_year = current_year
)

bonus_stn_by_stratum <- goa_allocation_bonus$ms_allocation$ms_allocation -
  goa_stn_allocation$ms_allocation$ms_allocation
names(x = bonus_stn_by_stratum) <- goa_stn_allocation$ms_allocation$stratum
bonus_stn_by_stratum <- bonus_stn_by_stratum[bonus_stn_by_stratum > 0]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Format drawn stations for final output
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_drawn_stations <- rbind(
  data.frame(stn_allocation[, c("GRIDID", "STRATUM", "STATION", "TRAWLABLE",
                                "VESSEL", "LONGITUDE", "LATITUDE")],
             STATION_TYPE = "prescribed"),
  data.frame(GRIDID = NA,
             STRATUM = rep(as.numeric(names(x = bonus_stn_by_stratum)),
                           bonus_stn_by_stratum),
             STATION = NA,
             TRAWLABLE = NA,
             VESSEL = NA,
             LONGITUDE = NA,
             LATITUDE = NA,
             STATION_TYPE = "bonus_stn")
)

goa_drawn_stations <- goa_drawn_stations[order(goa_drawn_stations$STRATUM), ]

## Create excel file and append station allocation
xlsx::write.xlsx(
  x = goa_drawn_stations,
  file = paste0(output_dir, "goa_2025_station_allocation_", total_n, ".xlsx"),
  sheetName = "Station Allocation",
  row.names = FALSE,
  showNA = FALSE
)

## Append the station allocation per stratum and vessel
xlsx::write.xlsx(
  x = cbind(
    STRATUM = as.numeric(rownames(table(goa_drawn_stations$STRATUM,
                                        goa_drawn_stations$VESSEL))),
    table(goa_drawn_stations$STRATUM, goa_drawn_stations$VESSEL)
  ),
  file = paste0(output_dir, "goa_2025_station_allocation_", total_n, ".xlsx"),
  sheetName = "Stratum Allocation",
  row.names = FALSE, showNA = FALSE, append = TRUE
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Append priority strata from which to drop stations to output
##  Loop through an internal of slighly lower total effort and rerun
##  station allocation. Calculate the difference in station effort for
##  each strata between the original allocation and the lowered allocation
##  and append to the contingency_table.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contingency_table <- subset(x = goa_stn_allocation$ms_allocation,
                            select = c(stratum, nmfs_area))

for (ieffort in seq(from = total_n - 10, to = total_n - 50, by = -10)) {
  temp_allocation <- StationAllocationAIGOA::goa_allocate_stations(
    n = ieffort,
    min_n_per_stratum = 4,
    survey_year = current_year)

  contingency_table[, paste0("ms_allocation_", ieffort)] <-
    goa_stn_allocation$ms_allocation$ms_allocation -
    temp_allocation$ms_allocation$ms_allocation
}

xlsx::write.xlsx(
  x = contingency_table,
  file = paste0(output_dir, "goa_2025_station_allocation_", total_n, ".xlsx"),
  sheetName = "Priority Strata for Stn Drops",
  row.names = FALSE, showNA = FALSE, append = TRUE
)

