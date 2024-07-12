##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Inventory stations by trawlability status
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

## Import Packages
library(terra)

## Import station trawlability info. This shapefile contains both the
## historical and 2025 stations
station_trawl_info <-
  terra::vect(x = "data/GOA/processed_shapefiles/station_trawl_info.shp")

## How many unique stations are in the historical vs 2025 stations
tapply(X = station_trawl_info$GOAGRID_ID,
       INDEX = station_trawl_info$DESIGN_YEA,
       FUN = function(x) length(unique(x)))

## Tabulate trawlability status for historical stations
with(as.data.frame(station_trawl_info[station_trawl_info$DESIGN_YEA == 1984, ]),
     table(TRAWLABLE))

stations_full_overlap <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 2025),
     as.numeric(x = names(x = which(x = table(GOAGRID_ID) == 1))))
length(stations_full_overlap)
# 20543 of the 24333 2025 stations have full overlap in trawlability status

table(subset(x = as.data.frame(x = station_trawl_info),
       subset = GOAGRID_ID %in% stations_full_overlap,
       select = "TRAWLABLE"))

stations_full_overlap <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 2025),
       as.numeric(x = names(x = which(x = table(GOAGRID_ID) == 1))))
length(stations_full_overlap)

stations_partial_overlap <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 2025),
       as.numeric(x = names(x = which(x = table(GOAGRID_ID) > 1))))
length(stations_partial_overlap)
# There are 3790 stations that are characterized by more than one trawlability
# status. Stations that are characterized by unknown and trawlable areas should
# not be a problem.

stations_partial_overlap_UT <-
  subset(x = as.data.frame(x = station_trawl_info),
         subset = GOAGRID_ID %in% stations_partial_overlap &
           TRAWLABLE == "N")$GOAGRID_ID
# stations_partial_overlap_UT <-

cells_partial_overlap_UT <-
  unique(subset(x = as.data.frame(x = station_trawl_info),
         subset = GOAGRID_ID %in% stations_partial_overlap &
           TRAWLABLE == "N")$STATIONID)

cells_partial_overlap_UT_polygon <-
  subset(x = as.data.frame(x = station_trawl_info),
         subset = STATIONID %in% cells_partial_overlap_UT &
           DESIGN_YEA == 2025)
cells_partial_overlap_UT_polygon$TRAWLABLE <-
  ifelse(test = cells_partial_overlap_UT_polygon$TRAWLABLE == "N",
         yes = "N", no = "UNK")

cells_partial_overlap_UT_polygon <-
  stats::aggregate(AREA_KM2 ~ TRAWLABLE + GOAGRID_ID + STATIONID,
                   data = cells_partial_overlap_UT_polygon,
                   FUN = sum)

stations_partial_overlap_UT_gt_5km2 <-
  cells_partial_overlap_UT_polygon$GOAGRID_ID[which(cells_partial_overlap_UT_polygon$TRAWLABLE == "UNK" & cells_partial_overlap_UT_polygon$AREA_KM2 > 5)]

stations_partial_overlap_UT_lt_5km2 <-
  cells_partial_overlap_UT_polygon$GOAGRID_ID[which(cells_partial_overlap_UT_polygon$TRAWLABLE == "UNK" & cells_partial_overlap_UT_polygon$AREA_KM2 < 5)]

subset(x = as.data.frame(station_trawl_info), GOAGRID_ID == 45094)
plot(station_trawl_info[station_trawl_info$STATIONID == "435-77" &
                          station_trawl_info$DESIGN_YEA == 2025, ])
subset(x = as.data.frame(station_trawl_info), STATIONID == "435-77" & DESIGN_YEA == 2025)
plot(station_trawl_info[station_trawl_info$GOAGRID_ID == 45094 &
                          station_trawl_info$TRAWLABLE == "N", ],
     add = T, col = "red")
plot(station_trawl_info[station_trawl_info$GOAGRID_ID == 45094 &
                          station_trawl_info$TRAWLABLE == "UNK", ],
     add = T, col = "lightgrey")

plot(station_trawl_info[station_trawl_info$GOAGRID_ID == 45095 &
                          station_trawl_info$TRAWLABLE == "N", ],
     add = T, col = "orange")
plot(station_trawl_info[station_trawl_info$GOAGRID_ID == 45095 &
                          station_trawl_info$TRAWLABLE == "UNK", ],
     add = T, col = "darkgrey")

## Scenario 1/2: historical UT station encompasses entire grid
##               new station also encompasses entire grid
##               NO CHANGE

## Of the cells in UT_cells, which historical stations encompassed a full grid cell?
one_station_grid21 <- names(which(table(stations_21$ID) == 1))
full_UT_cells <- UT_cells[UT_cells %in% one_station_grid21]
## Of the 1663 grid cells with UT area, 720 are fully untrawlable.

## Which new stations are contained within full_UT_cells?
nrow(stations_23[stations_23$STATIONID %in% full_UT_cells])
table(table(stations_23[stations_23$STATIONID %in% full_UT_cells]$STATIONID))
## Of the 727 fully untrawlable grid cells, there are 1070 new stations,
## 431 of which are a one-station-per-grid station, the rest are > 1 stations
## within the cell.

## Scenario 3/4: historical UT station only partially encompasses grid cell
##               new station either:
##               1) encompasses entire grid cell: partial trawlable area
##               2) partially encompass grid cell: partial overlap with
##                  trawlable and untrawlable area within grid cell

## Which UT stations only partially encompass grid cell
partial_UT_cells <-
  UT_cells[!UT_cells %in% one_station_grid21]
## Of the 1663 grid cells with UT area, 886 are partially untrawlable.

partial_UT_stations_23 <-
  stations_23[stations_23$STATIONID %in% partial_UT_cells, #&
                # stations_23$TRAWLABLE == "N",
                ]
nrow(partial_UT_stations_23)
## THere are 2045 new stations in these partially UT cells to account for

scenario_b <- partial_UT_stations_23[
  round(partial_UT_stations_23$TRAWLABLE_, 4) == 0
  ]$GOAGRID_ID
length(scenario_b)
## partial_a: of the 2045 new stations in partially UT cells, 320 stations
## fully overlap with the UT part of the grid cell

scenario_c <- partial_UT_stations_23[
  is.na(partial_UT_stations_23$TRAWLABLE)
  ]$GOAGRID_ID
length(scenario_c)
## partial_b: of the 2045 new stations in partially UT cells, 483 stations
## fully overlap with the trawlable part of the grid cell


scenario_d <- with(as.data.frame(partial_UT_stations_23[!partial_UT_stations_23$GOAGRID_ID %in% c(scenario_b, scenario_c)]),
     GOAGRID_ID[TRAWLABLE_ >= 5]
)
length(scenario_d)
scenario_e <- with(as.data.frame(partial_UT_stations_23[!partial_UT_stations_23$GOAGRID_ID %in% c(scenario_b, scenario_c)]),
                  GOAGRID_ID[TRAWLABLE_ < 5]
)
length(scenario_e)
## of the 2045 new stations in partially UT cells, 642 stations are < 5 km^2
## so even if they were not untrawlable, the station allocation protocol would
## still exclude these stations. 600 stations are >= 5 km^2 and could still be
## included in the staitonal allocation.

## Plots
for (ichunk in rev(c("scenario_b", "scenario_d", "scenario_e"))) {
  pdf(file = paste0("C:/Users/zack.oyafuso/Desktop/", ichunk, ".pdf"),
      width = 8, height = 10, onefile = TRUE)
  par(mfrow = c(5, 5), oma = c(2, 2, 2, 2))
  for (istation in get(ichunk)) {
    icell <- stations_23$STATIONID[stations_23$GOAGRID_ID == istation]
    plot(stations_23[stations_23$STATIONID == icell],
         mar = c(0,0,1,0),
         main = icell, axes = F)
    plot(stations_23[stations_23$GOAGRID_ID == istation],
         col = "red", add = TRUE)
    plot(stations_21UT[stations_21UT$ID == icell],
         add = TRUE, col = rgb(0,0,0,0.5))
  }
  dev.off()
}
