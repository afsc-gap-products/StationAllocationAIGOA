library(terra)

stations_23 <- terra::vect("data/GOA/processed_shapefiles/goa_stations_2023.shp")
stations_21UT <- terra::vect("data/GOA/processed_shapefiles/GOA_untrawl_2021.shp")
stations_21 <- terra::vect("data/GOA/shapefiles_from_GDrive/goagrid.shp")

## Take out stratum = 0 for stations_21
nrow(stations_21) #21588 stations in historical design
nrow(stations_21UT) #1695 stations in historical design are untrawlable

## Which grid cells have untrawlable areas
UT_cells <- sort(unique(stations_21UT$ID))
nrow(stations_23[stations_23$STATIONID %in% UT_cells])
## We need to account for 1613 grid cells and 3115 stations within these cells

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
