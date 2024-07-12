# ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ## Project:       Reformat Trawlability Infomration
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For use in the transition of trawlable information of the
##                historical 1984 stations (last updated after the 2023 GOA
##                survey) to the 2025 stations.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(RODBC)
library(getPass)
library(usethis)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to oracle (remember to connect to VPN!)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
assign(x = "channel",
       value = RODBC::odbcConnect("AFSC",
                                  uid = getPass::getPass("uid"),
                                  pwd = getPass::getPass("pwd")),
       envir = .GlobalEnv)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query historical station data from Oracle with trawlable informaiton
##   updated since 2023
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid_2023 <-
  RODBC::sqlQuery(channel = channel, query = "select * from GOA.GOAGRID_GIS")
attributes(goa_grid_2023)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import polygons of the historical stations and merge the trawlable
##   information from goa_grid_2023 onto the shapefile
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid_2023_shp <-
  terra::vect(x = "data/GOA/shapefiles_from_GDrive/goagrid.shp")

goa_grid_2023_shp <- merge(x = goa_grid_2023_shp,
                           y = goa_grid_2023[, c("GOAGRID_ID", "TRAWLABLE")],
                           by = "GOAGRID_ID")
goa_grid_2023_shp$TRAWLABLE <-
  ifelse(test = is.na(x = goa_grid_2023_shp$TRAWLABLE),
         yes = "UNK",
         no = goa_grid_2023_shp$TRAWLABLE)
goa_grid_2023_shp$STRATUM <- as.character(x = goa_grid_2023_shp$STRATUM)
names(x = goa_grid_2023_shp) <- c("GOAGRID_ID", "AREA_M2", "PERIMETER_M",
                                  "GOAGRID_", "AREA_KM2", "PERIMETER_KM",
                                  "STRATUM", "ID", "TRAWLABLE")

goa_grid_2023_shp$DESIGN_YEAR <- 1984

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import station polygons to be used in 2025
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_stations_2025 <-
  terra::vect(x = "data/GOA/processed_shapefiles/goa_stations_agg_2025.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Intersect the historical station polygons (w/trawlable info) onto the
##   2025 station polygons. Turn "NA" trawlable value to "UNK".
##   Aggregate multiple polygons within a goagrid_id/trawlable combination.
##   Calculate polygon areas and perimeters of the new station/trawlable polys.
##   Denote the design year that the stations belong to (2025)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_stations_2025_overlap_trawl <-
  terra::intersect(x = goa_stations_2025,
                   y = goa_grid_2023_shp[, c("ID", "TRAWLABLE")])

goa_stations_2025_overlap_trawl <-
  terra::aggregate(x = goa_stations_2025_overlap_trawl,
                   by = c("STRATUM", "GOAGRID_ID", "TRAWLABLE"),
                   count = FALSE)

goa_stations_2025_overlap_trawl$AREA_KM2 <-
  round(terra::expanse(x = goa_stations_2025_overlap_trawl) / 1e6, 7)

goa_stations_2025_overlap_trawl$PERIMETER_KM <-
  terra::perim(x = goa_stations_2025_overlap_trawl) / 1000

goa_stations_2025_overlap_trawl$DESIGN_YEAR <- 2025

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge 1984 and 2025 station information. Add a time stamp. Save.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(x = goa_grid_2023_shp)[names(x = goa_grid_2023_shp) == "ID"] <-
  "STATIONID"

station_trawl_info <-
  rbind(goa_grid_2023_shp[, c("GOAGRID_ID", "STATIONID", "STRATUM",
                              "TRAWLABLE", "AREA_KM2", "PERIMETER_KM",
                              "DESIGN_YEAR")],
        goa_stations_2025_overlap_trawl[, c("GOAGRID_ID", "STATIONID", "STRATUM", "TRAWLABLE",
                                       "AREA_KM2", "PERIMETER_KM", "DESIGN_YEAR")])

station_trawl_info <- station_trawl_info[order(station_trawl_info$GOAGRID_ID),]
attributes(station_trawl_info)$timestamp <- Sys.Date()

terra::writeVector(x = station_trawl_info,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "station_trawl_info.shp"),
                   overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Subset 2025 station information and save internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_stations_2025 <- subset(x = as.data.frame(x = station_trawl_info),
                            subset = DESIGN_YEAR == 2025)
usethis::use_data(goa_stations_2025, overwrite = TRUE, internal = F)
