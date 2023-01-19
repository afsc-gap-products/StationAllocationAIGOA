##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Reformat Trawlability Infomration
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For use in the transition of trawlable information of the
##                historical 1984 stations (last updated after the 2021 GOA
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
##   updated since 2021
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_q <- "select * from GOA.GOAGRID_GIS"
goa_grid_2021 <- RODBC::sqlQuery(channel = channel, query = grid_q)
attributes(goa_grid_2021)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import polygons of the historical stations and attach the trawlable
##   information from goa_grid_2021 onto the shapefile
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid_2021_shp <- terra::vect("data/GOA/shapefiles_from_GDrive/goagrid.shp")

goa_grid_2021_shp <- merge(x = goa_grid_2021_shp,
                           y = goa_grid_2021[, c("GOAGRID_ID", "TRAWLABLE")],
                           by = "GOAGRID_ID")
goa_grid_2021_shp$TRAWLABLE <-
  ifelse(test = is.na(goa_grid_2021_shp$TRAWLABLE),
         yes = "UNK",
         no = goa_grid_2021_shp$TRAWLABLE)
goa_grid_2021_shp$STRATUM <- as.character(goa_grid_2021_shp$STRATUM)
names(goa_grid_2021_shp) <- c("GOAGRID_ID", "AREA_M2", "PERIMETER_M",
                              "GOAGRID_", "AREA_KM2", "PERIMETER_KM",
                              "STRATUM", "ID", "TRAWLABLE")

goa_grid_2021_shp$DESIGN <- 1984

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import station polygons to be used in 2025
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_stations_2025 <-
  terra::vect(x = "data/GOA/processed_shapefiles/goa_stations_2023.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Intersect the historical station polygons (w/trawlable info) onto the
##   2025 station polygons. Turn "NA" trawlable value to "UNK".
##   Aggregate multiple polygons within a goagrid_id/trawlable.
##   Calculate polygon areas and perimeters of the new station/trawlable polys.
##   Denote the design year that the stations belong to (2025)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overlap_with_trawl_polygon <-
  terra::intersect(x = goa_stations_2025,
                   y = goa_grid_2021_shp)

overlap_with_trawl_polygon$TRAWLABLE <-
  ifelse(test = is.na(overlap_with_trawl_polygon$TRAWLABLE),
         yes = "UNK",
         no = overlap_with_trawl_polygon$TRAWLABLE)

overlap_with_trawl_polygon <-
  terra::aggregate(x = overlap_with_trawl_polygon,
                   by = c("GOAGRID_ID", "TRAWLABLE"))

overlap_with_trawl_polygon$AREA_KM2 <-
  round(terra::expanse(x = overlap_with_trawl_polygon) / 1e6, 7)

overlap_with_trawl_polygon$PERIMETER_KM <-
  terra::perim(x = overlap_with_trawl_polygon) / 1000

overlap_with_trawl_polygon$DESIGN <- 2025

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge 1984 and 2025 station information. Add a time stamp. Save.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
station_trawl_info <-
  rbind(goa_grid_2021_shp[, c("GOAGRID_ID", "STRATUM", "TRAWLABLE",
                              "AREA_KM2", "PERIMETER_KM", "DESIGN")],
        overlap_with_trawl_polygon[, c("GOAGRID_ID", "STRATUM", "TRAWLABLE",
                                       "AREA_KM2", "PERIMETER_KM", "DESIGN")])

station_trawl_info <- station_trawl_info[order(station_trawl_info$GOAGRID_ID),]
attributes(station_trawl_info)$timestamp <- Sys.Date()

terra::writeVector(x = station_trawl_info,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "station_trawl_info.shp"),
                   overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Subset 2025 station information and save internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_stations_2025 <- subset(x = as.data.frame(station_trawl_info),
                            subset = DESIGN == 2025)
usethis::use_data(goa_stations_2025, overwrite = TRUE, internal = F)
