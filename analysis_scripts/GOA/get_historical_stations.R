##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(RODBC)
library(getPass)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to VPN!
##   Connect to oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
assign(x = "channel",
       value = RODBC::odbcConnect("AFSC",
                                  uid = getPass::getPass("uid"),
                                  pwd = getPass::getPass("pwd")),
       envir = .GlobalEnv)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query Station data from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_q <- "select * from GOA.GOAGRID_GIS"
goa_grid_2021 <-  RODBC::sqlQuery(channel = channel, query = grid_q)
attributes(goa_grid_2021)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query Strata data from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_q <- "select * from GOA.GOA_STRATA where SURVEY = 'GOA'"
goa_strata_2021 <-  RODBC::sqlQuery(channel = channel, query = strata_q)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import shapefile of untrawlable areas as of 2019, Land, and sandman reef
##   that potentially coincide with the goa grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawl_2019_shp <- terra::vect(x = paste0("data/GOA/shapefiles_from_GDrive/",
                                           "goagrid2019_landuntrawlsndmn.shp"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Identify untrawlable stations as of 2021
##   Create a shapefile of untrawlable stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawlable_stations <-
  na.omit(goa_grid_2021$GOAGRID_ID[goa_grid_2021$TRAWLABLE == "N"])

goa_grid_2021_shp <- terra::vect("data/GOA/shapefiles_from_GDrive/goagrid.shp")
untrawl_2021_shp <-
  goa_grid_2021_shp[goa_grid_2021_shp$GOAGRID_ID %in% untrawlable_stations]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge 2019 shapefile of untrawlable areas and
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawl_2021_shp <- terra::union(untrawl_2019_shp, untrawl_2021_shp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Upload dfs to repo
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(goa_grid_2021, file = "data/GOA/grid_goa_2021.rds")
saveRDS(goa_strata_2021, file = "data/GOA/grid_strata_2021.rds")

if (!dir.exists("data/GOA/processed_shapefiles/"))
  dir.create("data/GOA/processed_shapefiles/")
terra::writeVector(x = untrawl_2021_shp,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "GOA_untrawl_2021.shp"),
                   overwrite=TRUE)
