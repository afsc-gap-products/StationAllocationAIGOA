##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Take 2021 stations polygons and assign trawlable information
##                based on information in GOA.GOAGRID_GIS from oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
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
##   Query Strata data from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_q <- "select * from GOA.GOA_STRATA where SURVEY = 'GOA'"
goa_strata_2021 <-  RODBC::sqlQuery(channel = channel, query = strata_q)
attributes(goa_strata_2021)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import shapefile of historical stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid_2021_shp <- terra::vect("data/GOA/shapefiles_from_GDrive/goagrid.shp")
goa_grid_df <- as.data.frame(goa_grid_2021_shp)
attributes(goa_grid_df)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Upload dfs to repo
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(as.data.frame(goa_grid_2021_shp), file = "data/GOA/grid_goa_2021.rds")
saveRDS(goa_strata_2021, file = "data/GOA/grid_strata_2021.rds")
