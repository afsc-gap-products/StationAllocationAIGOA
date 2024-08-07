##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a
##                Gulf of Alaska groundfish VAST spatiotemporal
##                operating single-species models
##
##                Set up other constants used in downstream processes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Library
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load the true density, and spatial domain dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
optim_df <- read.csv(file = "data/GOA/vast_data/goa_interpolation_grid.csv")

vast_data <- read.csv(file = "data/GOA/vast_data/goa_data_geostat.csv")
D_gct <- readRDS(file = "data/GOA/vast_data/VAST_fit_D_gct.RDS")
dimnames(x = D_gct)[[3]] <- sort(x = unique(x = vast_data$Year))

goa_stratum_boundaries <- as.data.frame(
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_strata_2025.gpkg")
)
goa_stratum_boundaries$USED <- goa_stratum_boundaries$DEPTH_MIN_M != 701

goa_stations <- as.data.frame(
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_stations_2025.gpkg")
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used throughout all scripts
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"
n_years <- dim(x = D_gct)[3]
n_spp <- dim(x = D_gct)[2]
n_cells <- dim(x = D_gct)[1]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usethis::use_data(D_gct, overwrite = TRUE)
usethis::use_data(optim_df, overwrite = TRUE)
usethis::use_data(goa_stations, overwrite = TRUE)
usethis::use_data(goa_stratum_boundaries, overwrite = TRUE)
