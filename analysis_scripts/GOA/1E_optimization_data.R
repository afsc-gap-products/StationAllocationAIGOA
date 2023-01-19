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
grid_goa <- read.csv(file = "data/GOA/vast_grid_goa.csv")
D_gct <- readRDS("data/GOA/VAST_fit_D_gct.RDS")

updated_goa_strata <-
  terra::vect(x = "data/GOA/processed_shapefiles/goa_strata_2023.shp")
depth_mods <- read.csv("data/GOA/strata_boundaries/depth_modifications_2023.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used throughout all scripts
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"
n_years <- dim(D_gct)[3]
n_spp <- dim(D_gct)[2]
n_cells <- dim(D_gct)[1]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Assign grid points to the new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_goa_sp <- terra::vect(x = cbind(ID = 1:nrow(grid_goa), grid_goa),
                           geom = c("Lon", "Lat"),
                           crs = lonlat_crs)
grid_goa_sp <- terra::project(x = grid_goa_sp,
                              y = updated_goa_strata)

grid_goa_sp <-
  terra::intersect(x = grid_goa_sp,
                   y = updated_goa_strata[, c("INPFC_AREA", "STRATUM")])
grid_goa_sp <- grid_goa_sp[grid_goa_sp$STRATUM %in%
                             depth_mods$stratum[depth_mods$used], ]

## remove duplicates
rm_idx <- c()

for (idx in names(which(table(grid_goa_sp$ID) == 2))) {
  rm_idx <- c(rm_idx, which(grid_goa_sp$ID == idx)[1])
}

grid_goa_sp <- grid_goa_sp[-rm_idx, ]

removed_cells <- (1:n_cells)[-grid_goa_sp$ID]
D_gct <- D_gct[-c(removed_cells), , ]
n_cells <- dim(D_gct)[1]

optim_df <- cbind(as.data.frame(grid_goa_sp),
                  terra::geom(grid_goa_sp)[, c("x", "y")])
attributes(optim_df)$crs <- terra::crs(grid_goa_sp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usethis::use_data(D_gct, overwrite = TRUE)
usethis::use_data(optim_df, overwrite = TRUE)
