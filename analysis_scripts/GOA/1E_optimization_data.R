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
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(terra)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load the true density, true index, and spatial domain dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load(file = "data/GOA/prednll_VAST_models.RData")
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
                   y = updated_goa_strata[, c("STRATUM")])
grid_goa_sp <- grid_goa_sp[grid_goa_sp$stratum %in%
                             depth_mods$stratum[depth_mods$used], ]


## remove duplicates
rm_idx <- c()

for (idx in names(which(table(grid_goa_sp$ID) == 2))) {
  rm_idx <- c(rm_idx, which(grid_goa_sp$ID == idx)[1])
}

grid_goa_sp <- grid_goa_sp[-rm_idx, ]

removed_cells <- (1:n_cells)[-grid_goa_sp$ID]
D_gct <- D_gct[-c(removed_cells, rm_idx), , ]
n_cells <- dim(D_gct)[1]

##################################################
####   Our df will have fields for:
####   domain: only one domain so the value is just 1
####   id: unique ID for each sampling cell
####   X1: strata variable 2: depth of cell (m)
####   X2: strata variable 1: longitude in eastings (km). Because the
####       optimization does not read in negative values, I shift the
####       values so that the lowest value is 0
####
####   Variables used to more efficiently calcualte stratum variance
####
####   WEIGHT: number of observed years
####   Y1, Y2, ... : density for a given cell summed across observed years
####   Y1_SQ_SUM, Y2_SQ_SUM, ... : density-squared for a given cell,
####           summed across observed years
##################################################
frame <- cbind(
  data.frame(domainvalue = 1,
             id = 1:n_cells,
             depth = grid_goa_sp$DEPTH_EFH[-c(removed_cells, rm_idx)],
             eastings_m = geom(grid_goa_sp)[-c(removed_cells, rm_idx), "x"],
             northings_m = geom(grid_goa_sp)[-c(removed_cells, rm_idx), "y"],
             stratum = grid_goa_sp$STRATUM[-c(removed_cells, rm_idx)],
             WEIGHT = n_years),

  matrix(data = apply(X = D_gct,
                      MARGIN = 1:2,
                      FUN = sum),
         ncol = n_spp,
         dimnames = list(NULL, paste0("Y", 1:n_spp))),

  matrix(data = apply(X = D_gct,
                      MARGIN = 1:2,
                      FUN = function(x) sum(x^2)),
         ncol = n_spp,
         dimnames = list(NULL, paste0("Y", 1:n_spp, "_SQ_SUM")))
)

attributes(frame)$spp_name <- dimnames(D_gct)[[2]]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usethis::use_data(frame, overwrite = TRUE)
