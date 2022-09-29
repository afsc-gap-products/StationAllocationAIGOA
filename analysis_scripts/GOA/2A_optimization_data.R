##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a
##                Gulf of Alaska groundfish VAST spatiotemporal
##                operating single-species models
##
##                Set up other constants used in downstream processes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Library
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(rgdal)
library(raster)
library(sp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load the true density, true index, and spatial domain dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load(file = "data/GOA/prednll_VAST_models.RData")
grid_goa <- read.csv(file = "data/GOA/grid_goa.csv")
D_gct <- readRDS("data/GOA/VAST_fit_D_gct.RDS")

## Think about ways to put this in package!!
updated_goa_strata <- rgdal::readOGR(dsn = "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/products/updated_goa_strata/updated_goa_strata.shp")
depth_mods <- read.csv("C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/products/depth_modifications.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used throughout all scripts
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"
n_years <- dim(D_gct)[3]
n_spp <- dim(D_gct)[2]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Assign grid points to the new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_goa_sp <- sp::SpatialPointsDataFrame(
  coords = grid_goa[, c("Lon", "Lat")],
  proj4string = sp::CRS(lonlat_crs),
  data = data.frame(ID = 1:nrow(grid_goa)))
grid_goa_sp <- sp::spTransform(x = grid_goa_sp,
                               CRSobj = crs(updated_goa_strata))

grid_goa_sp <- raster::intersect(x = grid_goa_sp,
                                 y = updated_goa_strata)
grid_goa_sp <- subset(x = grid_goa_sp,
                      subset = STRATUM %in% depth_mods$stratum[depth_mods$used])

grid_goa_sp <- sp::remove.duplicates(grid_goa_sp)

removed_cells <- (1:n_cells)[-grid_goa_sp$ID]
D_gct <- D_gct[-removed_cells, , ]
n_cells <- dim(D_gct)[1]

grid_goa_sp <- grid_goa_sp[, 1:3]

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
usethis::use_data(grid_goa_sp, overwrite = TRUE)

usethis::use_build_ignore("data/GOA")
