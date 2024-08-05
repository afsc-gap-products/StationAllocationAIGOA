###############################################################################
## Project:     VAST covariates across goa grid
## Author:      Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description: Assign bathymetry value from the EFH data lyaer to each grid in
##              the Gulf of Alaska grid.
##              Uses R version 4.3.x to incorporate updates from akgfmaps
###############################################################################
rm(list = ls())

##################################################
####   Import packages
##################################################
library(sf); library(stars); library(akgfmaps)

##################################################
####   CRSs used
##################################################
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##################################################
####   Merge together bathymetry rasters
##################################################
goa_bathy <-
  terra::rast(x = "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/goa_bathy")

##################################################
####   Import Extrapolation grid (goa_grid)
####   Import NMFS areas (nmfs)
####   Import CPUE data (goa_data_geostat)
####   Import stratum boundaries (stratum_boundaries)
##################################################

goa_grid <-
  read.csv(file = "data/GOA/extrapolation_grid/GOAThorsonGrid_Less700m.csv")
goa_grid <- goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2
names(x = goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

nmfs <-
  terra::vect(x = akgfmaps::get_nmfs_areas(set.crs = terra::crs(x = goa_bathy)))
nmfs <-
  nmfs[nmfs$REP_AREA %in% c(519, 610, 620, 630, 640, 650, 659), "REP_AREA"]

goa_data_geostat = read.csv(file = "data/GOA/vast_data/goa_data_geostat.csv")

stratum_boundaries <- read.csv(file = "data/GOA/strata_boundaries/depth_modifications_2025.csv")
stratum_boundaries <- rbind(stratum_boundaries,
                            data.frame(NMFS_AREA = c("Southeast Inside",
                                                     "NMFS519"),
                                       REP_AREA = c(659, 519),
                                       STRATUM = c(52, 16),
                                       DEPTH_MIN_M = 1,
                                       DEPTH_MAX_M = 1000))

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = terra::vect(x = goa_grid[, c("Lon", "Lat", "Area_km2")],
                         geom = c("Lon", "Lat"), keepgeom = TRUE,
                         crs = lonlat_crs)
grid_shape_aea <- terra::project(x = grid_shape, terra::crs(x = goa_bathy))
grid_shape_aea[, c("Eastings", "Northings")] <- terra::crds(x = grid_shape_aea)

grid_shape_aea$Depth_m <-
  terra::extract(x = goa_bathy, y = grid_shape_aea)$GOA_bathy

##################################################
####   Classify cells to NMFS reporting areas
##################################################
grid_shape_aea <- terra::intersect(x = grid_shape_aea, y = nmfs)

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
grid_shape_aea <- grid_shape_aea[
  (grid_shape_aea$Depth_m >= min(x = goa_data_geostat$Depth_m) &
     grid_shape_aea$Depth_m <= 700),
]

##################################################
#### Reclassify the stratum of the interpolation cells based on extracted depth
##################################################
grid_shape_aea$STRATUM <- NA

for (inmfs in nmfs$REP_AREA) {
  temp_boundaries <- subset(x = stratum_boundaries,
                            subset = REP_AREA == inmfs)
  for (istratum in 1:nrow(x = temp_boundaries)) {
    grid_shape_aea[
      grid_shape_aea$REP_AREA == inmfs &
        round(x = grid_shape_aea$Depth_m) >= temp_boundaries$DEPTH_MIN_M[istratum] &
        round(x = grid_shape_aea$Depth_m) <= temp_boundaries$DEPTH_MAX_M[istratum]
      ]$STRATUM <- temp_boundaries$STRATUM[istratum]
  }
}

##################################################
####   scale grid bathymetry values to standard normal, using the mean and sd
####   of the BTS data
##################################################
BTS_mean <- mean(x = log10(x = goa_data_geostat$Depth_m))
BTS_sd   <-  sd(x = log10(x = goa_data_geostat$Depth_m))

grid_shape_aea$LOG10_DEPTH_M <- log10(x = grid_shape_aea$Depth_m)
grid_shape_aea$LOG10_DEPTH_M_CEN <-
  (grid_shape_aea$LOG10_DEPTH_M - BTS_mean) / BTS_sd

##################################################
####   Save
##################################################
if (!dir.exists(paths = "data/GOA/vast_data/"))
  dir.create(path = "data/GOA/vast_data/")
write.csv(as.data.frame(x = grid_shape_aea),
          row.names = F,
          file = "data/GOA/vast_data/goa_interpolation_grid.csv")

saveRDS(object = as.data.frame(x = grid_shape_aea),
        file = "data/GOA/vast_data/goa_interpolation_grid.RDS")
