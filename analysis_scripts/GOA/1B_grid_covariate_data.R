##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       VAST covariates across goa grid
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Assign bathymetry value from the EFH data layer to each grid
##                in the Gulf of Alaska grid.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   CRSs used
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge bathy rasters ----
##   These rasters are really dense, so they are stored in parts and then
##   "puzzled" togethered using terra::merge()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_bathy <- list()

n_split_rasters <- length(grep(x = dir("data/GOA/processed_rasters/"),
                               pattern = "aigoa"))
for (i in 1:n_split_rasters) {
  split_bathy[[i]] <- terra::rast(x = paste0("data/GOA/processed_rasters/",
                                             "aigoa_", i, ".tif"))
}

bathy <- do.call(what = terra::merge, args = split_bathy)
rm(split_bathy, i, n_split_rasters)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import spatial objects:
##   Historical station polygons (historical stations)
##   Create mask by dissovling inner boundaries of the historical stations
##   Extrapolation grid used for VAST, 2 nmi resolution (goa_grid)
##   CPUE data (data)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
historical_stations <- terra::vect(x = "data/GOA/shapefiles_from_GDrive//goagrid.shp")
historical_stations <- terra::project(x = historical_stations, y = bathy)

historical_mask <- terra::aggregate(x = historical_stations)

vast_goa_grid <- read.csv("data/GOA/extrapolation_grid/GOAThorsonGrid.csv")
vast_goa_grid <- vast_goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
vast_goa_grid$Shape_Area <- vast_goa_grid$Shape_Area / 1e6 #Convert to km2
names(vast_goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

data = read.csv("data/GOA/goa_vast_data_input.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Transform extrapolation grid to aea, extract bathymetry values onto grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vast_goa_grid_shape <- terra::vect(x = vast_goa_grid[, c("Lon", "Lat")],
                                   geom = c("Lon", "Lat"),
                                   keepgeom = TRUE,
                                   crs = lonlat_crs)

vast_goa_grid_shape_aea = terra::project(x = vast_goa_grid_shape, y = bathy)
vast_goa_grid_shape_aea$Area_km2 <- vast_goa_grid$Area_km2
vast_goa_grid_shape_aea$DEPTH_EFH <-
  terra::extract(x = bathy, y = vast_goa_grid_shape_aea)$AIGOA_ba

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Remove cells not already in the goa grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vast_goa_grid_shape_aea <- terra::intersect(x = vast_goa_grid_shape_aea,
                                            y = historical_mask)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Remove cells with depths outside the observed range to the range
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vast_goa_grid_shape_aea <-
  vast_goa_grid_shape_aea[vast_goa_grid_shape_aea$DEPTH_EFH >= min(data$DEPTH_EFH) &
                            vast_goa_grid_shape_aea$DEPTH_EFH <= max(data$DEPTH_EFH),]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/")) dir.create("data/GOA/")
write.csv(x = as.data.frame(vast_goa_grid_shape_aea),
          row.names = F,
          file = "data/GOA/vast_grid_goa.csv")
