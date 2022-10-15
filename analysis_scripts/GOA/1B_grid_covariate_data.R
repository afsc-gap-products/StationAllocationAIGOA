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
utm_crs <- "+proj=utm +zone=5N +units=km"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge together bathymetry rasters
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_bathy <- list()
n_split_rasters <- length(dir("data/GOA/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) {
  split_bathy[[i]] <- terra::rast(x = paste0("data/GOA/",
                                             "split_goa_bathy_ras",
                                             "/goa_bathy_processed_",
                                             i, ".grd"))
}

bathy <- do.call(what = terra::merge, args = split_bathy)
rm(split_bathy, i)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import spatial objects:
##   Current Strata (current_survey_strata)
##   Spatial domain outline mask (current_survey_mask)
##   Extrapolation grid (goa_grid)
##   CPUE data (data)
##   Untrawlable areas (goa_grid_nountrawl)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
current_survey_strata <- terra::vect(x = "data/GOA/shapefiles/goa_strata.shp")

current_survey_mask <-
  terra::vect(x = "data/GOA/shapefiles/goagrid_polygon.shp")
current_survey_mask <- terra::project(x = current_survey_mask,
                                      y = bathy)
current_survey_mask <- terra::aggregate(x = current_survey_mask)
# current_survey_mask <- rgeos::gUnaryUnion(spgeom = current_survey_mask)

goa_grid <- read.csv("data/GOA/extrapolation_grid/GOAThorsonGrid.csv")
goa_grid <- goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2
names(goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

data = read.csv("data/GOA/goa_vast_data_input.csv")
goa_grid_nountrawl <- read.csv(paste0("data/GOA/extrapolation_grid/",
                                      "GOA_ALL_nountrawl.csv"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Transform extrapolation grid to aea, extract bathymetry values onto grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_shape <- terra::vect(x = goa_grid[, c("Lon", "Lat")],
                          geom = c("Lon", "Lat"),
                          keepgeom = TRUE,
                          crs = lonlat_crs)

grid_shape_aea = terra::project(x = grid_shape, y = bathy)
grid_shape_aea$DEPTH_EFH <- terra::extract(x = bathy,
                                           y = grid_shape_aea)$dbl

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Remove cells not already in the goa grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_shape_aea <- terra::intersect(x = grid_shape_aea,
                                   y = current_survey_mask)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Remove cells with depths outside the observed range to the range
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_shape_aea <-
  grid_shape_aea[grid_shape_aea$DEPTH_EFH >= min(data$DEPTH_EFH) &
                   grid_shape_aea$DEPTH_EFH <= max(data$DEPTH_EFH),]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/")) dir.create("data/GOA/")
write.csv(x = as.data.frame(grid_shape_aea),
          row.names = F,
          file = "data/GOA/grid_goa.csv")

test <- readRDS("data/frame.RDS")
