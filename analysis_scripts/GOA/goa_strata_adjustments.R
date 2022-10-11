##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Restratification Adjustments
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Workflow to take the new stratum boundaries and create new
##                strata polygons and intersect its boundaries with the 5 km
##                GOA grid. Then within each 5km grid cell, calculate the
##                total area and perimeter of the stratum component.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rgdal)
library(raster)
library(sp)
library(rgeos)
library(stars)
library(SpaDES)
library(RColorBrewer)
library(spatialEco)
library(lwgeom)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load depth modifications ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
depth_mods <- read.csv("data/GOA/strata_boundaries/depth_modifications.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge bathy rasters ----
##   These rasters are really dense, so they are stored in parts and then
##   "puzzled" togethered using the SpaDES.tools::mergeRaster() function.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_bathy <- list()
n_split_rasters <- length(dir("data/GOA/split_goa_bathy_ras/")) / 2

for (i in 1:n_split_rasters) { ## Loop over subraster -- start
  ## Import raster
  temp_raster <- raster::raster(paste0("data/GOA/split_goa_bathy_ras/",
                                       "goa_bathy_processed_", i, ".grd"))

  ## Append raster to the split_bathy list
  split_bathy[[i]] <- temp_raster
} ## Loop over subraster -- end

## Merge rasters together
bathy <- SpaDES.tools::mergeRaster(split_bathy)
rm(split_bathy, i, n_split_rasters, temp_raster)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Shapefiles ----
##   goa_domain is a mask of the survey footprint (lat-lon projection)
##   ak_land is AK land
##   goa_grid_untrawl is a polygon of untrawlable areas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- rgdal::readOGR("data/GOA/shapefiles/goagrid.shp")
goa_domain <- rgeos::gUnaryUnion(spgeom = goa_grid)
latlon_crs <- raster::crs(rgdal::readOGR("data/GOA/shapefiles/GOAdissolved.shp"))
goa_domain_latlon <- sp::spTransform(x = goa_domain, CRSobj = latlon_crs)

# data(Station)

ak_land <- rgdal::readOGR(dsn = "data/GOA/shapefiles/alaska_dcw.shp")
goa_grid_untrawl <- rgdal::readOGR("data/GOA/shapefiles/goagrid2019_landuntrawlsndmn.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Management areas ----
##   Create masks for each management area in aea projection and transform
##   goa_domain to aea projection.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Shumagin_shape <- raster::crop(x = goa_domain_latlon,
                               y = raster::extent(c(-176, -159, 50, 65)))
Shumagin_shape <- sp::spTransform(x = Shumagin_shape,
                                  CRSobj = raster::crs(ak_land))

Chirikof_shape <- raster::crop(x = goa_domain_latlon,
                               y = raster::extent(c(-159, -154, 50, 59)))
Chirikof_shape <- sp::spTransform(x = Chirikof_shape,
                                  CRSobj = raster::crs(ak_land))

Kodiak_shape <- raster::crop(x = goa_domain_latlon,
                             y = raster::extent(c(-154, -147, 50, 65)))
Kodiak_shape <- sp::spTransform(x = Kodiak_shape,
                                CRSobj = raster::crs(ak_land))
W_Cook_Inlet <- raster::crop(x = goa_domain_latlon,
                             y = raster::extent(c(-154.5, -154, 59.1,59.5)))
W_Cook_Inlet <- sp::spTransform(x = W_Cook_Inlet,
                                CRSobj = raster::crs(ak_land))
Kodiak_shape <- raster::bind(x = Kodiak_shape, y = W_Cook_Inlet)

Yakutat_shape <- raster::crop(x = goa_domain_latlon,
                              y = raster::extent(c(-147, -140, 50, 65)))
Yakutat_shape <- sp::spTransform(x = Yakutat_shape,
                                 CRSobj = raster::crs(ak_land))

Southeast_shape <- raster::crop(x = goa_domain_latlon,
                                y = raster::extent(c(-140, -132, 50, 65)))
Southeast_shape <- sp::spTransform(x = Southeast_shape,
                                   CRSobj = raster::crs(ak_land))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create strata polygons ----
##   For each management area, create new strata based on depth specifications
##   and append to strata_list
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

strata_list <- list()
for (idistrict in unique(depth_mods$manage_area)) {

  ## Crop bathymetry raster to just the management area
  district_outline <- get(paste0(idistrict, "_shape"))
  district_bathy <- raster::crop(x = bathy,
                                 y = district_outline)

  ## Define modified stratum depth boundaries
  depth_boundary <- subset(x = depth_mods,
                           subset = manage_area == idistrict,
                           select = c("lower_depth_m", "upper_depth_m"))


  ## Define each raster cell based on the stratum depth boundaries
  values(district_bathy) <-
    as.integer(as.character(cut(x = values(district_bathy),
                                breaks = c(0, depth_boundary$upper_depth_m),
                                labels = 1:nrow(depth_boundary)) ))

  ## crop out any part of district_bathy that is outside the survey footprint
  district_bathy <- raster::mask(x = district_bathy, mask = goa_domain)

  ## Convert raster to polygon
  strata_poly <- stars::st_as_stars(district_bathy) %>%
    sf::st_as_sf(merge = TRUE) ## this is the raster to polygons part

  ## Convert back to spatial object
  strata_poly <- sf::as_Spatial(strata_poly)

  ## Create dataframe of stratum information
  strata_poly@data <-
    data.frame(manage_area = idistrict,
               stratum = paste0(idistrict, "_", strata_poly@data$dblbnd))

  strata_poly <- raster::aggregate(x = strata_poly,
                                   by = c("manage_area", "stratum"))

  strata_poly$AREA_KM2 <- raster::area(strata_poly) / 1e6
  strata_poly$PER_KM <-
    as.numeric(lwgeom::st_perimeter(stars::st_as_stars(strata_poly)) / 1000)

  ## Append to strata_list
  strata_list <- c(strata_list, list(strata_poly))

  print(paste("Finished with the", idistrict, "region"))
}
rm(idistrict, strata_poly, district_outline, district_bathy)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge strata into one object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_list <- raster::bind(strata_list)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Take the current goa grid and merge all historical stations wihtin 5x5 km
##   cell. This essentially resets the stratum data in each grid cell.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_full_grid <- raster::aggregate(x = goa_grid,
                                   by = "ID",
                                   sums = list(list(sum, "AREA_KM2")))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   New stations ----
##
##   The station ID is a pair of numbers that indicate the position of a grid
##   cell. The first number indicates the
##
##   Intersect the new strata polygons with the 5x5 km grid. This function
##   takes a few minutes. Then calculate the area and perimeter of each station
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations <- raster::intersect(x = goa_full_grid, y = strata_list)
stations <- raster::aggregate(x = stations,
                              by = c("ID", "manage_area", "stratum"))
stations$PER_KM <- spatialEco::polyPerimeter(x = stations) / 1e3
stations$AREA_KM2 <- rgeos::gArea(spgeom = stations, byid = T) / 1e6
sp::proj4string(obj = stations) <-raster::crs(goa_full_grid)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Untrawlable area calculation ----
##   Intersect the untrawlable areas with the new stations and calculate
##   the area of the untrawlable area in the stations where there are
##   untrawlable areas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overlap_with_trawl_polygon <- raster::intersect(stations, goa_grid_untrawl)
overlap_with_trawl_polygon@data <- overlap_with_trawl_polygon@data[, 1:3]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Update stations shapefile with untrawlable information
##  Create a trawlable field in the stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations$trawlable <- TRUE

for (irow in 1:length(overlap_with_trawl_polygon)) {
  temp_id <- overlap_with_trawl_polygon$ID[irow]
  temp_str <- overlap_with_trawl_polygon$stratum[irow]

  temp_idx <- which(stations$ID == temp_id & stations$stratum == temp_str)
  temp_untrawl <- subset(x = overlap_with_trawl_polygon,
                         subset = ID == temp_id & stratum == temp_str)

  stations$trawlable[temp_idx] <- FALSE
}
rm(temp_id, temp_str, temp_idx)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data.frame()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save results ----
##  Create export directory and export
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/updated_goa_strata_2023/"))
  dir.create("data/GOA/updated_goa_strata_2023/")

names(strata_list) <- c("MGT_AREA", "STRATUM", "AREA_KM2", "PER_KM")
writeOGR(obj = strata_list,
         dsn = "data/GOA/updated_goa_strata_2023/updated_goa_strata.shp",
         layer = "updated_goa_strata",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

names(stations) <- c("ID", "MGT_AREA", "STRATUM", "PER_KM" ,
                     "AREA_KM2", "TRAWL", "UT_AR_KM2")
writeOGR(obj = stations,
         dsn = "data/GOA/updated_goa_strata_2023/updated_stations.shp",
         layer = "updated_stations",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

names(overlap_with_trawl_polygon) <- c("ID", "MGT_AREA", "STRATUM", "AREA_KM2")
writeOGR(obj = overlap_with_trawl_polygon,
         dsn = paste0("data/GOA/updated_goa_strata_2023/UT_areas.shp"),
         layer = "overlap_with_trawl_polygon",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

usethis::use_data(name = stations, overwrite = TRUE)
usethis::use_data(name = strata_list, overwrite = TRUE)
usethis::use_data(name = ak_land, overwrite = TRUE)
usethis::use_data(name = depth_mods, overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = "data/GOA/updated_goa_strata_2023/updated_strata.pdf",
    width = 6, height = 6, onefile = TRUE)
for (iarea in unique(depth_mods$manage_area)[]) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  n_strata <- with(depth_mods, table(manage_area))[iarea]
  temp_area <- subset(stations, MGT_AREA == iarea)
  temp_area <- temp_area[order(temp_area$STRATUM), ]

  plot(temp_area, axes = F, col = "white", border = F)

  for (istratum in 1:n_strata) {
    plot(subset(stations,
                STRATUM == unique(temp_area$STRATUM)[istratum]),
         col =   c(RColorBrewer::brewer.pal("Spectral",
                                            n = n_strata - 1),
                   "gray" )[istratum],
         border = F, add = TRUE)
  }

  plot(crop(x = overlap_with_trawl_polygon,
            y = get(paste0(iarea, "_shape"))),
       col = "black", add = TRUE, border = FALSE)
  # plot(ak_land, add = TRUE, col = "tan", border = F)

  plot(subset(strata_list, MGT_AREA == iarea),
       lwd = 0.1, add = TRUE)

  ## Legend
  legend_labels <- with(subset(depth_mods, manage_area == iarea),
                        paste0(lower_depth_m, " - ", upper_depth_m, " m"))
  legend_labels[length(legend_labels)] <- "> 700 m"
  legend_labels <- c(legend_labels, "Untrawlable")

  legend(c("Shumagin" = "bottomleft", "Chirikof" = "left",
           "Kodiak" = "bottomright", "Yakutat" = "bottomleft",
           "Southeast" = "bottomleft")[iarea],
         legend = legend_labels, title = "Stratum Legend",
         fill = c(RColorBrewer::brewer.pal("Spectral", n = n_strata - 1),
                  "gray", "black"))
  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}
dev.off()
