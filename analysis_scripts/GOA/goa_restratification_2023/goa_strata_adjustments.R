##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Restratification Adjustments
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Workflow to take the new stratum boundaries and create new
##                strata polygons and intersect its boundaries with the 5 km
##                GOA grid. Then within each 5km grid cell, calculate the
##                total area and perimeter of the stratum component.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(RColorBrewer)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load depth modifications ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
depth_mods <-
  read.csv("data/GOA/strata_boundaries/depth_modifications_2023.csv")

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
##   Import Data Polygons ----
##   goa_grid: historical survey stations
##   goa_domain: is a mask of the survey footprint made by dissolving the
##               inner boundaries of the goa_grid
##   goa_domain_latlon: goa_domain projected onto lon/lat
##
##   goa_grid_2021: historical stations from Oracle (GOA.GOAGRID_GIS)
##   goa_strata_2021: historical strata from Oracle (GOA.GOA_STRATA)
##
##   ak_land and ca_land are land polygons for Alaska and Canada
##   goa_grid_untrawl is a polygon of untrawlable areas
##
##   nmfs are polygons of the five management areas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/goagrid.shp")
goa_domain <- terra::aggregate(x = goa_grid)

latlon_crs <-
  terra::crs(x = terra::vect(x = paste0("data/GOA/shapefiles_from_GDrive/",
                                        "GOAdissolved.shp")))
goa_domain_latlon <- terra::project(x = goa_domain, y = latlon_crs)

goa_grid_2021 <- readRDS(file = "data/GOA/grid_goa_2021.rds")
goa_grid_2021 <- goa_grid_2021[order(goa_grid_2021$GOAGRID_ID), ]
goa_strata_2021 <- readRDS(file = "data/GOA/grid_strata_2021.rds")

ak_land <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/AKland.shp")
ak_land <- terra::project(x = ak_land, y = bathy)
ak_land <- terra::aggregate(x = ak_land)

ca_land <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/canada_dcw.shp")
ca_land <- ca_land[ca_land$POPYADMIN %in% c("BRITISH COLUMBIA",
                                            "YUKON TERRITORY",
                                            "NORTHWEST TERRITORIES") ]
ca_land <- terra::project(x = ca_land, y = bathy)
ca_land <- terra::aggregate(x = ca_land)

goa_grid_untrawl <-
  terra::vect(x = "data/GOA/processed_shapefiles/GOA_untrawl_2021.shp")

nmfs <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/GOA_Shapes.shp")
nmfs <- terra::project(x = nmfs, y = bathy)
nmfs$area_name <- c("Southeastern", "Southeastern", "Shumagin", "Chirikof",
               "Kodiak", "Yakutat", NA)
nmfs <- terra::aggregate(x = nmfs, by = "area_name")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create strata polygons ----
##   For each management area, create new strata based on depth specifications
##   and append to strata_list
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_list <- list()
for (idistrict in unique(depth_mods$manage_area)) { ## Loop over district --st.

  ## Crop bathymetry raster to just the management area
  district_outline <- nmfs[nmfs$area_name == idistrict]

  district_bathy <- terra::crop(x = bathy,
                                y = district_outline)
  district_bathy <- terra::mask(x = district_bathy,
                                mask = goa_domain)

  ## Define modified stratum depth boundaries
  depth_boundary <- subset(x = depth_mods,
                           subset = manage_area == idistrict,
                           select = c("lower_depth_m", "upper_depth_m"))

  ## Define each raster cell based on the stratum depth boundaries
  values(district_bathy) <-
    as.integer(as.character(cut(x = values(district_bathy),
                                breaks = c(0, depth_boundary$upper_depth_m),
                                labels = 1:nrow(depth_boundary)) ))

  ## Convert raster to polygon
  strata_poly <- terra::as.polygons(district_bathy)

  ## Create dataframe of stratum information
  strata_poly[, names(depth_mods)] <-
      subset(x = depth_mods, subset = manage_area == idistrict)

  strata_poly$AREA_KM2 <- terra::expanse(x = strata_poly, unit = "km")
  strata_poly$PER_KM <- terra::perim(x = strata_poly) / 1000

  ## Append to strata_list
  strata_list <- c(strata_list, list(strata_poly))

  print(paste("Finished with the", idistrict, "region"))
} ## Loop over district -- end
rm(idistrict, strata_poly, district_outline, district_bathy)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge strata into one object
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_list <- do.call(what = rbind, args = strata_list)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Take the current goa grid and merge all historical stations wihtin 5x5 km
##   cell. This essentially resets the stratum data in each grid cell.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_full_grid <- terra::aggregate(x = goa_grid,
                                  by = "ID",
                                  fun = sum)
goa_full_grid <- goa_full_grid[, c("ID", "agg_AREA_KM2")]
names(goa_full_grid) <- c("ID", "AREA_KM2")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   New stations ----
##
##   The station ID is a pair of numbers that indicate the position of a grid
##   cell. The first number indicates the
##
##   Intersect the new strata polygons with the 5x5 km grid. This function
##   takes a few minutes. Then calculate the area and perimeter of each station
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations <- terra::intersect(x = goa_full_grid, y = strata_list)
stations$PER_KM <- terra::perim(x = stations) / 1e3
stations$AREA_KM2 <- terra::expanse(x = stations) / 1e6
stations <- terra::project(x = stations, y = bathy)

stations[, c("E_m", "N_m")] <-
  terra::geom(terra::centroids(x = stations, inside = TRUE))[, c("x", "y")]
stations[, c("Lon", "Lat")] <-
  terra::geom(terra::centroids(x = terra::project(x = stations,
                                                  y = latlon_crs),
                               inside = TRUE))[, c("x", "y")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Untrawlable area calculation ----
##   Intersect the untrawlable areas with the new stations and calculate
##   the area of the untrawlable area in the stations where there are
##   untrawlable areas
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overlap_with_trawl_polygon <- terra::intersect(stations, goa_grid_untrawl)
overlap_with_trawl_polygon <- overlap_with_trawl_polygon[, c(1:3, 5:6)]
overlap_with_trawl_polygon$untrawl_area_km2 <-
  terra::expanse(x = overlap_with_trawl_polygon) / 1e6

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Update stations shapefile with untrawlable information
##  Create a trawlable field in the stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations$trawlable <- TRUE
stations$untrawl_area_km2 <- 0

temp_overlap_df <- as.data.frame(overlap_with_trawl_polygon)
temp_stations_df <- as.data.frame(stations)

temp_idx <- sapply(X = 1:length(overlap_with_trawl_polygon),
                   FUN = function(x) {
                     temp_id <- temp_overlap_df$ID[x]
                     temp_str <-  temp_overlap_df$stratum[x]
                     temp_idx <- which(temp_stations_df$ID == temp_id &
                                         temp_stations_df$stratum == temp_str)
                     return(temp_idx)
                   })

stations$trawlable[temp_idx] <- FALSE
stations$untrawl_area_km2[temp_idx] <- temp_overlap_df$untrawl_area_km2
stations$trawl_area_km2 <-
  round(x = stations$AREA_KM2 - stations$untrawl_area_km2, digits = 6)

rm(temp_overlap_df, temp_stations_df, temp_idx)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goagrid_ids <-
  (max(goa_grid_2021$GOAGRID_ID) + 1):
  (max(goa_grid_2021$GOAGRID_ID) + nrow(stations))

stations_2023 <- data.frame("GOAGRID#" = goagrid_ids + 1,
                            "GOAGRID_ID" = goagrid_ids,
                            "TRAWLABLE" = stations$trawlable,
                            "TRAWLABLE_AREA_KM2" = stations$trawl_area_km2,
                            "AREA_KM2" = round(stations$AREA_KM2, 6),
                            "PERIMETER_KM" = stations$PER_KM,
                            "STRATUM" = stations$stratum,
                            "STATIONID" = stations$ID,
                            "CENTER_LAT" = stations$Lat,
                            "CENTER_LONG" = stations$Lon,
                            check.names = FALSE)

stations_2023[, c("NORTH_LAT", "SOUTH_LAT", "EAST_LONG", "WEST_LONG")] <-
  goa_grid_2021[match(stations_2023$STATIONID, goa_grid_2021$STATIONID),
                c("NORTH_LAT", "SOUTH_LAT", "EAST_LONG", "WEST_LONG")]

stations[, names(stations_2023)] <- stations_2023
stations[, !names(stations) %in% names(stations_2023)] <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Which parts of the goa grid are not supported by the bathymetry raster
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
leftovers <- terra::erase(x = goa_grid, y = stations)
leftovers[, "unk_area_km2"] <- terra::expanse(x = leftovers) / 1e6
leftovers[, !names(leftovers) %in% c("ID", "STRATUM", "unk_area_km2")] <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_strata_2023 <-
  data.frame("SURVEY" = "GOA",
             "STRATUM" = strata_list$stratum,
             "AREA" = strata_list$AREA_KM2,
             "PERIMETER" = strata_list$PER_KM,
             "INPFC_AREA" = strata_list$manage_area,
             "MIN_DEPTH" = strata_list$lower_depth_m,
             "MAX_DEPTH" = strata_list$upper_depth_m,
             "DESCRIPTION" = NA,
             "SUMMARY_AREA" = NA,
             "SUMMARY_DEPTH"= NA,
             "SUMMARY_AREA_DEPTH" = NA,
             "REGULATORY_AREA_NAME" =
               sapply(X = strata_list$manage_area,
                      FUN = function(x) switch(x,
                                               "Shumagin" = "WESTERN GOA",
                                               "Chirikof" = "CENTRAL GOA",
                                               "Kodiak" = "CENTRAL GOA",
                                               "Yakutat" = "EASTERN GOA",
                                               "Southeastern" = "EASTERN GOA")),
             "STRATUM_TYPE" = NA)

strata_list[, names(goa_strata_2023)] <- goa_strata_2023
strata_list[, !names(strata_list) %in% names(goa_strata_2023)] <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save results ----
##  Create export directory and export
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/processed_shapefiles/"))
  dir.create("data/GOA/processed_shapefiles/")

terra::writeVector(x = strata_list,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "goa_strata_2023.shp"),
                   overwrite = TRUE)
terra::writeVector(x = stations,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "goa_stations_2023.shp"),
                   overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usethis::use_data(stations_2023, overwrite = TRUE, internal = F)
aea_crs <- terra::crs(bathy)
usethis::use_data(aea_crs, overwrite = TRUE)
depth_mods_2023 <- depth_mods
usethis::use_data(depth_mods_2023, overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = "data/GOA/processed_shapefiles/updated_strata.pdf",
    width = 8, height = 6, onefile = TRUE)
for (iarea in unique(depth_mods$manage_area)) { ## Loop over area -- start

  ## temporary objects
  n_strata <- with(depth_mods, table(manage_area))[iarea]
  temp_strata <- depth_mods$stratum[depth_mods$manage_area == iarea]
  temp_area <- terra::mask(x = stations, mask = nmfs[nmfs$area_name == iarea] )

  ## Base layer
  plot(temp_area, axes = F, col = "white", border = F)

  ## Plot stations color coded by strata
  for (istratum in 1:n_strata) {
    plot(stations[stations$STRATUM == temp_strata[istratum]],
         col =   c(RColorBrewer::brewer.pal("Spectral",
                                            n = n_strata - 1),
                   "purple")[istratum],
         border = F, add = TRUE)
  }

  ## Plot areas not covered by bathymery
  plot(terra::mask(x = leftovers,
                   mask = nmfs[nmfs$area_name == iarea]),
       col = "black", add = TRUE, border = F)

  ## Untrawlable areas
  plot(terra::mask(x = goa_grid_untrawl,
                   mask = nmfs[nmfs$area_name == iarea]),
       col = rgb(0, 0, 0, 0.5), add = TRUE, border = FALSE)

  ## Land
  plot(ak_land, add = TRUE, col = "tan", border = T, lwd = 0.5)
  plot(ca_land, add = TRUE, col = "tan", border = T, lwd = 0.5)

  ## Strata
  plot(strata_list[strata_list$INPFC_AREA == iarea], lwd = 0.05, add = TRUE)

  ## Legend
  legend_labels <- with(subset(depth_mods, manage_area == iarea),
                        paste0(lower_depth_m, " - ", upper_depth_m, " m"))
  legend_labels <- c(legend_labels, "Untrawlable", "Not Defined")

  legend(c("Shumagin" = "topleft", "Chirikof" = "topleft",
           "Kodiak" = "bottomright", "Yakutat" = "bottom",
           "Southeastern" = "bottomleft")[iarea],
         legend = legend_labels,
         title = "Stratum Legend",
         fill = c(RColorBrewer::brewer.pal("Spectral", n = n_strata - 1),
                  "purple", rgb(0, 0, 0, 0.5), "black"))
  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}  ## Loop over area -- end
dev.off()
