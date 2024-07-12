##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Restratification Adjustments
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Workflow to take the new stratum boundaries and create new
##                strata polygons and intersect its boundaries with the 5 km
##                GOA grid. Then within each 5km grid cell, calculate the
##                total area and perimeter of the stratum component.
##
##                When importing and manipulating spatial objects, make sure
##                the projection is aligned with the bathymetry raster for
##                consistency.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(RColorBrewer)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import stratum depth boundaries ----
##   Import depth
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
depth_mods <-
  read.csv(file = "data/GOA/strata_boundaries/depth_modifications_2025.csv")
bathy <-
  terra::rast("//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/goa_bathy")
# terra::rast("C:/Users/zack.oyafuso/Desktop/goa_bathy/")

old_goa_strata <-
  terra::vect(x = "data/GOA/shapefiles_from_GDrive/goa_strata.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import 2021 stations ----
##   First make sure the `goa_grid` is in the same projection as the `bathy`
##   raster. Next, dissolve the inner bondaries of the `goa_grid` polygons
##   to get an outline of the GOA domain (`goa_domain`). Then, reproject
##   `goa_domain` to latlon.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/goagrid.shp")
goa_grid <- terra::project(x = goa_grid, y = bathy)
goa_domain <- terra::aggregate(x = goa_grid)

latlon_crs <-
  terra::crs(x = terra::vect(x = paste0("data/GOA/shapefiles_from_GDrive/",
                                        "GOAdissolved.shp")))
goa_domain_latlon <- terra::project(x = goa_domain, y = latlon_crs)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Land Polys ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ak_land <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/AKland.shp")
ak_land <- terra::project(x = ak_land, y = bathy)
ak_land <- terra::aggregate(x = ak_land)

ca_land <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/canada_dcw.shp")
ca_land <- ca_land[ca_land$POPYADMIN %in% c("BRITISH COLUMBIA",
                                            "YUKON TERRITORY",
                                            "NORTHWEST TERRITORIES") ]
ca_land <- terra::project(x = ca_land, y = bathy)
ca_land <- terra::aggregate(x = ca_land)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import NMFS Areas ----
##   Management area is also a stratum variable. Reproject to the same
##   projection as the `bathy` raster. NMFS statistical areas 650 and 659 are
##   part of the Southeastern section of the goa domain, so we aggregate those
##   two areas together.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nmfs <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/GOA_Shapes.shp")
nmfs <- terra::project(x = nmfs, y = bathy)
nmfs$AREA_NAME <- c("Southeast Outside", "Southeast Outside", "Shumagin",
                    "Chirikof", "Kodiak", "West Yakutat", NA)
nmfs$REP_AREA <- c(650, 650, 610, 620, 630, 640, 649, NA)

nmfs <- terra::aggregate(x = nmfs, by = "REP_AREA")
nmfs <- subset(x = nmfs,
               select = c("AREA_NAME", "REP_AREA"),
               subset = !is.na(x = nmfs$AREA_NAME))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create stratum polygons ----
##   For each management area, create new strata based on depth specifications
##   and append to strata_list
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_list <- strata_agg_list <- list()
for (idistrict in unique(x = depth_mods$manage_area)) { ## Loop over area --st.

  ## Mask bathymetry raster to just the management area and goa_domain
  district_outline <- nmfs[nmfs$AREA_NAME  == idistrict]

  district_bathy <- terra::mask(x = bathy,
                                mask = goa_domain)
  district_bathy <- terra::mask(x = district_bathy,
                                mask = district_outline)
  district_bathy <- terra::crop(x = district_bathy,
                                y = district_outline)

  ## Define modified stratum depth boundaries
  depth_boundary <- subset(x = depth_mods,
                           subset = manage_area == idistrict,
                           select = c("lower_depth_m", "upper_depth_m"))

  ## Discretize the `bathy`` raster: Define each raster cell based on the
  ## defined stratum depth boundaries in `depth_mods` and create an arbitrary
  ## integer label for each stratum.
  terra::values(district_bathy) <-
    as.integer(as.character(cut(x = terra::values(x = district_bathy),
                                breaks = c(0, depth_boundary$upper_depth_m),
                                labels = 1:nrow(x = depth_boundary)) ))

  ## Convert discretized raster to polygon based on those discrete values
  strata_poly <- terra::as.polygons(x = district_bathy)
  strata_poly_disagg <- terra::disagg(x = strata_poly)
  strata_poly_disagg$area <- terra::expanse(x = strata_poly_disagg) / 1e6

  ## For each polygon, calculate the adjacent polygons. argument type == "rook"
  ## excludes polygons that touch at a single node.
  nearest_poly <- terra::adjacent(x = strata_poly_disagg, type = "rook")

  for (i in which(strata_poly_disagg$area < 5)) {
    temp_speck <- strata_poly_disagg[i, ]
    adj_polys <- nearest_poly[nearest_poly[, 2] == i, 1]

    if (length(x = adj_polys) != 0) {
      adj_poly <- adj_polys[which.max(x = strata_poly_disagg$area[adj_polys])]
      strata_poly_disagg$GOA_bathy[i] <- strata_poly_disagg$GOA_bathy[adj_poly]
    }
  }

  strata_poly_agg <- terra::aggregate(x = strata_poly_disagg,
                                      by = "GOA_bathy",
                                      fun = "sum",
                                      count = F)

  strata_poly_disagg <- terra::disagg(x = strata_poly_agg)
  strata_poly_disagg$area <- terra::expanse(x = strata_poly_disagg,
                                            unit = "km")

  nearest_poly <- terra::adjacent(x = strata_poly_disagg, type = "intersects")

  for (i in which(strata_poly_disagg$area < 5)) {
    temp_speck <- strata_poly_disagg[i, ]
    adj_polys <- nearest_poly[nearest_poly[, 2] == i, 1]

    if (length(adj_polys) != 0) {
      adj_poly <- adj_polys[which.max(x = strata_poly_disagg$area[adj_polys])]
      strata_poly_disagg$GOA_bathy[i] <- strata_poly_disagg$GOA_bathy[adj_poly]
    }
  }

  strata_poly_agg <- terra::aggregate(x = strata_poly_disagg,
                                      by = "GOA_bathy",
                                      fun = "sum",
                                      count = F)

  ## Create dataframe of stratum information
  strata_poly[, names(depth_mods)] <- strata_poly_agg[, names(depth_mods)] <-
    subset(x = depth_mods, subset = manage_area == idistrict)

  ## Calculate the total area and perimeter of the strata.
  strata_poly$AREA_KM2 <- terra::expanse(x = strata_poly, unit = "km")
  strata_poly$PER_KM <- terra::perim(x = strata_poly) / 1000

  strata_poly_agg$AREA_KM2 <- terra::expanse(x = strata_poly_agg, unit = "km")
  strata_poly_agg$PER_KM <- terra::perim(x = strata_poly_agg) / 1000

  ## Append to strata_list
  strata_list <- c(strata_list, list(strata_poly))
  strata_agg_list <- c(strata_agg_list, list(strata_poly_agg))

  print(paste("Finished with the", idistrict, "region"))
} ## Loop over district -- end
rm(idistrict, strata_poly, strata_poly_agg, district_outline, district_bathy)

##   Merge strata into one object
strata_list <- do.call(what = rbind, args = strata_list)
strata_agg_list <- do.call(what = rbind, args = strata_agg_list)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create 5x5 km survey grid ----
##   Take the goa grid and merge all historical stations wihtin 5x5 km
##   cell. This essentially resets the stratum data in each grid cell.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_full_grid <- terra::aggregate(x = goa_grid,
                                  by = "ID",
                                  fun = sum)
goa_full_grid <- goa_full_grid[, c("ID", "agg_AREA_KM2")]
names(goa_full_grid) <- c("ID", "AREA_KM2")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create new stations ----
##   Intersect the new stratum polygons with the 5x5 km grid to create new
##   station polygons. This function takes a few minutes. Then calculate the
##   area, perimeter, and location of each new station.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations <- terra::intersect(x = goa_full_grid,
                             y = strata_list[, c("manage_area", "stratum")])

stations$PER_KM <- terra::perim(x = stations) / 1e3
stations$AREA_KM2 <- terra::expanse(x = stations) / 1e6

stations[, c("E_m", "N_m")] <-
  terra::geom(terra::centroids(x = stations, inside = TRUE))[, c("x", "y")]
stations[, c("Lon", "Lat")] <-
  terra::geom(terra::centroids(x = terra::project(x = stations,
                                                  y = latlon_crs),
                               inside = TRUE))[, c("x", "y")]

stations_agg <- terra::intersect(x = goa_full_grid,
                                 y = strata_agg_list[, c("manage_area", "stratum")])

stations_agg$PER_KM <- terra::perim(x = stations_agg) / 1e3
stations_agg$AREA_KM2 <- terra::expanse(x = stations_agg) / 1e6

stations_agg[, c("E_m", "N_m")] <-
  terra::geom(terra::centroids(x = stations_agg, inside = TRUE))[, c("x", "y")]
stations_agg[, c("Lon", "Lat")] <-
  terra::geom(terra::centroids(x = terra::project(x = stations_agg,
                                                  y = latlon_crs),
                               inside = TRUE))[, c("x", "y")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format stations ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goagrid_ids <-
  (max(goa_grid$GOAGRID_ID) + 1):
  (max(goa_grid$GOAGRID_ID) +  nrow(stations))

stations_2025 <- data.frame("GOAGRID#" = goagrid_ids,
                            "GOAGRID_ID" = goagrid_ids,
                            "AREA_KM2" = round(stations$AREA_KM2, 6),
                            "PERIMETER_KM" = stations$PER_KM,
                            "STRATUM" = stations$stratum,
                            "STATIONID" = stations$ID,
                            "CENTER_LAT" = stations$Lat,
                            "CENTER_LONG" = stations$Lon,
                            check.names = FALSE)

stations[, names(stations_2025)] <- stations_2025
stations[, !names(stations) %in% names(stations_2025)] <- NULL

goagrid_ids <-
  (max(goa_grid$GOAGRID_ID) + 1):
  (max(goa_grid$GOAGRID_ID) +  nrow(stations_agg))

stations_agg_2025 <- data.frame("GOAGRID#" = goagrid_ids,
                                "GOAGRID_ID" = goagrid_ids,
                                "AREA_KM2" = round(stations_agg$AREA_KM2, 6),
                                "PERIMETER_KM" = stations_agg$PER_KM,
                                "STRATUM" = stations_agg$stratum,
                                "STATIONID" = stations_agg$ID,
                                "CENTER_LAT" = stations_agg$Lat,
                                "CENTER_LONG" = stations_agg$Lon,
                                check.names = FALSE)

stations_agg[, names(stations_agg_2025)] <- stations_agg_2025
stations_agg[, !names(stations_agg) %in% names(stations_agg_2025)] <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Areas not defined by bathy ----
##   Which parts of the goa grid are not supported by the bathymetry raster?
##   This is mostly for plotting purposes, still need to figure out how to
##   fill out these parts of the goa domain.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# leftovers <- terra::erase(x = goa_grid, y = stations)
# leftovers[, "unk_area_km2"] <- terra::expanse(x = leftovers) / 1e6
# leftovers[, !names(leftovers) %in% c("ID", "STRATUM", "unk_area_km2")] <- NULL

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format strata ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_strata_2025 <-
  data.frame("SURVEY" = "GOA",
             "SURVEY_DEFINITION_ID" = 47,
             "DESIGN_YEAR" = 2025,
             "AREA_ID" = strata_list$stratum,
             "AREA_TYPE" = "STRATUM",
             "AREA_NAME" = paste0(strata_list$manage_area, ", ",
                                  strata_list$lower_depth_m, "-",
                                  strata_list$upper_depth_m, " m"),
             "DESCRIPTION" = paste0(strata_list$manage_area, ", ",
                                  strata_list$lower_depth_m, "-",
                                  strata_list$upper_depth_m, " m"),
             "AREA_KM2" = strata_list$AREA_KM2,
             "MIN_DEPTH" = strata_list$lower_depth_m,
             "MAX_DEPTH" = strata_list$upper_depth_m,
             CRS = NA)

goa_strata_agg_2025 <-
  data.frame("SURVEY" = "GOA",
             "SURVEY_DEFINITION_ID" = 47,
             "DESIGN_YEAR" = 2025,
             "AREA_ID" = strata_agg_list$stratum,
             "AREA_TYPE" = "STRATUM",
             "AREA_NAME" = paste0(strata_agg_list$manage_area, ", ",
                                  strata_agg_list$lower_depth_m, "-",
                                  strata_agg_list$upper_depth_m, " m"),
             "DESCRIPTION" = paste0(strata_agg_list$manage_area, ", ",
                                    strata_agg_list$lower_depth_m, "-",
                                    strata_agg_list$upper_depth_m, " m"),
             "AREA_KM2" = strata_agg_list$AREA_KM2,
             "MIN_DEPTH" = strata_agg_list$lower_depth_m,
             "MAX_DEPTH" = strata_agg_list$upper_depth_m,
             CRS = NA)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result_areas <-
  merge(x = subset(x = as.data.frame(strata_list[, c("STRATUM", "AREA")]),
                   subset = STRATUM < 500),
        y = subset(x = as.data.frame(strata_agg_list[, c("STRATUM", "AREA")]),
                   subset = STRATUM < 500),
        suffixes = c("_ORIG", "_CLEAN"),
        by = "STRATUM")

result_areas$PERC_DIFF <-
  round(x = 100 * (result_areas$AREA_CLEAN - result_areas$AREA_ORIG) /
          result_areas$AREA_ORIG ,
        digits = 1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save shapefiles ----
##  Create export directory and export
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/processed_shapefiles/"))
  dir.create("data/GOA/processed_shapefiles/")

terra::writeVector(x = strata_agg_list,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "goa_strata_2025.shp"),
                   overwrite = TRUE)

terra::writeVector(x = stations_agg,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "goa_stations_2025.shp"),
                   overwrite = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot, finally ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = "data/GOA/processed_shapefiles/updated_strata.pdf",
    width = 8, height = 6, onefile = TRUE)
for (iarea in unique(depth_mods$manage_area)) { ## Loop over area -- start

  ## temporary objects
  n_strata <- with(depth_mods, table(manage_area))[iarea]
  temp_strata <- depth_mods$stratum[depth_mods$manage_area == iarea]
  temp_area <- terra::mask(x = stations, mask = nmfs[nmfs$AREA_NAME  == iarea] )

  ## Base layer
  plot(temp_area, axes = F, col = "white", border = F, mar = c(0,0,0,0))

  ## Plot stations color coded by strata
  for (istratum in 1:n_strata) {
    plot(stations[stations$STRATUM == temp_strata[istratum]],
         col =   c(RColorBrewer::brewer.pal("Set1",
                                            n = n_strata - 1),
                   "cyan")[istratum],
         border = F, add = TRUE)
  }

  ## Plot areas not covered by bathymery
  # plot(terra::mask(x = leftovers,
  #                  mask = nmfs[nmfs$AREA_NAME  == iarea]),
  #      col = "black", add = TRUE, border = F)

  ## Untrawlable areas
  # plot(terra::mask(x = UT_stations21[UT_stations21],
  #                  mask = nmfs[nmfs$AREA_NAME  == iarea]),
  #      col = rgb(0, 0, 0, 0.5), add = TRUE, border = FALSE)

  ## Land
  plot(ak_land, add = TRUE, col = "tan", border = T, lwd = 0.1)
  plot(ca_land, add = TRUE, col = "tan", border = T, lwd = 0.1)

  ## Strata
  plot(strata_list[strata_list$manage_area == iarea], lwd = 0.05, add = TRUE)

  ## Legend
  legend_labels <- with(subset(depth_mods, manage_area == iarea),
                        paste0(lower_depth_m, " - ", upper_depth_m, " m"))
  legend_labels <- c(legend_labels, "Untrawlable", "Not Defined")

  legend(c("Shumagin" = "topleft", "Chirikof" = "topleft",
           "Kodiak" = "bottomright", "West Yakutat" = "bottom",
           "Southeast Outside" = "bottomleft")[iarea],
         legend = legend_labels,
         title = "Stratum Legend",
         fill = c(RColorBrewer::brewer.pal("Set1", n = n_strata - 1),
                  "cyan", rgb(0, 0, 0, 0.5), "black"))
  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}  ## Loop over area -- end
dev.off()

## Open pdf file
if (file.exists("data/GOA/processed_shapefiles/updated_strata.pdf")) {
  system('open  data/GOA/processed_shapefiles/updated_strata.pdf')
}

pdf(file = "data/GOA/processed_shapefiles/updated_strata_agg.pdf",
    width = 8, height = 6, onefile = TRUE)
for (iarea in unique(depth_mods$manage_area)) { ## Loop over area -- start

  ## temporary objects
  n_strata <- with(depth_mods, table(manage_area))[iarea]
  temp_strata <- depth_mods$stratum[depth_mods$manage_area == iarea]
  temp_area <- terra::mask(x = stations_agg,
                           mask = nmfs[nmfs$AREA_NAME  == iarea] )

  ## Base layer
  plot(temp_area, axes = F, col = "white", border = F, mar = c(0,0,0,0))

  ## Plot stations color coded by strata
  for (istratum in 1:n_strata) {
    plot(stations_agg[stations_agg$STRATUM == temp_strata[istratum]],
         col =   c(RColorBrewer::brewer.pal("Set1",
                                            n = n_strata - 1),
                   "cyan")[istratum],
         border = F, add = TRUE)
  }

  ## Plot areas not covered by bathymery
  # plot(terra::mask(x = leftovers,
  #                  mask = nmfs[nmfs$AREA_NAME  == iarea]),
  #      col = "black", add = TRUE, border = F)

  ## Untrawlable areas
  # plot(terra::mask(x = UT_stations21[UT_stations21],
  #                  mask = nmfs[nmfs$AREA_NAME  == iarea]),
  #      col = rgb(0, 0, 0, 0.5), add = TRUE, border = FALSE)

  ## Land
  plot(ak_land, add = TRUE, col = "tan", border = T, lwd = 0.1)
  plot(ca_land, add = TRUE, col = "tan", border = T, lwd = 0.1)

  ## Strata
  plot(strata_agg_list[strata_agg_list$manage_area == iarea],
       lwd = 0.05, add = TRUE)

  ## Legend
  legend_labels <- with(subset(depth_mods, manage_area == iarea),
                        paste0(lower_depth_m, " - ", upper_depth_m, " m"))
  legend_labels <- c(legend_labels, "Untrawlable", "Not Defined")

  legend(c("Shumagin" = "topleft", "Chirikof" = "topleft",
           "Kodiak" = "bottomright", "West Yakutat" = "bottom",
           "Southeast Outside" = "bottomleft")[iarea],
         legend = legend_labels,
         title = "Stratum Legend",
         fill = c(RColorBrewer::brewer.pal("Set1", n = n_strata - 1),
                  "cyan", rgb(0, 0, 0, 0.5), "black"))
  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}  ## Loop over area -- end
dev.off()

## Open pdf file
if (file.exists("data/GOA/processed_shapefiles/updated_strata.pdf")) {
  system('open  data/GOA/processed_shapefiles/updated_strata.pdf')
}

## Open pdf file
if (file.exists("data/GOA/processed_shapefiles/updated_strata_agg.pdf")) {
  system('open  data/GOA/processed_shapefiles/updated_strata_agg.pdf')
}
