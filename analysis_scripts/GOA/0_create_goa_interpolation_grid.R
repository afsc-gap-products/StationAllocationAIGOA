##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA interpolation grid under 2025 GOA survey footprint
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Import libraries
library(akgfmaps); library(terra); library(sf)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create 2-nmi grid within the 2025 GOA footprint
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Import goa strata <= 700 m
goa_strata_2025 <-
  akgfmaps::get_base_layers(
    select.region = "goa",
    design.year = 2025,
    set.crs = "+proj=utm +zone=5 +units=km"
  )$survey.strata |>
  ## Remove strata in the 700 - 1000 m depth zone
  subset(subset = STRATUM < 500)

## Dissolve inner boundaries
goa_footprint_2025 <- goa_strata_2025 |> terra::vect() |> terra::aggregate()

## Create a 2-nautical-mile (3.704 km) rectangular grid that overlaps with the
## 2025 GOA survey footprint
goa_grid <-
  sf::st_make_grid(x = goa_footprint_2025,
                   square = T,
                   cellsize = 3.704) |>
  terra::vect() |>
  ## Intersect with the footprint
  terra::intersect(goa_footprint_2025)

## Pull the centroids of each gridcell
cell_centroids <- terra::centroids(x = goa_grid, inside = TRUE) |>
  terra::intersect(y = terra::vect(goa_strata_2025)[, "STRATUM"] )

## Attach the centroids (in UTMs)
goa_grid_df <- data.frame(stratum = cell_centroids$STRATUM,
                          cell_centroids |> terra::crds())
names(x = goa_grid_df) <- toupper(x = names(x = goa_grid_df))

## Project sampling unit centroids in lat/lon
goa_grid_df[, c("LON", "LAT")] <-
  terra::project(x = goa_grid, "+proj=longlat +datum=WGS84") |>
  terra::centroids() |>
  terra::crds()

## Calculate area of each sampling unit
goa_grid_df$AREA_KM2 <- terra::expanse(x = goa_grid) / 1e6

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Write to file
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(x = goa_grid_df,
          file = paste0("data/GOA/sdmtmb_data/",
                        "goa_2025_interpolation_grid.csv"),
          row.names = F)
