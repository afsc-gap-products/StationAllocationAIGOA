##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       EFH Bathymetry Layer
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create full bathymetry raster from patched rasters
##                in the RACE_EFH_Variables folder
##
##                Many of the rasters come from this RACE_EFH_Variables folder
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Absolute paths of rastesr used. You need access to the RACE_EFH_Variables
#    and VPN. The beginning of the path may be different across machines
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/goagrid.shp")
goa_domain <- terra::aggregate(x = goa_grid)

race_efh <- "Y:/RACE_EFH_variables/"
aigoa <- terra::rast(x = paste0(race_efh,
                                "Variables/ArcMap/AIGOAbathy_Combined/",
                                "aigoa_bathp1c/w001001.adf"))
aigoa <- terra::mask(terra::crop(x = aigoa, y = goa_domain),
                     mask = goa_domain)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   aigoa has some missing areas, mostly along the offshore edge
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(goa_domain, col = "red", border = F)
plot(aigoa, col = "black", add = TRUE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Regional bathys have more coverage of the domain
##   For each regional bathy, we crop and mask to the goa domain, trim the
##   empty spaces at the margins, and then call it region_patch
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ai_bathy <-
  terra::rast(x = paste0(race_efh, "ZimmBathymetry/AI/ai_bathy/w001001.adf"))
ai_bathy <- terra::mask(terra::crop(x = ai_bathy, y = goa_domain),
                        mask = goa_domain)
ai_patch <- terra::trim(ai_bathy)

wGOA <-
  terra::rast(x = paste0(race_efh, "ZimmBathymetry/wGOA/wgoa_bathy/w001001.adf"))
wGOA <- terra::mask(terra::crop(x = wGOA, y = goa_domain),
                    mask = goa_domain)
wGOA_patch <- terra::trim(wGOA)

cGOA <-
  terra::rast(x = paste0(race_efh, "ZimmBathymetry/CGOA/cgoa_bathy/w001001.adf"))
cGOA <- terra::mask(terra::crop(x = cGOA, y = goa_domain),
                    mask = goa_domain)
cGOA_patch <- terra::trim(cGOA)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Calculate the x-y locations of aigoa raster cells with no bathy info
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
missing_pts <- terra::extract(x = aigoa, y = goa_domain,
                              cells = TRUE, xy = TRUE)
missing_pts <-  missing_pts[is.na(missing_pts$AIGOA_ba), c("x", "y")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Extract raster values of the region_patch rasters at those missing points
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ai_patch_vals <- terra::extract(x = ai_patch,
                                y = missing_pts[, c("x", "y")],
                                cells = TRUE, xy = TRUE)
ai_patch_vals <- ai_patch_vals[!is.na(ai_patch_vals$cell), ]

wGOA_patch_vals <- terra::extract(x = wGOA_patch,
                                  y = missing_pts[, c("x", "y")],
                                  cells = TRUE, xy = TRUE)
wGOA_patch_vals <- wGOA_patch_vals[!is.na(wGOA_patch_vals$cell), ]

cGOA_patch_vals <- terra::extract(x = cGOA_patch,
                                  y = missing_pts[, c("x", "y")],
                                  cells = TRUE, xy = TRUE)
cGOA_patch_vals <- cGOA_patch_vals[!is.na(cGOA_patch_vals$cell), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set all values in the region_patch raster that overlaps with the
##   aigoa main bathy layer AND any depths < 100 to NA. The second part makes
##   sure that we are only, for now, dealing with the outer edge of the domain
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
values(ai_patch)[values(ai_patch) < 100] <- NA
values(ai_patch)[-ai_patch_vals$cell] <- NA

values(wGOA_patch)[values(wGOA_patch) < 100] <- NA
values(wGOA_patch)[-wGOA_patch_vals$cell] <- NA

values(cGOA_patch)[values(cGOA_patch) < 100] <- NA
values(cGOA_patch)[-cGOA_patch_vals$cell] <- NA

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   merge (or resample) the region_patch layers with the mian aigoa layer
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
aigoa_patched <- do.call(what = terra::merge,
                         args = list(aigoa, ai_patch, wGOA_patch, cGOA_patch))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(goa_domain, col = "red", border = F)
plot(aigoa_patched, add = TRUE, col = "black")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Split Bathy layer into smaller saveable pieces
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!dir.exists("data/GOA/processed_rasters/"))
  dir.create("data/GOA/processed_rasters/")

bathy_split <- terra::makeTiles(
  x = aigoa_patched,
  y = terra::rast(x = matrix(data = 1:9, nrow = 3, ncol = 3),
                  extent = terra::ext(aigoa_patched)),
  filename = "data/GOA/processed_rasters/aigoa_.tif")
