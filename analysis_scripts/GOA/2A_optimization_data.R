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
##   Load the true density, true index, and spatial domain dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(file = "data/GOA/prednll_VAST_models.RData")
grid_goa <- read.csv(file = "data/GOA/grid_goa.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used throughout all scripts
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Years to use
year_set <- 1996:2021
years_included <- c(1, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)
n_years <- length(years_included)

## Number of sampling grids
n_cells <- nrow(grid_goa)

## Scientific and common names used in optimization
common_names_all <- pred_jnll$spp_name

ns_all <- length(common_names_all)

spp_idx_opt <- c(25, 14, #cods
                 1, 7, 18, 12, 24, 5, 15, #flatfishes
                 16, 4, 23, 6, 13, 22 #rockfish types
)
common_names_opt <- common_names_all[spp_idx_opt]
ns_opt <- length(common_names_opt)

## Scientific and common names not used in optimization, but evaluated
## when simulating surveys
spp_idx_eval <- (1:ns_all)[-spp_idx_opt]
common_names_eval <- common_names_all[spp_idx_eval]
ns_eval <- length(common_names_eval)

## Specify Management Districts
districts <- data.frame("reg_area" = c("WRA", "CRA", "CRA", "ERA", "ERA"),
                        "district" = c("West", "Chirikof", "Kodiak",
                                       "Yakutat", "Southeast"),
                        "domainvalue" = 1:5,
                        "W_lon" = c(-170, -159, -154, -147, -140),
                        "E_lon" = c(-159, -154, -147, -140, -132))

district_vals <- cut(x = grid_goa$Lon,
                     breaks = c(-170, -159, -154, -147, -140, -132),
                     labels = 1:5)
districts[, c("W_UTM", "E_UTM")] <-
  do.call(rbind,tapply(X = grid_goa$E_km,
                       INDEX = district_vals,
                       FUN = range) )

n_districts <- nrow(districts)

## ranges of the spatial domain for plotting
x_range <- diff(range(grid_goa$E_km))
y_range <- diff(range(grid_goa$N_km))

## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save Data, Species densities separately
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(list = c("districts", "district_vals", "n_districts",
              "ns_all", "ns_eval", "ns_opt", "lonlat_crs", "utm_crs",
              "x_range", "y_range", "n_cells",
              "common_names_all", "common_names_eval", "common_names_opt",
              "spp_idx_eval", "spp_idx_opt",
              "year_set", "years_included", "n_years"),
     file = "data/GOA/optimization_data.RData")
